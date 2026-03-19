# Plan: TUV-x C++ Solver Library Rewrite

**TL;DR**: Create a new standalone C++ library (new repo, e.g., `tuv-x-cpp`) that replaces the Fortran-based TUV-x with a high-performance photolysis rate constant calculator. The library provides a pure programmatic C++ API — no configuration files — with composable lambda-based transforms for cross-sections, quantum yields, and dose rates. Data structures use template policies to support both CPU SIMD (multi-column vectorization) and GPU (CUDA/HIP) execution. Delta Eddington solver ships first; Discrete Ordinate follows. MUSICA wraps this library and handles YAML/JSON configuration and Python/JS/Julia bindings.

## Phase 0: Project Scaffolding

### SI Units Policy

The TUV-x C++ library performs **no unit conversions** whatsoever. All input data — whether provided programmatically or read from data files — must already be in SI units before use:

- **Wavelength**: meters (m), not nanometers
- **Distance/altitude**: meters (m), not kilometers
- **Angles**: radians (rad), not degrees
- **Temperature**: Kelvin (K)
- **Pressure**: Pascals (Pa)
- **Time**: seconds (s)
- **Cross-sections**: m² (not cm²)
- **Number density**: molecules/m³ (not molecules/cm³)

All outputs are in SI units. Any data files (e.g., legacy NetCDF files with wavelengths in nm) must be pre-converted to SI units before being consumed by TUV-x. The responsibility for unit conversion lies with the caller or the data preparation pipeline — not with the library itself.

- [ ] 1. **Create new repository** (`tuv-x-cpp` or chosen name) with CMake build system (CMake 3.21+ to match MUSICA). Configure languages `CXX` with optional `CUDA`/`HIP`. Set up GitHub Actions CI with strict quality gates: >95% unit test coverage enforced on PRs, valgrind memcheck on all tests, multi-compiler matrix (GCC, Clang, MSVC, Intel, NVHPC) across Linux/macOS/Windows, clang-tidy static analysis as a PR gate, and auto-generated formatting PRs via clang-format. Enforce legible naming conventions (no cryptic abbreviations).

- [ ] 2. **Migrate reusable C++ code from current repo.** Port the following already-complete implementations from `include/tuvx/` into the new repo's `include/` and `src/`:
   - `Array1D<T>` — new implementation (see below), modeled after `Array2D`/`Array3D`
   - `Array2D<T>`, `Array3D<T>` — `include/tuvx/util/array2d.hpp`, `include/tuvx/util/array3d.hpp`
   - `Grid<ArrayPolicy>` — `include/tuvx/grid.hpp`
   - `Profile<ArrayPolicy>` — `include/tuvx/profile.hpp`
   - `RadiatorState<ArrayPolicy>` — `include/tuvx/radiative_transfer/radiator.hpp`
   - `RadiationField`, `RadiationFieldComponents` — `include/tuvx/radiative_transfer/radiation_field.hpp`
   - `TridiagonalMatrix<T>`, `Solve()` — `include/tuvx/linear_algebra/linear_algebra.hpp`
   - Existing GTest tests for all of the above
   - Existing Google Benchmark for tridiagonal solver — `benchmark/benchmark_tridiagonal_solver.cpp`

- [ ] 3. **Define template policy system for CPU/GPU portability.** Extend the existing `ArrayPolicy` pattern into a coherent execution/memory policy. The solver algorithms are written once, templated on a single `ArrayPolicy` — from the algorithm's perspective, array access is identical regardless of backing storage:

   - `ArrayPolicy` controls: memory allocation (where), data layout (how), and element access (syntax)
   - `HostArrayPolicy` — `std::vector`-backed; data layout chosen for compiler auto-vectorization (SIMD)
   - `DeviceArrayPolicy` — CUDA/HIP device memory; data layout chosen for coalesced access; same element access syntax via `__host__ __device__` accessors
   - Solver functions are templated on `ArrayPolicy` and written exactly once — no `#ifdef` branching, no runtime dispatch
   - All core types (`Grid`, `Profile`, `RadiatorState`, `RadiationField`, `TridiagonalMatrix`) are already templated on `ArrayPolicy` in the existing codebase; this pattern is preserved and extended
   - **`Array1D<T>`**: Create a new 1D array type following the same pattern as `Array2D<T>` and `Array3D<T>`. This replaces all uses of `std::vector` in APIs that interact with `Grid`/`Profile`/policy-backed data, ensuring consistency when data is device-resident. `Array1D` must support the same `ArrayPolicy`-controlled memory allocation and access as the 2D/3D types. **No `std::vector` should appear in any function signature that also takes policy-backed types** — use `Array1D` instead

### Bulk Operations API (MICM-Consistent)

Raw element access (`array[i][j]`, `array[i][j][k]`) is available for setup and debugging, but **solver hot paths must use the bulk operations API** for layout-independent parallelization. This API follows the same design as MICM's `Matrix`/`VectorMatrix` pattern (both are MUSICA components and should present a consistent interface to users).

The bulk operations API provides three core abstractions:

#### Element access — square-bracket syntax

All array types use `operator[]` chaining for logical element access, matching MICM's convention:
```cpp
Array1D<double> a(num_columns);
a[col] = 1.0;

Array2D<double> b(num_sections, num_columns);
b[section][col] = 2.0;

Array3D<double> c(num_wavelengths, num_layers, num_columns);
c[wl][layer][col] = 3.0;
```
The `operator[]` returns a proxy/view object whose own `operator[]` resolves to the next dimension. The underlying index computation is policy-defined — the square-bracket chain provides logical access without prescribing physical layout.

#### `ForEachRow` — layout-agnostic bulk iteration

Operates on all elements along the innermost logical dimension (columns) simultaneously. The policy implementation decides how to vectorize/parallelize:
```cpp
// Given a 2D array [sections × columns]:
array.ForEachRow(
    [](const double& a, double& b) { b = a * 2.0; },
    array.GetConstColumnView(0),   // read from column 0
    array.GetColumnView(1)         // write to column 1
);
```
`ForEachRow` takes a lambda plus any number of column views (or `Array1D` vectors). The lambda receives one element per argument per row. The policy controls iteration order and vectorization.

Column views provide layout-abstracted access to a single column's data across all rows:
- `GetColumnView(col)` — mutable view
- `GetConstColumnView(col)` — read-only view
- `GetRowVariable()` — temporary per-row storage (for intermediate calculations)

#### `ArrayPolicy::Function` — reusable policy-compiled operations

For performance-critical operations that are called repeatedly (e.g., every solver timestep), create a reusable function object. This allows the policy to pre-compute layout metadata and generate optimized iteration patterns:
```cpp
auto func = Array3D<double>::Function(
    [](auto&& arr) {
        auto tmp = arr.GetRowVariable();
        arr.ForEachRow(
            [](const double& a, const double& b, double& t)
            { t = a + b; },
            arr.GetConstColumnView(0),
            arr.GetConstColumnView(1),
            tmp);
        arr.ForEachRow(
            [](const double& t, double& c) { c = t * 2.0; },
            tmp,
            arr.GetColumnView(2));
    }, prototype_array);  // prototype provides type/dimension info

// Reuse with different arrays of the same column count:
func(array_a);  // OK — applies optimized operation
func(array_b);  // OK — different row count is fine
```
Column counts are validated at invocation time; row counts may differ from the prototype. This matches MICM's `MatrixPolicy::Function` factory exactly.

### Multi-Column Architecture

The solver always operates on a **set of columns simultaneously**. Every `Array1D`, `Array2D`, and `Array3D` in the solver pipeline carries a column dimension:

| Array type | Logical dimensions | Column semantics | Examples |
|-----------|-------------------|------------------|----------|
| `Array1D` | `(column)` | One scalar per column | solar zenith angle, surface albedo |
| `Array2D` | `(X, column)` | 1D data per column | grid mid-points/edges, profile values |
| `Array3D` | `(X, Y, column)` | 2D data per column | optical depth `(wavelength, layer, column)`, cross-section output `(wavelength, height, column)` |

These are *logical* dimensions — the underlying memory layout and dimension ordering is determined entirely by the concrete `ArrayPolicy` implementation. Solver code accesses elements by logical index using square-bracket syntax (e.g., `array[wavelength_idx][layer_idx][column_idx]`) and the policy maps that to the optimal physical layout.

Solver hot paths use the **bulk operations API** (`ForEachRow`, `ColumnView`, `Function`) rather than element-by-element access. Solver code is written as bulk operations on `ArrayXD` data — the same operation applies to every column. The concrete `ArrayPolicy` implementation decides how to parallelize across the column dimension:
- `HostArrayPolicy`: layout optimized for compiler auto-vectorization (SIMD) across columns
- `DeviceArrayPolicy`: layout optimized for coalesced GPU memory access

This means solver functions never loop over columns explicitly — they express element-wise or reduction operations via `ForEachRow` on the full `ArrayXD`, and the policy handles parallelization. Performance-critical operations (e.g., tridiagonal solve, delta-scaling) should use `ArrayPolicy::Function` to create reusable, policy-optimized function objects.

**Exception**: `DataReader` returns (`Array1D` for wavelengths, `Array2D` for parameter tables) do *not* have a column dimension — these are static reference data shared across all columns. The column dimension is introduced when this data is broadcast into per-column arrays during transform evaluation.


## Phase 1: Delta Eddington Solver

- [ ] 4. **Implement spherical geometry utilities.** Port the pseudo-spherical correction from `src/spherical_geometry.F90` and the slant optical depth calculation from `src/radiative_transfer/solver.F90`. These compute the effective solar path through curved atmospheric layers. Pure functions, no state.

- [ ] 5. **Implement Delta Eddington solver.** Replace the placeholder in `include/tuvx/radiative_transfer/solvers/delta_eddington.inl` with the actual algorithm from `src/radiative_transfer/solvers/delta_eddington.F90`. The implementation follows 6 steps documented in [issue #64](https://github.com/NCAR/tuv-x/issues/64):
   - **Step 1**: Delta-scale optical properties; compute gamma coefficients ($\gamma_1$ through $\gamma_4$, $\lambda$, $\Gamma$) — per layer, per wavelength, per column
   - **Step 2**: Compute solar source functions $C^+(\tau)$ and $C^-(\tau)$
   - **Step 3**: Assemble tridiagonal system coefficients ($e_1$ through $e_4$; $A_n$, $B_n$, $D_n$ matrix diagonals)
   - **Step 4**: Assemble right-hand side $E_n$ with surface albedo boundary condition
   - **Step 5**: Solve tridiagonal system using existing `TridiagonalMatrix::Solve()` (Thomas algorithm)
   - **Step 6**: Back-substitute $Y$ to compute fluxes and actinic flux components
   - All operations process batches of columns; the `ArrayPolicy` determines the optimal data layout for parallelization
   - All units in SI (meters, radians, seconds) — no unit conversions in solver code; all input data must be provided in SI units

- [ ] 6. **Regression test against Fortran solver.** Adapt the existing test harness in `test/regression/solvers/` to validate the new C++ solver against pre-computed Fortran reference outputs stored as binary/NetCDF files. This decouples testing from the Fortran build.

## Phase 2: Composable Transform System

- [ ] 7. **Design the transform type system.** Define a `Transform<Policy>` as a callable that maps atmospheric state to a 3D array `[wavelength × height × column]`:
   ```cpp
   using TransformFunc = std::function<void(
       const Grid<Policy>& wavelength_grid,
       const Grid<Policy>& altitude_grid,
       const Profile<Policy>& temperature,
       const Profile<Policy>& pressure,
       const Profile<Policy>& air_density,
       Array3D<T>& output  // [wavelength × height × column]
   )>;
   ```
   This is the fundamental signature for cross-sections, quantum yields, and spectral weights. The output includes the column dimension because temperature, pressure, and density vary per column, so the computed values differ per column.

- [ ] 8. **Implement composable transform primitives.** Each returns a `TransformFunc` (or equivalent callable). These are the building blocks:

   | Primitive | Replaces | Description |
   |-----------|----------|-------------|
   | `from_data(reader, interpolator)` | Base `cross_section_t` | Load tabulated data via pluggable reader; interpolate to model grid |
   | `temperature_interpolation(reader)` | `tint` | Linear interpolation between tabulated T-dependent data sets |
   | `polynomial_scaling(coeffs, T_ref)` | CCl4, acetone, ClONO2 | $\sigma_0 \cdot P(T - T_{ref}, \lambda)$ |
   | `exponential_scaling(coeffs, T_ref)` | CFC-11, RONO2, N2O5 | $\sigma_0 \cdot \exp(f(T, \lambda))$ |
   | `linear_correction(slope_data)` | CH2O | $\sigma_0 + \beta \cdot \Delta T$ |
   | `analytic(lambda_expr)` | Rayleigh, nitroxy compounds | User-supplied lambda for analytic formulas |
   | `stern_volmer(phi0, k)` | CH2O quantum yield | Quenching: $1/(1/\phi_0 + k \cdot M)$ |
   | `parameterized(type, params)` | Taylor/Burkholder/Harwood | Pluggable parameterization strategies |

- [ ] 9. **Implement transform combinators.** These compose primitives:

   | Combinator | Description |
   |------------|-------------|
   | `in_region(lambda_min, lambda_max, transform)` | Apply transform only in a wavelength sub-range; zero outside |
   | `compose(base, modifier)` | Apply modifier to base output (multiply, add) |
   | `piecewise({region1: t1, region2: t2, ...})` | Different transforms in different wavelength regions |
   | `clamp(transform, min, max)` | Clamp output values |
   | `override_band(band_name, value, transform)` | Replace values in named spectral bands (Lyman-alpha, Schumann-Runge) |
   | `scale(factor, transform)` | Multiply by constant factor |

- [ ] 10. **Port existing cross-section algorithms as composed transforms.** For each of the 27 cross-section types in `src/cross_sections/`, express as a composition of primitives. Example for CCl4:
    ```cpp
    auto ccl4_xs = compose(
        from_data(netcdf_reader("CCl4.nc"), conserving_interpolator()),
        in_region(194e-9, 250e-9,
            polynomial_scaling({b0, b1, b2, b3, b4}, 295.0))
    );
    ```
    For complex cases like O3 (4 regions, refraction correction) or H2O2 (blended polynomials via Boltzmann $\chi$), use `piecewise()` and `analytic()` with custom lambdas.

- [ ] 11. **Port quantum yield algorithms** (~20 types in `src/quantum_yields/`) using the same primitive/combinator system. Most reuse the same primitives; `stern_volmer`, `analytic`, and `temperature_interpolation` cover the majority.

- [ ] 12. **Port spectral weight algorithms** (~12 types in `src/spectral_weights/`) — these are simpler wavelength-only transforms (`gaussian`, `notch_filter`, `exponential_decay`, etc.).

## Phase 3: Rate Calculation Engine

- [ ] 13. **Implement photolysis rate calculator.** Port `src/photolysis_rates.F90`. The core operation is:
    $$J_i(z) = \sum_\lambda F(\lambda, z) \cdot \sigma_i(\lambda, z) \cdot \phi_i(\lambda, z) \cdot \Delta\lambda$$
    This is a batched dot product per height layer per reaction — highly parallelizable across columns and reactions. The API takes a solved `RadiationField`, a list of `(cross_section_transform, quantum_yield_transform)` pairs, and returns an `Array3D [reaction × height × column]`.

- [ ] 14. **Implement dose rate calculator.** Port `src/dose_rates.F90`. Similar structure: spectral irradiance × spectral weight → dose rate per height.

- [ ] 15. **Implement heating rate calculator.** Port the heating rate logic from `src/heating_rates.F90`.

## Phase 4: Radiator and Atmospheric State

- [ ] 16. **Implement radiator accumulation.** Port the radiator system from `src/radiative_transfer/radiators/`. In the new API, radiators are objects that provide optical depth, single-scattering albedo, and asymmetry parameter arrays. The `RadiatorState::Accumulate()` method (currently a placeholder) sums contributions from multiple radiators. Host models provide radiator data directly via the C++ API.

- [ ] 17. **Implement Lyman-Alpha and Schumann-Runge band parameterization.** Port `src/la_sr_bands.F90` — this handles special UV absorption bands of O₂ that require separate treatment from the main solver.

- [ ] 18. **Implement interpolation utilities.** Port the 4 interpolator types from `src/interpolate.F90` (linear, conserving, fractional_source, fractional_target) as standalone functions/objects used by the data reader system.

## Phase 5: Data Reader Abstraction

- [ ] 19. **Define pluggable data reader interface.** Abstract interface for loading cross-section/quantum yield parameter data:
    ```cpp
    class DataReader {
        virtual Array2D<double> read_parameters(const std::string& variable_prefix) = 0;
        virtual Array1D<double> read_wavelengths() = 0;
        virtual Array1D<double> read_temperatures() = 0;  // optional
    };
    ```
    Implement `NetCDFReader` using NetCDF-C (no Fortran dependency). This reads the existing ~160 `.nc` files unchanged. Additional readers (HDF5, CSV, binary) can be added later.

## Phase 6: Public API

- [ ] 20. **Define the top-level solver API.** This replaces the Fortran `core_t`. The API is purely programmatic — no config files:
    ```cpp
    template<typename ArrayPolicy>
    class Solver {
        // Configure solver with template policy for CPU or GPU
        void solve(
            const AtmosphericState<ArrayPolicy>& state,
            const RadiatorState<ArrayPolicy>& radiator_state,
            RadiationField<ArrayPolicy>& radiation_field   // output
        );
    };

    // Full pipeline: solve radiation field then compute rates
    template<typename ArrayPolicy>
    class PhotolysisCalculator {
        void add_reaction(std::string name, TransformFunc<ArrayPolicy> cross_section, TransformFunc<ArrayPolicy> quantum_yield);
        void calculate(const RadiationField<ArrayPolicy>& field, const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& rates);
    };
    ```

- [ ] 21. **C API wrapper.** Provide a `extern "C"` API for FFI consumption by MUSICA's Fortran interface, Python (via ctypes as fallback), and other languages. This mirrors how the current regression test in `test/regression/solvers/delta_eddington.hpp` defines C-interop structs.

## Phase 7: Discrete Ordinate Solver (Later)

- [ ] 22. **Port the DISORT solver.** The 3,700-line Fortran DISORT implementation in `src/radiative_transfer/solvers/discrete_ordinant_util.F90` is a well-known community code. Port the `PSNDO` subroutine and supporting functions to C++ with the same template policy system. This supports 2–32 even streams and is computationally heavier — GPU acceleration here will be most impactful.

## Phase 8: MUSICA Integration

- [ ] 23. **Update MUSICA to consume the new library.** MUSICA's CMake system (`TUVX_GIT_REPOSITORY` / `TUVX_GIT_TAG`) uses FetchContent for tuv-x. Point it at the new repo. Build the YAML/JSON configuration layer in MUSICA that maps config files to the pure C++ API, translating `"type": "tint"` into the corresponding composed lambda transforms. MUSICA's PyBind11 / JS / Julia bindings wrap the new C API.

## Verification

- **Unit tests** (GTest): Each transform primitive, combinator, interpolator, data reader, and solver step gets individual tests. Migrate existing C++ tests from the current repo.
- **Regression tests**: Pre-compute reference outputs from the current Fortran solver for several standard configurations (`examples/tuv_5_4.json`, `examples/ts1_tsmlt.json`) and store as binary fixtures. The new C++ solver must match within floating-point tolerance.
- **Cross-validation of transforms**: For each of the 27 cross-section types and 20 quantum yield types, compare the new composed-transform output against Fortran reference values at multiple temperature/pressure points.
- **Performance benchmarks**: Google Benchmark tests for multi-column solver throughput (1, 10, 100, 1000 columns) on CPU. CUDA/HIP benchmarks once GPU policy is implemented.
- **MUSICA integration test**: End-to-end test through MUSICA's config-driven pipeline producing identical photolysis rates.

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Language | C++ | Builds on existing partial port; native MUSICA integration; standard HPC ecosystem |
| Transform API | Lambda/functional composables | Removes ~3,000 lines of factory/config dispatch code; composable and testable; config layer lives in MUSICA |
| GPU portability | Template policies | Lighter weight than Kokkos/SYCL; consistent with existing codebase; can wrap Kokkos later if needed |
| Migration strategy | New standalone repo | Clean separation; old Fortran repo stays as reference/validation; MUSICA switches dependency when ready |
| Data readers | Pluggable interface, NetCDF-C default | Future-proofs against format changes; no Fortran dependency |
| Solver order | Delta Eddington first | Simpler two-stream method; validates full pipeline before tackling DISORT's 3,700 lines |
| Configuration | No config files in solver library | Clean separation of concerns; MUSICA handles all configuration mapping |
| Units | SI only, no conversions | All inputs/outputs in SI units (m, rad, s, K, Pa); no unit conversion code in the library; callers and data pipelines are responsible for providing SI data |
| MICM consistency | Align API conventions with MICM | Both TUV-x and MICM are MUSICA components; consistent syntax (`operator[]`, `ForEachRow`, `ColumnView`, `Function` factory, `ArrayPolicy` template pattern) reduces cognitive load for users working across both libraries. Diverge only when TUV-x requirements demand it (e.g., `Array1D` instead of `std::vector`) |
