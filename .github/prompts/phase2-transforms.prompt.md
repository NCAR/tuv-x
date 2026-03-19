# Phase 2: Composable Transform System (Steps 7–12)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The current Fortran codebase has 27 cross-section types, 20 quantum yield types, and 12 spectral weight types, each implemented as a separate class with a factory pattern. Analysis shows these decompose into ~8 recurring mathematical patterns. The new library defines an open `TransformFunc` API — any lambda or callable matching the signature is a valid transform. A convenience library of pre-built factory functions covers the common patterns, and combinators allow composing transforms. Users can mix library-provided factories, combinators, and raw user-written lambdas freely. No config files — MUSICA handles configuration.

## Step 7: Design the transform type system

Define `TransformFunc` as a callable mapping atmospheric state → 3D output array `[wavelength × height × column]`:

```cpp
template<typename ArrayPolicy>
using TransformFunc = std::function<void(
    const Grid<ArrayPolicy>& wavelength_grid,
    const Grid<ArrayPolicy>& altitude_grid,
    const Profile<ArrayPolicy>& temperature,
    const Profile<ArrayPolicy>& pressure,
    const Profile<ArrayPolicy>& air_density,
    Array3D<typename ArrayPolicy::value_type>& output  // [wavelength × height × column]
)>;
```

**`TransformFunc` is the open, user-facing API.** Any callable matching this signature is a valid transform. Users can write raw lambdas, function objects, or free functions directly — no factory or wrapper is required:

```cpp
// User-written transform — just a lambda matching the TransformFunc signature:
TransformFunc<ArrayPolicy> rayleigh_xs = [](const auto& wl_grid, const auto& alt_grid,
    const auto& temperature, const auto& pressure, const auto& air_density, auto& output) {
    // Rayleigh scattering formula using ForEachRow, ColumnView, etc.
    // This IS the transform — no wrapper needed.
};
```

The output includes the column dimension because temperature, pressure, and density vary per column, so the computed values differ per column.

Transform implementations should use the **bulk operations API** (`ForEachRow`, `ColumnView`, `Function`) for all per-element computation. For transforms computed from static reference data (loaded via `DataReader`), the reference data is broadcast into the per-column output using `ForEachRow`. Performance-critical transforms (called every timestep) should return pre-compiled `ArrayPolicy::Function` objects.

## Step 8: Implement convenience transform library

Provide a library of pre-built factory functions that return `TransformFunc` callables. These are **convenience utilities built on the same public API a user would use** — they capture parameters and return a lambda matching the `TransformFunc` signature. Users can use these off-the-shelf, compose them with combinators (Step 9), or ignore them entirely and write raw lambdas.

### Factory functions (each returns a `TransformFunc`)

| Factory function | Math | Useful for |
|-----------|------|---------------------------|
| `from_data(reader, interpolator)` | Tabulated data → model grid interpolation | Base cross-sections, quantum yields |
| `temperature_interpolation(reader)` | $\sigma(\lambda,T) = \sigma_i + \frac{T-T_i}{T_{i+1}-T_i}(\sigma_{i+1}-\sigma_i)$ | T-dependent tint types |
| `polynomial_scaling(coeffs, T_ref)` | $\sigma_0 \cdot P(T-T_{ref}, \lambda)$ | CCl4, acetone, ClONO2, HCFC |
| `exponential_scaling(coeffs, T_ref)` | $\sigma_0 \cdot \exp(f(T,\lambda))$ | CFC-11, RONO2, N2O5, CHBr3, CH3ONO2 |
| `linear_correction(slope_data)` | $\sigma_0 + \beta \cdot \Delta T$ | CH2O |
| `stern_volmer(phi0, k)` | $\phi = 1/(1/\phi_0 + k \cdot M)$ | Quenching quantum yields |
| `parameterized(type, params)` | Taylor series, Burkholder, Harwood strategies | T-parameterized cross-sections |

### Signature adapter

| Adapter | Purpose |
|---------|---------|
| `wrap_analytic(f)` | Takes a simple `f(λ)→double` or `f(λ,T)→double` and wraps it into the full `TransformFunc` signature, handling grid iteration and column broadcasting internally. Useful for quick one-liner formulas. |

`wrap_analytic` is a convenience for simple formulas — it is **not** the primary way to create custom transforms. For anything beyond a trivial formula, users should write a full `TransformFunc` lambda directly.

## Step 9: Implement transform combinators

Combinators are higher-order functions: they take one or more `TransformFunc` values and return a new `TransformFunc`. They work with **any** transform — library-provided factories or user-written lambdas:

| Combinator | Semantics |
|------------|-----------|
| `in_region(λ_min, λ_max, t)` | Apply `t` only where λ_min ≤ λ ≤ λ_max; zero outside |
| `multiply(a, b)` | Run both transforms, multiply their outputs element-wise |
| `add(a, b)` | Run both transforms, add their outputs element-wise |
| `piecewise({regions...})` | Different transforms per wavelength region |
| `clamp(t, min, max)` | Clamp output values |
| `override_band(name, value, t)` | Replace values in named bands (Lyman-α, Schumann-Runge) |
| `scale(factor, t)` | Multiply by constant |

**Why combinators instead of sequential application?** A `TransformFunc` *writes* to the output array — it doesn't modify an existing value. Applying two transforms in sequence would just overwrite the first result with the second. Combinators like `multiply(a, b)` internally run both transforms (one into the output, one into a temporary), then combine the results element-wise.

## Step 10: Port cross-section algorithms

Express each of the 27 Fortran cross-section types as either a library factory, a composition of factories via combinators, or a user-written `TransformFunc` lambda. Reference files in `src/cross_sections/`:

### Simple cases (direct factory mapping)
- `base` → `from_data(reader, conserving_interpolator())`
- `tint` → `temperature_interpolation(reader)`
- `CCl4` → `multiply(from_data(...), in_region(194e-9, 250e-9, polynomial_scaling(...)))`
- `CFC-11` → `multiply(from_data(...), exponential_scaling(...))`
- `CH2O` → `add(from_data(...), linear_correction(slope_data))`
- `Rayleigh` → `wrap_analytic([λ]{ ... })` or a direct `TransformFunc` lambda — Rayleigh scattering formula with λ in meters

### Complex cases (user-written lambdas or mixed compositions)
- `O3` → `piecewise({<185e-9: data1, 185e-9–195e-9: data2, 195e-9–345e-9: tint, >345e-9: data4})` with refraction correction — the refraction correction can be a user-written `TransformFunc` lambda
- `H2O2` → blended polynomial via Boltzmann χ — write as a direct `TransformFunc` lambda
- `acetone` → Horner polynomial in temperature — write as a direct `TransformFunc` lambda or use `polynomial_scaling`

### Fortran source files to reference
All in `src/cross_sections/`: `tint.F90`, `temperature_based.F90`, `no2_tint.F90`, `o3_tint.F90`, `rayliegh.F90`, `ccl4.F90`, `cfc-11.F90`, `ch2o.F90`, `h2o2-oh_oh.F90`, `hno3-oh_no2.F90`, `n2o5-no2_no3.F90`, `clono2.F90`, `hobr-oh_br.F90`, `hcfc.F90`, `chbr3.F90`, `chcl3.F90`, `cl2-cl_cl.F90`, `rono2.F90`, `t_butyl_nitrate.F90`, `acetone-ch3co_ch3.F90`, `ch3ono2-ch3o_no2.F90`, `nitroxy_acetone.F90`, `nitroxy_ethanol.F90`, `n2o-n2_o1d.F90`, `oclo.F90`

Temperature parameterization utilities in `src/cross_sections/util/`: `temperature_parameterization.F90`, `temperature_parameterization_taylor_series.F90`, `temperature_parameterization_burkholder.F90`, `temperature_parameterization_harwood.F90`

## Step 11: Port quantum yield algorithms

~20 types in `src/quantum_yields/`. Most map to library factories, combinators, or direct `TransformFunc` lambdas:

| Fortran type | Transform |
|--------------|-------------------|
| `base` | `from_data(...)` or constant value |
| `tint`, `NO2 tint` | `temperature_interpolation(...)` |
| `CH2O` | `stern_volmer(phi0, k)` in region 330e-9–360e-9 m |
| `taylor series` | direct `TransformFunc` lambda with clamp(0,1) |
| `O3→O(1D)` | direct `TransformFunc` lambda — multi-region Gaussians (Matsumi formula) |
| `H2SO4 Mills` | direct `TransformFunc` lambda — mean free path formula |
| `acetone` | direct `TransformFunc` lambda — Blitz parameterization |

## Step 12: Port spectral weight algorithms

~12 types in `src/spectral_weights/`. Simpler — wavelength-only (no T/P dependence):

| Fortran type | Transform |
|-------------|-----------|
| `base` | `from_data(reader, conserving_interpolator())` |
| `gaussian` | `wrap_analytic([](λ, μ) { return exp(-ln2 * 0.04 * (λ-μ)²); })` |
| `notch_filter` | `in_region(λ_min, λ_max, constant(1.0))` |
| `exp_decay` | `wrap_analytic([](λ) { return pow(10, (300e-9-λ)/14e-9); })` |
| `PAR` | `in_region(400e-9, 700e-9, constant(1.0))` |
| Others | Literature-specific analytic formulas |
