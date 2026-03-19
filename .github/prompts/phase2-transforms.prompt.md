# Phase 2: Composable Transform System (Steps 7–12)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The current Fortran codebase has 27 cross-section types, 20 quantum yield types, and 12 spectral weight types, each implemented as a separate class with a factory pattern. Analysis shows these decompose into ~8 recurring mathematical patterns. The new library replaces this with composable lambda-based transforms configured via a programmatic C++ API (no config files — MUSICA handles configuration).

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

This is the universal signature for cross-sections, quantum yields, and spectral weights. The output includes the column dimension because temperature, pressure, and density vary per column, so the computed values differ per column.

Transform implementations should use the **bulk operations API** (`ForEachRow`, `ColumnView`, `Function`) for all per-element computation. For transforms computed from static reference data (loaded via `DataReader`), the reference data is broadcast into the per-column output using `ForEachRow`. Performance-critical transforms (called every timestep) should return pre-compiled `ArrayPolicy::Function` objects.

## Step 8: Implement composable transform primitives

Each primitive is a factory function returning a `TransformFunc`:

| Primitive | Math | Replaces (Fortran types) |
|-----------|------|--------------------------|
| `from_data(reader, interpolator)` | Tabulated data → model grid interpolation | Base `cross_section_t`, base `quantum_yield_t` |
| `temperature_interpolation(reader)` | $\sigma(\lambda,T) = \sigma_i + \frac{T-T_i}{T_{i+1}-T_i}(\sigma_{i+1}-\sigma_i)$ | `tint` xs, `tint` qy, `NO2 tint`, `O3 tint` |
| `polynomial_scaling(coeffs, T_ref)` | $\sigma_0 \cdot P(T-T_{ref}, \lambda)$ | CCl4, acetone, ClONO2, HCFC |
| `exponential_scaling(coeffs, T_ref)` | $\sigma_0 \cdot \exp(f(T,\lambda))$ | CFC-11, RONO2, N2O5, CHBr3, CH3ONO2 |
| `linear_correction(slope_data)` | $\sigma_0 + \beta \cdot \Delta T$ | CH2O |
| `analytic(lambda_expr)` | User-supplied lambda | Rayleigh, nitroxy acetone, nitroxy ethanol |
| `stern_volmer(phi0, k)` | $\phi = 1/(1/\phi_0 + k \cdot M)$ | CH2O quantum yield |
| `parameterized(type, params)` | Taylor series, Burkholder, Harwood strategies | `temperature_based` xs |

## Step 9: Implement transform combinators

| Combinator | Semantics |
|------------|-----------|
| `in_region(λ_min, λ_max, t)` | Apply `t` only where λ_min ≤ λ ≤ λ_max; zero outside |
| `compose(base, modifier)` | `output = base * modifier` or `base + modifier` |
| `piecewise({regions...})` | Different transforms per wavelength region |
| `clamp(t, min, max)` | Clamp output values |
| `override_band(name, value, t)` | Replace values in named bands (Lyman-α, Schumann-Runge) |
| `scale(factor, t)` | Multiply by constant |

## Step 10: Port cross-section algorithms

Express each of the 27 Fortran cross-section types as compositions. Reference files in `src/cross_sections/`:

### Simple cases (direct primitive mapping)
- `base` → `from_data(reader, conserving_interpolator())`
- `tint` → `temperature_interpolation(reader)`
- `CCl4` → `compose(from_data(...), in_region(194e-9, 250e-9, polynomial_scaling(...)))`
- `CFC-11` → `compose(from_data(...), exponential_scaling(...))`
- `CH2O` → `compose(from_data(...), linear_correction(slope_data))`
- `Rayleigh` → `analytic([](λ) { ... })` — Rayleigh scattering formula with λ in meters

### Complex cases (multi-region / custom logic)
- `O3` → `piecewise({<185e-9: data1, 185e-9–195e-9: data2, 195e-9–345e-9: tint, >345e-9: data4})` with refraction correction
- `H2O2` → blended polynomial via Boltzmann χ — use `analytic()` with custom lambda
- `acetone` → Horner polynomial in temperature

### Fortran source files to reference
All in `src/cross_sections/`: `tint.F90`, `temperature_based.F90`, `no2_tint.F90`, `o3_tint.F90`, `rayliegh.F90`, `ccl4.F90`, `cfc-11.F90`, `ch2o.F90`, `h2o2-oh_oh.F90`, `hno3-oh_no2.F90`, `n2o5-no2_no3.F90`, `clono2.F90`, `hobr-oh_br.F90`, `hcfc.F90`, `chbr3.F90`, `chcl3.F90`, `cl2-cl_cl.F90`, `rono2.F90`, `t_butyl_nitrate.F90`, `acetone-ch3co_ch3.F90`, `ch3ono2-ch3o_no2.F90`, `nitroxy_acetone.F90`, `nitroxy_ethanol.F90`, `n2o-n2_o1d.F90`, `oclo.F90`

Temperature parameterization utilities in `src/cross_sections/util/`: `temperature_parameterization.F90`, `temperature_parameterization_taylor_series.F90`, `temperature_parameterization_burkholder.F90`, `temperature_parameterization_harwood.F90`

## Step 11: Port quantum yield algorithms

~20 types in `src/quantum_yields/`. Most map directly to existing primitives:

| Fortran type | Composed transform |
|--------------|-------------------|
| `base` | `from_data(...)` or constant value |
| `tint`, `NO2 tint` | `temperature_interpolation(...)` |
| `CH2O` | `stern_volmer(phi0, k)` in region 330e-9–360e-9 m |
| `taylor series` | `analytic(polynomial_lambda)` with clamp(0,1) |
| `O3→O(1D)` | `analytic(matsumi_formula)` — multi-region Gaussians |
| `H2SO4 Mills` | `analytic(mean_free_path_formula)` |
| `acetone` | `analytic(blitz_parameterization)` |

## Step 12: Port spectral weight algorithms

~12 types in `src/spectral_weights/`. Simpler — wavelength-only (no T/P dependence):

| Fortran type | Transform |
|-------------|-----------|
| `base` | `from_data(reader, conserving_interpolator())` |
| `gaussian` | `analytic([](λ, μ) { return exp(-ln2 * 0.04 * (λ-μ)²); })` |
| `notch_filter` | `in_region(λ_min, λ_max, analytic([](λ) { return 1.0; }))` |
| `exp_decay` | `analytic([](λ) { return pow(10, (300e-9-λ)/14e-9); })` |
| `PAR` | `in_region(400e-9, 700e-9, constant(1.0))` |
| Others | Literature-specific analytic formulas |
