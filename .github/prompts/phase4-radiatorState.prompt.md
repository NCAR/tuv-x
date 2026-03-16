# Phase 4: Radiator and Atmospheric State (Steps 16–18)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

Before the solver can compute the radiation field, the optical properties of the atmosphere must be assembled. Radiators contribute optical depth, single-scattering albedo, and asymmetry parameter arrays. The Lyman-Alpha and Schumann-Runge bands require special parameterization. Interpolation utilities support data loading.

## Step 16: Implement radiator accumulation

Port the radiator system from `src/radiative_transfer/radiators/`.

In the new API, radiators are objects providing optical depth, single-scattering albedo, and asymmetry parameter arrays per wavelength per layer per column. The `RadiatorState::Accumulate()` method (currently a placeholder in the C++ code) sums contributions from multiple radiators using the correct weighting:

- Optical depth: additive ($\tau_{total} = \sum \tau_i$)
- Single-scattering albedo: weighted by optical depth ($\omega_{total} = \sum \omega_i \tau_i / \tau_{total}$)
- Asymmetry parameter: weighted by SSA × optical depth ($g_{total} = \sum g_i \omega_i \tau_i / (\omega_{total} \tau_{total})$)

Host models provide radiator data directly via the C++ API — no config-file-based radiator construction.

### Fortran reference
- `src/radiative_transfer/radiators/aerosol.F90`
- `src/radiative_transfer/radiators/from_host.F90`
- `src/radiative_transfer/radiators/from_netcdf_file.F90`
- `src/radiative_transfer/radiative_transfer.F90` — orchestrator that iterates over radiators

## Step 17: Implement Lyman-Alpha and Schumann-Runge band parameterization

Port `src/la_sr_bands.F90`. This handles special UV absorption bands of O₂ (121.4–175.4 nm) that require separate treatment:

- Lyman-Alpha band: strong O₂ absorption line at 121.6 nm
- Schumann-Runge bands: structured O₂ absorption between 175–205 nm
- These bands use lookup tables and parameterized transmission functions rather than the standard column-by-column solver

### Fortran reference
- `src/la_sr_bands.F90` — `la_sr_bands_t` class with LUT-based calculations

## Step 18: Implement interpolation utilities

Port the 4 interpolator types from `src/interpolate.F90`:

| Interpolator | Use case | Method |
|-------------|----------|--------|
| `linear` | Point-to-point (temperature profiles) | Linear interpolation |
| `conserving` | Cross-sections, quantum yields (default) | Area-preserving: trapezoidal area over source bins ÷ target bin width |
| `fractional_source` | Bin-to-bin relative to source width | For extraterrestrial flux |
| `fractional_target` | Bin-to-bin relative to target width | For bin-mapped data |

The conserving interpolator is the most important — it ensures spectral quantities are correctly mapped between wavelength grids while preserving integrated values. Supports `fold_in` for handling data that extends beyond the model grid.

### Fortran reference
- `src/interpolate.F90` — `interpolator_linear_t`, `interpolator_conserving_t`, `interpolator_fractional_source_t`, `interpolator_fractional_target_t`
