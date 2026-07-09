# Spectral-Weight Regression Reference Data

## Analytic spectral weights (PR 1)

The CSV files in `reference/` contain **unitless, wavelength-only** dose-rate
action-spectrum weights. They were produced by a standalone program evaluating
the same formulas as the Fortran `src/spectral_weights/` sources (no TUV-x or
NetCDF dependency), transcribed from the `main` branch at commit `c2a8f33`.

Each file has columns: `wavelength_m, weight`. Values are dimensionless (no
cm^2 -> m^2 conversion). Wavelength grid: 13 mid-points bracketing each
breakpoint (298 / 328 / 400 nm for erythema; 400 / 700 nm for PAR) *without*
landing exactly on one, since the m<->nm reconstruction could otherwise flip a
strict comparison by one ULP.

| File | Fortran source | Formula (`w` = wavelength nm) |
|------|----------------|-------------------------------|
| `standard_human_erythema.csv` | `standard_human_erythema.F90` | CIE erythema `sw_fery(w)` (piecewise, see below) |
| `uv_index.csv` | `UV_Index.F90` | `40 * sw_fery(w)` |
| `scup_mice.csv` | `scup_mice.F90` | `sw_futr(w) / sw_futr(300)` |
| `exp_decay.csv` | `exp_decay.F90` | `10^((300 - w)/14)` |
| `par.csv` | `par.F90` | `8.36e-3 * w` for `400 < w < 700`, else 0 |

`sw_fery(w)`:
`1` for `w <= 298`; `10^(0.094*(298-w))` for `298 < w <= 328`;
`10^(0.015*(140-w))` for `328 < w <= 400`; `1e-36` otherwise.

`sw_futr(w) = exp(P(w))`, where `P` is the 4th-order Lagrange interpolating
polynomial (deGruijl et al. 1993) through nodes
`w = 270, 302, 334, 367, 400 nm` with values
`-10.91, -0.86, -8.60, -9.36, -13.15`.

The C++ side builds each from `wrap_analytic` (no new library code). `sw_fery`
and `sw_futr` are shared helpers in `fixed_configuration.hpp`.

### Deferred

- `gaussian` and `eppley` normalize by a sum / integral over the wavelength grid
  (a reduction the element-wise forms don't provide); they wait for a `normalize`
  combinator, and `eppley` additionally needs tabulated data.
- `plant_damage`, `plant_damage_flint_caldwell(_ext)`, and `phytoplankton_boucher`
  are analytic and will follow in a second batch.
- `notch_filter` is simply `in_region(begin, end, constant(1))` — a composition,
  not a dedicated spectrum.
