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

## Analytic spectral weights (PR 2)

`gaussian.csv`, `plant_damage.csv`, `plant_damage_flint_caldwell.csv`,
`plant_damage_flint_caldwell_ext.csv`, and `phytoplankton_boucher.csv` complete
the analytic spectral weights (all but `eppley`, which needs tabulated data).
Same commit, dimensionless, 2-column format. Shared grid (13 mid-points,
285-405 nm) brackets the cutoffs 313 / 366 / 390 nm and the phytoplankton region
edges 290 / 400 nm without landing on them.

| File | Fortran source | Formula (`w` = wl nm) |
|------|----------------|------------------------|
| `gaussian.csv` | `gaussian.F90` | `exp(-ln2*0.04*(w-305)^2)`, normalized to unit sum over the grid |
| `plant_damage.csv` | `plant_damage.F90` | `570.25 - 4.70144 w + 0.01274 w^2 - 1.13118e-5 w^3`, zeroed if <0 or w>313 |
| `plant_damage_flint_caldwell.csv` | `plant_damage_flint_caldwell.F90` | Flint-Caldwell (below), cutoff 366 |
| `plant_damage_flint_caldwell_ext.csv` | `plant_damage_flint_caldwell_ext.F90` | same, cutoff 390 |
| `phytoplankton_boucher.csv` | `phytoplankton_boucher.F90` | `max(0, -3.17e-6 + exp(112.5 - 0.6223 w + 7.67e-4 w^2))` for 290<w<400 |

Flint-Caldwell: `exp( 4.688272*exp(-exp(0.1703411*(w-307.867)/1.15)) + ((390-w)/121.7557 - 4.183832) ) * w/300`,
zeroed if <0 or above the cutoff.

The C++ side builds each from `wrap_analytic`; `gaussian` composes it with the
new `normalize` combinator (divide by the wavelength-sum). The gaussian
reference is grid-dependent (the normalization sums over these exact
wavelengths), so it is reproduced on the same grid.

Only `eppley` remains from Step 12 -- it reads tabulated data and normalizes by
the band-width integral, so it waits for the NetCDF reader (Phases 4-5).
