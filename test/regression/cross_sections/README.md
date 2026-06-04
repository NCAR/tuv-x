# Cross-Section Regression Reference Data

## Analytic cross-sections (PR A)

The five CSV files in `reference/` contain cross-section values for the
temperature- and height-independent ("analytic") cross-sections. They were
produced by a standalone program that evaluates the formulas below — the same
formulas implemented by the Fortran source files, with no TUV-x or NetCDF
dependency.

These values are the ground truth the C++ implementation is regression-tested
against (relative tolerance 1e-10). The C++ side assembles each species from
the general analytic transform forms (`tuvx/transforms/analytic_forms.hpp`) in
the example header `tuvx/fixed_configuration.hpp`.

### Reference source

All formulas were transcribed from the Fortran on the `main` branch at commit:

```
c2a8f3385c45df08dd65a0d2503665182332b783
```

| File | Fortran source (at that commit) |
|------|---------------------------------|
| `rayleigh.csv` | `src/cross_sections/rayliegh.F90` |
| `hobr.csv` | `src/cross_sections/hobr-oh_br.F90` |
| `t_butyl_nitrate.csv` | `src/cross_sections/t_butyl_nitrate.F90` |
| `nitroxy_acetone.csv` | `src/cross_sections/nitroxy_acetone.F90` |
| `nitroxy_ethanol.csv` | `src/cross_sections/nitroxy_ethanol.F90` |

To inspect the original Fortran without checking the branch out:
`git show c2a8f33:src/cross_sections/rayliegh.F90` (etc.).

## How to recreate the reference data

The data are reproducible in any language (the original generator was a small
standalone Fortran program; the C++ regression test would produce identical
values). To regenerate:

1. Build the **wavelength grid**: 26 mid-points, 200 nm to 450 nm inclusive in
   10 nm steps (200, 210, …, 450). This range deliberately straddles every
   species' active band so the zero-weight regions are exercised too.
2. For each mid-point, evaluate the species formula below. The Fortran works in
   **nm** and returns **cm²**.
3. Convert to SI for the CSV: `wavelength_m = wavelength_nm × 1e-9` and
   `cross_section_m2 = cross_section_cm2 × 1e-4`.
4. Write one CSV per species with header `wavelength_m,cross_section_m2`, one
   row per wavelength, using ≥17 significant digits so the values round-trip
   exactly as IEEE-754 doubles.

Let `wl` be the wavelength mid-point in **nm** in all formulas below.

### rayleigh

No wavelength-range restriction (applies across the whole grid).

```
mu  = wl / 1000                                     # micrometers
pwr = 3.6772 + 0.389*mu + 0.09426/mu    if mu <= 0.55
    = 4.04                              if mu >  0.55
sigma_cm2 = 4.02e-28 / mu**pwr
```

### hobr

Active for `250 <= wl <= 550`; `sigma_cm2 = 0` outside.

```
sigma_cm2 = ( 24.77 * exp(-109.80 * (ln(284.01/wl))**2)
            + 12.22 * exp( -93.63 * (ln(350.57/wl))**2)
            + 2.283 * exp(-242.40 * (ln(457.38/wl))**2) ) * 1e-20
```

(Note: `hobr-oh_br.F90` also declares unused `a/b/c` constants identical to
`nitroxy_ethanol`'s; they are dead code and play no part in the HOBr value.)

### t_butyl_nitrate, nitroxy_acetone, nitroxy_ethanol

All three share the form, differing only in coefficients and active range.
Active for `wl_min <= wl <= wl_max`; `sigma_cm2 = 0` outside.

```
sigma_cm2 = exp( c + wl * (b + a*wl) )
```

| Species | a | b | c | wl_min (nm) | wl_max (nm) |
|---------|---|---|---|-------------|-------------|
| t_butyl_nitrate | -0.993e-3 | 0.5307 | -115.5 | 270 | 330 |
| nitroxy_acetone | -1.365e-3 | 0.7834 | -156.8 | 284 | 335 |
| nitroxy_ethanol | -2.359e-3 | 1.2478 | -210.4 | 270 | 306 |
