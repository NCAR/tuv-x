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

## Temperature-dependent cross-sections (PR B)

`n2o.csv`, `cl2.csv`, `h2o2.csv`, and `chbr3.csv` cover the temperature-dependent
analytic cross-sections. Same `main` commit (`c2a8f33`) and same SI conventions
(`wavelength_m = wl_nm * 1e-9`, `cross_section_m2 = cross_section_cm2 * 1e-4`),
but each file now carries a **temperature axis**:

```
header:  wavelength_m,temperature_K,cross_section_m2
```

Rows are the Cartesian product of the per-species wavelength grid and
temperature samples below. `wl` is the wavelength mid-point in **nm**, `T` in K.

| File | Fortran source |
|------|----------------|
| `n2o.csv` | `src/cross_sections/n2o-n2_o1d.F90` |
| `cl2.csv` | `src/cross_sections/cl2-cl_cl.F90` |
| `h2o2.csv` | `src/cross_sections/h2o2-oh_oh.F90` |
| `chbr3.csv` | `src/cross_sections/chbr3.F90` |

### n2o

Active for `173 <= wl <= 240`; `0` outside. `Tadj = clamp(T, 194, 320)`.

```
A = ((((A4*wl + A3)*wl + A2)*wl + A1)*wl + A0)        # Horner, deg 4
B = (((B3*wl + B2)*wl + B1)*wl + B0)                  # Horner, deg 3
sigma_cm2 = exp( A + (Tadj - 300) * exp(B) )
```

A0..A4 = 68.21023, -4.071805, 4.301146e-2, -1.777846e-4, 2.520672e-7
B0..B3 = 123.4014, -2.116255, 1.111572e-2, -1.881058e-5

- wavelengths (nm): 160,170,180,190,200,210,220,230,238,250  (238 rather than the
  240 band edge: a grid point exactly on the edge can flip inclusion by one ULP
  under the m<->nm conversion, so the regression avoids sitting on it)
- temperatures (K): 180, 200, 250, 300, 330  (180/330 exercise the clamp; 300 zeroes the B term)

### cl2

No wavelength restriction.

```
beta  = 402.7 / T
alpha = (exp(2*beta) - 1) / (exp(2*beta) + 1)         # = tanh(beta)
sigma_cm2 = 1e-20 * sqrt(alpha) *
            ( 27.3  * exp(-99.0 * alpha * (ln(329.5/wl))**2)
            + 0.932 * exp(-91.5 * alpha * (ln(406.5/wl))**2) )
```

- wavelengths (nm): 300, 329.5, 360, 406.5, 450  (329.5 / 406.5 are the band centres)
- temperatures (K): 200, 250, 300, 350

### h2o2 (hybrid)

In-band `260 <= wl < 350`: analytic; below the band: tabulated. `T` clamp `[200,400]`.

```
sumA = deg-7 Horner poly in wl    # A0..A7 = 6.4761e4, -9.2170972e2, 4.535649,
                                  #   -4.4589016e-3, -4.035101e-5, 1.6878206e-7,
                                  #   -2.652014e-10, 1.5534675e-13
sumB = deg-4 Horner poly in wl    # B0..B4 = 6.8123e3, -5.1351e1, 1.1522e-1,
                                  #   -3.0493e-5, -1.0924e-7
chi  = 1 / (1 + exp(-1265 / T))
sigma_cm2 = (chi*sumA + (1-chi)*sumB) * 1e-21
```

Note: the Fortran **source** uses `A0 = 6.4761E+04`. The bundled `.nc` header text
says `6.4761E-04`; that is a documentation typo — the source value is authoritative.

Tabulated (below-band) values are from `cross_section.h2o2-oh_oh.nc` (cm^2):

| wl (nm) | 190 | 195 | 200 | 205 | 210 |
|---------|-----|-----|-----|-----|-----|
| sigma   | 6.72e-19 | 5.64e-19 | 4.75e-19 | 4.08e-19 | 3.57e-19 |

- analytic wavelengths (nm): 270, 290, 310, 330, 345
- temperatures (K): 180, 250, 298, 350, 420  (180/420 exercise the clamp)

### chbr3 (hybrid)

Analytic only where `290 < wl < 340` **and** `210 < T < 300`; otherwise tabulated.

```
sigma_cm2 = exp( (C0 - C1*wl)*(T0 - T) - (C2 + C3*wl) )
```

C0,C1,C2,C3,T0 = 0.06183, 0.000241, 2.376, 0.14757, 273

Tabulated (below-band) values are from `cross_section.chbr3.nc` (cm^2):

| wl (nm) | 190 | 192 | 194 | 196 | 198 | 200 | 202 | 204 | 206 | 208 | 210 |
|---------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| sigma (1e-18) | 3.99 | 3.60 | 3.51 | 3.66 | 3.93 | 4.16 | 4.33 | 4.40 | 4.45 | 4.51 | 4.68 |

- analytic wavelengths (nm): 300, 310, 320, 330
- temperatures (K): 220, 250, 295  (all inside (210,300))

### Limitation of the hybrid reference grids

The hybrid species (`h2o2`, `chbr3`) fall back to data tabulated in NetCDF
outside their analytic band. In the production Fortran that data is resampled
onto the model grid by the *conserving interpolator*, which is not yet ported to
C++. To keep this reference faithful **without** that interpolator, the hybrid
grids are pinned so every tabulated bin sits exactly on a NetCDF native
wavelength (resampling is then identity) and no bin lands in a region that would
require interpolation or extrapolation:

- no bins in the 210–260 nm gap (h2o2) or above the tabulated range,
- CHBr3 temperatures kept inside (210,300) so in-band bins always take the
  analytic branch.

The arbitrary-grid resampling path and CHBr3's temperature-gated fallback are
therefore **not** exercised here; they will be covered once the conserving
interpolator and NetCDF reader land (Phases 4–5).
