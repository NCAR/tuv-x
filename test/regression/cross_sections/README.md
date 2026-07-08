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

The degree-7 polynomial has large alternating coefficients whose terms nearly
cancel, so the in-band result is ill-conditioned: different platforms' libm and
floating-point contraction evaluate it slightly differently (worst observed
~1.4e-10 relative). The H2O2 regression therefore uses a relative tolerance of
1e-8 (still ~8 significant figures) rather than the 1e-10 used for the others.

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

## Tabulated base x temperature correction (PR C)

`hno3.csv`, `rono2.csv`, `ch3ono2.csv`, `ch2o.csv`, and `cfc11.csv` cover
cross-sections that are a **tabulated base** cross-section combined with a
temperature correction. Same `main` commit (`c2a8f33`), SI conventions, and
temperature-axis CSV format (`wavelength_m,temperature_K,cross_section_m2`) as
the PR B files.

Because the base data comes from NetCDF and the conserving interpolator is not
yet ported, each grid is **pinned to the file's native wavelengths** (identity
resampling); the hardcoded base arrays in `fixed_configuration.hpp` are the same
values, so the comparison is exact. Temperatures sampled: 250, 280, 298, 330 K
(298 K = reference, zeroing every correction).

| File | Fortran source | Base data (`.nc`) | Formula (`l` = wl nm, `dT = T - 298`) |
|------|----------------|-------------------|----------------------------------------|
| `hno3.csv` | `hno3-oh_no2.F90` | `cross_section.hno3.nc` | `sigma0(l) * exp(sigma1(l)*dT)` |
| `rono2.csv` | `rono2.F90` | `cross_section.rono2.nc` | `sigma0(l) * exp(sigma1(l)*dT)` |
| `ch3ono2.csv` | `ch3ono2-ch3o_no2.F90` | `cross_section.ch3ono2_ch3o_no2.nc` | `sigma0(l) * exp(sigma1(l)*dT)` |
| `ch2o.csv` | `ch2o.F90` | `cross_section.ch2o.nc` | `sigma0(l) + sigma1(l)*dT` |
| `cfc11.csv` | `cfc-11.F90` | `cross_section.cfc-11.nc` | `sigma0(l) * exp((l - 184.9)*1e-4*dT)` |

`sigma0` (and, for CH2O, `sigma1`) are cross-sections in cm^2, converted to m^2
by 1e-4. For HNO3/RONO2/CH3ONO2 the exponent coefficient `sigma1` (1/K) is
tabulated; for CFC-11 only the base is tabulated and the temperature dependence
is a closed-form function of wavelength. The C++ side composes these from the
existing factory library: `multiply(from_data(sigma0), exponential_scaling(...))`,
`linear_correction(...)`, and `multiply(from_data(sigma0), bounded_analytic(...))`.

## Temperature-table interpolation (PR D)

`no2.csv` covers a cross-section tabulated at several temperatures and linearly
interpolated to the model temperature (the Fortran "tint" family). Same commit,
SI conventions, and temperature-axis format as PR B/C.

| File | Fortran source | Base data (`.nc`) | Operation |
|------|----------------|-------------------|-----------|
| `no2.csv` | `no2_tint.F90` | `cross_section.no2_tint.nc` | clamped linear interpolation in temperature |

NO2 is tabulated at 220 K and 294 K (wavelengths 300.7685 nm, 305.361 nm). The
weight at temperature `T` is the linear interpolation between the two bracketing
reference temperatures, with `T` clamped to `[220, 294]` (constant extrapolation
outside). Sampled temperatures 210/257/294/300 K exercise clamp-low, midpoint
interpolation, an exact node, and clamp-high.

The C++ side uses the `temperature_table(reference_temperatures, cross_sections)`
form. Because the `no2_tint` stub only varies mildly across temperature, the
interpolation math itself is validated separately by unit tests in
`test/unit/transforms/test_factories.cpp` (per-wavelength interpolation, bracket
selection, and end-clamping with non-degenerate synthetic data).

The other "tint"-type species (`tint`, `o3_tint`) and `oclo` are deferred:
`oclo`'s test exercises **wavelength** extrapolation modes (constant/boundary),
which belong with the conserving-interpolator work (Phases 4-5) rather than
temperature interpolation.
