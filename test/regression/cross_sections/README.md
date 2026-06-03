# Cross-Section Regression Reference Data

## Analytic cross-sections (PR A)

The five CSV files in `reference/` contain cross-section values computed by the
standalone Fortran program `generate_analytic.F90`, which implements the same
formulas as the Fortran source files on the `main` branch without any TUV-x or
NetCDF dependencies.

| File | Source formula | Active range |
|------|---------------|--------------|
| `rayleigh.csv` | `main:src/cross_sections/rayliegh.F90` | All wavelengths |
| `hobr.csv` | `main:src/cross_sections/hobr-oh_br.F90` | 250-550 nm |
| `t_butyl_nitrate.csv` | `main:src/cross_sections/t_butyl_nitrate.F90` | 270-330 nm |
| `nitroxy_acetone.csv` | `main:src/cross_sections/nitroxy_acetone.F90` | 284-335 nm |
| `nitroxy_ethanol.csv` | `main:src/cross_sections/nitroxy_ethanol.F90` | 270-306 nm |

Each file has columns: `wavelength_m, cross_section_m2`

Wavelength grid: 26 mid-points from 200 nm to 450 nm in 10 nm steps.
Values are in SI units (m²); the Fortran formulas output cm² which the
generator converts with a factor of 1e-4.

## Regenerating the reference data

```
cd test/regression/cross_sections
gfortran -O2 -o generate_analytic generate_analytic.F90
./generate_analytic
```

The five CSVs are written to the current directory; move them to `reference/`.
