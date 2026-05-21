# Delta-Eddington Regression Reference Data

The CSV files in `reference/` contain the expected output of the delta-Eddington
two-stream solver for a synthetic 3-layer atmosphere:

| Parameter | Value |
|-----------|-------|
| Layers | 3 |
| Optical depth per layer (τ) | 0.5 |
| Single-scattering albedo (ω) | 0.9 |
| Asymmetry parameter (g) | 0.85 |
| Surface albedo | 0.1 |
| Altitude edges | 0, 1, 2, 3 km |

One file is produced per solar zenith angle (SZA): 0.0, 0.5, and 1.0 radians.

Each file has columns:
`level, direct_irradiance, upwelling_irradiance, downwelling_irradiance,
direct_actinic_flux, upwelling_actinic_flux, downwelling_actinic_flux`

Level 0 is the ground; level `n_layers` is the top of atmosphere.

## How the data were generated

The reference data were computed by a standalone Fortran program that extracted
the core algorithm — with no framework dependencies — from three source files on
the `main` branch at commit `c2a8f3385c45df08dd65a0d2503665182332b783`:

| File | Role |
|------|------|
| `src/radiative_transfer/solvers/delta_eddington.F90` | Two-stream solver kernel |
| `src/spherical_geometry.F90` | Spherical-geometry slant-path (`set_parameters`) |
| `src/linear_algebras/linpack.F90` | Tridiagonal solver (`tridiag`) |

## Regenerating the reference data

1. Check out the `main` branch (or the commit above):

   ```
   git checkout c2a8f3385c45df08dd65a0d2503665182332b783
   ```

2. Copy the three source files listed above to a working directory.

3. Write a standalone Fortran program that:
   - Hardcodes the test atmosphere described in the table above.
   - Calls `set_parameters` (from `spherical_geometry.F90`) with the altitude
     grid and each SZA.
   - Calls the delta-Eddington inner loop (from `delta_eddington.F90`) with the
     delta-scaled optical properties and the spherical-geometry output.
   - Solves the tridiagonal system with `tridiag` (from `linpack.F90`).
   - Reverses the output arrays from top-to-bottom to bottom-to-top.
   - Writes one CSV per SZA named `delta_eddington_sza_X.XXXX.csv`.

4. Compile and run:

   ```
   gfortran -O2 -o generate_reference generate_reference.F90
   ./generate_reference
   ```

5. Move the generated CSVs into `reference/` and commit.
