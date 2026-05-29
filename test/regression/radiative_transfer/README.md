# Radiative Transfer Regression Reference Data

## Delta-Eddington

The `delta_eddington_sza_*.csv` files contain the expected output of the
delta-Eddington two-stream solver for a synthetic 3-layer atmosphere:

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

---

## Spherical Geometry

The `spherical_geometry_sza_*.csv` files contain the expected output of
`SphericalGeometry::SetParameters()` and `SlantOpticalDepth()` for the same
3-layer atmosphere:

| Parameter | Value |
|-----------|-------|
| Layers | 3 |
| Altitude edges | 0, 1, 2, 3 km |
| Layer optical depths (top-to-bottom) | 0.1, 0.2, 0.3 |

One file is produced per SZA: 0.0, 0.5, 1.0, 1.5, and 1.6 radians.
SZA > π/2 (≈ 1.5708 rad) places lower levels below the tangent height; those
rows have `nid = -1` and `slant_od = Inf`.

Each file has columns:
`level, nid, dsdh_0, dsdh_1, dsdh_2, slant_od`

- `nid`: number of layers crossed by the direct beam (`-1` = below tangent height)
- `dsdh_j`: slant/vertical-depth ratio for layer `j` (top-to-bottom, 0-based)
- `slant_od`: total slant-path optical depth using the layer optical depths above

Level 0 is the top of atmosphere (TOA); level `n_layers` is the ground.

### How the data were generated

The reference CSVs were generated from the C++ implementation on the
`phase1-delta-eddington` branch at commit `f940b0aa5e2d3db8470acb9efb217a3cee69dd4e`.
The implementation had passed the full delta-Eddington Fortran regression tests at
that point, establishing correctness of the spherical geometry.

### Regenerating the reference data

1. Check out the commit above:

   ```
   git checkout f940b0aa5e2d3db8470acb9efb217a3cee69dd4e
   ```

2. Write a small C++ program that:
   - Constructs `tuvx::SphericalGeometry` and calls `SetParameters` with
     `alt_edges = {0, 1000, 2000, 3000}` m and each SZA in radians.
   - For each level, records `nid_[level]` (-1 for `nullopt`), all `dsdh_[level][j]`
     values, and `SlantOpticalDepth(level, nid_[level], dsdh_[level], taun)` with
     `taun = {0.1, 0.2, 0.3}`.
   - Writes one CSV per SZA named `spherical_geometry_sza_X.XXXX.csv` (17 significant
     digits for floating-point, `Inf` for infinity).

3. Compile (C++17):

   ```
   g++ -std=c++17 -O2 -I include generate_spherical_reference.cpp -o gen_ref
   ./gen_ref
   ```

4. Move the generated CSVs into `reference/` and commit.
