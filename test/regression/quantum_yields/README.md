# Quantum-Yield Regression Reference Data

## Pure-wavelength quantum yields (PR 1)

The CSV files in `reference/` contain **unitless, wavelength-only**
photodissociation branching fractions. They were produced by a standalone
program evaluating the same formulas as the Fortran `src/quantum_yields/`
sources, transcribed from `main` at commit `c2a8f33`.

Each file has columns: `wavelength_m, quantum_yield`. Wavelength grid: 11
mid-points bracketing the thresholds (248 / 263.4 / 308 / 364 nm) without
landing exactly on one (the m<->nm round-trip could flip a strict comparison).

| File | Fortran source | Formula (`w` = wavelength nm) |
|------|----------------|-------------------------------|
| `clo_cl_o1d.csv` | `clo-cl_o1d.F90` | 1 if `w < 263.4`, else 0 |
| `clo_cl_o3p.csv` | `clo-cl_o3p.F90` | 0 if `w < 263.4`, else 1 |
| `clono2_cl_no3.csv` | `clono2-cl_no3.F90` | 0.6 (`w<308`); `7.143e-3 w - 1.6` (`308<=w<=364`); 1 (`w>364`) |
| `clono2_clo_no2.csv` | `clono2-clo_no2.F90` | `1 -` the Cl+NO3 branch |
| `ho2_oh_o.csv` | `ho2-oh_o.F90` | 1 (`w>=248`); else `max(0, (1 + 14(w-193)/55)/15)` |

All are pure functions of wavelength (no temperature / pressure / air-density
dependence) and build directly from `wrap_analytic` -- no new library code.

### Deferred

- Temperature- or air-density-dependent quantum yields (`no3_aq`, `mvk`,
  `ch3cocho`, `ch3coch2ch3`, `ch2chcho`, `o3-o2_o1d`/`o3-o2_o3p`, `acetone`, ...)
  follow in later batches (some use `parameterized` for air density; O3 and
  acetone are the complex Matsumi / Blitz forms).
- `tint`, `no2_tint`, `taylor_series`, `h2so4_mills`, and the external-array
  yields (`c2h5cho`, `ch3cho`, `ch2o`) read tabulated data -> wait on the NetCDF
  reader (Phases 4-5).
