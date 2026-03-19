# Phase 1: Delta Eddington Solver (Steps 4–6)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The Delta Eddington solver is the primary two-stream radiative transfer solver used for photolysis rate calculations. The Fortran implementation is ~468 lines. The C++ stub exists but `Solve()` is a placeholder returning dummy values. The algorithm follows Toon et al. (1989).

## Step 4: Implement spherical geometry utilities

Port from Fortran sources:
- `src/spherical_geometry.F90` — pseudo-spherical geometry corrections
- `src/radiative_transfer/solver.F90` — `slant_optical_depth` function

These compute the effective solar path through curved atmospheric layers accounting for Earth's curvature. Implement as pure functions (no state). All inputs/outputs in SI units (meters, radians).

## Step 5: Implement Delta Eddington solver

Replace the placeholder in `include/tuvx/radiative_transfer/solvers/delta_eddington.inl` with the actual algorithm from `src/radiative_transfer/solvers/delta_eddington.F90`.

### Algorithm (6 steps from Toon et al. 1989)

**Step 1 — Delta-scale and compute gamma coefficients** (per layer, per wavelength, per column):

$$\gamma_1^n = \frac{7 - \omega_0(4 + 3g_n)}{4}, \quad \gamma_2^n = \frac{-1 - \omega_0(4 - 3g_n)}{4}$$
$$\gamma_3^n = \frac{2 - 3g_n\mu_0}{4}, \quad \gamma_4^n = 1 - \gamma_3^n$$
$$\lambda_n = \sqrt{\gamma_1^2 - \gamma_2^2}, \quad \Gamma_n = \frac{\gamma_1^n - \lambda_n}{\gamma_2^n}$$

**Step 2 — Solar source functions** $C^+(\tau)$ and $C^-(\tau)$:

$$C^+(\tau) = \omega_0 \pi F_s e^{-\tau/\mu_0} \frac{(\gamma_1 - 1/\mu_0)\gamma_3 + \gamma_4\gamma_2}{\lambda^2 - 1/\mu_0^2}$$

**Step 3 — Tridiagonal system coefficients** ($e_1$ through $e_4$; $A_n$, $B_n$, $D_n$ diagonals)

**Step 4 — RHS assembly** ($E_n$) with surface albedo boundary condition

**Step 5 — Solve tridiagonal system** using existing `TridiagonalMatrix::Solve()` (Thomas algorithm)

**Step 6 — Back-substitute** $Y$ to compute fluxes and actinic flux:
$$F_n^+ = Y_1^n e_1^n + Y_2^n e_2^n + C_n^+(\tau_n)$$
$$F_n^- = Y_1^n e_3^n + Y_2^n e_4^n + C_n^-(\tau_n)$$

### Implementation requirements

- Template on `ArrayPolicy` — same algorithm for CPU and GPU
- Use the **bulk operations API** (`ForEachRow`, `ColumnView`, `Function` factory) for all hot-path computation — no explicit loops over columns
- Create `ArrayPolicy::Function` objects for the per-timestep solver steps (gamma computation, source functions, tridiagonal assembly, back-substitution) so the policy can pre-compile optimized iteration patterns
- Use `GetColumnView()` / `GetConstColumnView()` for accessing per-column data within `ForEachRow` lambdas
- Use `GetRowVariable()` for intermediate per-row temporaries ($\lambda_n$, $\Gamma_n$, $C^+$, $C^-$)
- Process batches of columns; the `ArrayPolicy` determines the optimal data layout for parallelization across columns
- All units in SI (meters, radians, seconds)
- Use existing `TridiagonalMatrix::Solve()` for step 5
- Refer to [issue #64](https://github.com/NCAR/tuv-x/issues/64) sub-issues #102, #103, #104 for step-by-step breakdown

### Key Fortran files to reference

- `src/radiative_transfer/solvers/delta_eddington.F90` — the Fortran implementation (~468 lines)
- `src/radiative_transfer/solver.F90` — base solver class with `slant_optical_depth`
- `src/spherical_geometry.F90` — spherical geometry corrections
- `src/linear_algebras/linpack.F90` — LINPACK tridiagonal solver (replaced by C++ Thomas algorithm)

## Step 6: Regression test against Fortran solver

Adapt the existing test harness in `test/regression/solvers/`:
- `test/regression/solvers/delta_eddington.hpp` — C-interop structs (`SolverInput`/`SolverOutput`)
- `test/regression/solvers/delta_eddington.cpp` — C++ wrapper that creates objects and calls solver
- `test/regression/solvers/delta_eddington.F90` — Fortran test that runs full pipeline

For the new library: pre-compute Fortran reference outputs and store as binary/NetCDF fixtures. The C++ tests load these fixtures and compare solver output within floating-point tolerance, decoupling from the Fortran build entirely.

Test configurations should include the standard examples:
- `examples/tuv_5_4.json` — standard 156-wavelength, 120-layer US standard atmosphere
- Multiple solar zenith angles (the existing Fortran test uses 1.8294° and 28.199°; C++ tests must specify these in radians: ~0.03194 rad and ~0.4922 rad)
- Multiple columns (verify multi-column batching produces identical per-column results)

### Documentation

Document all public functions (spherical geometry utilities, `RadiationField`, `Solve()`) with Doxygen `///` comments as they are written. For each ported algorithm, reference the Fortran subroutine name and source file. Keep comments succinct — state what the function does and its preconditions, not a line-by-line translation of the algorithm.
