# Phase 7: Discrete Ordinate Solver (Step 22)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The Discrete Ordinate solver (DISORT) is a multi-stream radiative transfer method supporting 2–32 even streams. It is more accurate than the two-stream Delta Eddington but computationally heavier. The Fortran implementation is ~3,700 lines. This phase should only begin after the Delta Eddington solver and full pipeline are working end-to-end.

## Step 22: Port the DISORT solver

Port `src/radiative_transfer/solvers/discrete_ordinant_util.F90` — the `PSNDO` subroutine and supporting functions.

### Key characteristics
- Multi-stream method: supports 2, 4, 8, 16, 32 streams (even numbers)
- Much more computationally intensive than Delta Eddington — stronger candidate for GPU acceleration
- Well-known community code with extensive literature
- The Fortran wrapper is in `src/radiative_transfer/solvers/discrete_ordinant.F90` (~282 lines)

### Implementation approach
- Same `ArrayPolicy` template system as Delta Eddington
- Same `Solver` API — the solver type is selected at construction, not runtime
- GPU acceleration will have the most impact here due to the computational intensity
- Consider batching strategies: multiple columns can share the angular quadrature setup

### Fortran reference
- `src/radiative_transfer/solvers/discrete_ordinant.F90` — wrapper (~282 lines)
- `src/radiative_transfer/solvers/discrete_ordinant_util.F90` — full DISORT implementation (~3,700 lines, `PSNDO` subroutine)
- `src/radiative_transfer/solver_factory.F90` — dispatches on `"discrete ordinate"` type string

### Documentation

Document the DISORT solver's public interface and key internal steps with Doxygen `///` comments. Reference the Fortran `PSNDO` subroutine and relevant literature (Stamnes et al., 1988). Note how it differs from Delta Eddington in stream count and accuracy.
