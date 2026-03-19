# Phase 6: Public API (Steps 20–21)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The public API replaces the Fortran `core_t` class. It is purely programmatic — no config files. Configuration mapping (YAML/JSON → API calls) is MUSICA's responsibility. The API must support both C++ direct use and FFI consumption via a C wrapper.

## Step 20: Define the top-level solver API

### C++ API design

```cpp
namespace tuvx {

/// Atmospheric column state for a batch of columns
template<typename ArrayPolicy>
struct AtmosphericState {
    Grid<ArrayPolicy> altitude_grid;           // per-column vertical grid [m]
    Grid<ArrayPolicy> wavelength_grid;         // shared wavelength grid [m]
    Profile<ArrayPolicy> temperature;          // temperature profile [K]
    Profile<ArrayPolicy> pressure;             // pressure profile [Pa]
    Profile<ArrayPolicy> air_density;          // air number density [molec/m³]
    Array1D<typename ArrayPolicy::value_type> solar_zenith_angles;   // per-column [rad]
    Array1D<typename ArrayPolicy::value_type> earth_sun_distances;   // per-column [AU]
    typename ArrayPolicy::value_type surface_albedo;
};

/// Solver for radiative transfer
template<typename ArrayPolicy>
class Solver {
public:
    void solve(
        const AtmosphericState<ArrayPolicy>& state,
        const RadiatorState<ArrayPolicy>& radiator_state,
        RadiationField<ArrayPolicy>& radiation_field
    );
};

/// Full pipeline: radiation field → photolysis/dose/heating rates
template<typename ArrayPolicy>
class PhotolysisCalculator {
public:
    void add_reaction(
        std::string name,
        TransformFunc<ArrayPolicy> cross_section,
        TransformFunc<ArrayPolicy> quantum_yield,
        double scaling_factor = 1.0
    );
    void calculate(
        const RadiationField<ArrayPolicy>& field,
        const AtmosphericState<ArrayPolicy>& state,
        Array3D<typename ArrayPolicy::value_type>& rates  // [reaction × height × column]
    );
};

} // namespace tuvx
```

### Design principles
- Template on `ArrayPolicy` — same API for CPU and GPU
- Batch-oriented: always operates on collections of columns via `ForEachRow`/`Function` bulk operations
- No internal state beyond registered reactions — thread-safe after setup
- All units in SI
- **MICM-consistent API conventions**: square-bracket element access (`array[i][j]`), `ForEachRow`/`ColumnView`/`Function` factory, `NumRows()`/`NumColumns()` queries — consistent with MICM's `Matrix`/`VectorMatrix` types since both are MUSICA components

### Fortran reference
- `src/core.F90` — `core_t` class: the current top-level API with `run()`, `get_grid()`, `get_profile()`, `get_photolysis_cross_section()`, etc.

## Step 21: C API wrapper

Provide `extern "C"` functions for FFI consumption:

```c
// Opaque handle types
typedef struct tuvx_solver* tuvx_solver_t;
typedef struct tuvx_calculator* tuvx_calculator_t;
typedef struct tuvx_radiation_field* tuvx_radiation_field_t;

// Lifecycle
tuvx_solver_t tuvx_create_solver();
void tuvx_destroy_solver(tuvx_solver_t solver);

// Solve
void tuvx_solve(tuvx_solver_t solver, /* flat array inputs */);

// Results
void tuvx_get_radiation_field(tuvx_solver_t solver, tuvx_radiation_field_t field);
void tuvx_get_photolysis_rates(tuvx_calculator_t calc, double* rates, int* dims);
```

This enables:
- MUSICA's Fortran interface via `iso_c_binding`
- Python bindings via ctypes (fallback) or direct PyBind11 wrapping of C++ API
- JavaScript/WASM compilation
- Julia `ccall`

### Reference
- `test/regression/solvers/delta_eddington.hpp` — existing C-interop struct pattern (`SolverInput`, `SolverOutput`)
