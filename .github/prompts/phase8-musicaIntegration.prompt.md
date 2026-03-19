# Phase 8: MUSICA Integration (Step 23)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

MUSICA (Multi-Scale Infrastructure for Chemistry and Aerosols) is the parent library that combines TUV-x with other chemistry tools (MICM, CARMA, MechanismConfiguration) and exposes them through Python, JavaScript, Julia, and Fortran interfaces. The new tuv-x C++ library must integrate seamlessly.

### MICM API Consistency

TUV-x and MICM are both MUSICA components and should present a consistent interface to users working across both libraries. The following conventions are shared:

| Convention | MICM | TUV-x | Notes |
|-----------|------|-------|-------|
| Element access | `matrix[row][col]` | `array[i][j]`, `array[i][j][k]` | Square-bracket syntax via proxy objects |
| Bulk iteration | `matrix.ForEachRow(lambda, views...)` | `array.ForEachRow(lambda, views...)` | Policy controls vectorization |
| Column access | `GetColumnView(col)`, `GetConstColumnView(col)` | Same | Read-only and mutable views |
| Temporaries | `GetRowVariable()` | Same | Per-row intermediate storage |
| Function factory | `MatrixPolicy<T>::Function(lambda, args...)` | `ArrayType<T>::Function(lambda, args...)` | Policy-compiled reusable operations |
| Dimension queries | `NumRows()`, `NumColumns()` | `NumRows()`, `NumColumns()` | Matching names |
| Template policy | `template<class> class MatrixPolicy` | `typename ArrayPolicy` | Same pattern, different naming (TUV-x uses existing convention) |
| `std::vector` in APIs | Used alongside `MatrixPolicy` | **Replaced by `Array1D`** | TUV-x is stricter — no `std::vector` in policy-adjacent APIs; MICM may adopt this later |

Diverge from MICM conventions **only** when TUV-x requirements demand it.

## Step 23: Update MUSICA to consume the new library

### CMake integration

MUSICA uses FetchContent to pull tuv-x. Update the dependency variables:

```cmake
# In MUSICA's cmake/dependencies.cmake
FetchContent_Declare(
    tuvx
    GIT_REPOSITORY ${TUVX_GIT_REPOSITORY}  # point to new repo
    GIT_TAG ${TUVX_GIT_TAG}
)
```

The new library should export a `musica::tuvx` CMake target (matching current convention).

### Configuration layer (built in MUSICA, not in tuv-x)

MUSICA builds the bridge between YAML/JSON configuration and the pure C++ API:

```cpp
// In MUSICA's code, not in tuv-x
TransformFunc<HostArrayPolicy> build_cross_section(const YAML::Node& config) {
    std::string type = config["type"].as<std::string>();
    if (type == "tint") {
        auto reader = NetCDFReader(config["netcdf files"][0]["file path"]);
        return temperature_interpolation(reader);
    } else if (type == "CCl4+hv->Products") {
        // ... build transform using multiply/add/piecewise per config
    }
    // etc.
}
```

This maps the existing config types (`"tint"`, `"O3"`, `"temperature based"`, etc.) to transforms built with factories and combinators.

### Language bindings

- **Python** (PyBind11): Wrap the C++ `Solver` and `PhotolysisCalculator` classes. Transforms can be specified as Python callables or by name (MUSICA maps names to C++ lambdas).
- **JavaScript** (WASM): Compile C++ solver to WebAssembly via Emscripten. Only CPU path.
- **Julia**: `ccall` into the C API wrapper.
- **Fortran**: `iso_c_binding` into the C API wrapper (for legacy model integration).

### Migration path

1. New tuv-x library reaches feature parity with Fortran for Delta Eddington solver
2. MUSICA adds the new library as an alternative dependency alongside the old one
3. Regression tests verify identical output
4. MUSICA switches default to new library
5. Old Fortran tuv-x repo archived (kept for reference)

### Key MUSICA references
- [MUSICA repo](https://github.com/NCAR/musica) — `cmake/dependencies.cmake` for FetchContent setup
- `include/musica/` — C++ API headers
- `python/` — PyBind11 bindings
- `javascript/` — WASM/Emscripten build
- `julia/` — Julia wrapper
- `fortran/` — Fortran interface via `iso_c_binding`

### Existing TUV-x integration in MUSICA
- `tutorials/` — TUV-x tutorials already exist in MUSICA
- Current MUSICA builds tuv-x as `musica::tuvx` target
- MUSICA uses `TUVX_GIT_REPOSITORY` and `TUVX_GIT_TAG` CMake variables for version control
