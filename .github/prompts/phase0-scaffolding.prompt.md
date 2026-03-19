# Phase 0: Project Scaffolding (Steps 1–3)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

We are creating a new standalone C++ library (`tuv-x-cpp`) to replace the Fortran-based TUV-x photolysis rate calculator. This phase sets up the new repo, migrates existing complete C++ code, and establishes the template policy system for CPU/GPU portability.

## Step 1: Create new repository

Create a new repo with:
- CMake 3.21+ build system (matching MUSICA requirements)
- Languages: `CXX` with optional `CUDA`/`HIP` support
- Google Test and Google Benchmark as test/benchmark dependencies (FetchContent)
- NetCDF-C as an optional dependency (for data readers)
- Project name: `tuv-x` or `tuv-x-cpp`, aliased as `musica::tuvx`
- Apache-2.0 license (matching MUSICA)

Directory structure:
```
include/tuvx/          # Public headers
src/                    # Implementation files
test/unit/              # GTest unit tests
test/regression/        # Regression tests against Fortran reference data
benchmark/              # Google Benchmark files
cmake/                  # CMake modules (dependencies, test utilities, coverage)
data/                   # Reference data files (NetCDF)
docker/                 # Dockerfiles for coverage, memcheck, compiler variants
.github/workflows/      # CI workflow files
.clang-format           # Formatting rules
.clang-tidy             # Static analysis rules
```

### CI/CD Requirements (GitHub Actions)

All of the following must be in place before any feature work begins.

#### 1. Multi-compiler test matrix (PR gate — must pass to merge)

Tests must pass with all major compilers on all three platforms:

| OS | Compilers |
|----|-----------|
| Ubuntu (24.04) | GCC 12, 13, 14; Clang 17, 18 |
| macOS (latest) | Apple Clang (Xcode default); GCC 14, 15 |
| Windows (latest) | MSVC (latest via Visual Studio); MinGW GCC |
| Docker (Linux) | Intel icx/icpx (oneAPI); NVIDIA nvc/nvc++ (NVHPC) |

Build types: both `Debug` and `Release` for each compiler.

#### 2. Code coverage (PR gate — >95% required)

- Use lcov/gcov for GCC builds or llvm-cov for Clang builds
- Coverage runs in a Docker container (like current `Dockerfile.coverage` pattern)
- Upload to Codecov via `codecov/codecov-action`
- **PR check must fail if line coverage drops below 95%** — configure via `codecov.yml`:
  ```yaml
  coverage:
    status:
      project:
        default:
          target: 95%
          threshold: 1%
      patch:
        default:
          target: 95%
  ```
- CMake flag: `-D TUVX_ENABLE_COVERAGE:BOOL=TRUE -D CMAKE_BUILD_TYPE=COVERAGE`

#### 3. Memory checking with Valgrind (PR gate — must pass clean)

- All unit and regression tests run under Valgrind memcheck
- Runs in a Docker container (like current `Dockerfile.memcheck` pattern)
- CMake flag: `-D TUVX_ENABLE_MEMCHECK:BOOL=TRUE -D CMAKE_BUILD_TYPE=DEBUG`
- CTest wraps test executables with `valgrind --leak-check=full --error-exitcode=1`
- **Any memory error fails the PR**

#### 4. Static analysis with clang-tidy (PR gate — must pass clean)

- Runs on every PR (unlike current repo where it's disabled)
- Uses `.clang-tidy` config file at repo root
- Recommended checks to enable:
  ```yaml
  Checks: >
    bugprone-*,
    cert-*,
    clang-analyzer-*,
    cppcoreguidelines-*,
    misc-*,
    modernize-*,
    performance-*,
    readability-*,
    -modernize-use-trailing-return-type,
    -readability-magic-numbers,
    -cppcoreguidelines-avoid-magic-numbers
  ```
- Run against `include/` and `src/` using the compilation database (`compile_commands.json`)
- **Any clang-tidy warning fails the PR** — no warnings-as-errors bypass

#### 5. Auto-formatting PRs via clang-format (push to main only)

- Triggered on pushes to `main` (not on PRs — to avoid circular triggers)
- Runs `clang-format -i --style=file` on all `.hpp`, `.h`, `.cpp`, `.cc` files in `include/`, `src/`, `test/`
- If changes are detected:
  1. Commits with message `"Auto-format code using clang-format"`
  2. Pushes to a `main-formatting` branch
  3. Opens a PR via `peter-evans/create-pull-request@v6` back to `main`
- `.clang-format` config should enforce the project style (recommend: based on LLVM or Google style, with project-specific overrides)

#### 6. Concurrency and workflow organization

- All workflows use concurrency groups to cancel superseded runs:
  ```yaml
  concurrency:
    group: ${{ github.workflow }}-${{ github.ref }}
    cancel-in-progress: true
  ```
- Workflow files:
  - `ubuntu.yml` — GCC + Clang matrix on Ubuntu
  - `mac.yml` — Apple Clang + GCC on macOS
  - `windows.yml` — MSVC + MinGW on Windows
  - `docker.yml` — Intel, NVHPC, coverage, memcheck containers
  - `clang-format.yml` — auto-formatting PRs
  - `clang-tidy.yml` — static analysis PR gate

### Naming Conventions

The Fortran codebase is full of cryptic abbreviations (`edr_`, `eup_`, `fdr_`, `fdn_`, `sirrad`, `xsI`, `nTemp`, `vertNdx`, `rateNdx`, `Tadj`, `PSNDO`, `LINPACK`, etc.). The new codebase must prioritize legibility:

- **Variables**: descriptive snake_case — `direct_irradiance`, `upwelling_flux`, `temperature_index`, `optical_depth_scaled`
- **Types/Classes**: descriptive PascalCase — `RadiationField`, `TridiagonalMatrix`, `SphericalGeometry`
- **Methods**: descriptive snake_case or camelCase (pick one, enforce via clang-tidy `readability-identifier-naming`)
- **Template parameters**: PascalCase — `ArrayPolicy`, `ValueType`
- **Constants**: `kPascalCase` or `SCREAMING_SNAKE_CASE`
- **No single-letter variables** except loop indices (`i`, `j`, `k`) and well-known math symbols in formulas (e.g., `tau`, `omega`, `mu`, `gamma` — but spelled out, not `t`, `w`, `m`, `g`)
- **No abbreviations** unless universally understood in the domain (e.g., `SSA` for single-scattering albedo is acceptable; `xs` for cross-section is not — use `cross_section`)
- Configure `readability-identifier-naming` in `.clang-tidy` to enforce these conventions:
  ```yaml
  CheckOptions:
    - key: readability-identifier-naming.ClassCase
      value: CamelCase
    - key: readability-identifier-naming.FunctionCase
      value: camelBack
    - key: readability-identifier-naming.VariableCase
      value: lower_case
    - key: readability-identifier-naming.ConstantCase
      value: UPPER_CASE
    - key: readability-identifier-naming.TemplateParameterCase
      value: CamelCase
  ```

## Step 2: Migrate reusable C++ code

Port these already-complete implementations from the current `tuv-x` repo:

### Data structures
- **`Array1D<T>` — new implementation** to create in `include/tuvx/util/array1d.hpp`, following the same pattern as `Array2D`/`Array3D`. This is a policy-aware 1D array that replaces `std::vector` in all interfaces that interact with policy-backed types (e.g., `Grid`, `Profile`, `RadiatorState`). When backed by `DeviceArrayPolicy`, data resides in device memory — `std::vector` cannot serve this role. Must support the same element access syntax (`operator[]`), bulk operations (`ForEachRow`, `Function` factory), and dimension queries (`NumRows()`) as `Array2D`/`Array3D` via `ArrayPolicy`. Note: for `Array1D`, `ForEachRow` iterates element-wise (each element is a "row").
- `Array2D<T>` — `include/tuvx/util/array2d.hpp`
- `Array3D<T>` — `include/tuvx/util/array3d.hpp`
- `Grid<ArrayPolicy>` — `include/tuvx/grid.hpp`
- `Profile<ArrayPolicy>` — `include/tuvx/profile.hpp`
- `RadiatorState<ArrayPolicy>` — `include/tuvx/radiative_transfer/radiator.hpp`
- `RadiationField`, `RadiationFieldComponents` — `include/tuvx/radiative_transfer/radiation_field.hpp`

### Linear algebra
- `TridiagonalMatrix<T>`, `Solve()` — `include/tuvx/linear_algebra/linear_algebra.hpp`

### Tests and benchmarks
- GTest tests for Array1D (new), Array2D, Array3D, Grid, Profile, TridiagonalSolver, ErrorFunction
- Google Benchmark for tridiagonal solver — `benchmark/benchmark_tridiagonal_solver.cpp`

Ensure all migrated tests pass in the new repo before proceeding.

## Step 3: Define template policy system

Extend the existing `ArrayPolicy` pattern. The key principle: **solver algorithms are written once**, templated on a single `ArrayPolicy`. From the algorithm's perspective, array access is identical regardless of backing storage.

`ArrayPolicy` controls:
- **Memory allocation**: where data lives (host heap, device memory)
- **Data layout**: how elements are physically arranged in memory (implementation-defined; optimized per target)
- **Element access**: square-bracket syntax (`array[i][j]`, `array[i][j][k]`) for reading/writing by logical index (identical API for host and device)

Concrete policies to implement:
- `HostArrayPolicy` — `std::vector`-backed; layout optimized for SIMD auto-vectorization
- `DeviceArrayPolicy` — CUDA/HIP device memory; layout optimized for coalesced access; `__host__ __device__` accessors

All core types (`Array1D`, `Grid`, `Profile`, `RadiatorState`, `RadiationField`, `TridiagonalMatrix`) are templated on `ArrayPolicy` — preserve and extend this pattern. `Array1D<T>` is new and must be created following the same conventions as `Array2D<T>` and `Array3D<T>`. **No `std::vector` should appear in any function signature that also accepts policy-backed types** — use `Array1D` instead to ensure correctness with device-resident data.

### MICM-Consistent Bulk Operations API

TUV-x array types must provide the same bulk operations API as MICM's `Matrix`/`VectorMatrix` types, since both are MUSICA components and a consistent interface reduces user friction. The API consists of:

1. **Square-bracket element access**: `operator[]` chaining for logical indexing (e.g., `array[i][j]` for 2D). Returns proxy/view objects. Must match MICM's `matrix[row][col]` convention.

2. **`ForEachRow(lambda, views...)`**: Layout-agnostic bulk iteration over all elements in a row (the solver's parallelization dimension). The policy controls vectorization/threading. Arguments are column views (`GetColumnView(col)` / `GetConstColumnView(col)`) and/or `Array1D` vectors.

3. **`GetRowVariable()`**: Temporary per-row storage for intermediate calculations within `ForEachRow` chains.

4. **`ArrayType<T>::Function(lambda, prototype_args...)`**: Factory that creates a reusable, policy-compiled function object. Pre-computes layout metadata at creation time. Column counts are validated at invocation; row counts may differ from the prototype. Matches MICM's `MatrixPolicy<T>::Function()` factory.

5. **`NumRows()`, `NumColumns()`**: Dimension queries matching MICM naming.

Solver hot paths **must** use `ForEachRow`/`Function` rather than element-by-element `operator[]` access — the bulk API is how the policy achieves vectorization/coalescing.

## Existing C++ code to reference

Key files in the current tuv-x repo:
- `include/tuvx/util/array1d.hpp` — policy-aware 1D array (**new — to be created**)
- `include/tuvx/util/array2d.hpp` — 2D array with `operator[][]`, `ForEachRow`, `ColumnView`, `Function` factory (layout determined by policy)
- `include/tuvx/util/array3d.hpp` — 3D array with `operator[][][]`, `ForEachRow`, `ColumnView`, `Function` factory (layout determined by policy)
- `include/tuvx/grid.hpp` — multi-column grid with mid_points/edges
- `include/tuvx/profile.hpp` — mid-point/edge values on a grid
- `include/tuvx/radiative_transfer/radiator.hpp` — RadiatorState with optical_depth, SSA, asymmetry
- `include/tuvx/radiative_transfer/radiation_field.hpp` — 3D arrays for direct/up/down flux
- `include/tuvx/linear_algebra/linear_algebra.hpp` — Thomas algorithm tridiagonal solver
- `test/unit/` — existing GTest files for the above
