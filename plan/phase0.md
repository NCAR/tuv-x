# Phase 0: Project Scaffolding — Revised Implementation Plan

> **Branch**: `phase0-scaffolding` (off `cpp-rewrite`, off `main`)
> **Original plan**: `.github/prompts/phase0-scaffolding.prompt.md` (unchanged)
> **Goal**: Set up a clean C++ project structure with CMake, CI, and the existing C++ data structures — while preserving all Fortran source as a living reference.

---

## Terminology: Radiator → Constituent

The Fortran codebase uses "radiator" — borrowed from automotive engineering. In radiative transfer theory, the correct term is **constituent** (a species that absorbs or scatters radiation). All C++ code renames: `radiator.hpp` → `constituent.hpp`, `RadiatorState` → `ConstituentState`, test directories, etc. The rename happens incrementally as files are touched, not as a bulk pass.

---

## Fortran Preservation Policy

The original plan calls for removing Fortran source on the `cpp-rewrite` branch. **We are revising this.** Fortran source is the primary reference for translation work in Phases 1–8. Deleting it early means we'd constantly be doing `git show main:src/...` to read originals, which is awkward and error-prone.

**Revised approach:**
- The current `src/` directory is **renamed to `fortran/`** — a clear signal that it's legacy reference code.
- `fortran/` is **not built, not tested, not linted** — it's reference material only.
- Fortran files remain in the tree for reading, searching, and side-by-side comparison throughout all phases.
- A final cleanup phase (after Phase 8) removes `fortran/`, legacy configs, and anything else that's no longer needed.
- `src/` is reclaimed as a clean directory for new C++ implementation files.
- `include/tuvx/` remains the home for public C++ headers (unchanged location).

---

## Directory Structure (End State of Phase 0)

```
fortran/                               # RENAMED from src/ — Fortran reference code, NOT BUILT
  core.F90, tuvx.F90, constants.F90 ...  (13 top-level .F90 files)
  cross_sections/                      # 31 .F90 files
  quantum_yields/                      # 20 .F90 files
  spectral_weights/                    # 13 .F90 files
  radiative_transfer/                  # solvers/, constituents/ (renamed from radiators/)
  grids/                               # 4 .F90 files
  profiles/                            # 10 .F90 files
  linear_algebras/                     # 2 .F90 files
  util/                                # 11 .F90 files

include/tuvx/                          # Public C++ headers (existing, reorganized)
  util/
    array1d.hpp                        # NEW — policy-aware 1D array
    array2d.hpp                        # EXISTING
    array3d.hpp                        # EXISTING
  grid.hpp                             # EXISTING
  profile.hpp                          # EXISTING
  linear_algebra/
    linear_algebra.hpp                 # EXISTING
    linear_algebra.inl                 # EXISTING
  radiative_transfer/
    constituent.hpp                    # EXISTING (renamed from radiator.hpp)
    radiation_field.hpp                # EXISTING
    solvers/
      delta_eddington.hpp              # EXISTING (placeholder)
      delta_eddington.inl              # EXISTING (placeholder)

src/                                   # NEW — C++ implementation files (.cpp)
  (empty for now — most Phase 0 code is header-only)
  (all new .cpp files go here as phases progress)

test/unit/                             # C++ unit tests (existing + new)
  util/
    array1d.cpp                        # NEW
    array2d.cpp                        # EXISTING
    array3d.cpp                        # EXISTING
  grid.cpp                             # EXISTING
  profile.cpp                          # EXISTING
  linear_algebra/
    test_tridiagonal_solver.cpp        # EXISTING
    test_error_function.cpp            # EXISTING

test/regression/                       # Regression tests (existing, not run until Phase 1+)
  solvers/
    delta_eddington.cpp                # EXISTING
    delta_eddington.hpp                # EXISTING

benchmark/
  benchmark_tridiagonal_solver.cpp     # EXISTING

cmake/                                 # CMake modules (existing + modified)
  dependencies.cmake                   # MODIFIED — drop Fortran deps, keep GTest/Benchmark/NetCDF-C
  test_util.cmake                      # MODIFIED — C++ tests only
  CodeCoverage.cmake                   # EXISTING
  (others existing)

docker/                                # Dockerfiles (existing + modified for C++-only builds)
.github/workflows/                     # CI workflows (new/modified)
.github/prompts/                       # Original planning docs (UNCHANGED)

plan/                                  # Revised implementation plans (this file)
data/                                  # NetCDF reference data (UNCHANGED)
examples/                              # JSON configs (UNCHANGED — useful for regression later)
```

---

## Step 1: CMake Build System (C++ Only)

### 1.1 Modify top-level `CMakeLists.txt`

The current `CMakeLists.txt` enables languages `Fortran, CXX, C` and builds both Fortran and C++ targets. We need to:

- **Remove `Fortran` from `project()` languages.** Keep `CXX` and `C` (C needed for NetCDF-C).
- **Add optional `CUDA` and `HIP` language support** (guarded by options, not enabled by default).
- **Require CMake 3.21+** (current minimum is already 3.21 per MUSICA).
- **Set C++ standard to C++20** (`CMAKE_CXX_STANDARD 20`).
- **Keep project name `tuv-x`** with alias `musica::tuvx`.
- **Remove all Fortran-specific compiler flag logic** (NAG, Intel Fortran, gfortran flags).
- **Remove the Fortran library target** (`musica::tuvx` currently links Fortran objects). Replace with a C++ header-only or static library target.
- **Disable/remove the current `tuv-x` CLI target in Phase 0.** The existing executable is a Fortran/YAML entry point and should not remain in the C++-only scaffold.
- **Rename `src/` to `fortran/`** — preserves all Fortran source as inert reference material. Not added to the build.
- **Reclaim `src/`** as a clean directory for C++ implementation files (initially empty — most Phase 0 code is header-only).

### 1.2 Modify `cmake/dependencies.cmake`

Current dependencies and their disposition:

| Dependency | Current | Phase 0 |
|-----------|---------|---------|
| netcdf-fortran | Required | **Remove** — no Fortran |
| netcdf (C) | Required | **Keep as optional** — needed for data readers (Phase 5), not yet used |
| yaml-cpp | FetchContent | **Remove** — config files are MUSICA's responsibility, not the solver library |
| GoogleTest | FetchContent | **Keep** — primary test framework |
| Google Benchmark | FetchContent | **Keep** — performance benchmarking |
| LAPACK/LAPACKE | Optional (find_package) | **Keep as optional** — used in benchmark comparisons |
| Doxygen/Sphinx | Optional | **Keep** — documentation generation |
| OpenMP | Optional | **Remove for now** — parallelization is via ArrayPolicy, not OpenMP directives |

### 1.3 Modify `cmake/test_util.cmake`

- **Remove `create_standard_test()`** (Fortran test helper).
- **Keep `create_standard_cxx_test()`** — this is the C++ GTest helper.
- **Keep Valgrind/memcheck integration** (`TUVX_ENABLE_MEMCHECK`).

### 1.4 CMake options (revised set)

```cmake
option(TUVX_ENABLE_CUDA    "Enable CUDA support"         OFF)
option(TUVX_ENABLE_HIP     "Enable HIP support"          OFF)
option(TUVX_ENABLE_TESTS   "Build unit tests"            ON)
option(TUVX_ENABLE_BENCHMARK "Build benchmarks"          OFF)
option(TUVX_ENABLE_COVERAGE "Enable code coverage"       OFF)
option(TUVX_ENABLE_MEMCHECK "Enable Valgrind memcheck"   OFF)
option(TUVX_ENABLE_LAPACK  "Enable LAPACK for benchmarks" OFF)
option(TUVX_ENABLE_NETCDF  "Enable NetCDF-C support"     OFF)
option(TUVX_BUILD_DOCS     "Build documentation"         OFF)
```

### 1.5 Remove YAML/CLI Surfaces from the Phase 0 Build Graph

Removing `yaml-cpp` is not just a dependency cleanup. It changes what can remain in-scope for Phase 0:

- **Remove or disable `TUVX_BUILD_CLI`** for this phase. The current `tuv-x` executable depends on Fortran sources and the YAML bridge.
- **Remove example-based tests that execute `tuv-x`** (`json`/`yaml` parity checks included). Those are validating the legacy front end, not the new C++ library scaffold.
- **Remove packaging/install logic tied to the CLI, yaml-cpp, and Fortran module files.** The Phase 0 artifact should be the C++ library plus public headers only.
- **Update workflows and Dockerfiles accordingly** so CI is not still installing or expecting `netcdf-fortran`, `yaml-cpp`, or the legacy CLI entry point.

---

## Step 2: Verify and Reorganize Existing C++ Code

The existing C++ headers are already in the right locations. This step verifies they compile and their tests pass under the new C++-only build system.

### 2.1 Audit existing headers

Each file needs to be checked for:
- Fortran interop code that should be removed (e.g., C bindings for Fortran callers)
- Dependencies on yaml-cpp or other removed dependencies
- Any `#include` paths that need updating

**Files to audit:**

| File | Key Types | Known Issues |
|------|-----------|-------------|
| `include/tuvx/util/array2d.hpp` | `Array2D<T>` | Simple `std::vector`-backed, row-major. No policy template yet — uses raw `T` not `ArrayPolicy`. Needs policy refactor (Step 3). |
| `include/tuvx/util/array3d.hpp` | `Array3D<T>` | Same pattern as Array2D. No policy template yet. |
| `include/tuvx/grid.hpp` | `Grid<ArrayPolicy>` | Already templated on `ArrayPolicy`, defaults to `Array2D<double>`. |
| `include/tuvx/profile.hpp` | `Profile<ArrayPolicy>` | Already templated on `ArrayPolicy`, defaults to `Array2D<double>`. |
| `include/tuvx/radiative_transfer/constituent.hpp` | `ConstituentState<ArrayPolicy>` (renamed from `RadiatorState`) | `Accumulate()` is a placeholder — returns first state. |
| `include/tuvx/radiative_transfer/radiation_field.hpp` | `RadiationFieldComponents<ArrayPolicy>`, `RadiationField<ComponentPolicy>` | Fully implemented containers. |
| `include/tuvx/linear_algebra/linear_algebra.hpp` | `TridiagonalMatrix<T>`, `Solve()` | Fully implemented Thomas algorithm. Uses `std::vector<T>` internally. |
| `include/tuvx/radiative_transfer/solvers/delta_eddington.hpp` | `DeltaEddington` | Placeholder — fills output with incrementing values. |
| `include/tuvx/cross_section.hpp` | `TuvxCrossSection` | Empty stub — **remove**, this will be replaced by the transform system (Phase 2). |
| `include/tuvx/cross_section_params.hpp` | (empty) | Empty stub — **remove**. |
| `include/tuvx/util/config_yaml.h` | C YAML API | **Remove** — yaml-cpp is a MUSICA concern, not the solver library. |

### 2.2 Verify existing tests compile and pass

These tests must compile and pass under the new build before proceeding:

| Test File | Tests | Dependencies |
|-----------|-------|-------------|
| `test/unit/util/array2d.cpp` | Array2D construction, sizing, access | GTest |
| `test/unit/util/array3d.cpp` | Array3D construction, sizing, access | GTest |
| `test/unit/grid.cpp` | Grid construction, column/section counts | GTest |
| `test/unit/profile.cpp` | Profile construction, units | GTest |
| `test/unit/linear_algebra/test_tridiagonal_solver.cpp` | Thomas algorithm vs LAPACKE | GTest, optional LAPACKE |
| `test/unit/linear_algebra/test_error_function.cpp` | Error computation | GTest |

The Fortran-based unit tests (`grid_warehouse`, `heating_rates`, `la_sr_bands`, `spherical_geometry`, etc.) are **not compiled or run** — they depend on the Fortran build.

The mixed-language regression test under `test/regression/solvers/` is also **out of scope for Phase 0 as currently written**. It is built through the Fortran-oriented test helper and depends on a Fortran driver. Either:

- disable it for Phase 0 alongside the other Fortran-backed tests, or
- replace it with a pure-C++ regression harness before claiming Phase 0 validation coverage for the solver path.

### 2.3 Remove stubs and dead code

- Delete `include/tuvx/cross_section.hpp` (empty stub)
- Delete `include/tuvx/cross_section_params.hpp` (empty stub)
- Delete `include/tuvx/util/config_yaml.h` and any associated `.c`/`.cpp` implementation (yaml-cpp wrapper)
- Remove any CMake, test, packaging, workflow, and Docker references to these files and the YAML-backed CLI path

---

## Step 3: Create `Array1D<T>` and Establish Policy Pattern

### 3.1 Current state of Array2D/Array3D

The existing `Array2D<T>` and `Array3D<T>` are **not yet policy-templated** in the MICM sense. They are simple `std::vector`-backed containers with `operator()` access (not `operator[]` chaining). The plan calls for a bulk operations API (`ForEachRow`, `ColumnView`, `Function` factory) matching MICM conventions.

This is the most significant code-authoring work in Phase 0. Options:

**Option A: Full policy refactor now.** Rewrite Array2D/Array3D with the full bulk operations API and `ArrayPolicy` template. Create Array1D to match. This is what the plan describes.

**Option B: Incremental refactor.** Keep existing Array2D/Array3D working, create Array1D following the same simple pattern, and add the bulk operations API as a second pass. This lets us verify the build/CI pipeline before tackling the harder API work.

**Recommendation: Option B.** The build system, CI, and test infrastructure are higher priority than the policy refactor. We can verify everything compiles and passes first, then layer on the bulk operations API.

### 3.2 Create `Array1D<T>` (initial version)

New file: `include/tuvx/util/array1d.hpp`

Follow the same pattern as Array2D/Array3D:
- `std::vector<T>` backing storage
- `operator[]` element access (not `operator()` — start with the target API)
- `Size()` / `NumElements()` dimension query
- Constructor taking size, optional fill value
- `begin()`/`end()` for iteration
- `AsVector()` for interop

New test file: `test/unit/util/array1d.cpp`

### 3.3 Bulk Operations API (second pass within Phase 0)

After the build system is solid and all existing tests pass, add the MICM-consistent bulk operations API to all three array types:

**For Array1D:**
- `operator[]` returns element reference (already simple)
- `ForEachRow(lambda, views...)` — iterates element-wise
- `NumRows()` (= number of elements, since 1D)
- `NumColumns()` (= 1, semantic compatibility)

**For Array2D:**
- `operator[]` returns a proxy whose `operator[]` resolves the second dimension
- `ForEachRow(lambda, views...)` — iterates over rows
- `GetColumnView(col)` / `GetConstColumnView(col)`
- `GetRowVariable()` — temporary storage
- `NumRows()`, `NumColumns()`
- Static `Function(lambda, prototype)` factory

**For Array3D:**
- `operator[]` returns a proxy whose `operator[]` returns another proxy
- Same bulk ops pattern

If this API lands in Phase 0, it needs dedicated GTests. The current constructor/access tests only validate `operator()` and contiguous iteration; they do **not** validate proxy access, row iteration helpers, views, or factory helpers.

**Key decision: `operator()` vs `operator[]`.**
The existing code uses `operator()(i, j)` for Array2D. The plan calls for `operator[][j]` chaining. We need to support both during transition (the plan's target API uses `[]` chaining). Add `operator[]` returning proxy objects; keep `operator()` as deprecated until all call sites are migrated.

### 3.4 Define `HostArrayPolicy`

Do **not** assume a single `HostArrayPolicy` struct can be dropped into the current templates unchanged. The existing public types consume concrete container types by dimension:

- `Grid` and `Profile` currently expect a 2D container type
- `ConstituentState` (formerly `RadiatorState`) and `RadiationFieldComponents` currently expect a 3D container type

Phase 0 should therefore choose one of these approaches explicitly:

- **Defer policy naming** and keep concrete host container types (`Array1D`, `Array2D`, `Array3D`) as the template arguments for now, or
- **Introduce dimension-specific host aliases/concepts** that match the existing contracts, then migrate the public types once the array refactor is complete.

`DeviceArrayPolicy` (CUDA/HIP) remains **deferred**. Do not make Phase 0 acceptance depend on a single host-policy struct unless the surrounding templates are redesigned to consume it correctly.

---

## Step 4: CI/CD Setup

### 4.1 GitHub Actions Workflows

Create/modify workflow files. Since we're removing Fortran, the existing workflows need significant changes.

**Priority order** (get the most valuable gates first):

1. **`ubuntu.yml`** — GCC + Clang on Ubuntu 24.04. This is the primary development gate.
   - Matrix: GCC 12/13/14, Clang 17/18
   - Build types: Debug, Release
   - Run: cmake, build, ctest

2. **`mac.yml`** — Apple Clang + GCC on macOS-latest.
   - Matrix: Apple Clang (default), GCC 14
   - Build types: Debug, Release

3. **`clang-tidy.yml`** — Static analysis PR gate.
   - Runs on PRs only
   - Uses compilation database from CMake
   - Any warning fails the PR

4. **`clang-format.yml`** — Auto-formatting.
   - Runs on push to `cpp-rewrite` (not `main` — we're on a branch)
   - Opens PR if formatting changes detected

5. **`docker.yml`** — Coverage + Memcheck + Intel/NVHPC.
   - Coverage container: lcov, 95% threshold
   - Memcheck container: Valgrind
   - Intel: icx/icpx via oneAPI
   - NVHPC: nvc++

6. **`windows.yml`** — MSVC + MinGW on windows-latest.
   - Lower priority — scientific HPC code is rarely Windows-first
   - Add after Linux/macOS are solid

### 4.2 `.clang-tidy` update

The current `.clang-tidy` only checks `readability-identifier-naming`. Expand per the plan:

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

WarningsAsErrors: '*'

CheckOptions:
  - key: readability-identifier-naming.ClassCase
    value: CamelCase
  - key: readability-identifier-naming.FunctionCase
    value: CamelCase
  - key: readability-identifier-naming.VariableCase
    value: lower_case
  - key: readability-identifier-naming.ConstantCase
    value: UPPER_CASE
  - key: readability-identifier-naming.TemplateParameterCase
    value: CamelCase
  - key: readability-identifier-naming.MemberSuffix
    value: '_'
```

**Risk**: Expanding clang-tidy checks will flag existing code. We need to either fix all existing headers to pass, or apply the expanded checks only to new code paths initially.

**Resolved**: function naming stays `CamelCase` to match existing public API (`NumberOfColumns`, `Size1`, `AsVector`). The `.clang-tidy` config reflects this.

### 4.3 `.clang-format` review

Current config is Google-based with Allman braces, 125-column limit. This is reasonable. No changes needed unless the team has preferences.

### 4.4 Codecov configuration

Create `codecov.yml` at repo root:

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

### 4.5 Docker updates

Modify Dockerfiles to build C++ only:
- Remove `gcc-fortran` / `gfortran` packages
- Remove `netcdf-fortran` dependency
- Keep `netcdf` (C) as optional
- Remove yaml-cpp
- Keep GTest, Google Benchmark, lcov, Valgrind

---

## Step 5: Documentation Baseline

### 5.1 Doxygen comments

As headers are audited in Step 2, add `///` Doxygen comments to every public class and function. Keep them succinct per the documentation policy:
- What it does
- What parameters mean
- Preconditions
- No `@author`, `@date`, `@file`

### 5.2 README update

Update the top-level README to reflect the C++ rewrite:
- Build instructions for the C++ library
- Brief API overview pointing at the target README example from the master plan
- Link to Doxygen docs (once generated)
- Note that Fortran source is preserved as reference material

---

## Execution Order

The steps above have dependencies. Here's the recommended execution sequence:

```
Step 0.0  Rename src/ → fortran/, create empty src/ for C++
    │
    ▼
Step 1.1  Modify CMakeLists.txt (remove Fortran, C++ only)
Step 1.2  Modify dependencies.cmake
Step 1.3  Modify test_util.cmake
Step 1.4  Update CMake options
Step 1.5  Remove YAML/CLI/example/packaging surfaces from the Phase 0 build graph
    │
    ▼
Step 2.1  Audit existing headers (remove stubs, fix includes)
Step 2.3  Remove dead code (cross_section stubs, yaml wrapper)
    │
    ▼
Step 2.2  Verify existing tests compile and pass
    │
    ▼
Step 3.2  Create Array1D<T> (initial, simple version)
          Write Array1D tests
    │
    ▼
Step 4.1  Set up CI workflows (ubuntu.yml first, then expand)
Step 4.2  Expand .clang-tidy
Step 4.4  Add codecov.yml
    │
    ▼
Step 3.3  Add bulk operations API (ForEachRow, ColumnView, Function)
Step 3.4  Define HostArrayPolicy
    │
    ▼
Step 4.5  Update Dockerfiles
Step 5.1  Doxygen comments on all public headers
Step 5.2  Update README
    │
    ▼
[Phase 0 complete — merge phase0-scaffolding → cpp-rewrite]
```

---

## Validation Strategy

### Principles

- **Validate at every phase checkpoint** with whatever is testable at that point.
- **Reference data generated from Fortran** — build and run the Fortran code, capture outputs as CSV.
- **Reference CSV committed to git** in `test/reference/` so it's tracked and reproducible.
- **No bit-for-bit requirement** — numerical operations may be reordered. Tolerances are empirical, adjusted per test based on what makes sense.
- **Intermediate outputs, not just final rates** — capture grid values, optical properties, solver intermediates, transform weights, etc. so divergences can be localized.
- **Defer optimization** — get correctness first, tune later.

### Reference Data Layout

```
test/reference/
  phase1/                              # Delta Eddington solver outputs
  phase2/                              # Cross-section/quantum yield weight matrices
  phase3/                              # Photolysis rates, dose rates, heating rates
  phase4/                              # Constituent accumulation, LA/SR bands
  phase5/                              # Data reader outputs (interpolated data)
  ...
```

### Phase 0 Validation

Phase 0 is scaffolding — no new numerical code. Validation is:

- **all retained C++ GTest unit tests pass under the new C++ build system**
- **new tests cover any Phase 0 array API additions beyond the current constructor/access surface**
- **mixed-language Fortran/C++ regression tests are either explicitly disabled for Phase 0 or ported to pure C++ before being counted as validation**

The tridiagonal solver already has tests comparing against LAPACKE; those continue to pass.

Numerical reference data generation begins at Phase 1.

### Phase 1+ Validation

Each subsequent phase generates its own reference data before C++ implementation begins:

1. Build Fortran code on this machine (homebrew clang/flang)
2. Run with controlled inputs (documented in the reference directory)
3. Capture intermediate and final outputs as CSV
4. C++ tests (GTest) load the CSV and compare within per-test tolerances
5. If tolerance needs adjusting, document why in the test

---

## Acceptance Criteria

Phase 0 is done when:

- [ ] `phase0-scaffolding` branch builds C++ only (no Fortran compilation)
- [ ] Fortran source preserved in `fortran/` (renamed from `src/`, not built)
- [ ] `Array1D<T>` exists; any Phase 0 bulk-operations API added to `Array1D<T>`, `Array2D<T>`, or `Array3D<T>` is covered by dedicated unit tests
- [ ] `Grid`, `Profile`, `ConstituentState` (renamed from `RadiatorState`), and `RadiationField` continue to compile against host-side array container types; if a host policy abstraction is introduced, it matches the dimensional contracts of those types
- [ ] All existing C++ unit tests pass
- [ ] New `Array1D` unit tests pass
- [ ] Benchmark compiles (with LAPACK optional)
- [ ] CLI-, YAML-, and example-runner code paths are removed or explicitly disabled for Phase 0
- [ ] CI runs on at least Ubuntu (GCC + Clang) and macOS
- [ ] clang-tidy passes on all headers in `include/tuvx/`
- [ ] clang-format config enforced
- [ ] Codecov configured with 95% target
- [ ] Public headers have Doxygen comments
- [ ] README reflects C++ rewrite status
- [ ] All retained GTest unit tests pass under the new build (no new numerical validation needed)
- [ ] Mixed-language regression tests are either disabled for Phase 0 or replaced with pure-C++ equivalents

---

## Codex Review

This plan was reviewed against the current repository state. The following adjustments are now reflected in the Phase 0 scope:

- **Mixed-language regression tests are not assumed to survive the C++-only cutover.** The existing Delta-Eddington regression target is wired through the Fortran test helper, so it must be disabled or ported before it counts as Phase 0 coverage.
- **Removing `yaml-cpp` also removes the current CLI path from scope.** The existing executable, example-runner tests, packaging rules, workflows, and Dockerfiles all still assume the YAML-backed Fortran front end. Phase 0 must retire or disable those surfaces, not just delete the dependency.
- **`HostArrayPolicy` cannot be treated as a drop-in single type under the current public templates.** Existing types expect concrete array containers with different dimensionality; any host-policy abstraction must respect that or be deferred.
- **The array API refactor needs direct tests.** Existing array tests only cover constructors, iteration, and `operator()`. New proxy access, row helpers, views, and factories must be validated explicitly if they are introduced in Phase 0.
- **Any clang-tidy naming change must be explicit about public API churn.** The repository currently exports `CamelCase` method names; a new naming policy should either preserve those names for now or document the migration as a deliberate compatibility break.
