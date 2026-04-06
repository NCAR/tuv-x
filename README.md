<h1 align="center">
<img src="docs/source/_static/logo.svg" width="300">
</h1><br>

Tropospheric ultraviolet-extended (TUV-x): A photolysis rate calculator

[![License](https://img.shields.io/github/license/NCAR/tuv-x.svg)](https://github.com/NCAR/tuv-x/blob/main/LICENSE)
[![Ubuntu](https://github.com/NCAR/tuv-x/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/ubuntu.yml)
[![Mac](https://github.com/NCAR/tuv-x/actions/workflows/mac.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/mac.yml)
[![Windows](https://github.com/NCAR/tuv-x/actions/workflows/windows.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/windows.yml)
[![Docker](https://github.com/NCAR/tuv-x/actions/workflows/docker.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/docker.yml)
[![codecov](https://codecov.io/gh/NCAR/tuv-x/branch/main/graph/badge.svg?token=H46AAEAQF9)](https://codecov.io/gh/NCAR/tuv-x)
[![DOI](https://zenodo.org/badge/396946468.svg)](https://zenodo.org/badge/latestdoi/396946468)

Copyright (C) 2020-2026 University Corporation for Atmospheric Research

## C++ Rewrite (In Progress)

TUV-x is being rewritten from Fortran to C++ as a header-only library.
The original Fortran source is preserved in the `fortran/` directory for
reference during the translation.

### Building

Requirements: CMake 3.21+, a C++20 compiler (GCC 12+, Clang 17+, Apple Clang, MSVC 2022+).

```bash
git clone https://github.com/NCAR/tuv-x.git
cd tuv-x
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `TUVX_ENABLE_TESTS` | ON | Build unit tests (requires GTest, fetched automatically) |
| `TUVX_ENABLE_BENCHMARK` | OFF | Build benchmarks (requires LAPACK) |
| `TUVX_ENABLE_LAPACK` | OFF | Enable LAPACK for benchmarks |
| `TUVX_ENABLE_COVERAGE` | OFF | Enable code coverage (lcov) |
| `TUVX_ENABLE_MEMCHECK` | OFF | Enable Valgrind memcheck |
| `TUVX_ENABLE_NETCDF` | OFF | Enable NetCDF-C support |
| `TUVX_ENABLE_CUDA` | OFF | Enable CUDA support |
| `TUVX_ENABLE_HIP` | OFF | Enable HIP support |
| `TUVX_BUILD_DOCS` | OFF | Build documentation |

### Using TUV-x as a Library

TUV-x provides a header-only interface library via the `musica::tuvx` CMake target:

```cmake
find_package(tuvx REQUIRED)
target_link_libraries(my_app PRIVATE musica::tuvx)
```

```cpp
#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>
```

### Repository Layout

| Directory | Contents |
|-----------|----------|
| `include/tuvx/` | Public C++ headers |
| `src/` | C++ implementation files (currently header-only) |
| `test/` | C++ unit and regression tests |
| `benchmark/` | Google Benchmark files |
| `fortran/` | Original Fortran source (reference only, not built) |
| `data/` | NetCDF reference data |
| `plan/` | Implementation plans for the C++ rewrite |

# Citation

> Madronich, Sasha, and Siri Flocke (1999), The role of solar radiation in atmospheric chemistry, in Handbook of Environmental Chemistry, edited by P. Boule, pp. 1-26, Springer-Verlag, Heidelberg.

```
@incollection{madronich_role_1999,
	address = {Berlin, Heidelberg},
	series = {The {Handbook} of {Environmental} {Chemistry}},
	title = {The {Role} of {Solar} {Radiation} in {Atmospheric} {Chemistry}},
	isbn = {978-3-540-69044-3},
	url = {https://doi.org/10.1007/978-3-540-69044-3_1},
	language = {en},
	booktitle = {Environmental {Photochemistry}},
	publisher = {Springer},
	author = {Madronich, Sasha and Flocke, Siri},
	editor = {Boule, Pierre},
	year = {1999},
	doi = {10.1007/978-3-540-69044-3_1},
	keywords = {Earth-Sun geometry., photolysis rate coefficients, radiative transfer, solar radiation, spectral actinic flux},
	pages = {1--26},
}
```

# Community and contributions

We welcome contributions and feedback from anyone.

- [Contact](https://github.com/NCAR/tuv-x/discussions) — start a conversation on our [discussion board](https://github.com/NCAR/tuv-x/discussions) or email musica-info@ucar.edu.
- [Collaboration](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf) — read the [MUSICA software development plan](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf).
- [Contributor's guide](https://ncar.github.io/tuv-x/contributing/contributors_guide.html) — please read before submitting a PR.

# Documentation

See the [TUV-x documentation](https://ncar.github.io/tuv-x/) for detailed instructions.

# License

- [Apache 2.0](/LICENSE)
- Copyright (C) 2020-2026 University Corporation for Atmospheric Research
