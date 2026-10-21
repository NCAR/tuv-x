# TurboTUV

TurboTUV is an LLM experiment by the NCAR Agentic Coding Task Force, sponsored
by the [Community Software Facility](https://ncarprojects.ucar.edu/csf) and the
[NCAR Office of Strategic Projects AI Initiative](https://ncarprojects.ucar.edu/ai).
In it, LLM coding agents translate NCAR's TUV-x photolysis model to C, validate
the translation against the original Fortran, and make it fast. The scientific
software engineer (SSE) supervising the project is David Fillmore.

C11 translation of NCAR TUV-x, with compiled-in scientific data and a direct
Fortran comparison suite. The library covers
interpolation, tridiagonal systems, spherical geometry, solar position,
radiator mixing, both Delta-Eddington and discrete-ordinate radiation solvers, and 26 direct
cross-section calculators, 21 quantum-yield calculators, and 13 spectral
weighting functions, all four temperature parameterizations, and photolysis,
heating, and dose integration, including the O2 band corrections. Grid and
profile constructors use embedded atmospheric and solar spectra; host updates
share the same structures. JSON and YAML parsing use an immutable C interface
in [include/turbotuv_config.h](include/turbotuv_config.h), backed by bundled
[LibYAML](https://github.com/yaml/libyaml) sources with private symbol prefixes.
Prepared spectral models in [include/turbotuv_model.h](include/turbotuv_model.h)
construct cross sections, quantum yields, and weights from the original
configuration formats and compiled datasets. They own their tables and support
repeated evaluation with different atmospheric profiles. The same header
provides grid and profile constructors for the original configuration formats.
The complete model API assembles those components, applies the O2 corrections,
and returns all six radiation fields plus photolysis, heating, and dose rates.
Host-supplied radiation uses the same rate pipeline.
Host applications can also supply and update column grids, profiles, and optical
properties. Each model owns its mutable state; create one model per worker.

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
build/release/turbotuv examples/tuv_5_4.json --output results.json
```

The production library requires only a C compiler and the system math library.
The standalone runner accepts the original JSON/YAML examples and writes all
radiation fields and rates. Optional `TURBOTUV_NETCDF_OUTPUT=ON` adds NetCDF
output to the runner; the library remains independent of NetCDF.
All 190 TUV-x data files are embedded in the library. NetCDF values are stored
as exact hexadecimal doubles; text tables retain their original bytes. The
public API is in [include/turbotuv.h](include/turbotuv.h). Units, dimensions,
array ordering, ownership, and errors are documented there.

Reference testing uses the sibling `../TUV-x` checkout at revision
`ebe5006` (0.17.0). `tools/reference.py gcc` and `tools/reference.py llvm` build
and run its existing tests using the Homebrew compiler sets. `--prefix` and
`--netcdf` select compiler-matched NetCDF Fortran and NetCDF C installations.
The reference Python analyses need numpy, scipy, and netCDF4.

To enable direct numerical comparisons, configure with
`TURBOTUV_REFERENCE_BUILD`, `TURBOTUV_REFERENCE_PREFIX`,
`TURBOTUV_NETCDF_PREFIX`, and the matching `CMAKE_Fortran_COMPILER`.
The test-only bridge links the original Fortran library; production code does
not call Fortran. Fast math and floating-point contraction are disabled to
preserve the reference's arithmetic.

`tools/embed_data.py --check` checks regeneration, and
`tests/verify_data.py build/release/tests/dump_data` audits every compiled value
and text byte against the original scientific files. Only these development
tools need Python/netCDF4 and the original input files. Source hashes and
upstream revision are recorded in [data/manifest.json](data/manifest.json).

Models run in prepared execution by default: spectral values that cannot
change between columns are computed once when the model is built, equal
spectra are shared, rates are summed over precomputed products, and
Delta-Eddington transfer solves blocks of wavelengths in lock step. Every result
is bit-identical to the direct translation, which remains available as generic
execution and as the `baseline` branch. A column runs 4.3–7.6 times faster than
the baseline and 5.4–6.8 times faster than the original Fortran TUV-x, whose
results it reproduces bit for bit with the GNU compilers; see
[docs/PERFORMANCE.md](docs/PERFORMANCE.md).

Configurations are validated as TUV-x validates them: misspelled or
unsupported keys are errors. Numbers are read the same way in every locale.
`turbotuv_bench` (built by default, not installed) times a model evaluating
many columns whose temperature, density, ozone and solar zenith angle change
through the host API, and prints a hash of every output value.

See [docs/USAGE.md](docs/USAGE.md) for the runner, output format, installation,
and host-model API, and [docs/PORTING.md](docs/PORTING.md) for validation details.
The upstream Apache-2.0 license and UCAR attribution apply to translated code
and scientific data.

Validated with Homebrew GCC/GFortran 16.2.0 and Clang/Flang 23.1.3: all 39 C
tests and 88 original reference tests pass on each toolchain. The 71.5 million
numerical comparisons per toolchain stay within double-precision roundoff;
the complete-model comparisons are bit-identical with GCC and within one
relative machine epsilon with Clang. Sanitizer and concurrent-model checks
are recorded in the validation document.
