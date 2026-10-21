# Translation and validation record

Reference: `../TUV-x`, revision `ebe5006` (version 0.17.0).

The acceptance target is the original unit/regression/example coverage plus
new boundary, randomized, and analytic tests, compared directly with Fortran
at double-precision roundoff. Passing the upstream percent-level regression
tolerances alone does not establish acceptance.

## Completed implementation tasks

* C11 build, public numerical API, explicit error returns and caller ownership.
* All 190 upstream scientific/configuration data files embedded, including
  metadata and provenance. Compiled-data audit: 4,547,069 bytes compared exactly
  on both compiler sets; no external runtime data files are required.
* Linear, conserving, fractional-source and fractional-target interpolation,
  including the reference overhang convention; tridiagonal solver; spherical
  ray geometry and air columns; Gregorian date and solar position routines.
* 863,465 direct Fortran/analytic numerical comparisons per compiler set;
  largest error is 2 epsilon-scaled units, with a 32-unit acceptance ceiling.
  Solar-angle error is scaled by 90 degrees to account for cancellation near
  the horizon; other comparisons use relative error. Exact zeros must match.

Compiler versions tested: Homebrew GCC/gfortran 16.2.0 and Clang/Flang 23.1.3.
Both use `-ffp-contract=off`; fast math is not enabled.

* Delta-Eddington solver and optical-property mixing: 275 atmospheres and
  377,125 comparisons per compiler set, including 1/2/5/16/81 layers, 11 solar
  angles from overhead to night, vacuum, thin/thick layers, pure scattering,
  Rayleigh contributions, and surface albedos from zero to one. All six fields
  match the original Fortran within 11.1 epsilon-scaled units (64-unit ceiling,
  normalized incident-beam scale). Independent surface-reflection and direct
  beam boundary checks pass. The C-only invariants also pass Clang address and
  undefined-behavior sanitizers.

* All 26 direct cross-section kernels use a common grid and parameter-table
  representation. Tests exercise 25 temperatures and 1,464 wavelengths,
  including neighboring representable values at wavelength cutoffs. The
  951,643 comparisons per compiler set have a largest relative error of
  0.624 epsilon-scaled units (32-unit ceiling). Embedded column regridding
  also checks conservation against an independently calculated trapezoid
  integral. C-only spectral tests pass both sanitizers.

* All 21 quantum-yield calculators and 13 spectral weighting functions, with
  temperature and density sweeps, adjacent floating-point values at cutoffs,
  all acetone branches, and the original Taylor-series constructor. The
  1,295,363 comparisons per compiler set have a largest relative error of
  13.5 epsilon-scaled units (32-unit ceiling). Independent complementary-branch
  and normalization tests pass, including under both sanitizers.

* All four temperature parameterizations (exponential polynomial, Taylor,
  Burkholder, Harwood), with a shared range/slice representation and the full
  temperature-dependent cross-section regridding path. All 417,727 comparisons
  are bit-for-bit identical on both compiler sets. Tests include overlapping
  inclusive ranges, fixed temperatures, both logarithm bases, inverted
  offsets, interior slices, and zero padding. Sanitizer checks pass.

* Actinic-flux and irradiance conversion, photolysis, heating, and dose
  integration. The 53,821 comparisons against original rate constructors and
  calculators are bit-for-bit identical on both compilers. They cover up to
  1,023 bands and 81 levels, three wavelength/temperature-dependent reaction
  spectra, variable thresholds/scaling, negative-total clipping, and darkness.
  Integrations traverse contiguous spectra and preserve each level's original
  summation order. Independent single-band and sanitizer checks pass.

* O2 Lyman-alpha and Schumann-Runge cross sections and optical depths, using
  compiled Chebyshev coefficients. All 318,240 comparisons are bit-for-bit
  identical on both compilers, including the 680 coefficients parsed by the
  original constructor, column cutoffs, nearly equal columns, extreme slant
  columns, temperature variation, and untouched neighboring bands. Sanitizers
  pass. The C code defines the all-thin-column limit using the reference's
  tabulated values; the original indexes outside its array in that case.

* Dense and packed-band LU factorization, reciprocal-condition estimation,
  and normal/transposed solves share one matrix representation and elimination
  path. All 753,290 factor/pivot/condition/solution comparisons are bit-for-bit
  identical on both compilers. Tests cover varied bandwidths and pivoting,
  sizes 1 through 81, extreme scales, independent normwise backward-error
  checks, and singular/invalid inputs. Sanitizers pass.

* Discrete-ordinate real eigenproblem, Gauss quadrature, and Legendre helpers.
  All 100,460 comparisons are bit-for-bit identical on both compilers, with
  independent eigen-equation residuals and polynomial quadrature integrals.
  Eigen tests include symmetric, similarity-scaled asymmetric, triangular,
  and isolated-block matrices through order 32; quadrature runs through order
  256. Private Fortran helpers are extracted verbatim into a test-only module
  with source provenance, leaving the original library unchanged. Sanitizers
  pass; a zero 2x2 matrix receives defined identity eigenvectors in C.

* Complete pseudo-spherical discrete-ordinate solver, specialized to TUV-x's
  fixed flux-only, delta-M, Lambertian-surface configuration. All 16 even
  stream counts from 2 through 32 pass 4,400 atmosphere cases and 16,368,400
  comparisons per compiler, within 46 epsilon-scaled units (64-unit ceiling).
  Cases cover 1/2/5/16/81 layers, 11 solar angles, thin/thick/vacuum optics,
  nearly conservative scattering, Rayleigh mixing, surface reflection and
  absorption-depth truncation. Both sanitizers pass. Layer scratch storage is
  reused across wavelengths. The band solver also covers bandwidths larger
  than matrix order, as needed for single-layer boundary systems.

* Shared owned grid/profile structures, equal-interval and supplied grids,
  compiled text-table parsing, air/O2 logarithmic profiles, O3 extension and
  column rescaling, solar position/distance profiles, solar-spectrum merging,
  and host updates. All 697 cases and 541,850 comparisons pass on both compiler
  sets within 5.85 epsilon-scaled units (32-unit ceiling), including zero-cell
  profiles, 240-layer profiles, upper columns, host midpoint overrides, and
  adjacent floating-point values at all solar splice boundaries. Sanitizers
  pass. The original solar conversion's default-REAL literal and opposite
  time-zone signs for zenith/distance are preserved explicitly.

* Gas and aerosol optical-property construction, plus immutable embedded
  radiator tables. All 432 cases and 2,316,702 comparisons pass on both compiler
  sets within 0.556 relative epsilon-scaled units (32-unit ceiling); sanitizers
  pass. Cases cover 1–120 layers, aerosol source grids of 2–121 points, column
  normalization, folded upper material, empty layers and air/absorber behavior.
  All 36 original NetCDF-radiator assertions are reproduced in C. The original
  single-sample aerosol constructor uses unallocated arrays; C defines it as
  uniform layer depth.
* The 110 modern upstream test input datasets are also compiled, in a separate
  test-only library. Both compiled-data audits compare 1,884,240 bytes exactly.
  Legacy TUV data and the generic NetCDF I/O fixture are excluded from this
  test catalogue; the original Fortran suites still run those tests normally.

* JSON/YAML configuration parsing uses a bundled, namespaced C LibYAML parser
  behind an immutable node API. All 281 documents and 26,385 checks pass on
  both compilers, including exact numeric/array/string/boolean comparisons
  through the original Fortran configuration interface, Unicode, aliases,
  null/string distinctions, file input, malformed/truncated input and bounded
  nesting. Both sanitizers pass. Upstream source hashes and MIT licensing are
  recorded under `third_party/libyaml`; no runtime dependency is introduced.

* Complete spectral constructors are connected to the C configuration API.
  The prepared objects own their fixed wavelength grid, parameter tables and
  options and never open scientific files. An optional resolver selects other
  compiled datasets, including the separate test catalogue. Shared validators
  keep direct-kernel and prepared-model requirements consistent.
  All 970 cases and 7,189,693 comparisons pass with GCC/GFortran (bit-identical)
  and Clang/Flang (worst 0.981 relative epsilon-scaled units, 32-unit ceiling).
  Address and undefined-behavior sanitizers pass. Cases include every spectral
  constructor from both complete examples, all modern spectral unit-test
  inputs, the production reaction/weight catalogues, midpoint evaluation,
  merged files, all interpolators, cutoffs, overrides, point data, extrapolation,
  and all four temperature parameterizations. The case generator records
  original source hashes and normalizes only the stale type labels that the
  original unit tests bypass by directly calling a specific constructor.
  The ozone constructor's replacement of the warm Hartley/Huggins column and
  quantum-yield upper-boundary padding convention are preserved explicitly.
  The CL2 reference allocates an extra undefined row at midpoints; comparisons
  include its defined rows. HNO3 and H2O2 ignore the reference midpoint flag;
  their edge evaluations are compared. C evaluates precisely the caller's
  supplied temperature array without that ambiguity.

* Configuration-driven profile construction now resolves named grids and
  builds every profile kind from compiled inputs. The catalogue test covers
  237 entries: 226 valid cases and 11 invalid dormant upstream definitions,
  plus nine explicit malformed-input checks. All 87,105 comparisons of
  defined Fortran fields are bit-identical on both compiler sets. Additional
  cases vary model tops, O2/O3 scale heights, dates and time zones. Sanitizers
  pass. The solar-flux profile defines its otherwise-unused fields as zero.
  Reference factory tests run in isolated build directories with links to
  the original data, so error/diagnostic files cannot alter the source tree.

* Complete model construction and repeated execution own their grids, profiles,
  spectra, optical properties and scratch arrays. O2 corrections, both solvers,
  host radiation and all three rate collections use the same pipeline. The
  whole-model suite covers 27 configurations and 380 runs, including both
  examples, all radiator regression configurations, the production photolysis
  and dose catalogues, heating, nonuniform columns, bright surfaces, vacuum,
  negative host radiation, omitted irradiance, and 17 malformed configurations.
  All 30,117,173 comparisons pass: bit-identical with GCC/GFortran and within
  0.999 relative epsilon-scaled units with Clang/Flang (64-unit ceiling).
  Radiation comparisons use incident-flux normalization; rates use their own
  magnitudes. Both sanitizers pass. Rate labels/order are also compared.
  A nonphysical NetCDF reader fixture is tested for rejection at the solve.
  The reference mixer dereferences its first radiator for an empty collection;
  empty C atmospheres are compared with one zero-depth Fortran absorber, which
  has exactly the same optical-floor arithmetic. Source hashes and this test
  adjustment are recorded by the generator. The reference's direct scaling by
  the Earth-Sun distance argument is preserved, without reinterpretation.

* Host model construction accepts copied named grids, profiles and radiators;
  matching names cannot appear twice. Repeated updates preserve dimensions and
  view addresses, support optional midpoint/column/exosphere values, and replace
  radiator optics explicitly. The wavelength grid is fixed after construction
  because spectral tables are prepared on it. Four configurations, 48 changing
  columns and 528 daylight/twilight/night runs pass 8,769,552 comparisons against
  the original core and host updater interfaces: bit-identical with GCC/GFortran,
  within 0.866 relative epsilon-scaled units with Clang/Flang. Both sanitizers pass.
  Tests release all input storage after construction, retain borrowed views,
  exercise omitted fields and reject duplicate names, wrong shapes and NaNs.

* The standalone runner accepts the unchanged JSON/YAML examples and writes
  round-trippable JSON, or optional NetCDF-4 output through the NetCDF C library.
  Its core dependency remains C11/libm. The CLI tests run in fresh temporary
  directories without scientific data files and compare both formats, angle
  overrides, darkness, malformed arguments and atomic output replacement.
  The 980,172 numerical output comparisons per compiler include the original
  executable's example NetCDF files and are all bit-identical. JSON and YAML
  results are identical. Address/undefined-behavior checks cover both writers.
  CMake installation exports a usable package target and includes the four
  original examples plus a minimal C consumer.
  On this macOS installation, enabling LeakSanitizer for the optional
  NetCDF/HDF5 process stalls after the output file has been closed and renamed.
  That CLI check uses `ASAN_OPTIONS=detect_leaks=0`; address and undefined-behavior
  instrumentation remain enabled. Core and JSON-only checks pass with leak
  detection enabled.

* Four concurrent models share immutable configuration and compiled data while
  owning their calculations. Both production configurations reproduce 4,410,208
  serial values bit-for-bit with both compilers and under sanitizers. Clang's
  static analyzer reports no findings in model construction, host updates,
  execution, or either output writer. Separate installed-package consumers run
  with GCC/static and Clang/shared builds; the shared core links only libSystem.

## Reference portability

`tools/reference.py` creates an isolated source copy under `build/` and applies
test-harness portability corrections: legacy default-real promotion for Flang,
and portable temporary-file replacement for GNU-only `sed -i` commands. The
original checkout and numerical equations are unchanged.

The first unmodified GNU run built and passed 76/88 tests; remaining failures
were the BSD sed incompatibility and missing SciPy. Flang initially failed to
build legacy tests because real-kind promotion was omitted. After correcting
those issues, the new numerical comparisons pass on both compilers.
Both corrected reference builds now pass all 88 upstream tests. Two additional
test-harness corrections are needed: explicitly rewind the reopened legacy
CH2O input, and route the C++ placeholder test's large structs through pointer
arguments. Instrumentation isolated the Flang crash to the by-value C call on
Darwin/arm64. The pointer shim preserves all original inputs, outputs, and
assertions. No numerical equations or tolerances are changed.

Undefined reference behavior is not copied: the C air-column routine defines
the top vertical-column entry, which Fortran leaves unwritten; comparisons
exclude that entry. C validates grids and prevents out-of-bounds access. The
reference conserving interpolator assumes padded source grids, so randomized
reference comparisons use interior target bins.

## Review corrections — 2026-10-19

An independent review compared every translated routine with the Fortran again,
fuzzed the parsers and solvers, injected allocation failures and exercised
inputs the reference suites do not reach. It found no memory errors. The
following defects are corrected, each with a regression test:

* A missing aerosol "550 nm optical depth" now defaults to 0, as in TUV-x,
  which normalizes the aerosol column to zero. C used -1, which skipped the
  normalization and changed rates by up to 2% (verified against the original
  executable).
* O3 sample profiles append 1-km points by accumulating the coordinate until
  one lies above the grid top, as TUV-x does. C counted floor(top - last) + 1,
  one too few when the difference rounds just below an integer, which zeroed
  the top edge. Shipped data were unaffected; the new case matches the
  original constructor's values.
* Configuration mappings are validated as TUV-x's `config%validate` does:
  required keys must be present and other keys must be accepted by that
  constructor, except comments beginning with "__". Misspelled keys were
  previously ignored. Validation covers exactly the constructors TUV-x
  validates, so the original configurations are accepted unchanged. As in
  TUV-x, merged files and point data are accepted only by base cross sections,
  where C had merged temperature tables under the last file's temperatures.
* Numbers in configurations and data tables are read independently of the
  process locale; a host using a ',' decimal point could not build a model.
  The conversion rewrites a decimal number as digits and an exponent, which
  denotes the same value, and is verified bit for bit against the "C" locale
  for all 186,533 numeric tokens in the embedded data.
* Spectra validate the temperatures and air densities their formulas use;
  temperature-independent kinds no longer reject a 0 K level, matching TUV-x.
* Malformed resolver datasets, NaN or unordered table temperatures, a zero
  Burkholder B and a nonpositive O2 scale factor are rejected at construction
  instead of overreading memory or failing every run.
* The runner publishes files with the replaced file's mode (or the umask
  default) rather than mkstemp's 0600, and a failed run leaves its JSON
  unterminated instead of a valid document. Exhausting memory while checking
  for a second configuration document reports TT_MEMORY, and duplicate keys
  in large mappings are found in O(n log n).

Constructors that TUV-x leaves unvalidated (for example the specialized
quantum yields) accept any keys, as there. The C port still rejects inputs
where the reference would continue with NaN, such as nonpositive temperatures
in temperature-dependent spectra.

Prepared execution, the default since 2026-10-19, reproduces the generic
results bit for bit; the reference suites above exercise it, and
[PERFORMANCE.md](PERFORMANCE.md) records its design and validation.

## Final validation — 2026-10-18

| Check | GCC/GFortran 16.2.0 | Clang/Flang 23.1.3 |
| --- | --- | --- |
| Complete TurboTUV C/reference/CLI suite | 39/39 pass | 39/39 pass |
| Existing TUV-x suite with documented harness portability fixes | 88/88 pass | 88/88 pass |
| Reference, analytic and output numerical comparisons | 71,501,781 pass | 71,501,781 pass |
| Complete-model comparisons, largest relative epsilon error | 0 | 0.999 |
| Original executable's example outputs | Bit-identical | Bit-identical |
| Production data audit | 190 datasets, exact | 190 datasets, exact |
| Test fixture data audit | 110 datasets, exact | 110 datasets, exact |
| Concurrent versus serial output | 4,410,208 values, exact | 4,410,208 values, exact |
| Installed CMake consumer | Static library, pass | Shared library, pass |

All 20 C/kernel/model/concurrency/locale/prepared-execution tests pass with Clang AddressSanitizer,
UndefinedBehaviorSanitizer and leak detection enabled. The optional NetCDF CLI
also passes address/undefined-behavior checks with the leak-checking exception
documented above; the JSON-only CLI passes with leak detection enabled.
The separate ThreadSanitizer build passes the concurrency test. Generated data,
factory cases and extracted reference helpers regenerate without changes.
The four shipped examples match their upstream bytes. The sibling TUV-x
checkout is unchanged.

The numerical limits are specified per test: relative scales for spectra/rates,
incident-flux normalization for transfer fields, and an angle scale near solar
cancellation. The largest transfer error is 46 epsilon-scaled units, within
the 64-unit limit; this remains double-precision roundoff. Undefined upstream
array entries and invalid dormant inputs are handled as described above.

For a C/reference run, first build the corresponding reference with
`tools/reference.py gcc` or `tools/reference.py llvm`. Configure TurboTUV with
the matching C/Fortran compilers and `TURBOTUV_REFERENCE_BUILD`,
`TURBOTUV_REFERENCE_PREFIX`, and `TURBOTUV_NETCDF_PREFIX`. Select a Python with
NumPy/netCDF4 using `Python3_EXECUTABLE` for the CLI reference-output tests.
Then build and run `ctest --output-on-failure`. The ordinary C-only build does
not need Fortran or scientific input files. See [USAGE.md](USAGE.md) for the
build, installation, output and host-interface commands.
