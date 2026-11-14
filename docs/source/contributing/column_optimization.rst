Column optimization
===================

The Fortran column path adapts the prepared execution in ``origin/astra_port``
and the sibling TurboTUV project (``src/spectrum_plan.c``, ``src/model_plan.c``
and ``src/radiation_block.c``):

* Spectrum plans classify constant wavelengths from each calculator's branches
  and table data. Constants are evaluated once; only varying wavelengths have
  height-resolved storage and are recalculated for each column. Prepared
  kernels hoist wavelength expressions and evaluate temperature/density
  expressions once per level, retaining the original clamps and formulas.
* Equal spectra share a plan across photolysis, heating and radiators, with
  distinct interface and midpoint plans. Sharing requires the same dynamic
  type and exact numeric configuration, including tables, parameterizations,
  branch settings and band overrides. Per-column epochs prevent reuse of
  atmospheric values in a later column.
* Temperature-parameterized cross sections prepare wavelength powers, the
  extended source grid and conserving interpolation segments once. Each column
  recalculates only target bins touching parameterized source points. Levels
  with identical effective temperatures in every active range share the same
  calculation, including fixed and overlapping ranges.
* Photolysis and heating plans precompute constant cross-section/yield products
  and traverse active wavelengths in ascending order while updating contiguous
  heights. Exactly zero constant terms are omitted only for finite flux.
  Heating retains the original energy-times-product grouping. O2 corrections
  are shared between rates within a column; the corrected bands remain
  height-dependent. Dose plans cache spectral weights and their active indices.
* Delta-Eddington processes blocks of 16 independent wavelengths with a
  contiguous lane dimension. Radiator mixing and delta scaling are fused;
  direct-beam exponentials are evaluated once per interface; matrix rows are
  assembled during tridiagonal elimination. Output is written bottom-up with
  radiation scaling fused into the final pass. Partial blocks process only
  valid lanes.

The ordinary ``calculate`` interfaces remain height-resolved. Prepared plans
materialize full arrays for radiator calculations and O2 corrections.
Photolysis diagnostics use the generic path; other diagnostic output keeps its
original dimensions and scaling. Discrete-ordinate transfer keeps its existing
solver and can use prepared spectra and rates.

Execution and cache lifetime
----------------------------

Prepared spectra and rates are the default. Select the generic path in the
top-level model configuration, or call ``core%set_prepared_execution(.false.)``:

.. code-block:: json

   { "column execution": "generic" }

The Delta-Eddington solver independently defaults to ``blocked``. Select its
retained scalar reference in the solver configuration:

.. code-block:: json

   { "type": "delta eddington", "execution": "scalar" }

Each core owns its cache; separate models do not share mutable plans. Standalone
photolysis/heating collections own a cache when no core cache is attached.
Temperature and air profiles are read afresh each column. Changes to wavelength
grid contents or dimensions invalidate the plans automatically. After changing
calculator tables, parameterization settings or spectral weights directly,
call ``core%invalidate_prepared()`` before the next run. MPI serialization
preserves the execution settings and scientific data; receiving models rebuild
their own transient plans without serializing cache pointers or values.

Analytic constant classifications are checked against generic evaluations on
three probe columns. A failed check or an unknown calculator subtype uses the
generic calculator. Configuration-only preparation and probes preserve the
caller's floating-point flags and halting modes. Invalid/nonfinite column
inputs, temperatures outside (0, 10000] K, or nonfinite flux select the generic
path. Unsupported radiator optical states use the scalar transfer reference.
Fast-math is unnecessary. Math-library evaluations and accumulation order are
retained where their rounding can be amplified, including chlorine's small
spectral tails and the transfer solve.

Validation and timing
---------------------

The ``column_optimization`` unit test compares all six radiation fields with
the retained scalar solver over 2,940 atmospheres: 1--81 layers, 1--33
wavelengths, zenith angles 0--180 degrees, surface albedos 0--1, vacuum,
thin/thick and purely scattering layers. It also checks compact spectra,
overrides, both exponential bases, inverted temperature offsets, overlapping
and fixed temperature ranges, and dispatch to a derived parameterization.
Radiation comparisons allow propagation of roundoff through the solve,
including cancellation at zero-flux boundaries. Existing scientific regression
tolerances are unchanged.

The ``prepared_spectra`` test checks 911 configurations generated from the
repository's scientific fixtures and examples, plus interpolation/range cases.
It compares 30,534,996 values over six changing columns, checks that equal
configurations share plans, and exercises guards, unknown subtypes, profile
units and epoch rollover. IEEE classifications are compared for deliberately
overflowing synthetic tables. Regenerate the catalogue with
``python tool/generate_prepared_cases.py``; ``--check`` verifies its provenance
and contents without writing it.

The ``prepared_model`` test compares all six radiation fields and every
photolysis, heating and dose rate against generic spectra/rates and scalar
transfer for both examples. It changes profiles, zenith angle and radiation
scaling, switches execution modes, invalidates caches, changes wavelength grid
contents, and packs/unpacks warmed models in MPI builds. It also exercises
out-of-range temperature and nonfinite-flux fallback.

Build the optional benchmark without requiring LAPACK or Google Benchmark:

.. code-block:: sh

   cmake -S . -B build -DTUVX_ENABLE_COLUMN_BENCHMARK=ON
   cmake --build build -j
   cd build
   ./tuvx_column_sweep ../examples/tuv_5_4.json 300 /tmp/optimized.bin

Run the same benchmark source linked against an unchanged ``main`` build,
then compare every output using NumPy:

.. code-block:: sh

   python ../tool/compare_column_outputs.py /tmp/main.bin /tmp/optimized.bin

Repeat for ``ts1_tsmlt.json``. The sweep changes temperature, air/O2/O3
profiles, solar zenith angle (0--100 degrees), and Earth--Sun distance for
each column. Only ``core%run`` is timed; serialization is outside the timed
region. The comparator checks dimensions, finiteness and every radiation
field and photolysis, heating and dose rate. Its default bound is 16 machine
epsilons, normalized by unit incident radiation or the individual reaction's
maximum magnitude over height in that column.

Run CTest serially when exercising the existing assertion-failure scripts:
some scripts share ``error.json`` in their working directory.

Local verification (GNU Fortran 15.2, Release ``-O3``) passed all 91 tests
in serial Release, checked Debug, and two-rank MPI builds. Comparing 300
changing columns for each example against ``main`` checked 63,960,600
values: the maximum normalized difference was 3.51 machine epsilons, and
all four diffuse radiation fields were bit-identical.

Three alternating timing rounds of 300 columns gave the following medians
on an Intel Core Ultra 7 255H with GNU Fortran 15.2 and ``-O3``. These use the
examples' diagnostic settings, including file output within ``core%run``.
These are local measurements; the benchmark above reproduces the workload.

.. list-table:: Median time per column
   :header-rows: 1

   * - Configuration
     - Main
     - Optimized Fortran
     - Speedup
   * - ``tuv_5_4``
     - 10.25 ms
     - 4.76 ms
     - 2.15x
   * - ``ts1_tsmlt``
     - 11.88 ms
     - 3.40 ms
     - 3.49x
