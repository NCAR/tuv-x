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
  Each kernel then loops over contiguous heights for one wavelength at a
  time, so compilers vectorize it as they do the generic array expressions.
  No element is evaluated speculatively.
* Equal spectra share a plan across photolysis, heating and radiators, with
  distinct interface and midpoint plans. Sharing requires the same dynamic
  type and exact numeric configuration, including tables, parameterizations,
  branch settings and band overrides. Per-column epochs prevent reuse of
  atmospheric values in a later column.
* Temperature-parameterized cross sections prepare wavelength powers, the
  extended source grid and conserving interpolation segments once. Each column
  recalculates only target bins touching parameterized source points. Levels
  with identical effective temperatures in every active range share the same
  calculation, including fixed and overlapping ranges. The remaining levels
  are interpolated together, each in the scalar segment order.
* Photolysis and heating plans precompute constant cross-section/yield products
  and traverse active wavelengths in ascending order while updating contiguous
  heights. Exactly zero constant terms are omitted only for finite flux.
  Heating retains the original energy-times-product grouping. O2 corrections
  are shared between rates within a column; the corrected bands remain
  height-dependent. Dose plans cache spectral weights and their active indices.
* Delta-Eddington processes blocks of 16 independent wavelengths with a
  contiguous lane dimension. Radiators are mixed in warehouse order over
  contiguous layers, in the same pass that checks their optical properties;
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
path. Radiator optical properties that are not finite, a negative optical
depth, a single scattering albedo outside 0--1 or an asymmetry factor of
magnitude 1 or more use the scalar transfer reference for every wavelength.

Fast-math is unnecessary. Every sum keeps the generic accumulation order, and
the transfer solve keeps scalar math-library calls. Spectrum kernels may use
the compiler's vector math library, as the generic array expressions already
do, so prepared and generic spectra agree to roundoff rather than bit for
bit. The wavelength logarithms of the chlorine spectrum use the generic
contiguous expression, because a rounding difference there is amplified in
the spectrum's small tails.

Parallel use
------------

A model computes one column per ``run`` call on one core; there is no
threading inside a column. Hosts parallelize over columns, with one model for
each MPI task or OpenMP thread, exactly as for generic execution. A model
unpacked from an MPI buffer builds its own plans on its first column. Plans
and caches are never shared between models, so separate models can run
concurrently.

Validation and timing
---------------------

The ``column_optimization`` unit test compares all six radiation fields with
the retained scalar solver over 2,940 atmospheres: 1--81 layers, 1--33
wavelengths, zenith angles 0--180 degrees, surface albedos 0--1, vacuum,
thin/thick and purely scattering layers. A further 210 atmospheres hold one
unsupported optical property (out of range, NaN or infinite) in the first,
an interior or the last block; their fields must equal the scalar solver's
exactly. It also checks compact spectra,
overrides, both exponential bases, inverted temperature offsets, overlapping
and fixed temperature ranges, and dispatch to a derived parameterization.
Radiation comparisons allow propagation of roundoff through the solve,
including cancellation at zero-flux boundaries. Existing scientific regression
tolerances are unchanged.

The ``prepared_spectra`` test checks 913 configurations generated from the
repository's scientific fixtures and examples, plus interpolation/range cases
and the JPL06 ozone cross section, which no fixture constructs.
It compares 61,596,064 values over eight changing columns, for plans at
interfaces and, for cross sections, at layer mid-points, whose profiles
differ. The columns include every branch temperature exactly and
temperatures from 1 K to the 10000 K limit. It checks that equal
configurations share plans, and exercises guards, unknown subtypes, profile
units and epoch rollover. IEEE classifications are compared for deliberately
overflowing synthetic tables and temperatures. Regenerate the catalogue with
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

Local verification (GNU Fortran 15.2) passed all 91 tests in serial Release
(``-O3``), in a checked Debug build (``-O0 -fcheck=all
-ffpe-trap=invalid,zero,overflow``) and in a two-rank MPI build with four
OpenMP threads. Comparing 300 changing columns for each example against
``main`` checked 63,960,600 values: the maximum normalized difference was
3.14 machine epsilons, and all four diffuse radiation fields were
bit-identical.

Alternating rounds of 60 columns gave the following medians on one
performance core of an Intel Core Ultra 7 255H with GNU Fortran 15.2 and
``-O3``: ten rounds for ``main`` and five for the optimized build. These use
the examples' diagnostic settings, including file output within
``core%run``; without diagnostics the optimized times are 1.34 ms and
1.11 ms. This processor slows by more than a factor of two under sustained
load, so the runs were short, pinned with ``taskset`` and separated by
pauses. ``main`` is the less repeatable of the two: its rounds ranged over
6.0--6.9 ms and 6.9--9.2 ms, while the optimized build stayed within
1 percent. These are local measurements; the benchmark above reproduces the
workload.

.. list-table:: Median time per column
   :header-rows: 1

   * - Configuration
     - Main
     - Optimized Fortran
     - Speedup
   * - ``tuv_5_4``
     - 6.56 ms
     - 1.41 ms
     - 4.7x
   * - ``ts1_tsmlt``
     - 8.49 ms
     - 1.13 ms
     - 7.5x

Separate models do not slow one another beyond the processor's own limits:
six models on the six performance cores computed 4.8 times as many
``tuv_5_4`` columns per second as one model.

Remaining cost
--------------

Stage timers in a scratch build (not part of the repository) attribute the
optimized column as follows, without diagnostics:

.. list-table:: Time per column by stage, in microseconds
   :header-rows: 1

   * - Stage
     - ``tuv_5_4``
     - ``ts1_tsmlt``
   * - Delta-Eddington solve
     - 535
     - 361
   * - Radiator states and O2 band optical depth
     - 137
     - 113
   * - Varying spectra: prepared kernels
     - 311
     - 137
   * - Varying spectra: temperature-parameterized cross sections
     - --
     - 297
   * - Rate integrals
     - 153
     - 90
   * - Dose rates
     - 102
     - --
   * - Spherical geometry
     - 25
     - 25
   * - Actinic flux, O2 band cross sections and plan lookup
     - 97
     - 94

Within the solve, the slant optical depths (138 microseconds for ``tuv_5_4``)
and the layer coefficients with their exponentials (136) are the largest
parts. Measured opportunities that remain:

* A default x86-64 build vectorizes two doubles at a time. The solver's lanes
  and the kernels' height loops widen with the target instruction set, for
  example AVX2, without source changes.
* The first loop of the layer coefficients contains conditional divisions and
  is not vectorized.
* The Schumann-Runge band cross sections are evaluated separately for the O2
  radiator and the O2 reactions.
* Radiator states are rebuilt from materialized cross sections for every
  column.
* Temperature-parameterized cross sections are bound by the power function of
  the parameterization (172 microseconds for ``ts1_tsmlt``).
