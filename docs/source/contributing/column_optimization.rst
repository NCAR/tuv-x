Column optimization
===================

The Fortran column path adapts optimizations from ``origin/astra_port``
(``astra_port/docs/PERFORMANCE.md`` and ``src/radiation_block.c``):

* Delta-Eddington coefficients, slant-depth sums and tridiagonal elimination
  process blocks of 16 independent wavelengths, with wavelength as the
  contiguous scratch-array dimension. Partial blocks process only valid lanes.
  Direct-beam exponentials are evaluated once per interface and reused.
* Base cross sections and quantum yields return one row to the rate integrals.
  Derived calculators retain dynamic dispatch and full height dependence.
  Band overrides are applied before broadcasting; diagnostics and O2-band
  corrections receive full arrays.
* Photolysis and heating integrals traverse wavelengths in ascending order
  and update contiguous heights, eliminating product and flux transposes.
* Exponential temperature parameterizations calculate wavelength powers once
  per column. All parameterizations extend their interpolation grid once per
  column, rather than reallocating it four times at every height.

The ordinary ``calculate`` interfaces remain height-resolved. The new
``calculate_column`` interfaces on cross sections and quantum yields return
one height row only for the exact base type; callers must broadcast that row.
No atmospheric values are cached between calls, so host updates cannot leave
stale column data. The C port's specialized plans for individual spectral
kernels, spectrum sharing and sparse rate terms are not included here.

The default Delta-Eddington execution is ``blocked``. Select the retained
scalar algorithm at runtime for comparison:

.. code-block:: json

   { "type": "delta eddington", "execution": "scalar" }

The execution setting is included in MPI serialization. Scalar exponential
calls preserve the reference math-library evaluations; changing summation
layout may still change results at numerical roundoff. Fast-math flags are
not required.

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

Local verification (GNU Fortran 15.2, Release ``-O3``) passed all 89 tests
in serial Release, checked Debug, and two-rank MPI builds. Comparing 300
changing columns for each example against ``main`` checked 63,960,600
values: the maximum normalized difference was 3.14 machine epsilons, and
all four diffuse radiation fields were bit-identical.

Three alternating timing rounds of 300 columns gave the following medians
on the development machine. These are local measurements, not portable
performance guarantees; the benchmark above reproduces the workload.

.. list-table:: Median time per column
   :header-rows: 1

   * - Configuration
     - Main
     - Optimized Fortran
     - Speedup
   * - ``tuv_5_4``
     - 11.00 ms
     - 8.11 ms
     - 1.36x
   * - ``ts1_tsmlt``
     - 12.03 ms
     - 7.04 ms
     - 1.71x
