# Prepared execution

TurboTUV runs a column in one of two ways, with bit-identical results:

* **Generic** (`TT_EXECUTION_GENERIC`): the direct translation of TUV-x. Every
  run calls the public kernels in sequence: each spectrum at every level and
  wavelength, mixing, the Delta-Eddington or discrete-ordinate solver, and each
  rate as a sum over every wavelength.
* **Prepared** (`TT_EXECUTION_PREPARED`, the default): work that cannot change
  between columns is done once, when the model is built, and the rest runs in
  kernels arranged for the hardware. Every output value is the generic value,
  bit for bit.

The baseline is the reviewed generic implementation, preserved as branch
`baseline` and tag `v0.1.0-baseline`. Generic execution in the current library
is that same code.

## Results

`bench/turbotuv_bench.c` builds one model and evaluates 300 columns. Each column
changes the temperature (±10 K offset plus a ±8 K wave), the air and O2 profiles
(±5%), the O3 profile (±30%), the solar zenith angle (0–95°) and the radiation
scale (±1.7%) through the host-update API, so no column-dependent work can be
reused. `tools/bench_compare.py` ran the baseline and prepared programs in five
rotating rounds on an idle Apple M5 Pro (macOS 26.7.1, AC power), single
threaded, CMake Release (`-O3 -DNDEBUG -ffp-contract=off`). Each run printed a
hash of every radiation field and rate of every column. Within each compiler,
the baseline and prepared hashes were equal: the results are bit-identical.

Median time per column (range over the five rounds):

| Configuration | Compiler | Baseline | Prepared | Speedup |
| --- | --- | ---: | ---: | ---: |
| tuv_5_4: 156 wavelengths, 121 levels, 109 photolysis + 28 dose | Clang 23.1.3 | 5456 µs (5426–5554) | 1054 µs (1048–1059) | **5.18x** |
| | GCC 16.2.0 | 8091 µs (8076–8097) | 1064 µs (1056–1071) | **7.61x** |
| ts1_tsmlt: 102 wavelengths, 121 levels, 73 photolysis + 4 heating | Clang 23.1.3 | 4089 µs (4068–4123) | 946 µs (945–966) | **4.32x** |
| | GCC 16.2.0 | 5170 µs (5162–5233) | 1097 µs (1090–1106) | **4.71x** |

Building a model takes 3.5 ms longer for tuv_5_4 (16.5 ms instead of 13.0 ms)
and 1.8 ms longer for ts1_tsmlt; this cost is paid once per model. Peak memory
of the benchmark grows from 12.5 to 17.1 MB (tuv_5_4) and from 10.5 to 13.6 MB
(ts1_tsmlt), mostly the per-run rows of the varying spectra.

Reproduce with:

```sh
git worktree add ../TurboTUV-baseline baseline
cmake -S ../TurboTUV-baseline -B ../TurboTUV-baseline/build -DCMAKE_BUILD_TYPE=Release
cmake --build ../TurboTUV-baseline/build -j
python3 tools/bench_compare.py baseline-llvm=../TurboTUV-baseline/build/turbotuv_bench \
    prepared-llvm=build/release/turbotuv_bench \
    --config examples/tuv_5_4.json --config examples/ts1_tsmlt.json --rounds 5
```

`turbotuv_bench --generic` runs the generic path of the current build.

## Against the original TUV-x

The redundant per-column work is part of TUV-x's own design, which evaluates
every calculator on demand. The C generic path reproduces those evaluations
on purpose, because bit-for-bit parity with TUV-x is what the port was
validated against. In TUV-x 0.17.0 (`../TUV-x`):

* `photolysis_rates.F90` (`get`) calls `calculate` for every reaction's cross
  section and quantum yield on every run, each producing a full height ×
  wavelength array. It then forms `transpose(cross_section * quantum_yield)`
  and takes one dot product per level.
* The base cross section (`cross_section.F90`, `calculate`) and quantum yield
  (`quantum_yield.F90`, `run`) copy their one column-independent wavelength
  table into every height level. They apply overrides level by level and
  transpose into a second new array.
* `grid_warehouse%get_grid` and `profile_warehouse%get_profile` return newly
  allocated copies (`allocate(grid, source = ...)`), which every `calculate`
  call makes.
* Temperature-parameterized cross sections (`temperature_based.F90`) recompute
  `(λ - λ0)**lp` for every term at every height level. They grow the spectrum
  through four `add_point` reallocations and re-interpolate all wavelength
  bins, including those temperature cannot affect.
* Delta-Eddington transfer calls the O(layers) `slant_optical_depth` for every
  layer of every wavelength and evaluates `exp(-τ_slant)` at both ends of each
  layer and again for the outputs.

`bench/tuvx_bench.F90` runs the original TUV-x through the column sweep of
`turbotuv_bench --sweep`. One model runs 300 columns at solar zenith angles
89 c / N degrees with unit radiation scale; only the runs are timed. Both
programs hash the six radiation fields and the photolysis, heating and dose
rates in the same order and memory layout. With a reference build configured
(see PORTING.md), the comparison below ran in five rotating rounds from the
reference runtime directory, where TUV-x's relative data paths resolve:

```sh
cd build/c-gcc/tests/reference-runtime
python3 ../../../../tools/bench_compare.py \
    tuvx-gnu=../tuvx_bench generic-gnu=../../turbotuv_bench:--sweep,--generic \
    prepared-gnu=../../turbotuv_bench:--sweep \
    tuvx-llvm=../../../c-llvm/tests/tuvx_bench \
    generic-llvm=../../../c-llvm/turbotuv_bench:--sweep,--generic \
    prepared-llvm=../../../c-llvm/turbotuv_bench:--sweep \
    --config ../../../../examples/tuv_5_4.json \
    --config ../../../../examples/ts1_tsmlt.json \
    --columns 300 --rounds 5 --allow-differences tuvx-llvm
```

Median time per column (range), with speedups over TUV-x built by the same
compiler family:

| Configuration | Program | GNU (GFortran / GCC 16.2.0) | LLVM (Flang / Clang 23.1.3) |
| --- | --- | ---: | ---: |
| tuv_5_4 | TUV-x | 6389 µs (6372–6793) | 5580 µs (5552–5598) |
| | C generic | 8076 µs, 0.79x | 5311 µs, 1.05x |
| | C prepared | 1038 µs, **6.15x** | 1028 µs, **5.43x** |
| ts1_tsmlt | TUV-x | 6551 µs (6526–6933) | 6242 µs (6172–6549) |
| | C generic | 5112 µs, 1.28x | 3967 µs, 1.57x |
| | C prepared | 1080 µs, **6.07x** | 923 µs, **6.76x** |

With the GNU compilers, all three programs printed the same hash for each
configuration: the prepared C results equal the original TUV-x's bit for bit
for every radiation field and rate of all 300 columns. Flang's results differ
from Clang's at roundoff, within the one relative epsilon recorded in
PORTING.md; the two C paths agree with each other exactly.

The C generic path is not uniformly faster than the Fortran: it avoids the
copies, transposes and reallocations but performs the same evaluations, and
GCC compiles the generic tuv_5_4 kernels less well than GFortran does.

## Where the time went

Sampled profiles of the baseline showed that radiative transfer was not the
main cost. For tuv_5_4, about 78% of a column went to evaluating cross sections
and quantum yields, including their `exp`, `log` and `pow` calls. Integrating
the rates took 6.5%, and Delta-Eddington transfer with its slant paths and
tridiagonal solves took 8.5%. For ts1_tsmlt, the 20 temperature-parameterized
cross sections alone took 42% of a column. Each one was re-evaluated at every
level, with a conserving regrid per level.

Yet most of that work could not change between columns. Only 5.7% (tuv_5_4) and
8.2% (ts1_tsmlt) of the cross-section and yield values vary with the
atmospheric state. 61 of 109 and 30 of 73 reactions have spectra that are
entirely independent of it.

Prepared execution, per stage (instrumented, Clang):

| Stage | tuv_5_4 | ts1_tsmlt |
| --- | ---: | ---: |
| Delta-Eddington transfer and mixing | 409 µs | 283 µs |
| Varying reaction spectra | 250 µs | 442 µs |
| Rate integrals | 263 µs | 103 µs |
| O2 band optics (radiator and reactions) | 77 µs | 77 µs |
| Geometry, gas optics, actinic flux | 38 µs | 29 µs |

## Design

The design adapts the TurboMICM reference solution (`../TurboMICM/fable`) to
the structure of TUV-x. That solution has four parts: vector-ordered data with
lanes processed in lock step, specialized kernels prepared from the mechanism,
a fingerprint so that stale kernels cannot run, and the generic kernels kept as
the fallback and test reference. Here they become:

* **Plans instead of generated code.** TUV-x's work is dense arrays over
  wavelength and level, not a sparse system, so literal indices in generated C
  would save little. Instead a plan is built from the model's own prepared
  spectra when the model is created. Because it is derived from the model
  itself, it cannot go stale, which is what the fingerprint guarded against.
  * `spectrum_plan.c` classifies every wavelength of every cross section and
    quantum yield. A wavelength is constant where the kernel's own branches and
    table data make the value independent of temperature and air density: base
    tables, out-of-band branches, zero bases under a bounded exponential,
    temperature tables whose columns are equal, and so on. The generic kernel
    evaluates constant wavelengths once.
  * Varying wavelengths are evaluated with the kernel's expressions, rearranged
    only by hoisting subexpressions of wavelength alone (such as `log(329.5/λ)`
    for Cl2 or `pow(λ - λ0, p)` for the exponential parameterization) and of
    level alone (such as the acetone yield's eleven `pow` calls). The
    remaining operations keep the generic order.
  * A guard re-evaluates the generic kernel on reference columns that span every
    clamp and branch. A wavelength classified as constant must be constant
    there. If the guard fails, that spectrum falls back to generic evaluation at
    every wavelength. The tests require that this never happens.
  * Temperature-parameterized cross sections recompute only the target bins
    that touch parameterized points, with the same conserving arithmetic.
* **Rates as precomputed terms.** `model_plan.c` turns each photolysis, heating
  and dose rate into terms in ascending wavelength order. Where both factors
  are constant, a term carries the product `σφ` exactly as the generic rate
  forms it. Terms whose product is zero are omitted: adding `+0` or `-0` to a
  sum that starts at `+0` cannot change it, provided the flux is finite, which
  every run checks. Equal spectra are evaluated once and shared, for example
  the O3 cross section used by two reactions. Sums are tiled over levels so
  that partial sums stay in registers; that helps GCC slightly and slows
  Clang, as measured below.
* **Lanes of wavelengths.** `radiation_block.c` solves Delta-Eddington transfer
  for blocks of 16 wavelengths in lock step. This is the analog of TurboMICM's
  blocks of cells, with each wavelength an independent problem sharing the ray
  geometry. Radiator mixing is fused into the first stage. The tridiagonal rows
  are assembled inside the elimination. `exp(-τ_slant)` is evaluated once per
  level instead of up to four times. Lanes past the end of the grid repeat the
  last wavelength. The compiler vectorizes across lanes from the data layout
  alone, with no intrinsics or pragmas.
* **The generic path stays.** `tt_model_set_execution` selects either path.
  When a prepared run meets an input outside its assumptions (invalid values,
  temperatures above 10⁴ K, or nonfinite radiation), it declines, and the
  generic sequence runs instead. Results and error statuses are therefore those
  of the generic path for every input. Discrete-ordinate transfer uses the
  generic solver; its spectra and rates are prepared.

TurboMICM also offered a "fast" mode, which changed results at roundoff
(reassociated sums and reciprocal multiplies). That was not needed here: every
transformation above preserves each floating-point operation and its order. No
fused multiply-add is introduced, and no vectorized math library is used. The
exponentials and powers are the platform's scalar libm calls, as in the
generic path.

## Speedup by optimization

Scratch builds of the library, not part of the repository, measured what each
optimization contributes. They time every stage of a column in both execution
modes, and can switch off one optimization at a time:

| Switch | What the prepared run does instead |
| --- | --- |
| no classification | evaluates every wavelength of every spectrum with the generic kernel, so every rate term reads both factors |
| no sharing | evaluates equal spectra separately |
| all bins | recomputes every bin of a parameterized cross section |
| zero terms kept | keeps rate terms whose product is zero |
| no tiling | sums each rate over all 121 levels in one pass |
| generic Delta-Eddington | mixes the radiators, calls `tt_delta_eddington` and scales, as the generic run does |
| one lane | runs the blocked solver one wavelength at a time |

Every variant printed the same output hash as the release build, so each
switch changes only speed. The stage timers change a column's time by about
1% or less, within the run-to-run variation. The programs ran the benchmark's
varied columns in five rotating rounds of 300 columns; times are medians.
Whole-column differences below about 3% are within the variation; the stage
times are steadier. The generic times come from the current library rather
than the baseline build, so they differ slightly from the results above.

Each optimization's saving is the time it removes from its stage between the
generic and the prepared run. Where several optimizations act on one stage,
the switches divide that stage's saving among them. Time saved per column,
with each optimization's share of the total saving:

| Optimization | tuv_5_4 Clang | tuv_5_4 GCC | ts1_tsmlt Clang | ts1_tsmlt GCC |
| --- | ---: | ---: | ---: | ---: |
| Spectrum plans: constant values once per model, varying values with hoisted subexpressions | 3,540 µs (83%) | 5,983 µs (86%) | 2,342 µs (76%) | 3,226 µs (79%) |
| Parameterized cross sections: only the affected bins | 1 µs (0%) | 6 µs (0%) | 219 µs (7%) | 204 µs (5%) |
| Equal spectra evaluated once | 32 µs (1%) | 35 µs (1%) | 69 µs (2%) | 71 µs (2%) |
| Rate terms, net | 299 µs (7%) | 410 µs (6%) | 126 µs (4%) | 168 µs (4%) |
| · zero products omitted | 325 µs | 170 µs | 114 µs | 54 µs |
| · constant products precomputed | 524 µs | 617 µs | 194 µs | 255 µs |
| · term loop in place of the generic sweep | −421 µs | −388 µs | −134 µs | −146 µs |
| · 16-level tiles | −129 µs | 11 µs | −48 µs | 5 µs |
| Blocked Delta-Eddington, net | 426 µs (10%) | 567 µs (8%) | 280 µs (9%) | 374 µs (9%) |
| · mixing and scale fused, one `exp` per level, rows built in the elimination | 166 µs | 162 µs | 124 µs | 129 µs |
| · 16 wavelengths in lock step | 260 µs | 405 µs | 156 µs | 245 µs |
| O2 band cross sections once per run | 0 µs (0%) | 1 µs (0%) | 34 µs (1%) | 53 µs (1%) |
| Flux finiteness check, which makes zero omission exact | −9 µs (−0.2%) | −14 µs (−0.2%) | −6 µs (−0.2%) | −10 µs (−0.2%) |
| **Generic → prepared** | 5,321 → 1,031 µs | 8,044 → 1,054 µs | 3,985 → 923 µs | 5,174 → 1,096 µs |

Applied one after another in this order, the optimizations multiply to the
full speedup:

| Step | tuv_5_4 Clang | tuv_5_4 GCC | ts1_tsmlt Clang | ts1_tsmlt GCC |
| --- | ---: | ---: | ---: | ---: |
| Generic | 5,321 µs | 8,044 µs | 3,985 µs | 5,174 µs |
| Spectrum plans, with parameterized bins | ×2.99 | ×3.92 | ×2.80 | ×2.96 |
| Equal spectra shared | ×1.02 | ×1.02 | ×1.05 | ×1.04 |
| Rate terms | ×1.21 | ×1.26 | ×1.10 | ×1.11 |
| Blocked Delta-Eddington | ×1.42 | ×1.54 | ×1.29 | ×1.33 |
| O2 cross sections once, flux check | ×0.99 | ×0.99 | ×1.03 | ×1.04 |
| **Prepared** | 1,031 µs: **5.16x** | 1,054 µs: **7.63x** | 923 µs: **4.32x** | 1,096 µs: **4.72x** |

Prepared time per column with one optimization switched off:

| Switched off | tuv_5_4 Clang | tuv_5_4 GCC | ts1_tsmlt Clang | ts1_tsmlt GCC |
| --- | ---: | ---: | ---: | ---: |
| nothing (prepared) | 1,031 µs | 1,054 µs | 923 µs | 1,096 µs |
| constant classification | 4,404 µs (×4.27) | 6,010 µs (×5.70) | 3,321 µs (×3.60) | 4,132 µs (×3.77) |
| sharing | 1,063 µs (×1.03) | 1,105 µs (×1.05) | 992 µs (×1.07) | 1,167 µs (×1.07) |
| parameterized bins | 1,036 µs (×1.00) | 1,083 µs (×1.03) | 1,195 µs (×1.29) | 1,332 µs (×1.22) |
| zero-term omission | 1,356 µs (×1.31) | 1,237 µs (×1.17) | 1,044 µs (×1.13) | 1,148 µs (×1.05) |
| 16-level tiles | 911 µs (×0.88) | 1,070 µs (×1.02) | 879 µs (×0.95) | 1,095 µs (×1.00) |
| blocked Delta-Eddington | 1,460 µs (×1.42) | 1,679 µs (×1.59) | 1,203 µs (×1.30) | 1,467 µs (×1.34) |
| 16 lanes | 1,292 µs (×1.25) | 1,465 µs (×1.39) | 1,080 µs (×1.17) | 1,334 µs (×1.22) |

* Spectrum plans are 76–86% of the saving. Alone, they make a column 2.8–3.9
  times faster. GCC gains most because its generic spectral kernels are the
  slowest.
* Sharing interacts with the plans. Applied to the generic evaluation,
  evaluating equal spectra once would save 1.1 ms (Clang) or 1.9 ms (GCC) per
  tuv_5_4 column, because many reactions share constant quantum yields. Once
  constant values are precomputed, sharing saves only 30–70 µs. The tables
  count it after the plans.
* Rate terms save 4–7% net. Omitting zero products and precomputing constant
  products together save 0.3–0.85 ms. The term loop gives some of that back,
  because per term it is slower than the generic sweep over contiguous arrays.
  With Clang, the 16-level tiles double the rate stage (261 µs against 132 µs
  for tuv_5_4); with GCC they save 5–11 µs.
* The blocked Delta-Eddington solver saves 8–10% of the total. Fusing the
  mixing and scale, evaluating each `exp` once and building rows in the
  elimination account for 30–45% of that; running 16 wavelengths in lock step
  for the rest.
* Recomputing only the affected bins of parameterized cross sections matters
  only for ts1_tsmlt, with its 20 such spectra: 5–7%. Computing the O2 band
  cross sections once saves one call for ts1_tsmlt's second O2 reaction. The
  flux finiteness check that makes zero omission exact costs 6–14 µs.

## Validation

* `tests/test_prepared.c` checks the prepared path against the generic path,
  bit for bit:
  * every spectrum plan against its generic kernel, for all 912 cross-section
    and quantum-yield configurations of the spectrum test catalogue, including
    every kind and all four parameterizations. Each is evaluated at 41 levels
    with randomized and branch-edge temperatures and densities;
  * the blocked Delta-Eddington solver against mixing plus the generic solver,
    on 550 atmospheres: 1–81 layers, zenith 0–180°, vacuum, thin and thick
    layers, pure scattering and surface albedos 0–1, with partial and multiple
    blocks;
  * all 26 runnable model configurations on 14 perturbed columns each. These
    include daylight, twilight, night, a level beyond the plans' temperature
    limit and an invalid level, and compare every field, rate and radiator view
    together with run statuses and error messages.
  In all, 124 million values are compared per toolchain.
* The reference suites now run the prepared path by default and still match
  the original Fortran: bit-identical with GCC/GFortran, within one relative
  epsilon with Clang/Flang, as before. The host-update suite runs every column
  both ways.
* All tests pass with GCC 16.2.0 and Clang 23.1.3 and under Clang's address,
  undefined-behavior and leak sanitizers. The concurrent-model test passes
  ThreadSanitizer.
* The benchmark hashes above compare the prepared path with the preserved
  baseline build over 300 columns of every output.

## Further opportunities

The prepared path is now bound by irreducible scalar `exp`, `pow` and `log`
calls, and by work that is inherently O(levels²): the slant optical depths of
the spherical geometry. Remaining measured opportunities, none needed for the
present speedups:

* summing each rate over all levels in one pass instead of 16-level tiles.
  With Clang, columns take 5–12% less time (911 instead of 1,031 µs for
  tuv_5_4); with GCC they change by at most 2% (see
  [Speedup by optimization](#speedup-by-optimization));
* register blocking of the rate integrals across reactions;
* computing the Schumann-Runge cross sections once for both the O2 radiator
  and the O2 photolysis reaction (77 µs per column today);
* sharing parameterized-cross-section levels that fall in the same fixed
  temperature range (most ts1 parameterizations fix T below 200–280 K);
* threads: models are independent, so a host can run one per thread, as
  `tests/test_parallel.c` does.
