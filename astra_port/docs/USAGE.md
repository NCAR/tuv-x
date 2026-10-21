# Running and embedding TurboTUV

Build with CMake and a C11 compiler. The core library uses the system math
library and bundled C parser sources. Scientific data are immutable arrays
compiled into the library; paths in a configuration select catalogue entries.
Only the configuration file is read at runtime.

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
build/release/turbotuv examples/tuv_5_4.json --output results.json
build/release/turbotuv examples/ts1_tsmlt.yml --zenith 30 --scale 1 > column.json
```

Without `--zenith`, the runner evaluates every edge of the configured time grid,
using the solar-zenith and Earth–Sun distance profiles. With `--zenith`, it runs
one column and defaults to unit radiation scaling. `--scale` overrides either
mode. The reference multiplies radiation by its distance argument directly;
TurboTUV preserves that convention. Output files are replaced only after a
successful calculation and close. `--help` and `--version` are available.

The four JSON/YAML examples are unchanged copies of TUV-x 0.17.0 examples. They
run from any working directory and require no `data/` directory. Test templates
containing `#PHOTO_RATES` or `#DOSE_RATES` still need those markers expanded into
configuration text, as in the original regression harness; the test generator
does that automatically.

Configurations are checked as TUV-x checks them: each section must contain
its required keys and no keys other than those it accepts, apart from
comments whose names begin with `__`. Numbers always use `.` as the decimal
point, whatever locale the host process selects.

## Benchmark

```sh
build/release/turbotuv_bench examples/tuv_5_4.json --columns 1000
```

The benchmark builds one model and evaluates `--columns` columns (default
1000), changing the temperature, air, O2 and O3 profiles, the solar zenith
angle and the radiation scale of each through the host-update API. It times
the updates and runs, after one untimed column, and prints the time per column
and an FNV-1a hash of every radiation field and rate of every column. Equal
hashes from two builds mean bit-identical results. `--repeat` repeats the
column set, `--seed` selects another, `--output` writes all values, and
`--generic` uses generic execution. `--sweep` holds the profiles fixed, runs
column c at a solar zenith angle of 89 c / N degrees, and hashes the radiation
fields and rates only. In a reference build, `tests/tuvx_bench` runs the same
sweep through the original TUV-x, so the two can be timed and hashed against
each other. `tools/bench_compare.py` times several
builds in rotating rounds and checks that their hashes agree; see
[PERFORMANCE.md](PERFORMANCE.md).

## Output

JSON uses the `TurboTUV-1` schema. Root fields contain altitude edges in km,
wavelength midpoints and edges in nm, and ordered photolysis/heating/dose labels.
Each `steps` entry contains time in hours, zenith in degrees, radiation scaling,
temperature in K, all six radiation fields, physical actinic flux and spectral
irradiance, and the three rate arrays. Doubles use 17 significant digits for
lossless round trips.

Radiation arrays have shape `[wavelength][level]`. Levels ascend in altitude.
The six solver fields are dimensionless, including the requested scaling.
Actinic flux is in photons cm⁻² s⁻¹ integrated over each wavelength bin;
spectral irradiance is in W m⁻² integrated over each bin. Rate arrays have shape
`[rate][level]`: photolysis in s⁻¹, heating in J s⁻¹, dose in W m⁻² with the
configured biological or instrumental weighting. Empty rate collections are
empty arrays. Labels map directly to the first index.

Optional NetCDF output uses the NetCDF C library, solely for writing results:

```sh
cmake -S . -B build/netcdf -DCMAKE_BUILD_TYPE=Release -DTURBOTUV_NETCDF_OUTPUT=ON
cmake --build build/netcdf -j
build/netcdf/turbotuv examples/ts1_tsmlt.json --netcdf results.nc
```

Use `PKG_CONFIG_PATH` if NetCDF is installed outside the default search path.
The core library remains independent of NetCDF. NetCDF-4 output contains root
coordinate/radiation variables and `photolysis`, `heating`, and `dose` groups.
Dimensions are `[time,wavelength,vertical_level]` for spectra and
`[time,vertical_level]` for rates. Rate variable names follow the reference's
slash-to-underscore convention and retain the exact label in an attribute.
Radiation names match the original output writer; their units are correctly
marked dimensionless. Heating and photolysis can use the same label because
they occupy separate groups. Time is the leading dimension in C output;
the original Fortran writer places it last.

## Library

Install and consume the exported CMake target:

```sh
cmake --install build/release --prefix /path/to/install
```

```cmake
find_package(TurboTUV CONFIG REQUIRED)
add_executable(column column.c)
target_link_libraries(column PRIVATE TurboTUV::turbotuv)
```

[examples/column.c](../examples/column.c) is a complete minimal caller. The
public headers document status returns, ownership, dimensions and units:
`turbotuv.h` provides the numerical kernels and shared structures,
`turbotuv_config.h` provides immutable configuration nodes, and
`turbotuv_model.h` provides prepared spectra and complete models.

A model owns all arrays and cached spectra. Destroying its input configuration
or host inputs is safe after construction. Result and lookup views remain owned
by the model; results are overwritten by its next run. Use separate model
instances for simultaneous calls, including separate MPI ranks or threads.
There is no mutable global model state or required MPI/OpenMP runtime.

`tt_model_create_with_inputs` accepts additional named host grids, profiles and
radiators. Omit their names from the corresponding configuration collections;
duplicates are errors. Grids are installed before profiles are constructed,
and profiles before dependent radiators/rates. All supplied data are copied.
The same model then supports repeated column updates:

* `tt_model_update_grid` changes altitude edges and optional midpoints.
* `tt_model_update_profile` changes edges and optional midpoints, columns and
  the extra column above the model top. Omitting columns preserves them.
* `tt_model_set_radiator` supplies optical depths, optional albedos and
  asymmetries; it replaces that radiator's calculated optical properties.
* `tt_model_set_radiation` supplies the six dimensionless fields for the
  `from host` solver. Flux components are required; omitted irradiances are zero.

Updates preserve dimensions and borrowed view addresses. A height change does
not resample profiles or static aerosol/radiator data: the host explicitly
updates those quantities. Construct a new model to change wavelength grids,
because spectra and weights are prepared on that grid. Host setters copy data
and reject incorrect dimensions and nonfinite input before changing state.

Runs use prepared execution by default. Column-independent spectral work is
done when the model is built, and each run evaluates only what can change.
`tt_model_set_execution(model, TT_EXECUTION_GENERIC)` selects the direct
translation instead. Both give bit-identical results and the same error status
for every input: a prepared run hands any input outside its assumptions to the
generic sequence. A prepared model holds its plans and the rows of its varying
spectra: 4.6 MB more peak memory for tuv_5_4 (17.1 MB in all) and 3.1 MB more
for ts1_tsmlt, and building it takes 2–4 ms longer.

For additional scientific datasets, regenerate the catalogue with
`tools/embed_data.py`, or provide immutable compiled arrays through
`tt_data_source`. Model factories never open scientific input files.
