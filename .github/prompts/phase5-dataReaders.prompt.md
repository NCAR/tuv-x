# Phase 5: Data Reader Abstraction (Step 19)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

The current Fortran code reads ~160 NetCDF files for cross-section and quantum yield parameter data at initialization time. The new library needs a pluggable data reader interface with a default NetCDF-C implementation (no Fortran dependency).

## Step 19: Define pluggable data reader interface

### Interface design

```cpp
class DataReader {
public:
    virtual ~DataReader() = default;

    /// Read parameter data as a 2D array [wavelength × parameter_index]
    virtual Array2D<double> read_parameters(const std::string& variable_prefix) = 0;

    /// Read the wavelength grid [m]
    virtual Array1D<double> read_wavelengths() = 0;

    /// Read temperature grid (optional — only for T-dependent data)
    virtual Array1D<double> read_temperatures() = 0;
};
```

**Note**: The `DataReader` returns `Array1D` (not `std::vector`) so that returned data can be directly used with policy-backed types (`Grid`, `Profile`, etc.) without copies. When the solver runs on GPU, these arrays may need to reside in device memory — `std::vector` cannot support this.

### NetCDF-C reader implementation

Implement `NetCDFReader : DataReader` using the NetCDF-C library (`libnetcdf`):
- Opens `.nc` files and reads `wavelength`, `temperature`, `{prefix}_parameters` variables
- The existing ~160 NetCDF files in `data/cross_sections/` and `data/quantum_yields/` use a consistent naming convention with variable prefixes `cross_section_`, `quantum_yield_`, `temperature_`
- No unit conversions: all data files must provide wavelengths in meters (m) and temperatures in Kelvin (K). If legacy data files use nm, they must be converted to SI units before use by TUV-x

### NetCDF variable structure (existing files)
- `wavelength` — 1D array of wavelengths [m]
- `temperature` — 1D array of temperatures [K] (optional)
- `{prefix}_parameters` — 2D array `[wavelength × parameter_type]`

### Future readers
The pluggable interface allows adding CSV, HDF5, binary, or in-memory readers without changing transform code. This is particularly useful for:
- Testing (in-memory data readers)
- Alternative data formats
- Embedded/compiled-in data for deployment without data files

### Fortran reference
- `src/netcdf.F90` — NetCDF reading utilities used by cross-section and quantum yield constructors
- `data/cross_sections/` — 145+ `.nc` files
- `data/quantum_yields/` — 14 `.nc` files

### Documentation

Document the `DataReader` interface and `NetCDFReader` implementation with Doxygen `///` comments. Document the expected file format, the SI-unit requirement, and how users implement custom readers.
