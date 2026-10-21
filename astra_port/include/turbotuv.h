/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef TURBOTUV_H
#define TURBOTUV_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* All physical quantities use the TUV-x units: nm, km, K, molecules/cm^3.
 * Arrays are contiguous, owned by the caller unless explicitly documented.
 * No routine terminates the process or holds mutable global state. */
typedef enum { TT_OK = 0, TT_INVALID, TT_RANGE, TT_MEMORY, TT_SINGULAR, TT_NOT_FOUND } tt_status;
const char *tt_status_string(tt_status status);

typedef enum {
    TT_LINEAR,
    TT_CONSERVING,
    TT_FRACTIONAL_SOURCE,
    TT_FRACTIONAL_TARGET
} tt_interpolation;

/* nx, nt count coordinates; all grids must be strictly increasing.
 * Linear interpolation writes nt values, the other methods nt-1 bin values.
 * Linear interpolation returns zero outside [x[0], x[nx-1]), including the
 * last source point, matching TUV-x. Conserving requires full source coverage.
 * fold applies the reference's historical overhang convention unchanged. */
tt_status tt_interpolate(tt_interpolation method, int fold, size_t nx, const double *x,
                         const double *y, size_t nt, const double *target, double *out);

/* Lower/upper diagonals have n entries; lower[0], upper[n-1] are unused.
 * work has n entries. Inputs and output/work must not overlap. */
tt_status tt_tridiagonal(size_t n, const double *lower, const double *diagonal, const double *upper,
                         const double *rhs, double *solution, double *work);

typedef enum { TT_MATRIX_DENSE, TT_MATRIX_BAND } tt_matrix_layout;
/* Column-major matrices. Dense A(i,j) is values[j*leading_dimension+i].
 * Band A(i,j) is values[j*leading_dimension+lower+upper+i-j], requiring
 * leading_dimension >= 2*lower+upper+1, including fill space. Lower/upper
 * are ignored for dense matrices. All indices and pivots are zero-based. */
typedef struct {
    tt_matrix_layout layout;
    size_t order, leading_dimension, lower, upper;
    double *values;
} tt_matrix;
/* Factor in place. If reciprocal_condition is non-NULL, work must have order
 * entries; the estimate follows LINPACK's original scaling and pivoting.
 * A singular matrix returns TT_SINGULAR, retaining factors and the estimate. */
tt_status tt_lu_factor(tt_matrix *matrix, size_t *pivots, double *reciprocal_condition,
                       double *work);
/* Solve using factors, replacing rhs with the solution. transpose != 0 solves
 * A^T x = rhs. Repeated solves do not change the factors or pivots. */
tt_status tt_lu_solve(const tt_matrix *factors, const size_t *pivots, int transpose, double *rhs);

/* The catalogue contains immutable, compiled-in copies of all upstream data.
 * Shape and storage follow NetCDF C order (last dimension contiguous).
 * Text entries are NUL terminated. Lookup accepts upstream relative paths. */
typedef struct {
    const char *name;
    size_t rank;
    size_t shape[4];
    size_t count;
    const double *values;
    const char *attributes_json;
} tt_data_variable;
typedef struct {
    const char *path;
    const char *sha256;
    size_t variable_count;
    const tt_data_variable *variables;
    const char *text;
    size_t text_length;
} tt_dataset;
size_t tt_dataset_count(void);
const tt_dataset *tt_dataset_at(size_t index);
const tt_dataset *tt_dataset_find(const char *path);
const tt_data_variable *tt_dataset_variable(const tt_dataset *data, const char *name);

/* Grid/profile constructors own their arrays. Initialize outputs to {0},
 * destroy before reuse, and do not free individual fields. A failed
 * constructor leaves the output unchanged. One edge (zero cells) is valid. */
typedef struct {
    size_t cells;
    double *edges, *midpoints, *widths;
} tt_grid;
tt_status tt_grid_from_edges(size_t count, const double *edges, tt_grid *out);
tt_status tt_grid_equal(double first, double last, double interval, tt_grid *out);
/* Embedded grid text has a header record followed by one coordinate per line. */
tt_status tt_grid_from_dataset(const tt_dataset *data, tt_grid *out);
/* Update existing storage; NULL midpoints computes the edge averages. */
tt_status tt_grid_update(tt_grid *grid, const double *edges, const double *midpoints);
void tt_grid_destroy(tt_grid *grid);

typedef struct {
    size_t cells;
    double *edges, *midpoints, *differences;
    double *columns, *extended_columns, *burden;
} tt_profile;
typedef enum { TT_PROFILE_LINEAR, TT_PROFILE_AIR, TT_PROFILE_O2, TT_PROFILE_O3 } tt_profile_kind;
typedef struct {
    tt_profile_kind kind;
    double scale_height, reference_column, length_conversion;
} tt_profile_options;
/* Default scale heights: air/O2 8.01 km, O3 4.5 km, linear 0. O3's default
 * reference column is 300 DU; 1 disables rescaling (the TUV-x convention).
 * length_conversion defaults to 1 for linear and 1e5 for the gas profiles;
 * set 1e5 explicitly for linear profiles on height grids in km. */
tt_profile_options tt_profile_options_default(tt_profile_kind kind);
/* Edges and extended_columns have cells+1 entries, other fields have cells.
 * columns includes the extra top column in its last cell; extended_columns
 * keeps that extra column separate. burden is the top-down column sum. */
tt_status tt_profile_from_edges(const tt_grid *grid, const double *values,
                                const tt_profile_options *options, tt_profile *out);
/* Air/O2 samples contain air densities: interpolate log(density), extend the
 * final coordinate by .001 km, then multiply O2 by .2095. O3 extends upward
 * in 1-km steps using scale_height before interpolating and rescaling. */
tt_status tt_profile_from_samples(const tt_grid *grid, size_t count, const double *coordinates,
                                  const double *values, const tt_profile_options *options,
                                  tt_profile *out);
tt_status tt_profile_from_dataset(const tt_grid *grid, const tt_dataset *data,
                                  const tt_profile_options *options, tt_profile *out);
void tt_profile_destroy(tt_profile *profile);
/* Create zero-filled host storage. update requires edges; other arguments
 * are optional. NULL columns preserves previous columns/burdens. A supplied
 * scale_height requires columns and excludes top_column; it uses km to cm.
 * Explicit top_column is added to the final column and stored separately. */
tt_status tt_profile_create(size_t cells, tt_profile *out);
tt_status tt_profile_update(tt_profile *profile, const double *edges, const double *midpoints,
                            const double *columns, const double *scale_height,
                            const double *top_column);
/* Both outputs are allocated. Preserve TUV-x's time-zone conventions:
 * zenith uses time+time_zone; Earth-Sun distance uses time-time_zone. */
tt_status tt_solar_profiles(const tt_grid *time, int year, int month, int day, double time_zone,
                            double latitude, double longitude, tt_profile *zenith,
                            tt_profile *distance);
/* Merge up to four embedded solar spectra on the historical boundaries
 * 150.01, 200.07, 1000.99 nm. All inputs are compiled text tables. The result
 * has grid->cells entries in photons cm^-2 s^-1, integrated over each bin. */
tt_status tt_extraterrestrial_flux(const tt_grid *grid, size_t count, const tt_dataset *const *data,
                                   const tt_interpolation *methods, double *out);

/* Geometry uses increasing height edges. nid has n entries. dsdh is an n by
 * (n-1) row-major array, with rows ordered from top to bottom of atmosphere.
 * This is the ray-path convention used by both TUV-x transfer solvers. */
tt_status tt_spherical_geometry(size_t n, const double *height, double zenith, int *nid,
                                double *dsdh);
tt_status tt_air_mass(size_t n, const double *air_column, const int *nid, const double *dsdh,
                      double *vertical, double *slant);

/* Calendar follows the reference's supported range, 1950 through 2050.
 * Solar position returns elevation in degrees; zenith = 90 - elevation. */
int tt_day_of_year(int year, int month, int day);
double tt_earth_sun_distance(int year, int day, double utc_hour);
double tt_solar_elevation(int year, int day, double utc_hour, double latitude, double longitude);

/* Spectral fields use [wavelength][layer/level] order, with height increasing.
 * Moments use [moment][wavelength][layer]; moment 0 is the first Legendre
 * coefficient. All input/output arrays are caller-owned. */
typedef struct {
    const double *optical_depth;
    const double *single_scattering_albedo;
    const double *asymmetry;
    int is_air;
} tt_radiator;
typedef struct {
    size_t layers, wavelengths, moments;
    const double *height;
    const double *optical_depth;
    const double *single_scattering_albedo;
    const double *legendre;
} tt_atmosphere;
typedef struct {
    double *direct_flux, *down_flux, *up_flux;
    double *direct_irradiance, *down_irradiance, *up_irradiance;
} tt_radiation;

/* Gas cross sections are evaluated at layer midpoints, [wavelength][layer].
 * All three outputs have layers*wavelengths entries. Air scatters; other
 * gases absorb, and both have zero asymmetry before moment mixing. */
tt_status tt_gas_optics(size_t layers, size_t wavelengths, const double *columns,
                        const double *cross_section, int is_air, double *optical_depth,
                        double *albedo, double *asymmetry);
/* Aerosol depth samples use the original 1-km grid beginning at zero.
 * Interpolation folds material above the model top into the last layer.
 * tau550 < -1e-35 disables column normalization; alpha defaults to 1 in TUV-x.
 * A single depth sample defines a uniform layer depth in C (undefined in the
 * reference constructor). Output ordering is [wavelength][layer]. */
tt_status tt_aerosol_optics(const tt_grid *height, const tt_grid *wavelength, size_t count,
                            const double *depth_samples, double albedo, double asymmetry,
                            double tau550, double alpha, double *optical_depth,
                            double *single_scattering_albedo, double *asymmetry_factor);
/* Return an immutable view of embedded optical_depth, single_scattering_albedo,
 * and asymmetry_factor variables; the dataset must outlive the returned view. */
tt_status tt_radiator_from_dataset(const tt_dataset *data, size_t layers, size_t wavelengths,
                                   tt_radiator *out);

tt_status tt_mix_radiators(size_t count, const tt_radiator *radiators, size_t layers,
                           size_t wavelengths, size_t moments, double *optical_depth,
                           double *albedo, double *legendre);
tt_status tt_delta_eddington(const tt_atmosphere *atmosphere, double zenith,
                             const double *surface_albedo, tt_radiation *out);
/* Even stream counts 2 through 32, with at least that many Legendre moments. */
tt_status tt_discrete_ordinates(const tt_atmosphere *atmosphere, size_t streams, double zenith,
                                const double *surface_albedo, tt_radiation *out);

typedef struct {
    size_t rows, columns, temperature_count;
    const double *values; /* [column][row] */
    const double *temperatures;
} tt_spectral_table;
typedef struct {
    size_t wavelengths, levels;
    const double *edges;       /* wavelengths + 1 */
    const double *midpoints;   /* wavelengths */
    const double *temperature; /* levels */
    const double *air_density; /* levels, required by quantum-yield calculations */
} tt_spectral_grid;
typedef enum {
    TT_CS_BASE,
    TT_CS_ACETONE,
    TT_CS_CCL4,
    TT_CS_CFC11,
    TT_CS_CH2O,
    TT_CS_CH3ONO2,
    TT_CS_CHBR3,
    TT_CS_CHCL3,
    TT_CS_CL2,
    TT_CS_CLONO2,
    TT_CS_H2O2,
    TT_CS_HCFC,
    TT_CS_HNO3,
    TT_CS_HOBR,
    TT_CS_N2O,
    TT_CS_N2O5,
    TT_CS_NITROXY_ACETONE,
    TT_CS_NITROXY_ETHANOL,
    TT_CS_NO2,
    TT_CS_O3_JPL06,
    TT_CS_O3,
    TT_CS_OCLO,
    TT_CS_RAYLEIGH,
    TT_CS_RONO2,
    TT_CS_T_BUTYL_NITRATE,
    TT_CS_TINT,
    TT_CS_COUNT
} tt_cross_section_kind;
const char *tt_cross_section_name(tt_cross_section_kind kind);
tt_status tt_cross_section(tt_cross_section_kind kind, const tt_spectral_grid *grid, size_t count,
                           const tt_spectral_table *tables, double *out);
double tt_air_refractive_index(double wavelength, double air_density);

typedef enum {
    TT_EXTRAPOLATE_ZERO,
    TT_EXTRAPOLATE_BOUNDARY,
    TT_EXTRAPOLATE_CONSTANT
} tt_extrapolation;
/* Prepare one parameter column using the reference's four padding points.
 * out has wavelengths entries; scientific data come from the embedded catalogue. */
tt_status tt_regrid_column(const tt_dataset *data, const char *variable, size_t column,
                           size_t wavelengths, const double *edges, tt_extrapolation lower,
                           double lower_value, tt_extrapolation upper, double upper_value,
                           double *out);

typedef enum {
    TT_TEMPERATURE_EXPONENTIAL,
    TT_TEMPERATURE_TAYLOR,
    TT_TEMPERATURE_BURKHOLDER,
    TT_TEMPERATURE_HARWOOD
} tt_temperature_kind;
typedef struct {
    double minimum, maximum, fixed_value;
    int is_fixed;
} tt_temperature_range;
/* The parameterized slice is [first, first + count) in the raw wavelength
 * grid. EXPONENTIAL uses terms entries of a,b,powers; the other kinds use
 * count entries of a,b (TAYLOR uses a[term][wavelength] and sigma[wavelength]).
 * Ranges have inclusive endpoints. Overlapping exponential ranges add;
 * the other kinds replace with the last matching range. No matching range
 * leaves the supplied cross section unchanged. A zero range_count uses T. */
typedef struct {
    tt_temperature_kind kind;
    size_t first, count, terms;
    const double *a, *b, *powers, *sigma;
    double base_temperature, base_wavelength;
    int base_10, invert_temperature;
    double burkholder_a, burkholder_b;
    size_t range_count;
    const tt_temperature_range *ranges;
} tt_temperature_model;
tt_status tt_temperature_apply(const tt_temperature_model *model, double temperature, size_t count,
                               const double *wavelengths, double *cross_section);
/* Apply at each level, pad the raw spectrum with zeroes, then conserve its
 * integral on the target grid. Output is [wavelength][level]. */
tt_status tt_temperature_cross_section(const tt_temperature_model *model, size_t count,
                                       const double *wavelengths, const double *cross_section,
                                       const tt_spectral_grid *grid, double *out);

typedef enum {
    TT_QY_BASE,
    TT_QY_TINT,
    TT_QY_NO2,
    TT_QY_O3_O1D,
    TT_QY_O3_O3P,
    TT_QY_HO2,
    TT_QY_NO3_AQ,
    TT_QY_CH2O,
    TT_QY_CH3CHO,
    TT_QY_C2H5CHO,
    TT_QY_CH2CHCHO,
    TT_QY_MVK,
    TT_QY_ACETONE,
    TT_QY_MEK,
    TT_QY_CH3COCHO,
    TT_QY_CLO_O1D,
    TT_QY_CLO_O3P,
    TT_QY_CLONO2_CL,
    TT_QY_CLONO2_CLO,
    TT_QY_H2SO4,
    TT_QY_TAYLOR,
    TT_QY_COUNT
} tt_quantum_yield_kind;
typedef struct {
    double low_wavelength_value, high_wavelength_value;
    double minimum_temperature, maximum_temperature;
    unsigned acetone_branches; /* 1 = CO, 2 = CH3CO, 3 = both */
    double minimum_wavelength, maximum_wavelength;
    size_t coefficient_count;
    const double *coefficients;                  /* Taylor coefficients, ascending powers */
    double molecular_diameter, molecular_weight; /* m, kg/mol */
    const double *collision_rates; /* wavelengths entries, s^-1; negative = use base table */
} tt_yield_options;
tt_yield_options tt_yield_options_default(void);
const char *tt_quantum_yield_name(tt_quantum_yield_kind kind);
tt_status tt_quantum_yield(tt_quantum_yield_kind kind, const tt_spectral_grid *grid, size_t count,
                           const tt_spectral_table *tables, const tt_yield_options *options,
                           double *out);

typedef enum {
    TT_WEIGHT_BASE,
    TT_WEIGHT_NOTCH,
    TT_WEIGHT_GAUSSIAN,
    TT_WEIGHT_EPPLEY,
    TT_WEIGHT_PAR,
    TT_WEIGHT_EXPONENTIAL,
    TT_WEIGHT_SCUP,
    TT_WEIGHT_ERYTHEMA,
    TT_WEIGHT_UV_INDEX,
    TT_WEIGHT_PHYTOPLANKTON,
    TT_WEIGHT_PLANT,
    TT_WEIGHT_PLANT_FLINT,
    TT_WEIGHT_PLANT_FLINT_EXT,
    TT_WEIGHT_COUNT
} tt_weight_kind;
/* a,b are notch boundaries; a is the Gaussian centroid. base is required
 * only for BASE/EPPLEY. Output has wavelengths entries. */
tt_status tt_spectral_weight(tt_weight_kind kind, size_t wavelengths, const double *edges,
                             const double *midpoints, const double *base, double a, double b,
                             double *out);

/* Convert normalized radiation to photon cm^-2 s^-1 or W m^-2 per band.
 * solar_flux is already integrated over each band, as in TUV-x. Negative
 * totals are clamped to zero. Fields and spectral outputs are [wavelength][level]. */
tt_status tt_actinic_flux(size_t wavelengths, size_t levels, const tt_radiation *radiation,
                          const double *solar_flux, double *out);
tt_status tt_spectral_irradiance(size_t wavelengths, size_t levels, const tt_radiation *radiation,
                                 const double *solar_flux, const double *midpoints, double *out);
/* Integrate one reaction/action spectrum. Rate outputs have levels entries.
 * Cross sections and yields are [wavelength][level]. Photolysis: s^-1;
 * heating: J s^-1 per molecule; dose: weighted W m^-2. The energy threshold
 * for heating is specified as its equivalent wavelength in nm. */
tt_status tt_photolysis_rate(size_t wavelengths, size_t levels, const double *actinic_flux,
                             const double *cross_section, const double *quantum_yield, double scale,
                             double *out);
tt_status tt_heating_rate(size_t wavelengths, size_t levels, const double *actinic_flux,
                          const double *cross_section, const double *quantum_yield,
                          const double *midpoints, double energy_wavelength, double scale,
                          double *out);
tt_status tt_dose_rate(size_t wavelengths, size_t levels, const double *irradiance,
                       const double *weights, double *out);

/* O2 corrections operate only on the prescribed Lyman-alpha and 17
 * Schumann-Runge bins. Absent bands use the index (size_t)-1. The detector
 * rejects incompatible edges when a grid spans either complete band. */
typedef struct {
    size_t lyman_alpha, schumann_runge;
} tt_oxygen_bands;
tt_status tt_oxygen_band_grid(size_t wavelengths, const double *edges, tt_oxygen_bands *bands);
/* Column arrays are in molecules/cm^2, ordered with increasing height.
 * vertical_air_column has levels-1 entries; slant columns have levels.
 * Pass scaled air columns or columns integrated from an explicit O2 profile.
 * Outputs are updated only in the selected bands: optical_depth has shape
 * [wavelength][levels-1], cross_section has [wavelength][levels]. Either
 * output may be NULL. Coefficients are compiled into the library. */
tt_status tt_oxygen_optics(const tt_oxygen_bands *bands, size_t wavelengths, size_t levels,
                           const double *temperature, const double *vertical_air_column,
                           const double *slant_air_column, const double *slant_oxygen_column,
                           double *optical_depth, double *cross_section);

#ifdef __cplusplus
}
#endif
#endif
