/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef TURBOTUV_MODEL_H
#define TURBOTUV_MODEL_H
#include "turbotuv_config.h"
#ifdef __cplusplus
extern "C" {
#endif

/* A resolver supplies immutable compiled datasets. NULL uses the built-in
 * catalogue. Scientific constructors never open data files. */
typedef const tt_dataset *(*tt_data_resolver)(void *context, const char *path);
typedef struct {
    tt_data_resolver resolve;
    void *context;
} tt_data_source;
tt_status tt_grid_from_config(tt_config_node config, const tt_data_source *source, tt_grid *out,
                              tt_error *error);

/* Named grid views let profiles choose their configured axis without a
 * mutable global warehouse. Names and units are matched exactly.
 * Extraterrestrial-flux profiles populate midpoints; their other fields are
 * defined as zero. All other profile constructors use the usual edge fields. */
typedef struct {
    const char *name, *units;
    const tt_grid *grid;
} tt_named_grid;
tt_status tt_profile_from_config(tt_config_node config, size_t grid_count,
                                 const tt_named_grid *grids, const tt_data_source *source,
                                 tt_profile *out, tt_error *error);

typedef enum {
    TT_SPECTRUM_CROSS_SECTION,
    TT_SPECTRUM_QUANTUM_YIELD,
    TT_SPECTRUM_WEIGHT
} tt_spectrum_role;
typedef struct tt_spectrum tt_spectrum;
/* Prepare tables once on a fixed wavelength grid, using upstream JSON/YAML
 * configuration keys. The result owns its grid, tables and options; neither
 * config nor source must remain alive. Output is unchanged on failure. */
tt_status tt_spectrum_create(tt_spectrum_role role, tt_config_node config,
                             const tt_grid *wavelength, const tt_data_source *source,
                             tt_spectrum **out, tt_error *error);
void tt_spectrum_destroy(tt_spectrum *spectrum);
/* Cross sections/yields write [wavelength][level]. Supply temperature and air
 * arrays for these roles. Weights write wavelength entries and ignore levels
 * and both profile arguments. Prepared objects support concurrent readers. */
tt_status tt_spectrum_calculate(const tt_spectrum *spectrum, size_t levels,
                                const double *temperature, const double *air, double *out);
/* Immutable table views, valid until destroy; useful for inspecting the
 * wavelength interpolation. Parameterized cross sections have no tables. */
size_t tt_spectrum_table_count(const tt_spectrum *spectrum);
const tt_spectral_table *tt_spectrum_table(const tt_spectrum *spectrum, size_t index);
/* A complete model owns its configuration-derived state and reusable working
 * arrays. Separate models may run concurrently; one model has one active run.
 * Views below remain valid until destroy and must not be modified by callers. */
typedef struct tt_model tt_model;
typedef enum { TT_RATE_PHOTOLYSIS, TT_RATE_HEATING, TT_RATE_DOSE } tt_rate_kind;
typedef struct {
    size_t levels, wavelengths, photolysis, heating, dose;
} tt_model_dimensions;
/* Actinic flux and spectral irradiance need an "extraterrestrial flux"
 * profile; without one (and so without rates) both are zero. */
typedef struct {
    tt_radiation radiation;
    const double *actinic_flux, *spectral_irradiance; /* [wavelength][level] */
    const double *photolysis, *heating, *dose;        /* [rate][level] */
} tt_model_result;
typedef struct {
    const char *name, *units;
    const tt_profile *profile;
} tt_named_profile;
typedef struct {
    const char *name;
    size_t layers, wavelengths;
    tt_radiator radiator;
} tt_named_radiator;
typedef struct {
    size_t grid_count;
    const tt_named_grid *grids;
    size_t profile_count;
    const tt_named_profile *profiles;
    size_t radiator_count;
    const tt_named_radiator *radiators;
} tt_model_inputs;
tt_status tt_model_create(tt_config_node config, const tt_data_source *source, tt_model **out,
                          tt_error *error);
/* Add host-owned inputs before resolving dependent components. Names must
 * be absent from the corresponding configuration collection. All data are
 * copied; input arrays may be freed after construction. */
tt_status tt_model_create_with_inputs(tt_config_node config, const tt_data_source *source,
                                      const tt_model_inputs *inputs, tt_model **out,
                                      tt_error *error);
void tt_model_destroy(tt_model *model);
tt_model_dimensions tt_model_shape(const tt_model *model);
const tt_grid *tt_model_grid(const tt_model *model, const char *name, const char *units);
const tt_profile *tt_model_profile(const tt_model *model, const char *name, const char *units);
const tt_radiator *tt_model_radiator(const tt_model *model, const char *name);
const char *tt_model_rate_name(const tt_model *model, tt_rate_kind kind, size_t index);
/* How runs are evaluated. PREPARED, the default, uses work done once at
 * construction: spectral values that cannot change between columns, spectra
 * shared by equal configurations, rates summed over precomputed products, and
 * Delta-Eddington transfer for blocks of wavelengths in lock step. GENERIC
 * calls the public kernels in sequence, as a direct translation of TUV-x. The
 * results are bit-identical; a prepared run falls back to the generic
 * sequence for inputs it does not cover, such as invalid values, so errors
 * are the same too. Returns TT_NOT_FOUND if PREPARED is selected but could
 * not be prepared (runs are then generic). */
typedef enum { TT_EXECUTION_PREPARED, TT_EXECUTION_GENERIC } tt_execution;
tt_status tt_model_set_execution(tt_model *model, tt_execution execution);
/* The reference scales all six fields by this argument directly. Preserve
 * that convention: pass the configured Earth-Sun distance, or 1 for unit
 * scaling. Result storage is replaced by the next run on the same model. */
tt_status tt_model_run(tt_model *model, double zenith, double radiation_scale,
                       const tt_model_result **out, tt_error *error);
/* Host solver input: three flux fields are mandatory. Missing irradiance
 * fields become zero. Shape is model wavelengths * levels. Inputs are copied. */
tt_status tt_model_set_radiation(tt_model *model, const tt_radiation *radiation);
/* Updates preserve dimensions and borrowed view addresses. Height updates
 * do not recompute profiles or static radiators: update those explicitly.
 * Wavelength grids are fixed after construction because spectra are cached.
 * Other grids accept cells+1 edges and optional cells midpoint values. */
tt_status tt_model_update_grid(tt_model *model, const char *name, const char *units, size_t cells,
                               const double *edges, const double *midpoints);
/* Same optional-field semantics as tt_profile_update. */
tt_status tt_model_update_profile(tt_model *model, const char *name, const char *units,
                                  size_t cells, const double *edges, const double *midpoints,
                                  const double *columns, const double *scale_height,
                                  const double *top_column);
/* A supplied radiator replaces its calculated optics for subsequent runs.
 * Optical depth is required; NULL SSA/asymmetry arrays become zero. */
tt_status tt_model_set_radiator(tt_model *model, const tt_named_radiator *radiator);

#ifdef __cplusplus
}
#endif
#endif
