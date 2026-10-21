#ifndef TT_MODEL_INTERNAL_H
#define TT_MODEL_INTERNAL_H
#include "build_internal.h"

typedef struct {
    char *name, *units;
    tt_grid value;
} tt_model_grid_entry;
typedef struct {
    char *name, *units;
    tt_profile value;
} tt_model_profile_entry;
typedef struct {
    char *name;
    tt_spectrum *spectrum;
} tt_model_spectrum_entry;
typedef struct {
    char *name;
    tt_radiator view;
    double *data;
    int gas;
    size_t profile, spectrum;
} tt_model_radiator_entry;
typedef struct {
    char *name;
    tt_spectrum *cross_section, *yield;
    int oxygen;
    size_t heating_index;
    double scale, energy_wavelength;
} tt_model_reaction;
typedef struct {
    char *name;
    double *weights;
} tt_model_dose;
struct tt_model {
    tt_model_dimensions shape;
    size_t grid_count, profile_count, spectrum_count, radiator_count;
    tt_model_grid_entry *grids;
    tt_model_profile_entry *profiles;
    tt_model_spectrum_entry *spectra;
    tt_model_radiator_entry *radiators;
    tt_radiator *radiator_views;
    tt_model_reaction *reactions;
    tt_model_dose *doses;
    size_t height, wavelength, temperature, air, oxygen, surface, solar, oxygen_radiator;
    size_t streams;
    int host_solver, scaled_oxygen;
    double oxygen_scale;
    tt_oxygen_bands bands;
    tt_atmosphere atmosphere;
    tt_model_result result;
    int *nid;
    double *path, *vertical, *slant_air, *slant_oxygen, *oxygen_vertical;
    double *work, *cross_section, *yield, *radiation_storage, *host_radiation;
    double *optical_storage, *rate_storage;
    struct tt_model_plan *plan; /* prepared runs; NULL if unavailable */
    int generic;                /* TT_EXECUTION_GENERIC selected */
};
size_t tt_model_grid_index(const tt_model *model, const char *name, const char *units);
size_t tt_model_profile_index(const tt_model *model, const char *name, const char *units);
tt_status tt_model_copy_profile(const tt_profile *source, tt_profile *out);
tt_status tt_model_copy_radiator(const tt_named_radiator *source, tt_model_radiator_entry *out,
                                 size_t layers, size_t wavelengths);
#endif
