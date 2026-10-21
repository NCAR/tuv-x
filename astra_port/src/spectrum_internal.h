#ifndef TT_SPECTRUM_INTERNAL_H
#define TT_SPECTRUM_INTERNAL_H
#include "build_internal.h"
typedef struct {
    size_t first, end;
    double value;
} tt_spectral_override;
struct tt_spectrum {
    tt_spectrum_role role;
    int kind;
    tt_grid wavelength;
    size_t table_count;
    tt_spectral_table *tables;
    size_t override_count;
    tt_spectral_override *overrides;
    tt_yield_options yield;
    double weight_a, weight_b;
    int parameterized;
    tt_temperature_model temperature;
    size_t raw_count;
    double *raw_wavelength, *raw_data;
};
void tt_spectrum_temperature_build(tt_builder *b, tt_config_node config, tt_spectrum *s);
#endif
