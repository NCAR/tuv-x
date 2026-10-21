/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "spectrum_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double *copy(tt_builder *b, size_t n, const double *values) {
    double *out = tt_build_alloc(b, n, sizeof(double));
    if (out && n)
        memcpy(out, values, n * sizeof(double));
    return out;
}
static const tt_data_variable *parameters(tt_builder *b, tt_config_node node, const char *path,
                                          const char *name, const tt_data_variable **wave) {
    const tt_dataset *data = tt_build_dataset(b, node, path);
    const tt_data_variable *p = tt_dataset_variable(data, name);
    *wave = tt_dataset_variable(data, "wavelength");
    if (!p || p->rank != 2 || !p->shape[0] || !p->shape[1] || !p->values ||
        p->shape[1] > SIZE_MAX / p->shape[0] || p->count != p->shape[0] * p->shape[1] || !*wave ||
        !(*wave)->values || (*wave)->count != p->shape[1]) {
        tt_build_fail(b, node, TT_INVALID, "Invalid parameterization dataset: %s", path);
        return NULL;
    }
    return p;
}
/* Keys TUV-x accepts for each parameterization (src/cross_sections/util/). */
static const char *const exponential_required[] = {
    "AA", "BB", "lp", "base temperature", "base wavelength", "logarithm", NULL};
static const char *const exponential_optional[] = {"minimum wavelength", "maximum wavelength",
                                                   "temperature ranges",
                                                   "invert temperature offset", NULL};
static const char *const taylor_required[] = {"netcdf file", "base temperature", NULL};
static const char *const taylor_optional[] = {"minimum wavelength", "maximum wavelength",
                                              "temperature ranges", "type", NULL};
static const char *const burkholder_required[] = {"netcdf file", "A", "B", NULL};
static const char *const burkholder_optional[] = {"type", "temperature ranges",
                                                  "minimum wavelength", "maximum wavelength", NULL};
static const char *const harwood_required[] = {
    "aa", "bb", "base temperature", "base wavelength", "logarithm", NULL};
static const char *const harwood_optional[] = {"type",
                                               "minimum wavelength",
                                               "maximum wavelength",
                                               "temperature ranges",
                                               "invert temperature offset",
                                               NULL};
static const char *const range_keys[] = {"minimum", "maximum", "fixed value", NULL};

static void ranges(tt_builder *b, tt_config_node config, tt_temperature_model *p) {
    tt_config_node list = tt_config_get(config, "temperature ranges");
    size_t n = tt_config_size(list);
    /* An explicitly empty list applies no parameterization at positive T. */
    p->range_count = n ? n : 1;
    tt_temperature_range *r = tt_build_alloc(b, p->range_count, sizeof *r);
    p->ranges = r;
    if (!r)
        return;
    if (!n)
        r[0].maximum = tt_build_has(config, "temperature ranges") ? 0 : DBL_MAX;
    for (size_t i = 0; i < n; ++i) {
        tt_config_node node = tt_config_at(list, i);
        tt_build_keys(b, node, NULL, range_keys);
        r[i].minimum = tt_build_number(b, node, "minimum", 0);
        r[i].maximum = tt_build_number(b, node, "maximum", DBL_MAX);
        r[i].fixed_value = tt_build_number(b, node, "fixed value", 0);
        r[i].is_fixed = tt_build_has(node, "fixed value");
        if (r[i].minimum > r[i].maximum)
            tt_build_fail(b, node, TT_RANGE, "Reversed temperature range");
    }
}
void tt_spectrum_temperature_build(tt_builder *b, tt_config_node config, tt_spectrum *s) {
    tt_config_node node = tt_config_get(config, "parameterization");
    const char *type = tt_build_string(b, node, "type", "");
    tt_temperature_model *p = &s->temperature;
    if (!strcmp(type, ""))
        p->kind = TT_TEMPERATURE_EXPONENTIAL;
    else if (!strcmp(type, "TAYLOR_SERIES"))
        p->kind = TT_TEMPERATURE_TAYLOR;
    else if (!strcmp(type, "BURKHOLDER"))
        p->kind = TT_TEMPERATURE_BURKHOLDER;
    else if (!strcmp(type, "HARWOOD"))
        p->kind = TT_TEMPERATURE_HARWOOD;
    else
        tt_build_fail(b, node, TT_INVALID, "Unknown temperature parameterization '%s'", type);
    switch (p->kind) {
    case TT_TEMPERATURE_EXPONENTIAL:
        tt_build_keys(b, node, exponential_required, exponential_optional);
        break;
    case TT_TEMPERATURE_TAYLOR:
        tt_build_keys(b, node, taylor_required, taylor_optional);
        break;
    case TT_TEMPERATURE_BURKHOLDER:
        tt_build_keys(b, node, burkholder_required, burkholder_optional);
        break;
    case TT_TEMPERATURE_HARWOOD:
        tt_build_keys(b, node, harwood_required, harwood_optional);
        break;
    }
    tt_grid custom = {0};
    const tt_grid *grid = &s->wavelength;
    if (tt_build_has(config, "parameterization wavelength grid")) {
        tt_config_node g = tt_config_get(config, "parameterization wavelength grid");
        const char *units = tt_build_string(b, g, "units", NULL);
        if (strcmp(units, "nm"))
            tt_build_fail(b, g, TT_INVALID, "Parameterization grid requires nm units");
        tt_error error = {0};
        tt_status status = tt_grid_from_config(g, &b->source, &custom, &error);
        if (status != TT_OK)
            tt_build_fail(b, g, status, "%s", error.message);
        grid = &custom;
    }
    double lo = tt_build_number(b, node, "minimum wavelength", 0);
    double hi = tt_build_number(b, node, "maximum wavelength", DBL_MAX);
    const double *coordinates = grid->midpoints;
    size_t n = grid->cells, first = 0, end = n;
    const tt_data_variable *v = NULL, *w = NULL;
    if (p->kind == TT_TEMPERATURE_TAYLOR || p->kind == TT_TEMPERATURE_BURKHOLDER) {
        tt_config_node file = tt_config_get(node, "netcdf file");
        const char *path = tt_build_string(b, file, "file path", NULL);
        v = parameters(b, file, path, "temperature_parameters", &w);
        if (v) {
            coordinates = w->values;
            n = end = w->count;
        }
        if (v && v->shape[0] < 2)
            tt_build_fail(b, file, TT_INVALID, "At least two parameter columns are required");
        if (v && p->kind == TT_TEMPERATURE_BURKHOLDER) {
            lo = fmax(lo, coordinates[0]);
            hi = fmin(hi, coordinates[n - 1]);
            if (lo >= hi)
                tt_build_fail(b, node, TT_RANGE, "Empty Burkholder wavelength range");
        }
    }
    if (p->kind != TT_TEMPERATURE_BURKHOLDER) {
        while (first < n && coordinates[first] < lo)
            ++first;
        while (end && coordinates[end - 1] > hi)
            --end;
    }
    if (first >= end || end > n)
        tt_build_fail(b, node, TT_RANGE, "Empty parameterization wavelength range");
    if (b->status != TT_OK) {
        tt_grid_destroy(&custom);
        return;
    }
    p->count = end - first;
    if (p->kind == TT_TEMPERATURE_TAYLOR) {
        p->terms = v->shape[0] - 1;
        p->base_temperature = tt_build_number(b, node, "base temperature", NAN);
        p->sigma = copy(b, p->count, v->values + first);
        double *a = tt_build_alloc(b, p->terms, p->count * sizeof(double));
        p->a = a;
        if (a)
            for (size_t j = 0; j < p->terms; ++j)
                memcpy(a + j * p->count, v->values + (j + 1) * n + first,
                       p->count * sizeof(double));
    } else if (p->kind == TT_TEMPERATURE_BURKHOLDER) {
        p->a = copy(b, n, v->values);
        p->b = copy(b, n, v->values + n);
        p->burkholder_a = tt_build_number(b, node, "A", NAN);
        p->burkholder_b = tt_build_number(b, node, "B", NAN);
        if (p->burkholder_b == 0)
            tt_build_fail(b, node, TT_INVALID, "Burkholder B must be nonzero");
    } else {
        size_t na = 0, nb = 0, np = 0;
        int exp = p->kind == TT_TEMPERATURE_EXPONENTIAL;
        p->a = tt_build_array(b, node, exp ? "AA" : "aa", &na);
        p->b = tt_build_array(b, node, exp ? "BB" : "bb", &nb);
        if (exp)
            p->powers = tt_build_array(b, node, "lp", &np);
        if (!na || na != nb || (exp ? na != np : na != p->count))
            tt_build_fail(b, node, TT_INVALID, "Parameter array sizes do not match");
        p->terms = na;
        p->base_temperature = tt_build_number(b, node, "base temperature", NAN);
        p->base_wavelength = tt_build_number(b, node, "base wavelength", NAN);
        const char *base = tt_build_string(b, node, "logarithm", NULL);
        p->base_10 = !strcmp(base, "base 10");
        if (!p->base_10 && strcmp(base, "natural"))
            tt_build_fail(b, node, TT_INVALID, "Unknown logarithm '%s'", base);
        p->invert_temperature = tt_build_boolean(b, node, "invert temperature offset", 0);
    }
    ranges(b, node, p);
    const tt_data_variable *file_data = NULL, *file_wave = NULL;
    size_t before = 0, after = 0, file_count = 0;
    if (tt_build_has(config, "netcdf file")) {
        const char *path = tt_build_string(b, config, "netcdf file", NULL);
        file_data = parameters(b, config, path, "cross_section_parameters", &file_wave);
        if (file_data && file_data->shape[0] != 1)
            tt_build_fail(b, config, TT_INVALID, "Base spectrum requires one column");
        if (file_wave) {
            file_count = file_wave->count;
            while (before < file_count && file_wave->values[before] < lo)
                ++before;
            after = before;
            while (after < file_count && file_wave->values[after] <= hi)
                ++after;
        }
    }
    s->raw_count = before + p->count + file_count - after;
    s->raw_wavelength = tt_build_alloc(b, s->raw_count, sizeof(double));
    s->raw_data = tt_build_alloc(b, s->raw_count, sizeof(double));
    p->first = before;
    if (b->status == TT_OK) {
        for (size_t i = 0; i < before; ++i) {
            s->raw_wavelength[i] = file_wave->values[i];
            s->raw_data[i] = file_data->values[i];
        }
        memcpy(s->raw_wavelength + before, coordinates + first, p->count * sizeof(double));
        for (size_t i = after; i < file_count; ++i) {
            size_t j = before + p->count + i - after;
            s->raw_wavelength[j] = file_wave->values[i];
            s->raw_data[j] = file_data->values[i];
        }
        for (size_t i = 0; i < s->raw_count; ++i)
            if (!(s->raw_wavelength[i] > 0) || !isfinite(s->raw_wavelength[i]) ||
                (i && !(s->raw_wavelength[i] > s->raw_wavelength[i - 1])))
                tt_build_fail(b, config, TT_RANGE, "Merged wavelength grid is not increasing");
    }
    tt_grid_destroy(&custom);
}
