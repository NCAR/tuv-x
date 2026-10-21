/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "spectrum_internal.h"
#include "spectral_validation.h"
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const char *const weight_names[] = {"base",
                                           "Notch Filter",
                                           "Gaussian",
                                           "Eppley UV Photometer",
                                           "PAR, 400-700 nm, umol m-2 s-1",
                                           "Exponential decay, 14 nm/10",
                                           "SCUP-mice (de Gruijl et al., 1993)",
                                           "Standard human erythema (Webb et al., 2011)",
                                           "UV index (WMO, 1994; Webb et al., 2011)",
                                           "Phytoplankton (Boucher et al., 1994)",
                                           "Plant damage (Caldwell, 1971)",
                                           "Plant damage,Flint&Caldwell,2003,orig.",
                                           "Plant damage,Flint&Caldwell,2003,ext390"};

/* Configuration keys TUV-x accepts: src/cross_section.F90, src/cross_sections/,
 * src/quantum_yield.F90, src/quantum_yields/, src/spectral_weight.F90 and
 * src/spectral_weights/. Kinds TUV-x does not validate are not validated. */
static const char *const type_key[] = {"type", NULL}, *const name_key[] = {"name", NULL};
static const char *const type_files[] = {"type", "netcdf files", NULL};
static const char *const extrapolated[] = {"lower extrapolation", "upper extrapolation", "name",
                                           NULL};
static const char *const base_cross_section[] = {
    "netcdf files", "name", "merge data", "override bands", "apply O2 bands", "data", NULL};
static const char *const jpl_files[] = {"type", "base netcdf file", "low-temperature netcdf file",
                                        "high-temperature netcdf file", NULL};
static const char *const parameterized_required[] = {"type", "parameterization", NULL};
static const char *const parameterized_optional[] = {"name", "parameterization wavelength grid",
                                                     "netcdf file", NULL};
static const char *const base_yield[] = {"netcdf files",
                                         "lower extrapolation",
                                         "upper extrapolation",
                                         "name",
                                         "constant value",
                                         "override bands",
                                         NULL};
static const char *const acetone_yield[] = {"name",
                                            "low wavelength value",
                                            "high wavelength value",
                                            "minimum temperature",
                                            "maximum temperature",
                                            "branch",
                                            NULL};
static const char *const h2so4_required[] = {"type",
                                             "parameterized wavelengths",
                                             "collision interval s",
                                             "molecular diameter m",
                                             "molecular weight kg mol-1",
                                             NULL};
static const char *const taylor_required[] = {"coefficients", NULL};
static const char *const taylor_optional[] = {
    "type",           "netcdf files",   "lower extrapolation", "upper extrapolation", "name",
    "constant value", "override bands", "minimum wavelength",  "maximum wavelength",  NULL};
static const char *const base_weight[] = {"netcdf files", "lower extrapolation",
                                          "upper extrapolation", "name", NULL};
static const char *const notch_required[] = {"type", "notch filter begin", "notch filter end",
                                             NULL};
static const char *const gaussian_required[] = {"type", "centroid", NULL};
static const char *const file_required[] = {"file path", NULL};
static const char *const file_optional[] = {
    "interpolator", "lower extrapolation", "upper extrapolation", "zero below", "zero above", NULL};
static const char *const value_key[] = {"value", NULL};
static const char *const point_data[] = {"default value", "point values", NULL};
static const char *const band_required[] = {"band", "value", NULL};
static const char *const band_range[] = {"minimum wavelength", "maximum wavelength", NULL};

static void spectrum_keys(tt_builder *b, tt_config_node config, const tt_spectrum *s) {
    const char *const *required = NULL, *const *optional = NULL;
    if (s->parameterized) {
        required = parameterized_required;
        optional = parameterized_optional;
    } else if (s->role == TT_SPECTRUM_CROSS_SECTION) {
        required = type_files;
        optional = extrapolated;
        switch (s->kind) {
        case TT_CS_BASE:
            required = type_key;
            optional = base_cross_section;
            break;
        case TT_CS_RAYLEIGH:
        case TT_CS_CL2:
        case TT_CS_HOBR:
        case TT_CS_N2O:
        case TT_CS_NITROXY_ACETONE:
        case TT_CS_NITROXY_ETHANOL:
        case TT_CS_T_BUTYL_NITRATE:
            required = type_key;
            optional = name_key;
            break;
        case TT_CS_HCFC:
            optional = name_key;
            break;
        case TT_CS_O3_JPL06:
            required = jpl_files;
            optional = name_key;
            break;
        default:
            break;
        }
    } else if (s->role == TT_SPECTRUM_QUANTUM_YIELD) {
        if (s->kind == TT_QY_BASE) {
            required = type_key;
            optional = base_yield;
        } else if (s->kind == TT_QY_ACETONE) {
            required = type_key;
            optional = acetone_yield;
        } else if (s->kind == TT_QY_H2SO4) {
            required = h2so4_required;
            optional = base_yield;
        } else if (s->kind == TT_QY_TAYLOR) {
            required = taylor_required;
            optional = taylor_optional;
        }
    } else if (s->kind != TT_WEIGHT_EPPLEY) {
        required = s->kind == TT_WEIGHT_NOTCH      ? notch_required
                   : s->kind == TT_WEIGHT_GAUSSIAN ? gaussian_required
                                                   : type_key;
        optional = s->kind == TT_WEIGHT_BASE ? base_weight : name_key;
    }
    if (required)
        tt_build_keys(b, config, required, optional);
}

static double extrapolation(tt_builder *b, tt_config_node config, const char *key,
                            double boundary) {
    tt_config_node node = tt_config_get(config, key);
    if (tt_config_type(node) == TT_CONFIG_MISSING)
        return 0;
    tt_build_keys(b, node, type_key, value_key);
    const char *type = tt_build_string(b, node, "type", NULL);
    if (!strcmp(type, "boundary"))
        return boundary;
    if (!strcmp(type, "constant"))
        return tt_build_number(b, node, "value", NAN);
    tt_build_fail(b, node, TT_INVALID, "Unknown extrapolation '%s'", type);
    return 0;
}

/* Base cross sections, temperature-interpolated spectra and quantum yields
 * have distinct configuration scopes and padding conventions in TUV-x. */
static void prepare_table(tt_builder *b, tt_spectrum *s, tt_config_node config, tt_config_node file,
                          tt_spectral_table *table) {
    int cs = s->role == TT_SPECTRUM_CROSS_SECTION;
    int qy = s->role == TT_SPECTRUM_QUANTUM_YIELD;
    int special = cs && (s->kind == TT_CS_HNO3 || s->kind == TT_CS_RONO2 || s->kind == TT_CS_TINT ||
                         s->kind == TT_CS_NO2 || s->kind == TT_CS_O3);
    int sort = cs ? (s->kind == TT_CS_TINT || s->kind == TT_CS_NO2 || s->kind == TT_CS_O3)
                  : qy && (s->kind == TT_QY_TINT || s->kind == TT_QY_NO2);
    /* TUV-x validates the file entries it reads with the base process_file:
     * the base-constructor kinds and O3 JPL06. The other special kinds read
     * their files themselves. */
    if (cs && !special)
        tt_build_keys(b, file, file_required, file_optional);
    const char *path = cs ? tt_build_string(b, file, "file path", NULL) : tt_config_string(file);
    if (!path) {
        tt_build_fail(b, file, TT_INVALID, "Expected dataset path");
        return;
    }
    const tt_dataset *data = tt_build_dataset(b, file, path);
    if (!data)
        return;
    const tt_data_variable *values = tt_dataset_variable(data, cs   ? "cross_section_parameters"
                                                               : qy ? "quantum_yield_parameters"
                                                                    : "spectral_weight_parameters");
    const tt_data_variable *wave = tt_dataset_variable(data, "wavelength");
    const tt_data_variable *temp = tt_dataset_variable(data, "temperature");
    if (!values || values->rank != 2 || !values->shape[0] || !values->shape[1] || !values->values ||
        values->shape[1] > SIZE_MAX / sizeof(double) / values->shape[0] ||
        values->count != values->shape[0] * values->shape[1] ||
        (wave && (wave->count != values->shape[1] || !wave->values)) || (temp && !temp->values)) {
        tt_build_fail(b, file, TT_INVALID, "Invalid spectral dataset: %s", path);
        return;
    }
    size_t nr = values->shape[1], nc = values->shape[0], nw = s->wavelength.cells;
    size_t rows = wave ? nw : nr;
    if (nr > SIZE_MAX / sizeof(double) / 2 - 4 || rows > SIZE_MAX / sizeof(double)) {
        tt_build_fail(b, file, TT_MEMORY, "Spectral table is too large");
        return;
    }
    int reverse = 0;
    if (sort && (!temp || temp->count != nc)) {
        tt_build_fail(b, file, TT_INVALID, "Temperature/parameter count mismatch: %s", path);
        return;
    }
    if (sort && nc > 1) {
        reverse = temp->values[1] < temp->values[0];
        for (size_t j = 0; j < nc; ++j)
            if (!isfinite(temp->values[j]) ||
                (j && !(reverse ? temp->values[j] < temp->values[j - 1]
                                : temp->values[j] > temp->values[j - 1])))
                tt_build_fail(b, file, TT_INVALID, "Nonmonotone temperature table: %s", path);
    }
    if (table->values && (table->rows != rows || table->columns != nc)) {
        tt_build_fail(b, file, TT_INVALID, "Cannot merge different table shapes: %s", path);
        return;
    }
    if (!table->values) {
        table->rows = rows;
        table->columns = nc;
        table->values = tt_build_alloc(b, nc, rows * sizeof(double));
    }
    if (temp) {
        free((void *)table->temperatures);
        table->temperature_count = temp->count;
        double *t = tt_build_alloc(b, temp->count, sizeof(double));
        table->temperatures = t;
        if (t)
            for (size_t j = 0; j < temp->count; ++j)
                t[j] = temp->values[reverse ? temp->count - 1 - j : j];
    }
    if (b->status != TT_OK)
        return;
    double *target = (double *)table->values;
    if (!wave) {
        for (size_t j = 0; j < nc; ++j)
            memcpy(target + j * nr, values->values + (reverse ? nc - 1 - j : j) * nr,
                   nr * sizeof(double));
        return;
    }
    int fold = 0;
    tt_interpolation method = TT_CONSERVING;
    if (cs && s->kind != TT_CS_HNO3)
        method = tt_build_interpolator(b, tt_config_get(file, "interpolator"), &fold);
    tt_config_node padding = cs && !special ? file : config;
    size_t first = 0, end = nw;
    if (cs && !special) {
        double lo = tt_build_number(b, file, "zero below", -DBL_MAX);
        double hi = tt_build_number(b, file, "zero above", DBL_MAX);
        while (first < nw && s->wavelength.edges[first] < lo)
            ++first;
        while (end && s->wavelength.edges[end - 1] > hi)
            --end;
        if (first >= end)
            tt_build_fail(b, file, TT_RANGE, "Empty wavelength range in %s", path);
    }
    double *x = tt_build_alloc(b, nr + 4, 2 * sizeof(double));
    double *work = tt_build_alloc(b, nw + 1, sizeof(double));
    if (!x || !work) {
        free(x);
        free(work);
        return;
    }
    double *y = x + nr + 4;
    for (size_t i = 0; i < nr; ++i) {
        double w = wave->values[i];
        x[i + 2] = cs && s->kind == TT_CS_O3 ? tt_air_refractive_index(w, 2.45e19) * w : w;
    }
    x[0] = 0;
    x[1] = (1 - 1.e-5) * x[2];
    x[nr + 2] = (1 + 1.e-5) * x[nr + 1];
    x[nr + 3] = 1.e38;
    for (size_t j = 0; j < nc && b->status == TT_OK; ++j) {
        const double *column = values->values + (reverse ? nc - 1 - j : j) * nr;
        double lo = extrapolation(b, padding, "lower extrapolation", column[0]);
        /* Quantum-yield add_points queries its upper boundary after inserting
         * two lower points, while retaining the original row count. */
        double upper_boundary = qy ? (nr > 2 ? column[nr - 3] : lo) : column[nr - 1];
        double hi = extrapolation(b, padding, "upper extrapolation", upper_boundary);
        if (cs && (s->kind == TT_CS_HNO3 || s->kind == TT_CS_RONO2) && j == 1) {
            lo = column[0];
            hi = column[nr - 1];
        }
        memcpy(y + 2, column, nr * sizeof(double));
        y[0] = y[1] = lo;
        y[nr + 2] = y[nr + 3] = hi;
        tt_status status =
            tt_interpolate(method, fold, nr + 4, x, y, nw + 1, s->wavelength.edges, work);
        if (status != TT_OK)
            tt_build_fail(b, file, status, "Cannot interpolate dataset %s", path);
        if (b->status == TT_OK)
            for (size_t i = first; i < end; ++i)
                target[j * nw + i] = target[j * nw + i] + work[i];
    }
    free(x);
    free(work);
}

static void constant_table(tt_builder *b, tt_spectrum *s, double value) {
    s->tables = tt_build_alloc(b, 1, sizeof *s->tables);
    if (!s->tables)
        return;
    s->table_count = 1;
    size_t nw = s->wavelength.cells;
    double *v = tt_build_alloc(b, nw, sizeof(double));
    s->tables[0] = (tt_spectral_table){nw, 1, 0, v, NULL};
    if (v)
        for (size_t i = 0; i < nw; ++i)
            v[i] = value;
}

static void prepare_tables(tt_builder *b, tt_config_node config, tt_spectrum *s) {
    int cs = s->role == TT_SPECTRUM_CROSS_SECTION;
    tt_config_node files = tt_config_get(config, "netcdf files");
    size_t count = tt_config_size(files);
    int jpl = cs && s->kind == TT_CS_O3_JPL06;
    int merge = cs && tt_build_boolean(b, config, "merge data", 0);
    /* TUV-x accepts merged files and point data only for base cross sections;
     * merging temperature tables would mislabel their columns. */
    if ((merge || (cs && tt_build_has(config, "data"))) && s->kind != TT_CS_BASE) {
        tt_build_fail(b, config, TT_INVALID,
                      "Merged files and point data require a base cross section");
        return;
    }
    if (jpl)
        count = 3;
    if (count) {
        s->table_count = merge ? 1 : count;
        s->tables = tt_build_alloc(b, s->table_count, sizeof *s->tables);
        const char *const jpl_keys[] = {"base netcdf file", "low-temperature netcdf file",
                                        "high-temperature netcdf file"};
        for (size_t i = 0; i < count && b->status == TT_OK; ++i) {
            tt_config_node file = jpl ? tt_config_get(config, jpl_keys[i]) : tt_config_at(files, i);
            prepare_table(b, s, config, file, &s->tables[merge ? 0 : i]);
        }
    } else if (s->role == TT_SPECTRUM_QUANTUM_YIELD) {
        constant_table(b, s, tt_build_number(b, config, "constant value", 0));
    }
    if (b->status == TT_OK && cs && s->kind == TT_CS_O3) {
        if (s->table_count < 4 || s->tables[1].columns < 4 ||
            s->tables[0].rows != s->wavelength.cells || s->tables[1].rows != s->wavelength.cells)
            tt_build_fail(b, config, TT_INVALID, "O3 requires four compatible parameter tables");
        else
            /* The warm Hartley/Huggins column uses the separate base spectrum. */
            memcpy((double *)s->tables[1].values + 3 * s->wavelength.cells, s->tables[0].values,
                   s->wavelength.cells * sizeof(double));
    }
    if (cs && tt_build_has(config, "data")) {
        tt_config_node data = tt_config_get(config, "data");
        tt_build_keys(b, data, NULL, point_data);
        if (!s->table_count)
            constant_table(b, s, tt_build_number(b, data, "default value", 0));
        if (s->table_count != 1 || !s->tables || s->tables[0].rows != s->wavelength.cells) {
            tt_build_fail(b, data, TT_INVALID, "Point data require one wavelength table");
            return;
        }
        tt_config_node points = tt_config_get(data, "point values");
        for (size_t i = 0; i < tt_config_size(points) && b->status == TT_OK; ++i) {
            tt_config_node point = tt_config_at(points, i);
            double w = tt_build_number(b, point, "wavelength", NAN);
            double v = tt_build_number(b, point, "value", NAN);
            size_t j = 0;
            while (j < s->wavelength.cells && !tt_build_equal(w, s->wavelength.midpoints[j]))
                ++j;
            if (j == s->wavelength.cells)
                tt_build_fail(b, point, TT_RANGE, "Point wavelength %.17g is not on grid", w);
            else
                ((double *)s->tables[0].values)[j] = v;
        }
    }
}

static void prepare_overrides(tt_builder *b, tt_config_node config, tt_spectrum *s) {
    tt_config_node bands = tt_config_get(config, "override bands");
    s->override_count = tt_config_size(bands);
    s->overrides = tt_build_alloc(b, s->override_count, sizeof *s->overrides);
    for (size_t i = 0; i < s->override_count && b->status == TT_OK; ++i) {
        tt_config_node node = tt_config_at(bands, i);
        /* Cross-section bands take no wavelength range in TUV-x. */
        tt_build_keys(b, node, band_required,
                      s->role == TT_SPECTRUM_QUANTUM_YIELD ? band_range : NULL);
        const char *name = tt_build_string(b, node, "band", NULL);
        tt_spectral_override *o = s->overrides + i;
        o->value = tt_build_number(b, node, "value", NAN);
        size_t nw = s->wavelength.cells;
        if (!strcmp(name, "range")) {
            double lo = tt_build_number(b, node, "minimum wavelength", 0);
            double hi = tt_build_number(b, node, "maximum wavelength", DBL_MAX);
            o->end = nw;
            while (o->first < nw && s->wavelength.midpoints[o->first] < lo)
                ++o->first;
            while (o->end && s->wavelength.midpoints[o->end - 1] > hi)
                --o->end;
        } else {
            double lo = 0, hi = 0;
            if (!strcmp(name, "lyman-alpha")) {
                lo = 121.4;
                hi = 121.9;
            } else if (!strcmp(name, "schumann-runge")) {
                lo = 175.4;
                hi = 206.2;
            } else if (!strcmp(name, "schumann-runge continuum")) {
                lo = 121.9;
                hi = 175.4;
            } else
                tt_build_fail(b, node, TT_INVALID, "Unknown override band '%s'", name);
            while (o->first < nw && !tt_build_equal(s->wavelength.edges[o->first], lo))
                ++o->first;
            o->end = 0;
            while (o->end <= nw && !tt_build_equal(s->wavelength.edges[o->end], hi))
                ++o->end;
        }
        if (o->first >= o->end || o->end > nw)
            tt_build_fail(b, node, TT_RANGE, "Override band is outside the wavelength grid");
    }
}

static void prepare_options(tt_builder *b, tt_config_node config, tt_spectrum *s) {
    tt_yield_options *o = &s->yield;
    if (s->role == TT_SPECTRUM_QUANTUM_YIELD && s->kind == TT_QY_ACETONE) {
        o->low_wavelength_value = tt_build_number(b, config, "low wavelength value", .95);
        o->high_wavelength_value = tt_build_number(b, config, "high wavelength value", 0);
        o->minimum_temperature = tt_build_number(b, config, "minimum temperature", 0);
        o->maximum_temperature = tt_build_number(b, config, "maximum temperature", DBL_MAX);
        const char *branch = tt_build_string(b, config, "branch", "CH3CO");
        if (!strcmp(branch, "CO"))
            o->acetone_branches = 1;
        else if (!strcmp(branch, "CH3CO"))
            o->acetone_branches = 2;
        else if (!strcmp(branch, "CO+CH3CO"))
            o->acetone_branches = 3;
        else
            tt_build_fail(b, config, TT_INVALID, "Unknown acetone branch '%s'", branch);
    } else if (s->role == TT_SPECTRUM_QUANTUM_YIELD && s->kind == TT_QY_TAYLOR) {
        o->coefficients = tt_build_array(b, config, "coefficients", &o->coefficient_count);
        o->minimum_wavelength = tt_build_number(b, config, "minimum wavelength", 0);
        o->maximum_wavelength = tt_build_number(b, config, "maximum wavelength", DBL_MAX);
    } else if (s->role == TT_SPECTRUM_QUANTUM_YIELD && s->kind == TT_QY_H2SO4) {
        size_t n = 0, m = 0;
        double *w = tt_build_array(b, config, "parameterized wavelengths", &n);
        double *dt = tt_build_array(b, config, "collision interval s", &m);
        o->molecular_diameter = tt_build_number(b, config, "molecular diameter m", NAN);
        o->molecular_weight = tt_build_number(b, config, "molecular weight kg mol-1", NAN);
        if (n != m)
            tt_build_fail(b, config, TT_INVALID, "Collision intervals do not match wavelengths");
        double *rates = tt_build_alloc(b, s->wavelength.cells, sizeof(double));
        o->collision_rates = rates;
        if (rates)
            for (size_t j = 0; j < s->wavelength.cells; ++j)
                rates[j] = -1;
        for (size_t i = 0; i < n && b->status == TT_OK; ++i) {
            size_t j = 0;
            while (j < s->wavelength.cells && s->wavelength.midpoints[j] != w[i])
                ++j;
            if (j == s->wavelength.cells || dt[i] <= 0)
                tt_build_fail(b, config, TT_RANGE,
                              "Invalid H2SO4 collision interval or wavelength");
            else
                rates[j] = 1.0 / dt[i];
        }
        free(w);
        free(dt);
        /* H2SO4 calls the base yield before replacing its parameterized bins. */
        if (b->status == TT_OK && (!s->table_count || s->tables[0].rows != s->wavelength.cells))
            tt_build_fail(b, config, TT_INVALID, "H2SO4 requires a wavelength table");
        if (b->status == TT_OK)
            for (size_t i = 0; i < s->override_count; ++i)
                for (size_t j = s->overrides[i].first; j < s->overrides[i].end; ++j)
                    ((double *)s->tables[0].values)[j] = s->overrides[i].value;
    } else if (s->role == TT_SPECTRUM_WEIGHT && s->kind == TT_WEIGHT_NOTCH) {
        s->weight_a = tt_build_number(b, config, "notch filter begin", NAN);
        s->weight_b = tt_build_number(b, config, "notch filter end", NAN);
    } else if (s->role == TT_SPECTRUM_WEIGHT && s->kind == TT_WEIGHT_GAUSSIAN)
        s->weight_a = tt_build_number(b, config, "centroid", NAN);
}

tt_status tt_spectrum_create(tt_spectrum_role role, tt_config_node config,
                             const tt_grid *wavelength, const tt_data_source *source,
                             tt_spectrum **out, tt_error *error) {
    if (error)
        memset(error, 0, sizeof *error);
    tt_builder b = {source ? *source : (tt_data_source){0}, error, TT_OK};
    if (!out || !wavelength || !wavelength->cells || (unsigned)role > TT_SPECTRUM_WEIGHT) {
        tt_build_fail(&b, config, TT_INVALID, "Invalid spectrum arguments");
        return b.status;
    }
    tt_spectrum *s = tt_build_alloc(&b, 1, sizeof *s);
    if (!s)
        return b.status;
    s->role = role;
    s->kind = -1;
    s->yield = tt_yield_options_default();
    tt_status status = tt_grid_from_edges(wavelength->cells + 1, wavelength->edges, &s->wavelength);
    if (status != TT_OK)
        tt_build_fail(&b, config, status, "Invalid wavelength grid");
    else if (wavelength->midpoints)
        memcpy(s->wavelength.midpoints, wavelength->midpoints, wavelength->cells * sizeof(double));
    const char *type = tt_build_string(&b, config, "type", NULL);
    if (role == TT_SPECTRUM_CROSS_SECTION) {
        for (int k = 0; k < TT_CS_COUNT; ++k)
            if (!strcmp(type, tt_cross_section_name(k)))
                s->kind = k;
        if (!strcmp(type, "SO2"))
            s->kind = TT_CS_BASE;
        if (!strcmp(type, "temperature based"))
            s->parameterized = 1;
    } else if (role == TT_SPECTRUM_QUANTUM_YIELD) {
        for (int k = 0; k < TT_QY_COUNT; ++k)
            if (!strcmp(type, tt_quantum_yield_name(k)))
                s->kind = k;
    } else
        for (int k = 0; k < TT_WEIGHT_COUNT; ++k)
            if (!strcmp(type, weight_names[k]))
                s->kind = k;
    if (s->kind < 0 && !s->parameterized)
        tt_build_fail(&b, config, TT_INVALID, "Unknown spectrum type '%s'", type);
    spectrum_keys(&b, config, s);
    if (b.status == TT_OK && s->parameterized)
        tt_spectrum_temperature_build(&b, config, s);
    else if (b.status == TT_OK) {
        prepare_tables(&b, config, s);
        prepare_overrides(&b, config, s);
        prepare_options(&b, config, s);
    }
    if (b.status == TT_OK && !s->parameterized) {
        if (role == TT_SPECTRUM_CROSS_SECTION)
            status = tt_cross_section_validate(s->kind, wavelength->cells, s->table_count,
                                               s->tables, NULL);
        else if (role == TT_SPECTRUM_QUANTUM_YIELD)
            status = tt_quantum_yield_validate(s->kind, wavelength->cells, s->table_count,
                                               s->tables, &s->yield, NULL);
        else if (s->kind == TT_WEIGHT_BASE || s->kind == TT_WEIGHT_EPPLEY)
            status = s->table_count && s->tables[0].rows == wavelength->cells ? TT_OK : TT_INVALID;
        if (status != TT_OK)
            tt_build_fail(&b, config, status, "Incompatible tables or options for '%s'", type);
    }
    if (b.status != TT_OK) {
        tt_spectrum_destroy(s);
        return b.status;
    }
    *out = s;
    return TT_OK;
}
void tt_spectrum_destroy(tt_spectrum *s) {
    if (!s)
        return;
    tt_grid_destroy(&s->wavelength);
    for (size_t i = 0; s->tables && i < s->table_count; ++i) {
        free((void *)s->tables[i].values);
        free((void *)s->tables[i].temperatures);
    }
    free(s->tables);
    free(s->overrides);
    free((void *)s->yield.coefficients);
    free((void *)s->yield.collision_rates);
    free((void *)s->temperature.a);
    free((void *)s->temperature.b);
    free((void *)s->temperature.powers);
    free((void *)s->temperature.sigma);
    free((void *)s->temperature.ranges);
    free(s->raw_wavelength);
    free(s->raw_data);
    free(s);
}
size_t tt_spectrum_table_count(const tt_spectrum *s) {
    return s ? s->table_count : 0;
}
const tt_spectral_table *tt_spectrum_table(const tt_spectrum *s, size_t i) {
    return s && i < s->table_count ? s->tables + i : NULL;
}
tt_status tt_spectrum_calculate(const tt_spectrum *s, size_t levels, const double *temperature,
                                const double *air, double *out) {
    if (!s || !out)
        return TT_INVALID;
    tt_spectral_grid grid = {s->wavelength.cells,     levels,      s->wavelength.edges,
                             s->wavelength.midpoints, temperature, air};
    tt_status status;
    if (s->parameterized)
        return tt_temperature_cross_section(&s->temperature, s->raw_count, s->raw_wavelength,
                                            s->raw_data, &grid, out);
    if (s->role == TT_SPECTRUM_CROSS_SECTION)
        status = tt_cross_section(s->kind, &grid, s->table_count, s->tables, out);
    else if (s->role == TT_SPECTRUM_QUANTUM_YIELD)
        status = tt_quantum_yield(s->kind, &grid, s->table_count, s->tables, &s->yield, out);
    else
        return tt_spectral_weight(s->kind, grid.wavelengths, grid.edges, grid.midpoints,
                                  s->table_count ? s->tables[0].values : NULL, s->weight_a,
                                  s->weight_b, out);
    int apply = (s->role == TT_SPECTRUM_CROSS_SECTION && s->kind == TT_CS_BASE) ||
                (s->role == TT_SPECTRUM_QUANTUM_YIELD &&
                 (s->kind == TT_QY_BASE || s->kind == TT_QY_TAYLOR));
    if (status == TT_OK && apply)
        for (size_t i = 0; i < s->override_count; ++i)
            for (size_t w = s->overrides[i].first; w < s->overrides[i].end; ++w)
                for (size_t z = 0; z < levels; ++z)
                    out[w * levels + z] = s->overrides[i].value;
    return status;
}
