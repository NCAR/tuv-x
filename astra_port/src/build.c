/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "build_internal.h"
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void tt_build_fail(tt_builder *b, tt_config_node node, tt_status status, const char *format, ...) {
    if (b->status != TT_OK)
        return;
    b->status = status;
    if (b->error) {
        b->error->status = status;
        tt_config_location(node, &b->error->line, &b->error->column);
        va_list args;
        va_start(args, format);
        vsnprintf(b->error->message, sizeof b->error->message, format, args);
        va_end(args);
    }
}
int tt_build_has(tt_config_node node, const char *key) {
    return tt_config_type(tt_config_get(node, key)) != TT_CONFIG_MISSING;
}
const char *tt_build_string(tt_builder *b, tt_config_node node, const char *key,
                            const char *fallback) {
    tt_config_node value = tt_config_get(node, key);
    const char *s = tt_config_string(value);
    if (tt_config_type(value) == TT_CONFIG_MISSING && fallback)
        return fallback;
    if (!s)
        tt_build_fail(b, node, TT_INVALID, "Expected string '%s'", key);
    return s ? s : "";
}
double tt_build_number(tt_builder *b, tt_config_node node, const char *key, double fallback) {
    tt_config_node value = tt_config_get(node, key);
    if (tt_config_type(value) == TT_CONFIG_MISSING && !isnan(fallback))
        return fallback;
    double number = 0;
    if (tt_config_number(value, &number) != TT_OK)
        tt_build_fail(b, node, TT_INVALID, "Expected finite number '%s'", key);
    return number;
}
int tt_build_integer(tt_builder *b, tt_config_node node, const char *key) {
    int value = 0;
    if (tt_config_integer(tt_config_get(node, key), &value) != TT_OK)
        tt_build_fail(b, node, TT_INVALID, "Expected integer '%s'", key);
    return value;
}
int tt_build_boolean(tt_builder *b, tt_config_node node, const char *key, int fallback) {
    tt_config_node value = tt_config_get(node, key);
    if (tt_config_type(value) == TT_CONFIG_MISSING)
        return fallback;
    int result = 0;
    if (tt_config_boolean(value, &result) != TT_OK)
        tt_build_fail(b, node, TT_INVALID, "Expected boolean '%s'", key);
    return result;
}
void *tt_build_alloc(tt_builder *b, size_t count, size_t size) {
    if (b->status != TT_OK)
        return NULL;
    void *p = size && count <= SIZE_MAX / size ? calloc(count ? count : 1, size) : NULL;
    if (!p)
        tt_build_fail(b, (tt_config_node){0}, TT_MEMORY, "Allocation failed");
    return p;
}
double *tt_build_array(tt_builder *b, tt_config_node node, const char *key, size_t *count) {
    tt_config_node values = tt_config_get(node, key);
    *count = tt_config_size(values);
    if (tt_config_type(values) != TT_CONFIG_SEQUENCE)
        tt_build_fail(b, node, TT_INVALID, "Expected numeric array '%s'", key);
    double *out = tt_build_alloc(b, *count, sizeof(double));
    for (size_t i = 0; out && i < *count; ++i)
        if (tt_config_number(tt_config_at(values, i), out + i) != TT_OK)
            tt_build_fail(b, values, TT_INVALID, "Invalid number in '%s'", key);
    return out;
}
const tt_dataset *tt_build_dataset(tt_builder *b, tt_config_node node, const char *path) {
    if (b->status != TT_OK)
        return NULL;
    const tt_dataset *data =
        b->source.resolve ? b->source.resolve(b->source.context, path) : tt_dataset_find(path);
    if (!data)
        tt_build_fail(b, node, TT_NOT_FOUND, "Dataset is not compiled in: %s", path);
    return data;
}
tt_interpolation tt_build_interpolator(tt_builder *b, tt_config_node node, int *fold) {
    static const char *const required[] = {"type", NULL}, *const optional[] = {"fold in", NULL};
    tt_build_keys(b, node, required, optional);
    const char *type = tt_build_string(b, node, "type", "conserving");
    *fold = tt_build_boolean(b, node, "fold in", 0);
    static const char *const names[] = {"linear", "conserving", "fractional source",
                                        "fractional target"};
    for (int i = 0; i < 4; ++i)
        if (!strcmp(type, names[i])) {
            if (*fold && i < 2)
                tt_build_fail(b, node, TT_INVALID, "This interpolator cannot fold in");
            return (tt_interpolation)i;
        }
    tt_build_fail(b, node, TT_INVALID, "Unknown interpolator '%s'", type);
    return TT_CONSERVING;
}
static int listed(const char *key, const char *const *list) {
    for (; list && *list; ++list)
        if (!strcmp(key, *list))
            return 1;
    return 0;
}
void tt_build_keys(tt_builder *b, tt_config_node node, const char *const *required,
                   const char *const *optional) {
    tt_config_kind kind = tt_config_type(node);
    if (b->status != TT_OK || kind == TT_CONFIG_MISSING)
        return;
    if (kind != TT_CONFIG_MAPPING) {
        tt_build_fail(b, node, TT_INVALID, "Expected a configuration mapping");
        return;
    }
    for (size_t i = 0; i < tt_config_size(node); ++i) {
        const char *key = tt_config_key(node, i);
        if (strncmp(key, "__", 2) && !listed(key, required) && !listed(key, optional)) {
            tt_build_fail(b, tt_config_at(node, i), TT_INVALID, "Unknown configuration key '%s'",
                          key);
            return;
        }
    }
    for (; required && *required; ++required)
        if (!tt_build_has(node, *required)) {
            tt_build_fail(b, node, TT_INVALID, "Missing configuration key '%s'", *required);
            return;
        }
}
int tt_build_equal(double a, double b) {
    return a == b || 2.0 * fabs(a - b) / (fabs(a) + fabs(b)) < 1.e-10 || fabs(a - b) <= 1.e-30;
}
tt_status tt_grid_from_config(tt_config_node config, const tt_data_source *source, tt_grid *out,
                              tt_error *error) {
    if (error)
        memset(error, 0, sizeof *error);
    tt_builder b = {source ? *source : (tt_data_source){0}, error, TT_OK};
    const char *type = tt_build_string(&b, config, "type", NULL);
    tt_status status = TT_INVALID;
    /* Keys TUV-x accepts (src/grids/). */
    static const char *const equal[] = {"type",       "units", "begins at", "ends at",
                                        "cell delta", "name",  NULL};
    static const char *const values[] = {"type", "units", "values", "name", NULL};
    static const char *const file[] = {"type", "units", "file path", "name", NULL};
    const char *const *keys = !strcmp(type, "equal interval")     ? equal
                              : !strcmp(type, "from config file") ? values
                              : !strcmp(type, "from csv file")    ? file
                                                                  : NULL;
    if (keys)
        tt_build_keys(&b, config, keys, NULL);
    if (!out)
        tt_build_fail(&b, config, TT_INVALID, "Missing grid output");
    if (!strcmp(type, "equal interval")) {
        double first = tt_build_number(&b, config, "begins at", NAN);
        double last = tt_build_number(&b, config, "ends at", NAN);
        double delta = tt_build_number(&b, config, "cell delta", NAN);
        if (b.status == TT_OK)
            status = tt_grid_equal(first, last, delta, out);
    } else if (!strcmp(type, "from config file")) {
        size_t n = 0;
        double *v = tt_build_array(&b, config, "values", &n);
        if (b.status == TT_OK)
            status = tt_grid_from_edges(n, v, out);
        free(v);
    } else if (!strcmp(type, "from csv file")) {
        const char *path = tt_build_string(&b, config, "file path", NULL);
        const tt_dataset *data = tt_build_dataset(&b, config, path);
        if (b.status == TT_OK)
            status = tt_grid_from_dataset(data, out);
    } else
        tt_build_fail(&b, config, TT_INVALID, "Unknown grid type '%s'", type);
    if (status != TT_OK)
        tt_build_fail(&b, config, status, "Invalid grid: %s", tt_status_string(status));
    return b.status;
}
