/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "build_internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const tt_grid *find_grid(tt_builder *b, tt_config_node node, size_t count,
                                const tt_named_grid *grids, const char *name, const char *units) {
    for (size_t i = 0; i < count; ++i)
        if (grids[i].name && grids[i].units && !strcmp(name, grids[i].name) &&
            !strcmp(units, grids[i].units))
            return grids[i].grid;
    tt_build_fail(b, node, TT_NOT_FOUND, "Missing grid '%s' with units '%s'", name, units);
    return NULL;
}
static tt_status solar_flux(tt_builder *b, tt_config_node config, const tt_grid *grid,
                            tt_profile *out) {
    tt_config_node files = tt_config_get(config, "file path");
    size_t n = tt_config_size(files);
    if (tt_config_type(files) != TT_CONFIG_SEQUENCE || !n || n > 4) {
        tt_build_fail(b, config, TT_INVALID, "Solar flux requires one to four compiled spectra");
        return b->status;
    }
    const tt_dataset *data[4] = {0};
    tt_interpolation methods[4];
    tt_config_node interpolators = tt_config_get(config, "interpolator");
    if (tt_config_type(interpolators) != TT_CONFIG_MISSING && tt_config_size(interpolators) != n)
        tt_build_fail(b, config, TT_INVALID,
                      "Solar interpolator count does not match input spectra");
    for (size_t i = 0; i < n; ++i) {
        const char *path = tt_config_string(tt_config_at(files, i));
        if (!path) {
            tt_build_fail(b, files, TT_INVALID, "Expected solar spectrum path");
            break;
        }
        data[i] = tt_build_dataset(b, files, path);
        const char *name = tt_config_string(tt_config_at(interpolators, i));
        methods[i] = TT_CONSERVING;
        if (name && *name && strcmp(name, "conserving")) {
            if (!strcmp(name, "fractional source"))
                methods[i] = TT_FRACTIONAL_SOURCE;
            else if (!strcmp(name, "fractional target"))
                methods[i] = TT_FRACTIONAL_TARGET;
            else
                tt_build_fail(b, interpolators, TT_INVALID, "Unsupported solar interpolator '%s'",
                              name);
        }
    }
    if (b->status != TT_OK)
        return b->status;
    tt_status status = tt_profile_create(grid->cells, out);
    if (status == TT_OK)
        status = tt_extraterrestrial_flux(grid, n, data, methods, out->midpoints);
    return status;
}

tt_status tt_profile_from_config(tt_config_node config, size_t count, const tt_named_grid *grids,
                                 const tt_data_source *source, tt_profile *out, tt_error *error) {
    if (error)
        memset(error, 0, sizeof *error);
    tt_builder b = {source ? *source : (tt_data_source){0}, error, TT_OK};
    if (!out || (count && !grids)) {
        tt_build_fail(&b, config, TT_INVALID, "Invalid profile arguments");
        return b.status;
    }
    const char *type = tt_build_string(&b, config, "type", NULL);
    /* Keys TUV-x accepts for each profile type (src/profiles/). */
    static const char *const air[] = {"type", "units", "file path", "name", NULL};
    static const char *const scaled[] = {"interpolator", "scale height", NULL};
    static const char *const ozone[] = {"interpolator", "scale height", "reference column", NULL};
    static const char *const solar_angle[] = {"type",      "units",    "year", "month", "day",
                                              "longitude", "latitude", "name", NULL};
    static const char *const distance_keys[] = {"type", "units", "year", "month",
                                                "day",  "name",  NULL};
    static const char *const zone[] = {"time zone", NULL};
    static const char *const flux_keys[] = {"type",         "units", "file path",
                                            "interpolator", "name",  NULL};
    static const char *const diagnostics[] = {"enable diagnostics", NULL};
    static const char *const values_required[] = {"type", "units", "grid", "name", NULL};
    static const char *const values_optional[] = {"values", "uniform value", "scale height", NULL};
    static const char *const file_required[] = {"type", "units", "file path", "grid", "name", NULL};
    if (!strcmp(type, "air"))
        tt_build_keys(&b, config, air, NULL);
    else if (!strcmp(type, "O2"))
        tt_build_keys(&b, config, air, scaled);
    else if (!strcmp(type, "O3"))
        tt_build_keys(&b, config, air, ozone);
    else if (!strcmp(type, "solar zenith angle"))
        tt_build_keys(&b, config, solar_angle, zone);
    else if (!strcmp(type, "Earth-Sun distance"))
        tt_build_keys(&b, config, distance_keys, zone);
    else if (!strcmp(type, "extraterrestrial flux"))
        tt_build_keys(&b, config, flux_keys, diagnostics);
    else if (!strcmp(type, "from config file"))
        tt_build_keys(&b, config, values_required, values_optional);
    else if (!strcmp(type, "from csv file"))
        tt_build_keys(&b, config, file_required, scaled);
    tt_profile_kind kind = TT_PROFILE_LINEAR;
    int file = !strcmp(type, "from csv file"), values = !strcmp(type, "from config file");
    int zenith = !strcmp(type, "solar zenith angle"),
        distance = !strcmp(type, "Earth-Sun distance");
    int flux = !strcmp(type, "extraterrestrial flux");
    if (!strcmp(type, "air"))
        kind = TT_PROFILE_AIR;
    else if (!strcmp(type, "O2"))
        kind = TT_PROFILE_O2;
    else if (!strcmp(type, "O3"))
        kind = TT_PROFILE_O3;
    else if (!file && !values && !zenith && !distance && !flux)
        tt_build_fail(&b, config, TT_INVALID, "Unknown profile type '%s'", type);
    const char *name = "height", *units = "km";
    if (zenith || distance) {
        name = "time";
        units = "hours";
    } else if (flux) {
        name = "wavelength";
        units = "nm";
    } else if (file || values) {
        tt_config_node grid = tt_config_get(config, "grid");
        name = tt_build_string(&b, grid, "name", NULL);
        units = tt_build_string(&b, grid, "units", NULL);
    }
    const tt_grid *grid = find_grid(&b, config, count, grids, name, units);
    if (!grid)
        tt_build_fail(&b, config, TT_INVALID, "Missing profile grid storage");
    if (b.status != TT_OK)
        return b.status;
    tt_status status = TT_OK;
    tt_profile profile = {0};
    if (zenith || distance) {
        int year = tt_build_integer(&b, config, "year"),
            month = tt_build_integer(&b, config, "month");
        int day = tt_build_integer(&b, config, "day");
        double zone = tt_build_number(&b, config, "time zone", 0);
        double latitude = zenith ? tt_build_number(&b, config, "latitude", NAN) : 0;
        double longitude = zenith ? tt_build_number(&b, config, "longitude", NAN) : 0;
        tt_profile other = {0};
        if (b.status == TT_OK)
            status = tt_solar_profiles(grid, year, month, day, zone, latitude, longitude,
                                       zenith ? &profile : &other, distance ? &profile : &other);
        tt_profile_destroy(&other);
    } else if (flux)
        status = solar_flux(&b, config, grid, &profile);
    else {
        tt_profile_options options = tt_profile_options_default(kind);
        options.scale_height = tt_build_number(&b, config, "scale height", options.scale_height);
        options.reference_column =
            tt_build_number(&b, config, "reference column", options.reference_column);
        if (!strcmp(name, "height") && !strcmp(units, "km"))
            options.length_conversion = 1.e5;
        if (values) {
            size_t n = grid->cells + 1;
            double *v = NULL;
            if (tt_build_has(config, "values")) {
                v = tt_build_array(&b, config, "values", &n);
                if (n != grid->cells + 1)
                    tt_build_fail(&b, config, TT_INVALID, "Profile values do not match grid edges");
            } else {
                double value = tt_build_number(&b, config, "uniform value", NAN);
                v = tt_build_alloc(&b, n, sizeof(double));
                if (v)
                    for (size_t i = 0; i < n; ++i)
                        v[i] = value;
            }
            if (b.status == TT_OK)
                status = tt_profile_from_edges(grid, v, &options, &profile);
            free(v);
        } else {
            const char *interpolator = tt_build_string(&b, config, "interpolator", "linear");
            /* The other reference choices return bin values, then incorrectly
             * index them as edge values. Reject that undefined construction. */
            if (strcmp(interpolator, "linear"))
                tt_build_fail(&b, config, TT_INVALID,
                              "Sample profiles require linear interpolation");
            const char *path = tt_build_string(&b, config, "file path", NULL);
            const tt_dataset *data = tt_build_dataset(&b, config, path);
            if (b.status == TT_OK)
                status = tt_profile_from_dataset(grid, data, &options, &profile);
        }
    }
    if (status != TT_OK)
        tt_build_fail(&b, config, status, "Cannot construct profile '%s': %s", type,
                      tt_status_string(status));
    if (b.status != TT_OK) {
        tt_profile_destroy(&profile);
        return b.status;
    }
    *out = profile;
    return TT_OK;
}
