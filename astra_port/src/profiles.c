/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "table_internal.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void tt_profile_destroy(tt_profile *profile) {
    if (profile) {
        free(profile->edges);
        *profile = (tt_profile){0};
    }
}

tt_profile_options tt_profile_options_default(tt_profile_kind kind) {
    return (tt_profile_options){kind,
                                kind == TT_PROFILE_AIR || kind == TT_PROFILE_O2 ? 8.01
                                : kind == TT_PROFILE_O3                         ? 4.5
                                                                                : 0.,
                                300., kind == TT_PROFILE_LINEAR ? 1. : 1.e5};
}

tt_status tt_profile_create(size_t n, tt_profile *out) {
    if (!out || n > (SIZE_MAX / sizeof(double) - 2) / 6)
        return TT_INVALID;
    double *p = calloc(6 * n + 2, sizeof(double));
    if (!p)
        return TT_MEMORY;
    *out =
        (tt_profile){n, p, p + n + 1, p + 2 * n + 1, p + 3 * n + 1, p + 4 * n + 1, p + 5 * n + 2};
    return TT_OK;
}

tt_status tt_profile_update(tt_profile *p, const double *edges, const double *midpoints,
                            const double *columns, const double *scale_height,
                            const double *top_column) {
    if (!p || !p->edges || !edges || (scale_height && (!columns || top_column)) ||
        (scale_height && !isfinite(*scale_height)) || (top_column && !isfinite(*top_column)))
        return TT_INVALID;
    size_t n = p->cells;
    for (size_t i = 0; i <= n; ++i)
        if (!isfinite(edges[i]) || (i < n && midpoints && !isfinite(midpoints[i])) ||
            (i < n && columns && !isfinite(columns[i])))
            return TT_INVALID;
    memmove(p->edges, edges, (n + 1) * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        p->differences[i] = p->edges[i + 1] - p->edges[i];
        p->midpoints[i] = midpoints ? midpoints[i] : (p->edges[i] + p->edges[i + 1]) * .5;
    }
    if (columns) {
        memmove(p->columns, columns, n * sizeof(double));
        memcpy(p->extended_columns, p->columns, n * sizeof(double));
        p->extended_columns[n] =
            top_column ? *top_column : p->edges[n] * (scale_height ? *scale_height : 0) * 1.e5;
        if (n)
            p->columns[n - 1] += p->extended_columns[n];
        double accum = 0;
        for (size_t i = n; i-- > 0;) {
            accum += p->columns[i];
            p->burden[i] = accum;
        }
    }
    return TT_OK;
}

static void profile_columns(const tt_grid *grid, const tt_profile_options *options, tt_profile *p) {
    size_t n = grid->cells;
    for (size_t i = 0; i < n; ++i) {
        p->midpoints[i] = .5 * (p->edges[i] + p->edges[i + 1]);
        p->differences[i] = p->edges[i + 1] - p->edges[i];
        if (options->kind == TT_PROFILE_AIR || options->kind == TT_PROFILE_O2)
            p->columns[i] = grid->widths[i] * sqrt(p->edges[i]) * sqrt(p->edges[i + 1]) *
                            options->length_conversion;
        else if (options->kind == TT_PROFILE_O3)
            p->columns[i] = grid->widths[i] * p->midpoints[i] * options->length_conversion;
        else
            p->columns[i] = p->midpoints[i] * grid->widths[i] * options->length_conversion;
        p->extended_columns[i] = p->columns[i];
    }
    p->extended_columns[n] = p->edges[n] * options->scale_height * options->length_conversion;
    if (n)
        p->columns[n - 1] = p->columns[n - 1] + p->extended_columns[n];
    double accum = 0;
    for (size_t i = n; i-- > 0;) {
        accum = accum + p->columns[i];
        p->burden[i] = accum;
    }
}

tt_status tt_profile_from_edges(const tt_grid *grid, const double *values,
                                const tt_profile_options *options, tt_profile *out) {
    if (!grid || !grid->edges || !values || !options || !out || options->kind < TT_PROFILE_LINEAR ||
        options->kind > TT_PROFILE_O3 || !isfinite(options->scale_height) ||
        options->scale_height < 0 || !isfinite(options->length_conversion) ||
        options->length_conversion <= 0 || !isfinite(options->reference_column) ||
        options->reference_column < 0)
        return TT_INVALID;
    for (size_t i = 0; i <= grid->cells; ++i)
        if (!isfinite(values[i]) || (options->kind != TT_PROFILE_LINEAR && values[i] < 0))
            return TT_INVALID;
    tt_profile profile;
    tt_status status = tt_profile_create(grid->cells, &profile);
    if (status)
        return status;
    memcpy(profile.edges, values, (grid->cells + 1) * sizeof(double));
    profile_columns(grid, options, &profile);
    if (options->kind == TT_PROFILE_O3 && options->reference_column != 1) {
        double total = 0;
        for (size_t i = 0; i < grid->cells; ++i)
            total = total + profile.columns[i];
        double input_du = total / 2.687e16, factor = options->reference_column / input_du;
        if (!isfinite(factor)) {
            tt_profile_destroy(&profile);
            return TT_RANGE;
        }
        if (factor != 1) {
            for (size_t i = 0; i <= grid->cells; ++i)
                profile.edges[i] = factor * profile.edges[i];
            profile_columns(grid, options, &profile);
        }
    }
    *out = profile;
    return TT_OK;
}

tt_status tt_profile_from_samples(const tt_grid *grid, size_t count, const double *coordinates,
                                  const double *values, const tt_profile_options *options,
                                  tt_profile *out) {
    if (!grid || !grid->edges || count < 2 || count > SIZE_MAX / sizeof(double) / 2 ||
        !coordinates || !values || !options || !out)
        return TT_INVALID;
    size_t extra = 0;
    if (options->kind == TT_PROFILE_O3) {
        /* TUV-x appends 1-km points, accumulating the coordinate, until one
         * lies above the grid top. floor(top - last) + 1 can be one short when
         * the difference rounds below an integer, so count the same way. */
        double top = grid->edges[grid->cells], last = coordinates[count - 1];
        double estimate = floor(fmax(0, top - last)) + 2;
        if (!isfinite(estimate) || estimate > (double)(SIZE_MAX / sizeof(double) / 2 - count))
            return TT_RANGE;
        size_t limit = (size_t)estimate;
        do {
            last = last + 1;
            ++extra;
        } while (!(last > top) && extra < limit);
        if (!(last > top)) /* coordinates beyond double resolution */
            return TT_RANGE;
    }
    if (count > SIZE_MAX / sizeof(double) / 2 - extra ||
        grid->cells >= SIZE_MAX / sizeof(double) - 1)
        return TT_INVALID;
    size_t n = count + extra;
    double *memory = malloc(2 * n * sizeof(double));
    double *edges = malloc((grid->cells + 1) * sizeof(double));
    if (!memory || !edges) {
        free(memory);
        free(edges);
        return TT_MEMORY;
    }
    double *x = memory, *y = memory + n;
    memcpy(x, coordinates, count * sizeof(double));
    memcpy(y, values, count * sizeof(double));
    if (options->kind == TT_PROFILE_AIR || options->kind == TT_PROFILE_O2) {
        x[count - 1] = x[count - 1] + .001;
        for (size_t i = 0; i < count; ++i)
            y[i] = log(y[i]);
    } else if (options->kind == TT_PROFILE_O3) {
        double factor = options->scale_height != 0 ? exp(-1. / options->scale_height) : 0;
        for (size_t i = count; i < n; ++i) {
            x[i] = x[i - 1] + 1;
            y[i] = y[i - 1] * factor;
        }
    }
    tt_status status = tt_interpolate(TT_LINEAR, 0, n, x, y, grid->cells + 1, grid->edges, edges);
    if (!status) {
        if (options->kind == TT_PROFILE_AIR || options->kind == TT_PROFILE_O2)
            for (size_t i = 0; i <= grid->cells; ++i) {
                edges[i] = exp(edges[i]);
                if (options->kind == TT_PROFILE_O2)
                    edges[i] = .2095 * edges[i];
            }
        status = tt_profile_from_edges(grid, edges, options, out);
    }
    free(memory);
    free(edges);
    return status;
}

tt_status tt_profile_from_dataset(const tt_grid *grid, const tt_dataset *data,
                                  const tt_profile_options *options, tt_profile *out) {
    double *values;
    size_t rows;
    tt_status status = tt_text_columns(data, 2, 0, &rows, &values);
    if (status)
        return status;
    status = tt_profile_from_samples(grid, rows, values, values + rows, options, out);
    free(values);
    return status;
}

tt_status tt_solar_profiles(const tt_grid *time, int year, int month, int day, double time_zone,
                            double latitude, double longitude, tt_profile *zenith,
                            tt_profile *distance) {
    int day_of_year = tt_day_of_year(year, month, day);
    if (!time || !time->edges || !zenith || !distance || zenith == distance || day_of_year < 1 ||
        year < 1950 || year > 2050 || !isfinite(time_zone))
        return TT_INVALID;
    tt_profile z, d;
    tt_status status = tt_profile_create(time->cells, &z);
    if (status)
        return status;
    status = tt_profile_create(time->cells, &d);
    if (status) {
        tt_profile_destroy(&z);
        return status;
    }
    for (size_t i = 0; i <= time->cells; ++i) {
        /* Preserve the two distinct time-zone conventions of the reference. */
        z.edges[i] = 90 - tt_solar_elevation(year, day_of_year, time->edges[i] + time_zone,
                                             latitude, longitude);
        d.edges[i] = tt_earth_sun_distance(year, day_of_year, time->edges[i] - time_zone);
        if (!isfinite(z.edges[i]) || !isfinite(d.edges[i])) {
            tt_profile_destroy(&z);
            tt_profile_destroy(&d);
            return TT_INVALID;
        }
    }
    tt_profile_options options = tt_profile_options_default(TT_PROFILE_LINEAR);
    profile_columns(time, &options, &z);
    profile_columns(time, &options, &d);
    *zenith = z;
    *distance = d;
    return TT_OK;
}
