/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "table_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

tt_status tt_extraterrestrial_flux(const tt_grid *grid, size_t count, const tt_dataset *const *data,
                                   const tt_interpolation *methods, double *out) {
    if (!grid || !grid->cells || !grid->edges || !grid->midpoints || !grid->widths || !count ||
        count > 4 || !data || !out || grid->cells >= SIZE_MAX / sizeof(double) - 1)
        return TT_INVALID;
    static const double boundary[] = {0., 150.01, 200.07, 1000.99, DBL_MAX};
    const double hc = 6.626068e-34 * 2.99792458e8;
    double *interpolated = malloc((grid->cells + 1) * sizeof(double));
    if (!interpolated)
        return TT_MEMORY;
    memset(out, 0, grid->cells * sizeof(double));
    tt_status status = TT_OK;
    for (size_t file = 0; file < count; ++file) {
        double *raw = NULL, *memory = NULL;
        size_t rows;
        status = tt_text_columns(data[file], 2, 0, &rows, &raw);
        if (status)
            break;
        if (rows > SIZE_MAX / sizeof(double) / 2 - 4) {
            free(raw);
            status = TT_INVALID;
            break;
        }
        memory = malloc(2 * (rows + 4) * sizeof(double));
        if (!memory) {
            free(raw);
            status = TT_MEMORY;
            break;
        }
        double *x = memory, *y = memory + rows + 4;
        size_t n;
        int neckel = data[file]->path && strstr(data[file]->path, "neckel.flx");
        if (neckel) {
            for (size_t i = 0; i < rows; ++i) {
                x[i] = raw[i] - (raw[i] < 630 ? .5 : raw[i] < 870 ? 1. : 2.5);
                y[i] = 1.e13 * hc * raw[rows + i] / raw[i];
            }
            x[rows] = x[rows - 1] + 2.5;
            y[rows] = 0;
            n = rows + 1;
        } else {
            x[0] = 0;
            x[1] = (1 - 1.e-5) * raw[0];
            y[0] = y[1] = 0;
            memcpy(x + 2, raw, rows * sizeof(double));
            memcpy(y + 2, raw + rows, rows * sizeof(double));
            x[rows + 2] = (1 + 1.e-5) * raw[rows - 1];
            x[rows + 3] = 1.e38;
            y[rows + 2] = y[rows + 3] = 0;
            n = rows + 4;
        }
        tt_interpolation method = methods ? methods[file] : TT_CONSERVING;
        status = tt_interpolate(method, 0, n, x, y, grid->cells + 1, grid->edges, interpolated);
        free(raw);
        free(memory);
        if (status)
            break;
        for (size_t w = 0; w < grid->cells; ++w)
            if (boundary[file] <= grid->edges[w] && grid->edges[w] < boundary[file + 1])
                out[w] = interpolated[w];
    }
    if (!status)
        for (size_t w = 0; w < grid->cells; ++w)
            /* The original literal is default REAL, then promoted to double. */
            out[w] = (double)1.e-13f * out[w] * grid->midpoints[w] * grid->widths[w] / hc;
    free(interpolated);
    return status;
}
