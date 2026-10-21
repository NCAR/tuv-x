/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "table_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void tt_grid_destroy(tt_grid *grid) {
    if (grid) {
        free(grid->edges);
        *grid = (tt_grid){0};
    }
}

tt_status tt_grid_update(tt_grid *grid, const double *edges, const double *midpoints) {
    if (!grid || !grid->edges || !edges || (grid->cells && (!grid->midpoints || !grid->widths)))
        return TT_INVALID;
    size_t n = grid->cells;
    for (size_t i = 0; i <= n; ++i)
        if (!isfinite(edges[i]) || (i && edges[i] <= edges[i - 1]) ||
            (i < n && midpoints && !isfinite(midpoints[i])))
            return TT_INVALID;
    memmove(grid->edges, edges, (n + 1) * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        grid->widths[i] = grid->edges[i + 1] - grid->edges[i];
        grid->midpoints[i] = midpoints ? midpoints[i] : (grid->edges[i] + grid->edges[i + 1]) * .5;
    }
    return TT_OK;
}

tt_status tt_grid_from_edges(size_t count, const double *edges, tt_grid *out) {
    if (!count || !edges || !out || count > SIZE_MAX / sizeof(double) / 3)
        return TT_INVALID;
    for (size_t i = 0; i < count; ++i)
        if (!isfinite(edges[i]) || (i && edges[i] <= edges[i - 1]))
            return TT_INVALID;
    double *memory = malloc((3 * count - 2) * sizeof(double));
    if (!memory)
        return TT_MEMORY;
    tt_grid grid = {count - 1, memory, memory + count, memory + 2 * count - 1};
    memcpy(grid.edges, edges, count * sizeof(double));
    for (size_t i = 0; i < grid.cells; ++i) {
        grid.midpoints[i] = .5 * (edges[i] + edges[i + 1]);
        grid.widths[i] = edges[i + 1] - edges[i];
    }
    *out = grid;
    return TT_OK;
}

tt_status tt_grid_equal(double first, double last, double interval, tt_grid *out) {
    if (!out || !isfinite(first) || !isfinite(last) || !isfinite(interval) || last < first ||
        interval <= 0)
        return TT_INVALID;
    double cells = (last - first) / interval;
    if (!isfinite(cells) || cells >= (double)(SIZE_MAX / sizeof(double) / 3 - 2))
        return TT_RANGE;
    size_t n = (size_t)cells;
    double endpoint = first + (double)n * interval;
    /* Fortran SPACING uses TINY for zero and subnormal arguments. */
    double spacing =
        endpoint == 0 ? DBL_MIN : fmax(DBL_MIN, scalbn(1., ilogb(fabs(endpoint)) - 52));
    if (last - endpoint > spacing)
        ++n;
    double *edges = malloc((n + 1) * sizeof(double));
    if (!edges)
        return TT_MEMORY;
    for (size_t i = 0; i <= n; ++i)
        edges[i] = fmin(first + (double)i * interval, last);
    tt_status status = tt_grid_from_edges(n + 1, edges, out);
    free(edges);
    return status;
}

tt_status tt_grid_from_dataset(const tt_dataset *data, tt_grid *out) {
    double *edges;
    size_t count;
    tt_status status = tt_text_columns(data, 1, 1, &count, &edges);
    if (status)
        return status;
    status = tt_grid_from_edges(count, edges, out);
    free(edges);
    return status;
}
