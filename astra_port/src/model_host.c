/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "model_internal.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

tt_status tt_model_copy_profile(const tt_profile *source, tt_profile *out) {
    if (!source)
        return TT_INVALID;
    size_t n = source->cells;
    if (n > (SIZE_MAX / sizeof(double) - 2) / 6)
        return TT_INVALID;
    const double *fields[] = {source->edges,   source->midpoints,        source->differences,
                              source->columns, source->extended_columns, source->burden};
    for (size_t f = 0; f < 6; ++f) {
        size_t count = n + (f == 0 || f == 4);
        if (count && !fields[f])
            return TT_INVALID;
        for (size_t i = 0; i < count; ++i)
            if (!isfinite(fields[f][i]))
                return TT_INVALID;
    }
    tt_status status = tt_profile_create(n, out);
    if (status)
        return status;
    double *target[] = {out->edges,   out->midpoints,        out->differences,
                        out->columns, out->extended_columns, out->burden};
    for (size_t f = 0; f < 6; ++f) {
        size_t count = n + (f == 0 || f == 4);
        if (count)
            memcpy(target[f], fields[f], count * sizeof(double));
    }
    return TT_OK;
}
tt_status tt_model_copy_radiator(const tt_named_radiator *source, tt_model_radiator_entry *out,
                                 size_t nl, size_t nw) {
    if (!source || source->layers != nl || source->wavelengths != nw ||
        !source->radiator.optical_depth)
        return TT_INVALID;
    const tt_radiator *r = &source->radiator;
    size_t n = nl * nw;
    for (size_t i = 0; i < n; ++i) {
        double od = r->optical_depth[i],
               a = r->single_scattering_albedo ? r->single_scattering_albedo[i] : 0,
               g = r->asymmetry ? r->asymmetry[i] : 0;
        if (!isfinite(od) || od < 0 || !isfinite(a) || a < 0 || a > 1 || !isfinite(g) ||
            fabs(g) > 1)
            return TT_INVALID;
    }
    memmove(out->data, r->optical_depth, n * sizeof(double));
    if (r->single_scattering_albedo)
        memmove(out->data + n, r->single_scattering_albedo, n * sizeof(double));
    else
        memset(out->data + n, 0, n * sizeof(double));
    if (r->asymmetry)
        memmove(out->data + 2 * n, r->asymmetry, n * sizeof(double));
    else
        memset(out->data + 2 * n, 0, n * sizeof(double));
    out->view.is_air = r->is_air;
    out->gas = 0;
    return TT_OK;
}
tt_status tt_model_update_grid(tt_model *m, const char *name, const char *units, size_t cells,
                               const double *edges, const double *midpoints) {
    if (!m || !name || !units)
        return TT_INVALID;
    size_t i = tt_model_grid_index(m, name, units);
    if (i == SIZE_MAX)
        return TT_NOT_FOUND;
    if (i == m->wavelength || cells != m->grids[i].value.cells)
        return TT_INVALID;
    return tt_grid_update(&m->grids[i].value, edges, midpoints);
}
tt_status tt_model_update_profile(tt_model *m, const char *name, const char *units, size_t cells,
                                  const double *edges, const double *midpoints,
                                  const double *columns, const double *scale_height,
                                  const double *top_column) {
    if (!m || !name || !units)
        return TT_INVALID;
    size_t i = tt_model_profile_index(m, name, units);
    if (i == SIZE_MAX)
        return TT_NOT_FOUND;
    if (cells != m->profiles[i].value.cells)
        return TT_INVALID;
    return tt_profile_update(&m->profiles[i].value, edges, midpoints, columns, scale_height,
                             top_column);
}
tt_status tt_model_set_radiator(tt_model *m, const tt_named_radiator *r) {
    if (!m || !r || !r->name)
        return TT_INVALID;
    for (size_t i = 0; i < m->radiator_count; ++i)
        if (!strcmp(m->radiators[i].name, r->name)) {
            tt_status status = tt_model_copy_radiator(r, m->radiators + i, m->shape.levels - 1,
                                                      m->shape.wavelengths);
            if (!status)
                m->radiator_views[i] = m->radiators[i].view;
            return status;
        }
    return TT_NOT_FOUND;
}
