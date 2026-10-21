/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

tt_status tt_gas_optics(size_t layers, size_t wavelengths, const double *columns,
                        const double *cross_section, int is_air, double *optical_depth,
                        double *albedo, double *asymmetry) {
    if (!layers || !wavelengths || wavelengths > SIZE_MAX / layers || !columns || !cross_section ||
        !optical_depth || !albedo || !asymmetry)
        return TT_INVALID;
    for (size_t i = 0; i < layers; ++i)
        if (!isfinite(columns[i]) || columns[i] < 0)
            return TT_INVALID;
    for (size_t w = 0; w < wavelengths; ++w)
        for (size_t i = 0; i < layers; ++i) {
            size_t j = w * layers + i;
            if (!isfinite(cross_section[j]) || cross_section[j] < 0)
                return TT_INVALID;
            optical_depth[j] = columns[i] * cross_section[j];
            albedo[j] = is_air ? 1 : 0;
            asymmetry[j] = 0;
        }
    return TT_OK;
}

tt_status tt_aerosol_optics(const tt_grid *height, const tt_grid *wavelength, size_t count,
                            const double *depth_samples, double albedo, double asymmetry,
                            double tau550, double alpha, double *optical_depth,
                            double *single_scattering_albedo, double *asymmetry_factor) {
    if (!height || !height->cells || !height->edges || !wavelength || !wavelength->cells ||
        !wavelength->midpoints || !count || !depth_samples || !optical_depth ||
        !single_scattering_albedo || !asymmetry_factor || !isfinite(albedo) || albedo < 0 ||
        albedo > 1 || !isfinite(asymmetry) || asymmetry < -1 || asymmetry > 1 ||
        !isfinite(tau550) || !isfinite(alpha) || count > SIZE_MAX / sizeof(double) / 3 ||
        height->cells > SIZE_MAX / sizeof(double) / 3 - count ||
        wavelength->cells > SIZE_MAX / height->cells)
        return TT_INVALID;
    size_t n = height->cells;
    for (size_t i = 0; i < count; ++i)
        if (!isfinite(depth_samples[i]) || depth_samples[i] < 0)
            return TT_INVALID;
    double *memory = malloc(3 * (n + count) * sizeof(double));
    if (!memory)
        return TT_MEMORY;
    double *x = memory, *y = x + count, *weighted = y + count;
    double *depth = weighted + count, *ssa = depth + n, *g = ssa + n;
    tt_status status = TT_OK;
    if (count == 1) {
        for (size_t i = 0; i < n; ++i) {
            depth[i] = depth_samples[0];
            ssa[i] = depth[i] * albedo;
            g[i] = depth[i] * asymmetry;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            x[i] = (double)i;
            y[i] =
                i + 1 < count ? .5 * (depth_samples[i] + depth_samples[i + 1]) : depth_samples[i];
            weighted[i] = y[i] * albedo;
        }
        status = tt_interpolate(TT_FRACTIONAL_SOURCE, 1, count, x, y, n + 1, height->edges, depth);
        if (!status)
            status = tt_interpolate(TT_FRACTIONAL_SOURCE, 1, count, x, weighted, n + 1,
                                    height->edges, ssa);
        for (size_t i = 0; i < count; ++i)
            weighted[i] = y[i] * asymmetry;
        if (!status)
            status = tt_interpolate(TT_FRACTIONAL_SOURCE, 1, count, x, weighted, n + 1,
                                    height->edges, g);
    }
    if (!status) {
        double scale = 1;
        if (tau550 > -10. / 1.e36) {
            double column = 0;
            for (size_t i = 0; i < n; ++i)
                column = column + depth[i];
            column = fmax(column, 10. / 1.e36);
            scale = (tau550 / column) * pow(550. / 340., alpha);
        }
        for (size_t w = 0; w < wavelength->cells; ++w) {
            if (!isfinite(wavelength->midpoints[w]) || wavelength->midpoints[w] <= 0) {
                status = TT_INVALID;
                break;
            }
            double spectral_scale = pow(340. / wavelength->midpoints[w], alpha);
            for (size_t i = 0; i < n; ++i) {
                size_t j = w * n + i;
                optical_depth[j] = depth[i] * scale;
                optical_depth[j] = optical_depth[j] * spectral_scale;
                single_scattering_albedo[j] = depth[i] > 0 ? ssa[i] / depth[i] : 1;
                asymmetry_factor[j] = depth[i] > 0 ? g[i] / depth[i] : 0;
            }
        }
    }
    free(memory);
    return status;
}

tt_status tt_radiator_from_dataset(const tt_dataset *data, size_t layers, size_t wavelengths,
                                   tt_radiator *out) {
    if (!data || !layers || !wavelengths || !out || wavelengths > SIZE_MAX / layers)
        return TT_INVALID;
    static const char *names[] = {"optical_depth", "single_scattering_albedo", "asymmetry_factor"};
    const tt_data_variable *variables[3];
    for (size_t i = 0; i < 3; ++i) {
        variables[i] = tt_dataset_variable(data, names[i]);
        if (!variables[i])
            return TT_NOT_FOUND;
        if (variables[i]->rank != 2 || variables[i]->shape[0] != wavelengths ||
            variables[i]->shape[1] != layers || variables[i]->count != layers * wavelengths ||
            !variables[i]->values)
            return TT_INVALID;
    }
    *out = (tt_radiator){variables[0]->values, variables[1]->values, variables[2]->values, 0};
    return TT_OK;
}
