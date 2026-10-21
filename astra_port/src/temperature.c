/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double integer_power(double x, size_t n) {
    double result = 1;
    while (n) {
        if (n & 1)
            result *= x;
        n >>= 1;
        if (n)
            x *= x;
    }
    return result;
}

tt_status tt_temperature_apply(const tt_temperature_model *p, double temperature, size_t n,
                               const double *wavelengths, double *cross_section) {
    if (!p || !n || !wavelengths || !cross_section || !isfinite(temperature) || temperature <= 0 ||
        (unsigned)p->kind > TT_TEMPERATURE_HARWOOD || p->first > n || p->count > n - p->first ||
        !p->count || !p->a || (p->range_count && !p->ranges))
        return TT_INVALID;
    if (p->kind == TT_TEMPERATURE_TAYLOR) {
        if (!p->terms || !p->sigma || p->terms > SIZE_MAX / p->count)
            return TT_INVALID;
    } else if (!p->b)
        return TT_INVALID;
    if (p->kind == TT_TEMPERATURE_EXPONENTIAL && (!p->terms || !p->powers))
        return TT_INVALID;
    size_t nr = p->range_count ? p->range_count : 1;
    for (size_t r = 0; r < nr; ++r) {
        double t = temperature;
        if (p->range_count) {
            const tt_temperature_range *range = p->ranges + r;
            if (range->minimum > range->maximum)
                return TT_INVALID;
            if (t < range->minimum || t > range->maximum)
                continue;
            if (range->is_fixed)
                t = range->fixed_value;
        }
        t = t - p->base_temperature;
        if (p->kind == TT_TEMPERATURE_EXPONENTIAL && p->invert_temperature)
            t = -t;
        if ((p->kind == TT_TEMPERATURE_HARWOOD || p->kind == TT_TEMPERATURE_BURKHOLDER) && t == 0)
            return TT_RANGE;
        double q = 1;
        if (p->kind == TT_TEMPERATURE_BURKHOLDER) {
            if (p->burkholder_b == 0)
                return TT_INVALID;
            q = 1 + exp(p->burkholder_a / (p->burkholder_b * t));
        }
        for (size_t i = 0; i < p->count; ++i) {
            size_t w = p->first + i;
            double value = 0;
            switch (p->kind) {
            case TT_TEMPERATURE_EXPONENTIAL:
                for (size_t j = 0; j < p->terms; ++j)
                    value = value + (p->a[j] + t * p->b[j]) *
                                        pow(wavelengths[w] - p->base_wavelength, p->powers[j]);
                cross_section[w] = cross_section[w] + (p->base_10 ? pow(10, value) : exp(value));
                break;
            case TT_TEMPERATURE_TAYLOR:
                value = 1;
                for (size_t j = 0; j < p->terms; ++j)
                    value = value + p->a[j * p->count + i] * integer_power(t, j + 1);
                cross_section[w] = value * p->sigma[i];
                break;
            case TT_TEMPERATURE_BURKHOLDER:
                cross_section[w] = (p->a[i] / q + p->b[i] * (1 - 1 / q)) * 1.e-20;
                break;
            case TT_TEMPERATURE_HARWOOD:
                cross_section[w] = pow(10, p->a[i] + p->b[i] / t);
                break;
            }
        }
    }
    return TT_OK;
}

tt_status tt_temperature_cross_section(const tt_temperature_model *model, size_t n,
                                       const double *wavelengths, const double *cross_section,
                                       const tt_spectral_grid *grid, double *out) {
    if (!model || !n || !wavelengths || !cross_section || !grid || !grid->wavelengths ||
        !grid->levels || !grid->edges || !grid->temperature || !out ||
        n > SIZE_MAX / sizeof(double) / 4 - 4 ||
        grid->wavelengths > SIZE_MAX / sizeof(double) / 2 ||
        grid->levels > SIZE_MAX / grid->wavelengths)
        return TT_INVALID;
    size_t np = n + 4, nw = grid->wavelengths;
    double *work = calloc(2 * np + nw, sizeof(double));
    if (!work)
        return TT_MEMORY;
    double *x = work, *y = x + np, *bins = y + np;
    /* TUV-x declares this padding offset as a default-real literal. */
    double delta = (double)1.e-5f;
    x[0] = 0;
    x[1] = (1 - delta) * wavelengths[0];
    memcpy(x + 2, wavelengths, n * sizeof(double));
    x[n + 2] = (1 + delta) * wavelengths[n - 1];
    x[n + 3] = 1.e38;
    tt_status status = TT_OK;
    for (size_t z = 0; z < grid->levels; ++z) {
        memcpy(y + 2, cross_section, n * sizeof(double));
        status = tt_temperature_apply(model, grid->temperature[z], n, wavelengths, y + 2);
        if (status)
            break;
        status = tt_interpolate(TT_CONSERVING, 0, np, x, y, nw + 1, grid->edges, bins);
        if (status)
            break;
        for (size_t w = 0; w < nw; ++w)
            out[w * grid->levels + z] = bins[w];
    }
    free(work);
    return status;
}
