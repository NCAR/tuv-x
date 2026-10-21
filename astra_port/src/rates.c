/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>

/* Retain the physical constants used by the reference, including its value
 * of Planck's constant, so updating SI constants cannot silently change rates. */
static const double hc = 6.626068e-34 * 2.99792458e8;

static int shape_valid(size_t nw, size_t nz) {
    return nw && nz && nz <= SIZE_MAX / nw;
}

tt_status tt_actinic_flux(size_t nw, size_t nz, const tt_radiation *radiation,
                          const double *solar_flux, double *out) {
    if (!shape_valid(nw, nz) || !radiation || !solar_flux || !out || !radiation->direct_flux ||
        !radiation->up_flux || !radiation->down_flux)
        return TT_INVALID;
    for (size_t w = 0; w < nw; ++w)
        for (size_t z = 0; z < nz; ++z) {
            size_t i = w * nz + z;
            double f =
                (radiation->direct_flux[i] + radiation->up_flux[i] + radiation->down_flux[i]) *
                solar_flux[w];
            out[i] = f < 0 ? 0 : f;
        }
    return TT_OK;
}

tt_status tt_spectral_irradiance(size_t nw, size_t nz, const tt_radiation *radiation,
                                 const double *solar_flux, const double *midpoints, double *out) {
    if (!shape_valid(nw, nz) || !radiation || !solar_flux || !midpoints || !out ||
        !radiation->direct_irradiance || !radiation->up_irradiance || !radiation->down_irradiance)
        return TT_INVALID;
    for (size_t w = 0; w < nw; ++w) {
        if (!isfinite(midpoints[w]) || midpoints[w] <= 0)
            return TT_INVALID;
        for (size_t z = 0; z < nz; ++z) {
            size_t i = w * nz + z;
            double e = (radiation->direct_irradiance[i] + radiation->up_irradiance[i] +
                        radiation->down_irradiance[i]) *
                       solar_flux[w] * hc / (midpoints[w] * 1.e-13);
            out[i] = e < 0 ? 0 : e;
        }
    }
    return TT_OK;
}

tt_status tt_photolysis_rate(size_t nw, size_t nz, const double *flux, const double *cross,
                             const double *yield, double scale, double *out) {
    if (!shape_valid(nw, nz) || !flux || !cross || !yield || !out || !isfinite(scale))
        return TT_INVALID;
    for (size_t z = 0; z < nz; ++z)
        out[z] = 0;
    for (size_t w = 0; w < nw; ++w)
        for (size_t z = 0; z < nz; ++z) {
            size_t i = w * nz + z;
            out[z] = out[z] + flux[i] * (cross[i] * yield[i]);
        }
    for (size_t z = 0; z < nz; ++z)
        out[z] = out[z] * scale;
    return TT_OK;
}

tt_status tt_heating_rate(size_t nw, size_t nz, const double *flux, const double *cross,
                          const double *yield, const double *midpoints, double threshold,
                          double scale, double *out) {
    if (!shape_valid(nw, nz) || !flux || !cross || !yield || !midpoints || !out ||
        !isfinite(threshold) || threshold <= 0 || !isfinite(scale))
        return TT_INVALID;
    for (size_t w = 0; w < nw; ++w)
        if (!isfinite(midpoints[w]) || midpoints[w] <= 0)
            return TT_INVALID;
    for (size_t z = 0; z < nz; ++z)
        out[z] = 0;
    for (size_t w = 0; w < nw; ++w) {
        double energy =
            fmax(0, hc * 1.e9 * (threshold - midpoints[w]) / (threshold * midpoints[w]));
        for (size_t z = 0; z < nz; ++z) {
            size_t i = w * nz + z;
            out[z] = out[z] + flux[i] * (energy * (cross[i] * yield[i]));
        }
    }
    for (size_t z = 0; z < nz; ++z)
        out[z] = out[z] * scale;
    return TT_OK;
}

tt_status tt_dose_rate(size_t nw, size_t nz, const double *irradiance, const double *weights,
                       double *out) {
    if (!shape_valid(nw, nz) || !irradiance || !weights || !out)
        return TT_INVALID;
    for (size_t z = 0; z < nz; ++z)
        out[z] = 0;
    for (size_t w = 0; w < nw; ++w)
        for (size_t z = 0; z < nz; ++z)
            out[z] = out[z] + irradiance[w * nz + z] * weights[w];
    return TT_OK;
}
