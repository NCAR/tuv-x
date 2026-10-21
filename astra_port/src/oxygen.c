/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "generated/oxygen.inc"

static const double lyman_edges[] = {121.4, 121.9};
static const double schumann_edges[] = {175.4, 177,   178.6, 180.2, 181.8, 183.5,
                                        185.2, 186.9, 188.7, 190.5, 192.3, 194.2,
                                        196.1, 198,   200,   202,   204.1, 206.2};
static const double limiting_cross_section[17] = {
    6.2180730e-21, 5.8473627e-22, 5.6996334e-22, 4.5627094e-22, 1.7668250e-22, 1.1178808e-22,
    1.2040544e-22, 4.0994668e-23, 1.8450616e-23, 1.5639540e-23, 8.7961075e-24, 7.6475608e-24,
    7.6260556e-24, 7.5565696e-24, 7.6334338e-24, 7.4371992e-24, 7.3642966e-24};

static tt_status locate_band(size_t nw, const double *edges, size_t n, const double *band,
                             size_t *index) {
    *index = SIZE_MAX;
    if (edges[0] > band[0] || edges[nw] < band[n - 1])
        return TT_OK;
    for (size_t w = 0; w <= nw; ++w)
        if (fabs(edges[w] - band[0]) < 10 * 1.e-7) {
            if (n - 1 > nw - w)
                return TT_RANGE;
            for (size_t i = 1; i < n; ++i)
                if (fabs(edges[w + i] - band[i]) > 10 * 1.e-7)
                    return TT_RANGE;
            *index = w;
            return TT_OK;
        }
    return TT_RANGE;
}

tt_status tt_oxygen_band_grid(size_t nw, const double *edges, tt_oxygen_bands *bands) {
    if (!nw || !edges || !bands)
        return TT_INVALID;
    for (size_t w = 0; w <= nw; ++w)
        if (!isfinite(edges[w]) || edges[w] < 0 || (w && edges[w] <= edges[w - 1]))
            return TT_INVALID;
    tt_status status = locate_band(nw, edges, 2, lyman_edges, &bands->lyman_alpha);
    return status ? status : locate_band(nw, edges, 18, schumann_edges, &bands->schumann_runge);
}

static double chebyshev(const double *coefficients, double x) {
    double di = 0, di1 = 0, y = (2 * x - (38 + 56)) / (56 - 38), y2 = 2 * y;
    for (size_t i = 20; i-- > 1;) {
        double previous = di;
        di = y2 * di - di1 + coefficients[i];
        di1 = previous;
    }
    return y * di - di1 + .5 * coefficients[0];
}

static void lyman_alpha(size_t nz, const double *column, const double *secchi, double *rm,
                        double *depth, double *cross) {
    const double b[] = {6.8431e-1, 2.29841e-1, 8.65412e-2};
    const double c[] = {8.22114e-21, 1.77556e-20, 8.22112e-21};
    const double d[] = {6.0073e-21, 4.28569e-21, 1.28059e-20};
    const double e[] = {8.21666e-21, 1.63296e-20, 4.85121e-17};
    for (size_t z = 0; z < nz; ++z) {
        double ro2 = 0;
        rm[z] = 0;
        for (size_t i = 0; i < 3; ++i) {
            double sigma = c[i] * column[z];
            double tau = sigma < 100.e8 ? exp(-sigma) : 0;
            rm[z] = rm[z] + tau * b[i];
            if (cross) {
                sigma = e[i] * column[z];
                tau = sigma < 100.e8 ? exp(-sigma) : 0;
                ro2 = ro2 + tau * d[i];
            }
        }
        if (cross)
            cross[z] = rm[z] > 1.e-100 && (z + 1 == nz || ro2 > 1.e-100) ? ro2 / rm[z] : 1.e-20;
    }
    if (depth)
        for (size_t z = 0; z + 1 < nz; ++z)
            depth[z] = rm[z] > 1.e-100 && rm[z + 1] > 0
                           ? log(rm[z + 1]) / secchi[z + 1] - log(rm[z]) / secchi[z]
                           : 1000;
}

static void schumann_runge(size_t nz, const double *column, const double *temperature,
                           const double *secchi, double *clamped, double *cross, double *depth) {
    double colmin = exp(38), norm = 1 / (double)(nz - 1);
    size_t top = nz, bottom = 0;
    for (size_t z = 0; z < nz; ++z) {
        clamped[z] = fmax(column[z], colmin);
        if (column[z] < colmin) {
            if (z < top)
                top = z;
        } else if (log(clamped[z]) > 56)
            bottom = z + 1;
        for (size_t w = 0; w < 17; ++w) {
            cross[w * nz + z] = limiting_cross_section[w];
            if (column[z] >= colmin) {
                double x = log(clamped[z]);
                if (x <= 56)
                    cross[w * nz + z] = exp(chebyshev(oxygen_a[w], x) * (temperature[z] - 220) +
                                            chebyshev(oxygen_b[w], x));
            }
        }
    }
    if (bottom == nz)
        bottom = nz - 1;
    for (size_t w = 0; w < 17; ++w) {
        double *xs = cross + w * nz;
        for (size_t z = 0; z < bottom; ++z)
            xs[z] = xs[bottom];
        /* TUV-x indexes zero when every column is below exp(38). Use the
         * already initialized limiting cross section in that undefined case. */
        if (top)
            for (size_t z = top; z < nz; ++z)
                xs[z] = xs[top - 1];
        if (depth)
            for (size_t z = 0; z + 1 < nz; ++z) {
                double value;
                if (fabs(1 - clamped[z + 1] / clamped[z]) <= 2 * 1.e-7)
                    value = xs[z + 1] * clamped[z + 1] * norm;
                else {
                    value = fabs((xs[z + 1] * clamped[z + 1] - xs[z] * clamped[z]) /
                                 (1 + log(xs[z + 1] / xs[z]) / log(clamped[z + 1] / clamped[z])));
                    value = 2 * value / (secchi[z] + secchi[z + 1]);
                }
                depth[w * (nz - 1) + z] = value;
            }
    }
}

tt_status tt_oxygen_optics(const tt_oxygen_bands *bands, size_t nw, size_t nz,
                           const double *temperature, const double *vertical, const double *slant,
                           const double *oxygen, double *depth, double *cross) {
    if (!bands || !nw || nz < 2 || !temperature || !vertical || !slant || !oxygen ||
        (!depth && !cross) || nz > SIZE_MAX / sizeof(double) / 20 || nw > SIZE_MAX / nz)
        return TT_INVALID;
    size_t la = bands->lyman_alpha, sr = bands->schumann_runge;
    if ((la != SIZE_MAX && la >= nw) || (sr != SIZE_MAX && (sr >= nw || nw - sr < 17)))
        return TT_INVALID;
    if (la == SIZE_MAX && sr == SIZE_MAX)
        return TT_OK;
    for (size_t z = 0; z < nz; ++z)
        if (!isfinite(temperature[z]) || temperature[z] <= 0 || !isfinite(slant[z]) ||
            slant[z] < 0 || !isfinite(oxygen[z]) || oxygen[z] < 0 ||
            (z + 1 < nz && (!isfinite(vertical[z]) || vertical[z] <= 0)))
            return TT_INVALID;
    double *work = malloc(20 * nz * sizeof(double));
    if (!work)
        return TT_MEMORY;
    double *secchi = work, *buffer = secchi + nz, *xs = buffer + nz;
    for (size_t z = 0; z + 1 < nz; ++z)
        secchi[z] = slant[z] <= .1 * 1.e36 ? slant[z] / vertical[z] : 2;
    secchi[nz - 1] = secchi[nz - 2];
    if (depth)
        for (size_t z = 0; z < nz; ++z)
            if (secchi[z] <= 0) {
                free(work);
                return TT_RANGE;
            }
    if (la != SIZE_MAX)
        lyman_alpha(nz, oxygen, secchi, buffer, depth ? depth + la * (nz - 1) : NULL,
                    cross ? cross + la * nz : NULL);
    if (sr != SIZE_MAX) {
        double *target = cross ? cross + sr * nz : xs;
        schumann_runge(nz, oxygen, temperature, secchi, buffer, target,
                       depth ? depth + sr * (nz - 1) : NULL);
    }
    free(work);
    return TT_OK;
}
