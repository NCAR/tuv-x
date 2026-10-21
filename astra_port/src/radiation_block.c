/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Delta-Eddington transfer for blocks of wavelengths in lock step. Each lane
 * performs exactly the arithmetic of tt_mix_radiators (one moment) and
 * tt_delta_eddington for its wavelength; the lanes of a block share the ray
 * geometry and loop structure, so the compiler can vectorize across them.
 * The solver also evaluates each exp(-slant) once per level, where the
 * generic code evaluates the same value up to four times. */
#include "model_plan.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Wavelengths per block. Lanes past the grid repeat its last wavelength. */
enum { BLOCK = 16 };

typedef struct {
    double *g, *om, *taun, *slant, *ex, *e1, *e2, *e3, *e4, *cup, *cdn, *cuptn, *cdntn, *cp, *u;
} block_work;

size_t tt_block_delta_eddington_work(size_t n) {
    return (17 * n + 2) * BLOCK;
}

static void partition(size_t n, double *work, block_work *s) {
    double **field = (double **)s;
    /* g, om, taun: n; slant, ex: n + 1; e1..cdntn: n; cp, u: 2n. */
    size_t sizes[] = {n, n, n, n + 1, n + 1, n, n, n, n, n, n, n, n, 2 * n, 2 * n};
    for (size_t i = 0; i < sizeof sizes / sizeof *sizes; ++i) {
        field[i] = work;
        work += sizes[i] * BLOCK;
    }
}

/* Mix the radiators and apply the delta scaling of one block. */
static tt_status scale_layers(size_t nr, const tt_radiator *r, size_t n, size_t w0, size_t count,
                              block_work *s) {
    for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < BLOCK; ++c) {
            size_t k = (w0 + (c < count ? c : count - 1)) * n + n - 1 - i;
            double scatter = 0.0, absorb = 0.0, g = 0.0;
            for (size_t j = 0; j < nr; ++j) {
                double od = r[j].optical_depth[k], a = r[j].single_scattering_albedo[k],
                       h = r[j].asymmetry[k];
                if (!isfinite(od) || od < 0 || !isfinite(a) || a < 0 || a > 1 || !isfinite(h) ||
                    fabs(h) > 1)
                    return TT_INVALID;
                double ds = od * a;
                scatter = scatter + ds;
                absorb = absorb + od * (1.0 - a);
                g = g + h * ds;
            }
            scatter = fmax(scatter, 1.0 / 1.e36);
            absorb = fmax(absorb, 1.0 / 1.e36);
            double tau = scatter + absorb;
            double omu = scatter == 1.0 / 1.e36 ? 1.0 / 1.e36 : scatter / tau;
            double gu = g / scatter;
            if (!isfinite(tau) || tau < 0 || !isfinite(omu) || omu < 0 || omu > 1 ||
                !isfinite(gu) || fabs(gu) >= 1)
                return TT_INVALID;
            double f = gu * gu;
            s->g[i * BLOCK + c] = (gu - f) / (1.0 - f);
            s->om[i * BLOCK + c] = (1.0 - f) * omu / (1.0 - omu * f);
            s->taun[i * BLOCK + c] = (1.0 - omu * f) * tau;
        }
    return TT_OK;
}

/* Slant optical depth at each level, as slant_depth in radiation.c. */
static void slant_depths(size_t n, const int *nid, const double *path, double mu, block_work *s) {
    for (size_t l = 0; l <= n; ++l) {
        double *slant = s->slant + l * BLOCK;
        int crossed = nid[l];
        if (l == 0 && mu >= 0) {
            for (size_t c = 0; c < BLOCK; ++c)
                slant[c] = 0.0;
            continue;
        }
        if (crossed < 0) {
            for (size_t c = 0; c < BLOCK; ++c)
                slant[c] = (double)1.e36f; /* reference literal has default real kind */
            continue;
        }
        const double *p = path + l * n;
        size_t single = (size_t)crossed < l ? (size_t)crossed : l;
        double acc[BLOCK];
        for (size_t c = 0; c < BLOCK; ++c)
            acc[c] = 0.0;
        for (size_t j = 0; j < single; ++j) {
            const double *t = s->taun + j * BLOCK;
            for (size_t c = 0; c < BLOCK; ++c)
                acc[c] = acc[c] + t[c] * p[j];
        }
        for (size_t j = single; j < (size_t)crossed; ++j) {
            const double *t = s->taun + j * BLOCK;
            for (size_t c = 0; c < BLOCK; ++c)
                acc[c] = acc[c] + 2.0 * t[c] * p[j];
        }
        memcpy(slant, acc, sizeof acc);
    }
    for (size_t l = 0; l <= n; ++l)
        for (size_t c = 0; c < BLOCK; ++c)
            s->ex[l * BLOCK + c] = exp(-s->slant[l * BLOCK + c]);
}

/* Layer coefficients e1..e4 and the particular solutions, as in radiation.c. */
static void coefficients(size_t n, const int *nid, double mu, block_work *s) {
    double tauc[BLOCK];
    for (size_t c = 0; c < BLOCK; ++c)
        tauc[c] = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double *gp = s->g + i * BLOCK, *omp = s->om + i * BLOCK, *taun = s->taun + i * BLOCK;
        const double *s0 = s->slant + i * BLOCK, *s1 = s0 + BLOCK, *x0 = s->ex + i * BLOCK,
                     *x1 = x0 + BLOCK;
        double lam[BLOCK], arg[BLOCK], expon[BLOCK], gam1[BLOCK], gam2[BLOCK], gam3[BLOCK],
            gam4[BLOCK], omi[BLOCK], m[BLOCK];
        for (size_t c = 0; c < BLOCK; ++c) {
            double previous = tauc[c];
            tauc[c] = previous + taun[c];
            double gi = copysign(fmin(fabs(gp[c]), 1.0 - 1.e-7), gp[c]);
            omi[c] = fmin(omp[c], 1.0 - 1.e-7);
            double mm = 1.0 / sqrt(1.e36);
            if (nid[i + 1] >= 0) {
                if (s1[c] == s0[c])
                    mm = sqrt(1.e36);
                else {
                    mm = (tauc[c] - previous) / (s1[c] - s0[c]);
                    mm = copysign(fmax(fabs(mm), 1.0 / sqrt(1.e36)), mm);
                }
            }
            m[c] = mm;
            gam1[c] = (7.0 - omi[c] * (4.0 + 3.0 * gi)) / 4.0;
            gam2[c] = -(1.0 - omi[c] * (4.0 - 3.0 * gi)) / 4.0;
            gam3[c] = (2.0 - 3.0 * gi * mu) / 4.0;
            gam4[c] = 1.0 - gam3[c];
            lam[c] = sqrt(gam1[c] * gam1[c] - gam2[c] * gam2[c]);
            arg[c] = -lam[c] * taun[c];
        }
        for (size_t c = 0; c < BLOCK; ++c)
            expon[c] = exp(arg[c]);
        for (size_t c = 0; c < BLOCK; ++c) {
            size_t q = i * BLOCK + c;
            double bgam = gam2[c] != 0.0 ? (gam1[c] - lam[c]) / gam2[c] : 0.0;
            s->e1[q] = 1.0 + bgam * expon[c];
            s->e2[q] = 1.0 - bgam * expon[c];
            s->e3[q] = bgam + expon[c];
            s->e4[q] = bgam - expon[c];
            double divisor = lam[c] * lam[c] - 1.0 / (m[c] * m[c]);
            divisor = copysign(fmax(1.e-3, fabs(divisor)), divisor);
            double up = omi[c] * ((gam1[c] - 1.0 / m[c]) * gam3[c] + gam4[c] * gam2[c]) / divisor;
            double dn = omi[c] * ((gam1[c] + 1.0 / m[c]) * gam4[c] + gam2[c] * gam3[c]) / divisor;
            s->cup[q] = up * x0[c];
            s->cdn[q] = dn * x0[c];
            s->cuptn[q] = up * x1[c];
            s->cdntn[q] = dn * x1[c];
        }
    }
}

/* The tridiagonal system of tt_delta_eddington, assembled row by row inside
 * the elimination of tt_tridiagonal. */
static tt_status solve(size_t n, double mu, const double *albedo, block_work *s) {
    double *cp = s->cp, *u = s->u;
    for (size_t c = 0; c < BLOCK; ++c) {
        double b = s->e1[c], d = -s->e2[c], r = -s->cdn[c];
        if (b == 0.0)
            return TT_SINGULAR;
        cp[c] = d / b;
        u[c] = r / b;
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        const double *e1 = s->e1 + i * BLOCK, *e2 = s->e2 + i * BLOCK, *e3 = s->e3 + i * BLOCK,
                     *e4 = s->e4 + i * BLOCK, *n1 = e1 + BLOCK, *n2 = e2 + BLOCK, *n3 = e3 + BLOCK,
                     *n4 = e4 + BLOCK;
        const double *cup1 = s->cup + (i + 1) * BLOCK, *cdn1 = s->cdn + (i + 1) * BLOCK,
                     *cuptn = s->cuptn + i * BLOCK, *cdntn = s->cdntn + i * BLOCK;
        size_t odd = 2 * i + 1, even = 2 * i + 2;
        double pivots[BLOCK];
        for (size_t c = 0; c < BLOCK; ++c) {
            double a = n2[c] * e1[c] - e3[c] * n4[c];
            double b = e2[c] * n2[c] - e4[c] * n4[c];
            double d = n1[c] * n4[c] - n2[c] * n3[c];
            double r = (cup1[c] - cuptn[c]) * n2[c] - (cdn1[c] - cdntn[c]) * n4[c];
            double pivot = b - a * cp[(odd - 1) * BLOCK + c];
            pivots[c] = pivot;
            double denom = 1.0 / pivot;
            cp[odd * BLOCK + c] = d * denom;
            u[odd * BLOCK + c] = (r - a * u[(odd - 1) * BLOCK + c]) * denom;
        }
        for (size_t c = 0; c < BLOCK; ++c)
            if (pivots[c] == 0.0)
                return TT_SINGULAR;
        for (size_t c = 0; c < BLOCK; ++c) {
            double a = e2[c] * e3[c] - e4[c] * e1[c];
            double b = e1[c] * n1[c] - e3[c] * n3[c];
            double d = e3[c] * n4[c] - e1[c] * n2[c];
            double r = e3[c] * (cup1[c] - cuptn[c]) + e1[c] * (cdntn[c] - cdn1[c]);
            double pivot = b - a * cp[odd * BLOCK + c];
            pivots[c] = pivot;
            double denom = 1.0 / pivot;
            cp[even * BLOCK + c] = d * denom;
            u[even * BLOCK + c] = (r - a * u[odd * BLOCK + c]) * denom;
        }
        for (size_t c = 0; c < BLOCK; ++c)
            if (pivots[c] == 0.0)
                return TT_SINGULAR;
    }
    size_t last = 2 * n - 1, i = n - 1;
    for (size_t c = 0; c < BLOCK; ++c) {
        size_t q = i * BLOCK + c;
        double rs = albedo[c];
        double a = s->e1[q] - rs * s->e3[q];
        double b = s->e2[q] - rs * s->e4[q];
        double ssfc = rs * mu * s->ex[n * BLOCK + c];
        double r = ssfc - s->cuptn[q] + rs * s->cdntn[q];
        double pivot = b - a * cp[(last - 1) * BLOCK + c];
        if (pivot == 0.0)
            return TT_SINGULAR;
        double denom = 1.0 / pivot;
        cp[last * BLOCK + c] = 0.0 * denom;
        u[last * BLOCK + c] = (r - a * u[(last - 1) * BLOCK + c]) * denom;
    }
    for (size_t row = last; row > 0; --row)
        for (size_t c = 0; c < BLOCK; ++c)
            u[(row - 1) * BLOCK + c] =
                u[(row - 1) * BLOCK + c] - cp[(row - 1) * BLOCK + c] * u[row * BLOCK + c];
    return TT_OK;
}

static void output(size_t n, size_t w0, size_t count, double mu, double scale, const block_work *s,
                   tt_radiation *out) {
    const double *u = s->u;
    for (size_t lev = 0; lev <= n; ++lev)
        for (size_t c = 0; c < count; ++c) {
            size_t k = (w0 + c) * (n + 1) + n - lev;
            double up, down;
            if (!lev) {
                down = 0.0;
                up = u[c] * s->e3[c] - u[BLOCK + c] * s->e4[c] + s->cup[c];
            } else {
                size_t i = lev - 1, row = 2 * i, q = i * BLOCK + c;
                down = u[row * BLOCK + c] * s->e3[q] + u[(row + 1) * BLOCK + c] * s->e4[q] +
                       s->cdntn[q];
                up = u[row * BLOCK + c] * s->e1[q] + u[(row + 1) * BLOCK + c] * s->e2[q] +
                     s->cuptn[q];
            }
            double direct = s->ex[lev * BLOCK + c];
            out->direct_flux[k] = direct * scale;
            out->direct_irradiance[k] = mu * direct * scale;
            out->down_irradiance[k] = down * scale;
            out->up_irradiance[k] = up * scale;
            out->down_flux[k] = down / 0.5 * scale;
            out->up_flux[k] = up / 0.5 * scale;
        }
}

tt_status tt_block_delta_eddington(size_t nr, const tt_radiator *radiators, size_t n, size_t nw,
                                   const int *nid, const double *path, double zenith,
                                   const double *surface, double scale, double *work,
                                   tt_radiation *out) {
    if ((nr && !radiators) || !n || !nw || !nid || !path || !surface || !work || !out ||
        !isfinite(zenith))
        return TT_INVALID;
    block_work s;
    partition(n, work, &s);
    double mu = cos(zenith * (3.1415926535898 / 180.0));
    for (size_t w0 = 0; w0 < nw; w0 += BLOCK) {
        size_t count = nw - w0 < BLOCK ? nw - w0 : BLOCK;
        double albedo[BLOCK];
        for (size_t c = 0; c < BLOCK; ++c) {
            albedo[c] = surface[w0 + (c < count ? c : count - 1)];
            if (!isfinite(albedo[c]) || albedo[c] < 0 || albedo[c] > 1)
                return TT_INVALID;
        }
        tt_status status = scale_layers(nr, radiators, n, w0, count, &s);
        if (status)
            return status;
        slant_depths(n, nid, path, mu, &s);
        coefficients(n, nid, mu, &s);
        status = solve(n, mu, albedo, &s);
        if (status)
            return status;
        output(n, w0, count, mu, scale, &s, out);
    }
    return TT_OK;
}
