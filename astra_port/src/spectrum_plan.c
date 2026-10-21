/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Prepared spectra. Each kind below restates the generic kernel in spectra.c,
 * yields.c or temperature.c for one wavelength at many levels. Subexpressions
 * of wavelength alone are computed when the plan is built, those of level
 * alone once per evaluation, and the remaining operations keep the generic
 * order, so each value is identical to the generic kernel's. Wavelengths are
 * classified as constant only where the kernel's branches or data make the
 * value independent of temperature and air density; the generic kernel
 * evaluates those once, and a sweep of reference columns guards the
 * classification. */
#include "spectrum_plan.h"
#include "spectral_validation.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Constant classifications that rely on a finite exponential hold for all
 * temperatures up to this limit; hotter columns use the generic path. */
#define TEMPERATURE_LIMIT 1.0e4
/* exp and pow(10, x) stay finite below these arguments. */
#define EXP_LIMIT 700.0
#define POW10_LIMIT 300.0

struct tt_spectrum_plan {
    const tt_spectrum *s;
    int cs, qy, kind, generic;
    size_t nw, count, stride, tables, required;
    double *constant, *coef, maximum_temperature;
    size_t *index;
    /* Parameterized cross sections: padded raw grid and values, powers of
     * the exponential form, varying padded points and target bins. */
    size_t raw, first_bin, end_bin;
    double *x, *y, *powers;
};

static double clamp(double t, double lo, double hi) {
    return fmax(lo, fmin(t, hi));
}
static double sq(double x) {
    return x * x;
}
/* Both generic files use this repeated-squaring power. */
static double integer_power(double x, size_t n) {
    double y = 1;
    while (n) {
        if (n & 1)
            y *= x;
        n >>= 1;
        if (n)
            x *= x;
    }
    return y;
}
static double parameter(const tt_spectral_table *p, size_t w, size_t c) {
    return p->values[c * p->rows + w];
}
/* Temperature-table interpolation index and fraction, as in spectra.c. */
static void bracket(const tt_spectral_table *q, double t, size_t *j, double *f) {
    size_t k = 0;
    while (k + 2 < q->temperature_count && t > q->temperatures[k + 1])
        ++k;
    *j = k;
    *f = (t - q->temperatures[k]) / (q->temperatures[k + 1] - q->temperatures[k]);
}
/* True when every interpolation column at w holds the same finite value, so
 * interpolation returns that value for any fraction. */
static int flat(const tt_spectral_table *q, size_t w) {
    double v = parameter(q, w, 0);
    if (!isfinite(v))
        return 0;
    for (size_t c = 1; c < q->temperature_count; ++c)
        if (parameter(q, w, c) != v || !isfinite(parameter(q, w, c)))
            return 0;
    return 1;
}
/* Largest |T - offset| for temperatures in (0, TEMPERATURE_LIMIT]. */
static double offset_bound(double offset) {
    return fmax(fabs(offset), fabs(TEMPERATURE_LIMIT - offset));
}

/* ---- Cross sections ---------------------------------------------------- */

enum { CS_STRIDE = 4 };

/* Fill coefficients for wavelength w and return whether its value varies. */
static int cs_prepare(tt_spectrum_plan *p, size_t w, double *c) {
    const tt_spectrum *s = p->s;
    const tt_spectral_table *t = s->tables;
    double lambda = s->wavelength.midpoints[w], edge = s->wavelength.edges[w];
    double base = p->required && p->kind != TT_CS_HCFC ? parameter(t, w, 0) : 0, coeff;
    switch (p->kind) {
    case TT_CS_ACETONE:
    case TT_CS_CH2O:
    case TT_CS_CLONO2:
        c[0] = base;
        for (size_t k = 1; k < (p->kind == TT_CS_ACETONE ? 4u
                                : p->kind == TT_CS_CH2O  ? 2u
                                                         : 3u);
             ++k)
            c[k] = parameter(t, w, k);
        return 1;
    case TT_CS_CCL4:
    case TT_CS_CHCL3:
        if (p->kind == TT_CS_CCL4 ? !(lambda > 194 && lambda < 250)
                                  : !(lambda > 190 && lambda < 240))
            return 0;
        coeff =
            p->kind == TT_CS_CCL4
                ? 1.0739 +
                      lambda * (-1.6275e-2 +
                                lambda * (8.8141e-5 + lambda * (-1.9811e-7 + 1.5022e-10 * lambda)))
                : 3.7973 +
                      lambda * (-7.0913e-2 +
                                lambda * (4.9397e-4 + lambda * (-1.5226e-6 + lambda * 1.7555e-9)));
        /* t is clamped to [-85, 5]: a zero base stays zero. */
        if (base == 0 && fabs(coeff) * 85 < POW10_LIMIT)
            return 0;
        c[0] = base;
        c[1] = coeff;
        return 1;
    case TT_CS_CFC11:
        c[0] = base;
        c[1] = lambda - 184.9;
        return !(base == 0 && fabs(c[1]) * 1.e-4 * offset_bound(298) < EXP_LIMIT);
    case TT_CS_CH3ONO2:
    case TT_CS_HNO3:
    case TT_CS_RONO2:
        c[0] = base;
        c[1] = parameter(t, w, 1);
        if (c[1] == 0 || (base == 0 && fabs(c[1]) * offset_bound(298) < EXP_LIMIT))
            return 0;
        return 1;
    case TT_CS_CHBR3:
        if (!(lambda > 290 && lambda < 340))
            return 0;
        c[0] = .06183 - .000241 * lambda;
        c[1] = 2.376 + .14757 * lambda;
        c[2] = base;
        return 1;
    case TT_CS_CL2:
        c[0] = sq(log(329.5 / lambda));
        c[1] = sq(log(406.5 / lambda));
        return 1;
    case TT_CS_H2O2:
        if (!(edge >= 260 && edge < 350))
            return 0;
        c[0] = ((((((1.5534675e-13 * lambda - 2.652014e-10) * lambda + 1.6878206e-7) * lambda -
                   4.035101e-5) *
                      lambda -
                  4.4589016e-3) *
                     lambda +
                 4.535649) *
                    lambda -
                9.2170972e2) *
                   lambda +
               6.4761e4;
        c[1] = (((-1.0924e-7 * lambda - 3.0493e-5) * lambda + 1.1522e-1) * lambda - 5.1351e1) *
                   lambda +
               6.8123e3;
        return 1;
    case TT_CS_HCFC:
        if (!(lambda >= 190 && lambda <= t[0].temperatures[1]))
            return 0;
        for (size_t k = 0; k < t[0].rows; ++k)
            c[k] = integer_power(lambda - 206.214, k);
        return 1;
    case TT_CS_N2O:
        if (!(lambda >= 173 && lambda <= 240))
            return 0;
        c[0] = (((2.520672e-7 * lambda - 1.777846e-4) * lambda + 4.301146e-2) * lambda - 4.071805) *
                   lambda +
               68.21023;
        c[1] = exp(((-1.881058e-5 * lambda + 1.111572e-2) * lambda - 2.116255) * lambda + 123.4014);
        return 1;
    case TT_CS_N2O5:
        c[0] = base;
        c[1] = 1000 * parameter(t + 1, w, 0);
        /* (300 - t) / (300 t) lies in [0, 1/1043] for t clamped to [233, 300]. */
        return !(c[1] == 0 || (base == 0 && fabs(c[1]) / 1043 < POW10_LIMIT));
    case TT_CS_NO2:
    case TT_CS_TINT:
        for (size_t k = 0; k < p->tables; ++k)
            if (!flat(t + k, w))
                return 1;
        return 0;
    case TT_CS_O3_JPL06:
        if (!(edge >= 196.078 && s->wavelength.edges[w + 1] <= 342.5))
            return 0;
        c[0] = parameter(t + 1, w, 0);
        c[1] = (parameter(t + 2, w, 0) - c[0]) / (298 - 218);
        c[2] = parameter(t + 2, w, 0);
        return 1;
    case TT_CS_O3: {
        double v185 = 185 * tt_air_refractive_index(185, 2.45e19);
        double v195 = 195 * tt_air_refractive_index(195, 2.45e19);
        double v345 = 345 * tt_air_refractive_index(345, 2.45e19);
        size_t table = edge < v185 ? 2 : edge < v195 ? 3 : edge < v345 ? 1 : 0;
        c[0] = (double)table;
        return table && !flat(t + table, w);
    }
    case TT_CS_OCLO:
        return 1;
    default:
        return 0;
    }
}

/* Level quantities shared by all wavelengths, then rows for each varying one. */
static void cs_evaluate(const tt_spectrum_plan *p, size_t nz, const double *temp, double *z,
                        double *rows) {
    const tt_spectral_table *t = p->s->tables;
    size_t nt = p->tables;
    switch (p->kind) {
    case TT_CS_ACETONE:
        for (size_t i = 0; i < nz; ++i)
            z[i] = clamp(temp[i], 235, 298);
        break;
    case TT_CS_CCL4:
    case TT_CS_CHCL3:
        for (size_t i = 0; i < nz; ++i)
            z[i] = clamp(temp[i], 210, 300) - 295;
        break;
    case TT_CS_CFC11:
        for (size_t i = 0; i < nz; ++i)
            z[i] = 1.e-4 * (temp[i] - 298);
        break;
    case TT_CS_CH2O:
    case TT_CS_CH3ONO2:
    case TT_CS_HNO3:
    case TT_CS_RONO2:
        for (size_t i = 0; i < nz; ++i)
            z[i] = temp[i] - 298;
        break;
    case TT_CS_CHBR3:
        for (size_t i = 0; i < nz; ++i)
            z[i] = 273 - temp[i];
        break;
    case TT_CS_CL2:
        for (size_t i = 0; i < nz; ++i) {
            double bb = exp(402.7 / temp[i]), bbsq = bb * bb, alpha = (bbsq - 1) / (bbsq + 1);
            z[i] = -99.0 * alpha;
            z[nz + i] = -91.5 * alpha;
            z[2 * nz + i] = 1.e-20 * sqrt(alpha);
        }
        break;
    case TT_CS_CLONO2:
        for (size_t i = 0; i < nz; ++i)
            z[i] = temp[i] - 296;
        break;
    case TT_CS_H2O2:
        for (size_t i = 0; i < nz; ++i) {
            z[i] = 1 / (1 + exp(-1265 / clamp(temp[i], 200, 400)));
            z[nz + i] = 1 - z[i];
        }
        break;
    case TT_CS_HCFC:
        /* Per level, each polynomial coefficient of the exponent. */
        for (size_t i = 0; i < nz; ++i) {
            double u = clamp(temp[i], 203, 295) - t[0].temperatures[0];
            for (size_t k = 0; k < t[0].rows; ++k)
                z[k * nz + i] =
                    parameter(t, k, 0) + u * (parameter(t, k, 1) + u * parameter(t, k, 2));
        }
        break;
    case TT_CS_N2O:
        for (size_t i = 0; i < nz; ++i)
            z[i] = clamp(temp[i], 194, 320) - 300;
        break;
    case TT_CS_N2O5:
        for (size_t i = 0; i < nz; ++i) {
            double u = clamp(temp[i], 233, 300);
            z[i] = 300 - u;
            z[nz + i] = 300 * u;
        }
        break;
    case TT_CS_NO2:
    case TT_CS_TINT:
        for (size_t k = 0; k < nt; ++k)
            for (size_t i = 0; i < nz; ++i) {
                const tt_spectral_table *q = t + k;
                size_t j;
                bracket(
                    q,
                    clamp(temp[i], q->temperatures[0], q->temperatures[q->temperature_count - 1]),
                    &j, z + (2 * k + 1) * nz + i);
                z[2 * k * nz + i] = (double)j;
            }
        break;
    case TT_CS_O3:
        for (size_t k = 1; k < 4; ++k)
            for (size_t i = 0; i < nz; ++i) {
                const tt_spectral_table *q = t + k;
                size_t j;
                bracket(
                    q,
                    clamp(temp[i], q->temperatures[0], q->temperatures[q->temperature_count - 1]),
                    &j, z + (2 * k + 1) * nz + i);
                z[2 * k * nz + i] = (double)j;
            }
        break;
    case TT_CS_OCLO:
        for (size_t i = 0; i < nz; ++i) {
            double u = temp[i];
            size_t k = 0;
            if (u <= t[0].temperatures[0])
                z[i] = -1;
            else if (u >= t[nt - 1].temperatures[0])
                z[i] = -2;
            else {
                while (k + 2 < nt && t[k + 1].temperatures[0] <= u)
                    ++k;
                z[i] = (double)k;
                z[nz + i] =
                    (u - t[k].temperatures[0]) / (t[k + 1].temperatures[0] - t[k].temperatures[0]);
            }
        }
        break;
    default:
        break;
    }
    for (size_t k = 0; k < p->count; ++k) {
        size_t w = p->index[k];
        const double *c = p->coef + k * p->stride;
        double *out = rows + k * nz;
        switch (p->kind) {
        case TT_CS_ACETONE:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * (1.0 + z[i] * (c[1] + z[i] * (c[2] + z[i] * c[3])));
            break;
        case TT_CS_CCL4:
        case TT_CS_CHCL3:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * pow(10, c[1] * z[i]);
            break;
        case TT_CS_CFC11:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * exp(c[1] * z[i]);
            break;
        case TT_CS_CH2O:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] + c[1] * z[i];
            break;
        case TT_CS_CH3ONO2:
        case TT_CS_HNO3:
        case TT_CS_RONO2:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * exp(c[1] * z[i]);
            break;
        case TT_CS_CHBR3:
            for (size_t i = 0; i < nz; ++i)
                out[i] = temp[i] > 210 && temp[i] < 300 ? exp(c[0] * z[i] - c[1]) : c[2];
            break;
        case TT_CS_CL2:
            for (size_t i = 0; i < nz; ++i) {
                double ex1 = 27.3 * exp(z[i] * c[0]), ex2 = .932 * exp(z[nz + i] * c[1]);
                out[i] = z[2 * nz + i] * (ex1 + ex2);
            }
            break;
        case TT_CS_CLONO2:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * (1 + z[i] * (c[1] + z[i] * c[2]));
            break;
        case TT_CS_H2O2:
            for (size_t i = 0; i < nz; ++i)
                out[i] = (z[i] * c[0] + z[nz + i] * c[1]) * 1.e-21;
            break;
        case TT_CS_HCFC:
            for (size_t i = 0; i < nz; ++i) {
                double sigma = 0;
                for (size_t j = 0; j < t[0].rows; ++j)
                    sigma = sigma + z[j * nz + i] * c[j];
                out[i] = exp(sigma);
            }
            break;
        case TT_CS_N2O:
            for (size_t i = 0; i < nz; ++i)
                out[i] = exp(c[0] + z[i] * c[1]);
            break;
        case TT_CS_N2O5:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] * pow(10, c[1] * z[i] / z[nz + i]);
            break;
        case TT_CS_NO2:
        case TT_CS_TINT:
            for (size_t i = 0; i < nz; ++i) {
                double value = 0;
                for (size_t q = 0; q < nt; ++q) {
                    size_t j = (size_t)z[2 * q * nz + i];
                    double f = z[(2 * q + 1) * nz + i];
                    value = value + parameter(t + q, w, j) +
                            f * (parameter(t + q, w, j + 1) - parameter(t + q, w, j));
                }
                out[i] = value;
            }
            break;
        case TT_CS_O3_JPL06:
            for (size_t i = 0; i < nz; ++i)
                out[i] = temp[i] < 218    ? c[0]
                         : temp[i] <= 298 ? c[0] + c[1] * (temp[i] - 218)
                                          : c[2];
            break;
        case TT_CS_O3: {
            size_t q = (size_t)c[0];
            for (size_t i = 0; i < nz; ++i) {
                size_t j = (size_t)z[2 * q * nz + i];
                double f = z[(2 * q + 1) * nz + i];
                out[i] = parameter(t + q, w, j) +
                         f * (parameter(t + q, w, j + 1) - parameter(t + q, w, j));
            }
            break;
        }
        case TT_CS_OCLO:
            for (size_t i = 0; i < nz; ++i) {
                if (z[i] == -1)
                    out[i] = parameter(t, w, 0);
                else if (z[i] == -2)
                    out[i] = parameter(t + nt - 1, w, 0);
                else {
                    size_t q = (size_t)z[i];
                    out[i] = parameter(t + q, w, 0) +
                             z[nz + i] * (parameter(t + q + 1, w, 0) - parameter(t + q, w, 0));
                }
            }
            break;
        default:
            break;
        }
    }
}

/* ---- Quantum yields ---------------------------------------------------- */

enum { QY_STRIDE = 5 };

static int qy_prepare(tt_spectrum_plan *p, size_t w, double *c) {
    const tt_spectrum *s = p->s;
    const tt_spectral_table *t = s->tables;
    const tt_yield_options *o = &s->yield;
    double lambda = s->wavelength.midpoints[w];
    double base = p->tables ? parameter(t, w, 0) : 0;
    switch (p->kind) {
    case TT_QY_TINT:
    case TT_QY_NO2:
        for (size_t k = 0; k < p->tables; ++k)
            if (!flat(t + k, w))
                return 1;
        return 0;
    case TT_QY_O3_O1D:
    case TT_QY_O3_O3P:
        if (lambda <= 305 || lambda > 328)
            return 0;
        c[0] = exp(-sq(sq((304.225 - lambda) / 5.576)));
        c[1] = exp(-sq((314.957 - lambda) / 6.601));
        c[2] = exp(-sq((310.737 - lambda) / 2.187));
        return 1;
    case TT_QY_NO3_AQ:
    case TT_QY_CH2CHCHO:
    case TT_QY_MEK:
        return 1;
    case TT_QY_CH2O: {
        double q1 = base, q2 = parameter(t, w, 1), qt = 1 - q1;
        if (!(lambda >= 330 && lambda < 360 && q2 > 0))
            return 0;
        c[0] = (1 - (q1 + q2)) / (2.45e19 * q2 * qt);
        c[1] = .05 * (lambda - 329);
        c[2] = 1 / qt;
        return 1;
    }
    case TT_QY_CH3CHO: {
        double q1 = parameter(t, w, 1), q2 = 1 - base, q = q1 > 0 ? q2 / q1 - 1 : 0;
        c[0] = q1 * (1 + q);
        c[1] = q;
        /* With q zero the denominator is exactly one for any air density. */
        return q != 0;
    }
    case TT_QY_C2H5CHO:
        if (base < 10.0 / 1.e36)
            return 0;
        c[0] = 1 / base - 1;
        return c[0] != 0;
    case TT_QY_MVK:
        c[0] = exp(-.055 * (lambda - 308));
        return 1;
    case TT_QY_ACETONE:
        if (lambda < 279 || lambda > 327)
            return 0;
        c[0] = lambda - 248;
        c[1] = 1.e7 / lambda;
        c[2] = c[1] - 33113;
        c[3] = c[1] - 30488;
        c[4] = lambda < 302;
        return 1;
    case TT_QY_CH3COCHO: {
        double phi = clamp(1 - (lambda - 380) / 60, 0, 1);
        if (!(phi > 0 && lambda >= 380 && lambda <= 440))
            return 0;
        c[0] = phi;
        c[1] = 1.36e8 * exp(-8793 / lambda);
        return 1;
    }
    case TT_QY_H2SO4:
        if (!(o->collision_rates[w] >= 0))
            return 0;
        c[0] = o->collision_rates[w];
        return 1;
    default:
        return 0;
    }
}

static void qy_evaluate(const tt_spectrum_plan *p, size_t nz, const double *temp, const double *m,
                        double *z, double *rows) {
    const tt_spectral_table *t = p->s->tables;
    const tt_yield_options *o = &p->s->yield;
    size_t nt = p->tables;
    switch (p->kind) {
    case TT_QY_TINT:
    case TT_QY_NO2:
        for (size_t k = 0; k < nt; ++k)
            for (size_t i = 0; i < nz; ++i) {
                const tt_spectral_table *q = t + k;
                double u = p->kind == TT_QY_TINT ? clamp(temp[i], q->temperatures[0],
                                                         q->temperatures[q->temperature_count - 1])
                                                 : temp[i];
                size_t j;
                bracket(q, u, &j, z + (2 * k + 1) * nz + i);
                z[2 * k * nz + i] = (double)j;
            }
        break;
    case TT_QY_O3_O1D:
    case TT_QY_O3_O3P:
        for (size_t i = 0; i < nz; ++i) {
            double kt = .695 * temp[i], q2 = exp(-825.518 / kt), qfac1 = 1 / (1 + q2),
                   qfac2 = q2 / (1 + q2), t300 = temp[i] / 300;
            z[i] = .8036 * qfac1;
            z[nz + i] = 8.9061 * t300 * t300 * qfac2;
            z[2 * nz + i] = .1192 * pow(t300, 1.5);
        }
        break;
    case TT_QY_NO3_AQ:
        for (size_t i = 0; i < nz; ++i)
            z[i] = exp(-2400 / temp[i] + 3.6);
        break;
    case TT_QY_CH2O:
        for (size_t i = 0; i < nz; ++i)
            z[i] = (300 - temp[i]) / 80;
        break;
    case TT_QY_CH3CHO:
        for (size_t i = 0; i < nz; ++i)
            z[i] = m[i] / 2.465e19;
        break;
    case TT_QY_C2H5CHO:
        for (size_t i = 0; i < nz; ++i)
            z[i] = m[i] / 2.45e19;
        break;
    case TT_QY_CH2CHCHO:
        for (size_t i = 0; i < nz; ++i) {
            if (m[i] > 2.6e19)
                z[i] = .004;
            else {
                double phi = m[i] <= 8.e17 ? .086 + 1.613e-17 * 8.e17 : .086 + 1.613e-17 * m[i];
                z[i] = .004 + 1 / phi;
            }
        }
        break;
    case TT_QY_MVK:
        for (size_t i = 0; i < nz; ++i)
            z[i] = 5.5 + 9.2e-19 * m[i];
        break;
    case TT_QY_ACETONE:
        for (size_t i = 0; i < nz; ++i) {
            double u = clamp(temp[i], o->minimum_temperature, o->maximum_temperature) / 295;
            double a0 = .350 * pow(u, -1.28);
            z[i] = a0;
            z[nz + i] = 1 - a0;
            z[2 * nz + i] = .068 * pow(u, -2.65);
            z[3 * nz + i] = 1.600e-19 * pow(u, -2.38);
            z[4 * nz + i] = .55e-3 * pow(u, -3.19);
            z[5 * nz + i] = 1.62e-17 * pow(u, -10.03);
            z[6 * nz + i] = 1.79e-3 * pow(u, -1.364);
            z[7 * nz + i] = 26.29 * pow(u, -6.59);
            z[8 * nz + i] = 5.72e-7 * pow(u, -2.93);
            z[9 * nz + i] = 30006 * pow(u, -.064);
            z[10 * nz + i] = 1.67e-15 * pow(u, -7.25);
            z[11 * nz + i] = 2.08e-3 * pow(u, -1.16);
        }
        break;
    case TT_QY_MEK:
        for (size_t i = 0; i < nz; ++i) {
            double ptorr = 760 * m[i] / 2.69e19;
            z[i] = fmin(1 / (.96 + 2.22e-3 * ptorr), 1);
        }
        break;
    case TT_QY_CH3COCHO:
        for (size_t i = 0; i < nz; ++i)
            z[i] = m[i] * 760 / 2.456e19;
        break;
    case TT_QY_H2SO4: {
        double pi = 3.1415926535898;
        for (size_t i = 0; i < nz; ++i)
            z[i] = sqrt(16.0 * 8.31446261815324 * temp[i] / (pi * o->molecular_weight)) *
                   (pi * sq(o->molecular_diameter) * m[i] * 1.e6);
        break;
    }
    default:
        break;
    }
    for (size_t k = 0; k < p->count; ++k) {
        size_t w = p->index[k];
        const double *c = p->coef + k * p->stride;
        double *out = rows + k * nz;
        switch (p->kind) {
        case TT_QY_TINT:
        case TT_QY_NO2:
            for (size_t i = 0; i < nz; ++i) {
                double value = 0;
                for (size_t q = 0; q < nt; ++q) {
                    size_t j = (size_t)z[2 * q * nz + i];
                    double f = z[(2 * q + 1) * nz + i];
                    value = value + parameter(t + q, w, j) +
                            f * (parameter(t + q, w, j + 1) - parameter(t + q, w, j));
                }
                out[i] = p->kind == TT_QY_NO2 ? fmax(value, 0) : value;
            }
            break;
        case TT_QY_O3_O1D:
        case TT_QY_O3_O3P:
            for (size_t i = 0; i < nz; ++i) {
                double ozone = .0765 + z[i] * c[0] + z[nz + i] * c[1] + z[2 * nz + i] * c[2];
                out[i] = p->kind == TT_QY_O3_O1D ? ozone : 1 - ozone;
            }
            break;
        case TT_QY_NO3_AQ:
        case TT_QY_CH2CHCHO:
        case TT_QY_MEK:
            memcpy(out, z, nz * sizeof(double));
            break;
        case TT_QY_CH2O:
            for (size_t i = 0; i < nz; ++i) {
                double q = c[0] * (1 + c[1] * z[i]);
                out[i] = 1 / (c[2] + q * m[i]);
            }
            break;
        case TT_QY_CH3CHO:
            for (size_t i = 0; i < nz; ++i)
                out[i] = clamp(c[0] / (1 + c[1] * z[i]), 0, 1);
            break;
        case TT_QY_C2H5CHO:
            for (size_t i = 0; i < nz; ++i)
                out[i] = fmin(1, 1 / (1 + c[0] * z[i]));
            break;
        case TT_QY_MVK:
            for (size_t i = 0; i < nz; ++i)
                out[i] = fmin(c[0] / z[i], 1);
            break;
        case TT_QY_ACETONE:
            for (size_t i = 0; i < nz; ++i) {
                double arg = z[2 * nz + i] * c[0];
                double ca0 = arg > 80 ? 5.e34 : exp(arg) * z[i] / z[nz + i], fco = 1 / (1 + ca0),
                       fac;
                if (c[4] != 0) {
                    double ca1 = z[3 * nz + i] * exp(-z[4 * nz + i] * c[2]);
                    fac = (1 - fco) / (1 + ca1 * m[i]);
                } else {
                    double ca2 = z[5 * nz + i] * exp(-z[6 * nz + i] * c[3]);
                    double ca3 = z[7 * nz + i] * exp(-z[8 * nz + i] * sq(c[1] - z[9 * nz + i]));
                    double ca4 = z[10 * nz + i] * exp(-z[11 * nz + i] * c[3]);
                    fac = (1 - fco) * (1 + ca3 + ca4 * m[i]) /
                          ((1 + ca3 + ca2 * m[i]) * (1 + ca4 * m[i]));
                }
                double qy = 0;
                if (o->acetone_branches & 1)
                    qy = qy + fco;
                if (o->acetone_branches & 2)
                    qy = qy + fac;
                out[i] = qy;
            }
            break;
        case TT_QY_CH3COCHO:
            for (size_t i = 0; i < nz; ++i)
                out[i] = c[0] / (c[0] + c[1] * z[i]);
            break;
        case TT_QY_H2SO4:
            for (size_t i = 0; i < nz; ++i)
                out[i] = i + 1 == nz ? 1 : c[0] / (c[0] + z[i]);
            break;
        default:
            break;
        }
    }
}

/* ---- Temperature-parameterized cross sections -------------------------- */

/* The parameterization as in temperature.c, on the padded raw values y. */
static tt_status apply(const tt_spectrum_plan *p, double temperature, double *y) {
    const tt_temperature_model *m = &p->s->temperature;
    size_t nr = m->range_count ? m->range_count : 1;
    for (size_t r = 0; r < nr; ++r) {
        double t = temperature;
        if (m->range_count) {
            const tt_temperature_range *range = m->ranges + r;
            if (t < range->minimum || t > range->maximum)
                continue;
            if (range->is_fixed)
                t = range->fixed_value;
        }
        t = t - m->base_temperature;
        if (m->kind == TT_TEMPERATURE_EXPONENTIAL && m->invert_temperature)
            t = -t;
        if ((m->kind == TT_TEMPERATURE_HARWOOD || m->kind == TT_TEMPERATURE_BURKHOLDER) && t == 0)
            return TT_RANGE;
        double *cs = y + 2 + m->first;
        switch (m->kind) {
        case TT_TEMPERATURE_EXPONENTIAL: {
            double factor[64];
            for (size_t j = 0; j < m->terms; ++j)
                factor[j] = m->a[j] + t * m->b[j];
            for (size_t i = 0; i < m->count; ++i) {
                double value = 0;
                for (size_t j = 0; j < m->terms; ++j)
                    value = value + factor[j] * p->powers[j * m->count + i];
                cs[i] = cs[i] + (m->base_10 ? pow(10, value) : exp(value));
            }
            break;
        }
        case TT_TEMPERATURE_TAYLOR: {
            double tp[64];
            for (size_t j = 0; j < m->terms; ++j)
                tp[j] = integer_power(t, j + 1);
            for (size_t i = 0; i < m->count; ++i) {
                double value = 1;
                for (size_t j = 0; j < m->terms; ++j)
                    value = value + m->a[j * m->count + i] * tp[j];
                cs[i] = value * m->sigma[i];
            }
            break;
        }
        case TT_TEMPERATURE_BURKHOLDER: {
            if (m->burkholder_b == 0)
                return TT_INVALID;
            double q = 1 + exp(m->burkholder_a / (m->burkholder_b * t));
            for (size_t i = 0; i < m->count; ++i)
                cs[i] = (m->a[i] / q + m->b[i] * (1 - 1 / q)) * 1.e-20;
            break;
        }
        case TT_TEMPERATURE_HARWOOD:
            for (size_t i = 0; i < m->count; ++i)
                cs[i] = pow(10, m->a[i] + m->b[i] / t);
            break;
        }
    }
    return TT_OK;
}

/* Conserving interpolation of target bins [first, end), as in numerics.c. */
static void conserve(size_t nx, const double *x, const double *y, const double *target,
                     size_t first, size_t end, double *out, size_t stride) {
    size_t j = 0;
    for (size_t i = first; i < end; ++i) {
        double sum = 0.0;
        while (j + 1 < nx && x[j + 1] <= target[i])
            ++j;
        for (size_t k = j; k + 1 < nx && x[k] < target[i + 1]; ++k) {
            double a1 = fmax(x[k], target[i]), a2 = fmin(x[k + 1], target[i + 1]);
            double slope = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
            double b1 = y[k] + slope * (a1 - x[k]), b2 = y[k] + slope * (a2 - x[k]);
            sum = sum + 0.5 * (a2 - a1) * (b2 + b1);
        }
        out[(i - first) * stride] = sum / (target[i + 1] - target[i]);
    }
}

static tt_status parameterized_prepare(tt_spectrum_plan *p) {
    const tt_spectrum *s = p->s;
    const tt_temperature_model *m = &s->temperature;
    size_t n = s->raw_count, np = n + 4, nw = p->nw;
    if (m->terms > 64)
        return TT_INVALID;
    p->raw = np;
    p->x = malloc(2 * np * sizeof(double));
    if (!p->x)
        return TT_MEMORY;
    p->y = p->x + np;
    double delta = (double)1.e-5f;
    p->x[0] = 0;
    p->x[1] = (1 - delta) * s->raw_wavelength[0];
    memcpy(p->x + 2, s->raw_wavelength, n * sizeof(double));
    p->x[n + 2] = (1 + delta) * s->raw_wavelength[n - 1];
    p->x[n + 3] = 1.e38;
    p->y[0] = p->y[1] = p->y[n + 2] = p->y[n + 3] = 0;
    memcpy(p->y + 2, s->raw_data, n * sizeof(double));
    if (m->kind == TT_TEMPERATURE_EXPONENTIAL) {
        p->powers = malloc(m->terms * m->count * sizeof(double));
        if (!p->powers)
            return TT_MEMORY;
        for (size_t j = 0; j < m->terms; ++j)
            for (size_t i = 0; i < m->count; ++i)
                p->powers[j * m->count + i] =
                    pow(s->raw_wavelength[m->first + i] - m->base_wavelength, m->powers[j]);
    }
    /* A target bin varies when a segment touching a parameterized point overlaps it. */
    size_t lo = m->first + 2, hi = lo + m->count;
    const double *edges = s->wavelength.edges;
    p->first_bin = nw;
    p->end_bin = 0;
    for (size_t i = 0; i < nw; ++i)
        for (size_t k = lo ? lo - 1 : 0; k < hi && k + 1 < np; ++k)
            if (p->x[k + 1] > edges[i] && p->x[k] < edges[i + 1]) {
                if (i < p->first_bin)
                    p->first_bin = i;
                p->end_bin = i + 1;
            }
    if (p->first_bin > p->end_bin)
        p->first_bin = p->end_bin = 0;
    return TT_OK;
}

static tt_status parameterized_evaluate(const tt_spectrum_plan *p, size_t nz, const double *temp,
                                        double *y, double *rows) {
    for (size_t i = 0; i < nz; ++i) {
        memcpy(y, p->y, p->raw * sizeof(double));
        tt_status status = apply(p, temp[i], y);
        if (status)
            return status;
        conserve(p->raw, p->x, y, p->s->wavelength.edges, p->first_bin, p->end_bin, rows + i, nz);
    }
    return TT_OK;
}

/* ---- Construction and evaluation --------------------------------------- */

void tt_spectrum_plan_destroy(tt_spectrum_plan *p) {
    if (!p)
        return;
    free(p->constant);
    free(p->coef);
    free(p->index);
    free(p->x);
    free(p->powers);
    free(p);
}

/* Reference columns for the constant values and the classification guard:
 * temperatures and air densities spanning every branch of the kernels. */
enum { REFERENCE_LEVELS = 12 };
static const double reference_temperature[3][REFERENCE_LEVELS] = {
    {250, 150, 199.5, 205, 212.25, 234.5, 236, 272.73, 296.5, 299.75, 340, 1500},
    {251.37, 120, 180, 210, 233, 260, 280.25, 295, 298, 330, 450, 9000},
    {222.2, 100, 194, 203, 218, 235, 255.5, 290, 300, 305, 800, 3000}};
static const double reference_air[3][REFERENCE_LEVELS] = {
    {2.45e19, 0, 1.e10, 1.e15, 7.e17, 8.e17, 9.e17, 1.e18, 1.e19, 2.6e19, 2.7e19, 1.e20},
    {2.5e19, 3.e19, 1.e17, 2.e18, 5.e18, 2.465e19, 2.69e19, 4.e16, 6.e19, 0, 1.e21, 3.3e18},
    {1.e19, 5.e18, 2.45e19, 1.e14, 2.6e19, 8.e17, 1.e16, 3.e17, 1.e18, 4.e19, 0, 7.e19}};

tt_status tt_spectrum_plan_create(const tt_spectrum *s, tt_spectrum_plan **out) {
    if (!s || !out || s->role == TT_SPECTRUM_WEIGHT)
        return TT_INVALID;
    tt_spectrum_plan *p = calloc(1, sizeof *p);
    if (!p)
        return TT_MEMORY;
    p->s = s;
    p->cs = s->role == TT_SPECTRUM_CROSS_SECTION && !s->parameterized;
    p->qy = s->role == TT_SPECTRUM_QUANTUM_YIELD;
    p->kind = s->parameterized ? -1 : s->kind;
    p->nw = s->wavelength.cells;
    p->tables = s->table_count;
    p->maximum_temperature = TEMPERATURE_LIMIT;
    size_t nw = p->nw, nz = REFERENCE_LEVELS, usable = 0;
    tt_status status = TT_OK;
    if (p->cs)
        status = tt_cross_section_validate(p->kind, nw, p->tables, s->tables, &p->required);
    p->stride = p->cs && p->kind == TT_CS_HCFC ? s->tables[0].rows : p->cs ? CS_STRIDE : QY_STRIDE;
    if (!p->stride)
        p->stride = 1;
    double *reference = malloc(3 * nw * nz * sizeof(double));
    int valid[3] = {0};
    p->constant = malloc(nw * sizeof(double));
    p->coef = malloc(nw * p->stride * sizeof(double));
    p->index = malloc(nw * sizeof(size_t));
    if (!status && !(reference && p->constant && p->coef && p->index))
        status = TT_MEMORY;
    /* Generic values on the reference columns. A parameterization rejects a
     * temperature equal to its base; such a column is simply not used. */
    for (size_t r = 0; r < 3 && !status; ++r) {
        tt_status column = tt_spectrum_calculate(s, nz, reference_temperature[r], reference_air[r],
                                                 reference + r * nw * nz);
        if (column == TT_OK)
            valid[r] = 1, ++usable;
        else if (column != TT_RANGE)
            status = column;
    }
    if (!status && !usable)
        status = TT_RANGE;
    if (!status && s->parameterized)
        status = parameterized_prepare(p);
    size_t first = valid[0] ? 0 : valid[1] ? 1 : 2;
    for (size_t w = 0; w < nw && !status && !p->generic; ++w) {
        double *c = p->coef + p->count * p->stride;
        int varies = s->parameterized ? w >= p->first_bin && w < p->end_bin
                     : p->cs          ? cs_prepare(p, w, c)
                                      : qy_prepare(p, w, c);
        p->constant[w] = reference[first * nw * nz + w * nz];
        /* Guard: a constant wavelength must be constant on every reference.
         * A failure would be a classification error; this spectrum then
         * falls back to the generic kernel at every wavelength. */
        for (size_t r = 0; r < 3 && !varies; ++r)
            for (size_t z = 0; z < nz && valid[r] && !p->generic; ++z)
                if (memcmp(reference + r * nw * nz + w * nz + z, p->constant + w, sizeof(double)))
                    p->generic = 1;
        if (varies)
            p->index[p->count++] = w;
    }
    if (!status && p->generic) {
        p->count = nw;
        for (size_t w = 0; w < nw; ++w)
            p->index[w] = w;
    }
    free(reference);
    if (status) {
        tt_spectrum_plan_destroy(p);
        return status;
    }
    *out = p;
    return TT_OK;
}

const double *tt_spectrum_plan_constant(const tt_spectrum_plan *p) {
    return p->constant;
}
size_t tt_spectrum_plan_varying(const tt_spectrum_plan *p, const size_t **indices) {
    if (indices)
        *indices = p->index;
    return p->count;
}
int tt_spectrum_plan_generic(const tt_spectrum_plan *p) {
    return p->generic;
}
double tt_spectrum_plan_maximum_temperature(const tt_spectrum_plan *p) {
    return p->maximum_temperature;
}
size_t tt_spectrum_plan_scratch(const tt_spectrum_plan *p, size_t nz) {
    if (p->generic)
        return 0;
    if (p->kind < 0)
        return p->raw;
    if (p->cs && p->kind == TT_CS_HCFC)
        return p->s->tables[0].rows * nz;
    if (p->qy && p->kind == TT_QY_ACETONE)
        return 12 * nz;
    return 8 * nz + 2 * p->tables * nz;
}
tt_status tt_spectrum_plan_evaluate(const tt_spectrum_plan *p, size_t nz, const double *temperature,
                                    const double *air, double *scratch, double *rows) {
    if (!p->count)
        return TT_OK;
    if (p->generic) /* every wavelength varies: rows are the generic layout */
        return tt_spectrum_calculate(p->s, nz, temperature, air, rows);
    if (p->kind < 0)
        return parameterized_evaluate(p, nz, temperature, scratch, rows);
    if (p->cs)
        cs_evaluate(p, nz, temperature, scratch, rows);
    else
        qy_evaluate(p, nz, temperature, air, scratch, rows);
    return TT_OK;
}

static int same(const double *a, const double *b, size_t n) {
    return n == 0 || (a && b && !memcmp(a, b, n * sizeof(double)));
}
static int same_tables(const tt_spectrum *a, const tt_spectrum *b) {
    if (a->table_count != b->table_count)
        return 0;
    for (size_t i = 0; i < a->table_count; ++i) {
        const tt_spectral_table *x = a->tables + i, *y = b->tables + i;
        if (x->rows != y->rows || x->columns != y->columns ||
            x->temperature_count != y->temperature_count ||
            !same(x->values, y->values, x->rows * x->columns) ||
            !same(x->temperatures, y->temperatures, x->temperature_count))
            return 0;
    }
    return 1;
}
static int same_options(const tt_yield_options *a, const tt_yield_options *b, size_t nw) {
    return a->acetone_branches == b->acetone_branches &&
           same(&a->low_wavelength_value, &b->low_wavelength_value, 1) &&
           same(&a->high_wavelength_value, &b->high_wavelength_value, 1) &&
           same(&a->minimum_temperature, &b->minimum_temperature, 1) &&
           same(&a->maximum_temperature, &b->maximum_temperature, 1) &&
           same(&a->molecular_diameter, &b->molecular_diameter, 1) &&
           same(&a->molecular_weight, &b->molecular_weight, 1) &&
           !a->collision_rates == !b->collision_rates &&
           (!a->collision_rates || same(a->collision_rates, b->collision_rates, nw));
}
static int same_ranges(const tt_temperature_model *x, const tt_temperature_model *y) {
    for (size_t r = 0; r < x->range_count; ++r) {
        const tt_temperature_range *a = x->ranges + r, *b = y->ranges + r;
        if (a->is_fixed != b->is_fixed || !same(&a->minimum, &b->minimum, 1) ||
            !same(&a->maximum, &b->maximum, 1) || !same(&a->fixed_value, &b->fixed_value, 1))
            return 0;
    }
    return 1;
}
int tt_spectrum_plan_equal(const tt_spectrum_plan *a, const tt_spectrum_plan *b) {
    if (a == b)
        return 1;
    if (a->cs != b->cs || a->qy != b->qy || a->kind != b->kind || a->generic != b->generic ||
        a->nw != b->nw || a->count != b->count || a->stride != b->stride ||
        !same(a->constant, b->constant, a->nw) || !same(a->coef, b->coef, a->count * a->stride) ||
        (a->count && memcmp(a->index, b->index, a->count * sizeof(size_t))) ||
        !same(a->s->wavelength.edges, b->s->wavelength.edges, a->nw + 1) ||
        !same(a->s->wavelength.midpoints, b->s->wavelength.midpoints, a->nw))
        return 0;
    if (a->kind < 0) {
        const tt_temperature_model *x = &a->s->temperature, *y = &b->s->temperature;
        size_t na = x->kind == TT_TEMPERATURE_TAYLOR        ? x->terms * x->count
                    : x->kind == TT_TEMPERATURE_EXPONENTIAL ? x->terms
                                                            : x->count;
        size_t nb = x->kind == TT_TEMPERATURE_EXPONENTIAL ? x->terms : na;
        return x->kind == y->kind && x->first == y->first && x->count == y->count &&
               x->terms == y->terms && x->base_10 == y->base_10 &&
               x->invert_temperature == y->invert_temperature && x->range_count == y->range_count &&
               a->raw == b->raw && same(&x->base_temperature, &y->base_temperature, 1) &&
               same(&x->base_wavelength, &y->base_wavelength, 1) &&
               same(&x->burkholder_a, &y->burkholder_a, 1) &&
               same(&x->burkholder_b, &y->burkholder_b, 1) && same(x->a, y->a, na) &&
               (x->kind == TT_TEMPERATURE_TAYLOR || same(x->b, y->b, nb)) &&
               (x->kind != TT_TEMPERATURE_TAYLOR || same(x->sigma, y->sigma, x->count)) &&
               same(a->x, b->x, 2 * a->raw) &&
               (!a->powers || same(a->powers, b->powers, x->terms * x->count)) && same_ranges(x, y);
    }
    return same_tables(a->s, b->s) && (!a->qy || same_options(&a->s->yield, &b->s->yield, a->nw));
}
