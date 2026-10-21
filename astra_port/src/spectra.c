#include "spectral_validation.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double clamp(double t, double lo, double hi) {
    return fmax(lo, fmin(t, hi));
}
static double sq(double x) {
    return x * x;
}
static double power(double x, size_t n) {
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

const char *tt_cross_section_name(tt_cross_section_kind kind) {
    static const char *const names[] = {"base",
                                        "CH3COCH3+hv->CH3CO+CH3",
                                        "CCl4+hv->Products",
                                        "CCl3F+hv->Products",
                                        "CH2O",
                                        "CH3ONO2+hv->CH3O+NO2",
                                        "CHBr3+hv->Products",
                                        "CHCl3+hv->Products",
                                        "Cl2+hv->Cl+Cl",
                                        "ClONO2",
                                        "H2O2+hv->OH+OH",
                                        "HCFC+hv->Products",
                                        "HNO3+hv->OH+NO2",
                                        "HOBr+hv->OH+Br",
                                        "N2O+hv->N2+O(1D)",
                                        "N2O5+hv->NO2+NO3",
                                        "CH3COCH2(ONO2)+hv->CH3COCH2(O.)+NO2",
                                        "CH2(OH)CH2(ONO2)+hv->CH2(OH)CH2(O.)+NO2",
                                        "NO2 tint",
                                        "O3 JPL06",
                                        "O3",
                                        "OClO+hv->Products",
                                        "air",
                                        "RONO2",
                                        "t_butyl_nitrate+hv->Products",
                                        "tint"};
    return (unsigned)kind < TT_CS_COUNT ? names[kind] : NULL;
}

double tt_air_refractive_index(double wavelength, double density) {
    double sig = 1.e3 / clamp(wavelength, 200, 2000), sigsq = sig * sig;
    double value = 8342.13 + 2406030.0 / (130.0 - sigsq) + 15997.0 / (38.9 - sigsq);
    value = value * density / (2.69e19 * 273.15 / 288.15);
    return 1.0 + 1.e-8 * value;
}

tt_status tt_regrid_column(const tt_dataset *data, const char *name, size_t column, size_t nw,
                           const double *edges, tt_extrapolation lower, double lo,
                           tt_extrapolation upper, double hi, double *out) {
    const tt_data_variable *w = tt_dataset_variable(data, "wavelength");
    const tt_data_variable *v = tt_dataset_variable(data, name);
    if (!w || !v || !w->values || !v->values || !edges || !out || !nw || w->rank != 1 ||
        w->count < 2 || column >= v->count / w->count || v->count % w->count ||
        (unsigned)lower > TT_EXTRAPOLATE_CONSTANT || (unsigned)upper > TT_EXTRAPOLATE_CONSTANT)
        return TT_INVALID;
    size_t n = w->count;
    if (n > (SIZE_MAX / sizeof(double)) / 2 - 4)
        return TT_MEMORY;
    double *x = malloc(2 * (n + 4) * sizeof(double));
    if (!x)
        return TT_MEMORY;
    double *y = x + n + 4;
    if (lower == TT_EXTRAPOLATE_ZERO)
        lo = 0.0;
    if (lower == TT_EXTRAPOLATE_BOUNDARY)
        lo = v->values[column * n];
    if (upper == TT_EXTRAPOLATE_ZERO)
        hi = 0.0;
    if (upper == TT_EXTRAPOLATE_BOUNDARY)
        hi = v->values[(column + 1) * n - 1];
    x[0] = 0.0;
    x[1] = (1.0 - 1.e-5) * w->values[0];
    y[0] = y[1] = lo;
    memcpy(x + 2, w->values, n * sizeof(double));
    memcpy(y + 2, v->values + column * n, n * sizeof(double));
    x[n + 2] = (1.0 + 1.e-5) * w->values[n - 1];
    x[n + 3] = 1.e38;
    y[n + 2] = y[n + 3] = hi;
    tt_status status = tt_interpolate(TT_CONSERVING, 0, n + 4, x, y, nw + 1, edges, out);
    free(x);
    return status;
}

static double parameter(const tt_spectral_table *p, size_t w, size_t c) {
    return p->values[c * p->rows + w];
}
static double interpolate_temperature(const tt_spectral_table *p, size_t w, double temp) {
    size_t n = p->temperature_count, k = 0;
    temp = clamp(temp, p->temperatures[0], p->temperatures[n - 1]);
    while (k + 2 < n && temp > p->temperatures[k + 1])
        ++k;
    double f = (temp - p->temperatures[k]) / (p->temperatures[k + 1] - p->temperatures[k]);
    return parameter(p, w, k) + f * (parameter(p, w, k + 1) - parameter(p, w, k));
}

tt_status tt_cross_section_validate(tt_cross_section_kind kind, size_t nw, size_t count,
                                    const tt_spectral_table *p, size_t *required_out) {
    if ((unsigned)kind >= TT_CS_COUNT || !nw || (count && !p))
        return TT_INVALID;
    int analytic = kind == TT_CS_CL2 || kind == TT_CS_HOBR || kind == TT_CS_N2O ||
                   kind == TT_CS_NITROXY_ACETONE || kind == TT_CS_NITROXY_ETHANOL ||
                   kind == TT_CS_RAYLEIGH || kind == TT_CS_T_BUTYL_NITRATE;
    size_t required = analytic ? 0 : 1;
    if (kind == TT_CS_N2O5)
        required = 2;
    if (kind == TT_CS_O3_JPL06)
        required = 3;
    if (kind == TT_CS_O3)
        required = 4;
    if (count < required)
        return TT_INVALID;
    for (size_t i = 0; i < count; ++i)
        if (!p[i].values || !p[i].rows || !p[i].columns ||
            (kind != TT_CS_HCFC && p[i].rows != nw) ||
            (p[i].temperature_count && !p[i].temperatures))
            return TT_INVALID;
    size_t columns = 1;
    if (kind == TT_CS_ACETONE)
        columns = 4;
    if (kind == TT_CS_CLONO2 || kind == TT_CS_HCFC)
        columns = 3;
    if (kind == TT_CS_CH2O || kind == TT_CS_CH3ONO2 || kind == TT_CS_HNO3 || kind == TT_CS_RONO2)
        columns = 2;
    if (required && p[0].columns < columns)
        return TT_INVALID;
    if (kind == TT_CS_HCFC && (p[0].temperature_count < 2 || !isfinite(p[0].temperatures[0]) ||
                               !isfinite(p[0].temperatures[1])))
        return TT_INVALID;
    if (kind == TT_CS_OCLO) {
        for (size_t i = 0; i < count; ++i)
            if (p[i].temperature_count < 1 || !isfinite(p[i].temperatures[0]) ||
                (i && !(p[i].temperatures[0] > p[i - 1].temperatures[0])))
                return TT_INVALID;
    }
    if (kind == TT_CS_TINT || kind == TT_CS_NO2 || kind == TT_CS_O3)
        for (size_t i = kind == TT_CS_O3 ? 1 : 0; i < count; ++i) {
            if (p[i].temperature_count < 2 || p[i].columns < p[i].temperature_count)
                return TT_INVALID;
            for (size_t k = 0; k < p[i].temperature_count; ++k)
                if (!isfinite(p[i].temperatures[k]) ||
                    (k && !(p[i].temperatures[k] > p[i].temperatures[k - 1])))
                    return TT_INVALID;
        }
    if (required_out)
        *required_out = required;
    return TT_OK;
}

tt_status tt_cross_section(tt_cross_section_kind kind, const tt_spectral_grid *grid, size_t count,
                           const tt_spectral_table *p, double *out) {
    if ((unsigned)kind >= TT_CS_COUNT || !grid || !out || !grid->edges || !grid->midpoints ||
        !grid->temperature || !grid->wavelengths || !grid->levels ||
        grid->levels > SIZE_MAX / grid->wavelengths || (count && !p))
        return TT_INVALID;
    size_t nw = grid->wavelengths, nz = grid->levels;
    size_t required = 0;
    tt_status status = tt_cross_section_validate(kind, nw, count, p, &required);
    if (status != TT_OK)
        return status;
    /* Temperatures must be positive where a kernel uses them; as in TUV-x,
     * temperature-independent kinds ignore them. */
    int thermal = kind != TT_CS_BASE && kind != TT_CS_HOBR && kind != TT_CS_NITROXY_ACETONE &&
                  kind != TT_CS_NITROXY_ETHANOL && kind != TT_CS_RAYLEIGH &&
                  kind != TT_CS_T_BUTYL_NITRATE;
    for (size_t z = 0; z < nz && thermal; ++z)
        if (!isfinite(grid->temperature[z]) || grid->temperature[z] <= 0)
            return TT_INVALID;
    double v185 = 185 * tt_air_refractive_index(185, 2.45e19);
    double v195 = 195 * tt_air_refractive_index(195, 2.45e19);
    double v345 = 345 * tt_air_refractive_index(345, 2.45e19);
    for (size_t w = 0; w < nw; ++w) {
        double lambda = grid->midpoints[w], edge = grid->edges[w];
        if (!isfinite(lambda) || lambda <= 0 || !isfinite(edge) || grid->edges[w + 1] <= edge)
            return TT_INVALID;
        for (size_t z = 0; z < nz; ++z) {
            double temp = grid->temperature[z], value = 0, t, coeff;
            double base = required && kind != TT_CS_HCFC ? parameter(p, w, 0) : 0;
            switch (kind) {
            case TT_CS_BASE:
                value = base;
                break;
            case TT_CS_ACETONE:
                t = clamp(temp, 235, 298);
                value = base * (1.0 + t * (parameter(p, w, 1) +
                                           t * (parameter(p, w, 2) + t * parameter(p, w, 3))));
                break;
            case TT_CS_CCL4:
                value = base;
                if (lambda > 194 && lambda < 250) {
                    t = clamp(temp, 210, 300) - 295;
                    coeff =
                        1.0739 + lambda * (-1.6275e-2 +
                                           lambda * (8.8141e-5 +
                                                     lambda * (-1.9811e-7 + 1.5022e-10 * lambda)));
                    value = base * pow(10, coeff * t);
                }
                break;
            case TT_CS_CFC11:
                t = 1.e-4 * (temp - 298);
                value = base * exp((lambda - 184.9) * t);
                break;
            case TT_CS_CH2O:
                value = base + parameter(p, w, 1) * (temp - 298);
                break;
            case TT_CS_CH3ONO2:
            case TT_CS_HNO3:
            case TT_CS_RONO2:
                value = base * exp(parameter(p, w, 1) * (temp - 298));
                break;
            case TT_CS_CHBR3:
                value = lambda > 290 && lambda < 340 && temp > 210 && temp < 300
                            ? exp((.06183 - .000241 * lambda) * (273 - temp) -
                                  (2.376 + .14757 * lambda))
                            : base;
                break;
            case TT_CS_CHCL3:
                value = base;
                if (lambda > 190 && lambda < 240) {
                    t = clamp(temp, 210, 300) - 295;
                    coeff = 3.7973 +
                            lambda *
                                (-7.0913e-2 +
                                 lambda * (4.9397e-4 + lambda * (-1.5226e-6 + lambda * 1.7555e-9)));
                    value = base * pow(10, coeff * t);
                }
                break;
            case TT_CS_CL2: {
                double bb = exp(402.7 / temp), bbsq = bb * bb, alpha = (bbsq - 1) / (bbsq + 1);
                double ex1 = 27.3 * exp(-99.0 * alpha * sq(log(329.5 / lambda)));
                double ex2 = .932 * exp(-91.5 * alpha * sq(log(406.5 / lambda)));
                value = 1.e-20 * sqrt(alpha) * (ex1 + ex2);
                break;
            }
            case TT_CS_CLONO2:
                t = temp - 296;
                value = base * (1 + t * (parameter(p, w, 1) + t * parameter(p, w, 2)));
                break;
            case TT_CS_H2O2:
                value = base;
                if (edge >= 260 && edge < 350) {
                    double a =
                        ((((((1.5534675e-13 * lambda - 2.652014e-10) * lambda + 1.6878206e-7) *
                                lambda -
                            4.035101e-5) *
                               lambda -
                           4.4589016e-3) *
                              lambda +
                          4.535649) *
                             lambda -
                         9.2170972e2) *
                            lambda +
                        6.4761e4;
                    double b = (((-1.0924e-7 * lambda - 3.0493e-5) * lambda + 1.1522e-1) * lambda -
                                5.1351e1) *
                                   lambda +
                               6.8123e3;
                    double chi = 1 / (1 + exp(-1265 / clamp(temp, 200, 400)));
                    value = (chi * a + (1 - chi) * b) * 1.e-21;
                }
                break;
            case TT_CS_HCFC:
                if (lambda >= 190 && lambda <= p[0].temperatures[1]) {
                    t = clamp(temp, 203, 295) - p[0].temperatures[0];
                    double sigma = 0;
                    for (size_t k = 0; k < p[0].rows; ++k)
                        sigma = sigma + (parameter(p, k, 0) +
                                         t * (parameter(p, k, 1) + t * parameter(p, k, 2))) *
                                            power(lambda - 206.214, k);
                    value = exp(sigma);
                }
                break;
            case TT_CS_HOBR:
                if (lambda >= 250 && lambda <= 550) {
                    value = 24.77 * exp(-109.80 * sq(log(284.01 / lambda))) +
                            12.22 * exp(-93.63 * sq(log(350.57 / lambda))) +
                            2.283 * exp(-242.40 * sq(log(457.38 / lambda)));
                    value = value * 1.e-20;
                }
                break;
            case TT_CS_N2O:
                if (lambda >= 173 && lambda <= 240) {
                    double a =
                        (((2.520672e-7 * lambda - 1.777846e-4) * lambda + 4.301146e-2) * lambda -
                         4.071805) *
                            lambda +
                        68.21023;
                    double b =
                        ((-1.881058e-5 * lambda + 1.111572e-2) * lambda - 2.116255) * lambda +
                        123.4014;
                    b = (clamp(temp, 194, 320) - 300) * exp(b);
                    value = exp(a + b);
                }
                break;
            case TT_CS_N2O5:
                t = clamp(temp, 233, 300);
                coeff = 1000 * parameter(p + 1, w, 0) * (300 - t) / (300 * t);
                value = base * pow(10, coeff);
                break;
            case TT_CS_NITROXY_ACETONE:
                if (lambda >= 284 && lambda <= 335)
                    value = exp(-156.8 + lambda * (.7834 - 1.365e-3 * lambda));
                break;
            case TT_CS_NITROXY_ETHANOL:
                if (lambda >= 270 && lambda <= 306)
                    value = exp(-210.4 + lambda * (1.2478 - 2.359e-3 * lambda));
                break;
            case TT_CS_NO2:
            case TT_CS_TINT:
                /* Keep the reference's sequential addition order. */
                for (size_t k = 0; k < count; ++k) {
                    const tt_spectral_table *q = p + k;
                    size_t j = 0;
                    t = clamp(temp, q->temperatures[0], q->temperatures[q->temperature_count - 1]);
                    while (j + 2 < q->temperature_count && t > q->temperatures[j + 1])
                        ++j;
                    double f =
                        (t - q->temperatures[j]) / (q->temperatures[j + 1] - q->temperatures[j]);
                    value = value + parameter(q, w, j) +
                            f * (parameter(q, w, j + 1) - parameter(q, w, j));
                }
                break;
            case TT_CS_O3_JPL06:
                value = base;
                if (edge >= 196.078 && grid->edges[w + 1] <= 342.5) {
                    double low = parameter(p + 1, w, 0), high = parameter(p + 2, w, 0);
                    value = temp < 218    ? low
                            : temp <= 298 ? low + (high - low) / (298 - 218) * (temp - 218)
                                          : high;
                }
                break;
            case TT_CS_O3:
                value = edge < v185   ? interpolate_temperature(p + 2, w, temp)
                        : edge < v195 ? interpolate_temperature(p + 3, w, temp)
                        : edge < v345 ? interpolate_temperature(p + 1, w, temp)
                                      : base;
                break;
            case TT_CS_OCLO:
                if (temp <= p[0].temperatures[0])
                    value = base;
                else if (temp >= p[count - 1].temperatures[0])
                    value = parameter(p + count - 1, w, 0);
                else {
                    size_t k = 0;
                    while (k + 2 < count && p[k + 1].temperatures[0] <= temp)
                        ++k;
                    t = (temp - p[k].temperatures[0]) /
                        (p[k + 1].temperatures[0] - p[k].temperatures[0]);
                    value = parameter(p + k, w, 0) +
                            t * (parameter(p + k + 1, w, 0) - parameter(p + k, w, 0));
                }
                break;
            case TT_CS_RAYLEIGH:
                t = 1.e-3 * lambda;
                coeff = t <= .55 ? 3.6772 + .389 * t + .09426 / t : 4.04;
                value = 4.02e-28 / pow(t, coeff);
                break;
            case TT_CS_T_BUTYL_NITRATE:
                if (lambda >= 270 && lambda <= 330)
                    value = exp(-115.5 + lambda * (.5307 - .993e-3 * lambda));
                break;
            case TT_CS_COUNT:
                return TT_INVALID;
            }
            out[w * nz + z] = value;
        }
    }
    return TT_OK;
}
