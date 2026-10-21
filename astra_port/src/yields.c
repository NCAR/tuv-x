#include "spectral_validation.h"
#include <float.h>
#include <math.h>
#include <stdint.h>

static double clamp(double x, double a, double b) {
    return fmax(a, fmin(x, b));
}
static double sq(double x) {
    return x * x;
}
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
static double parameter(const tt_spectral_table *p, size_t w, size_t c) {
    return p->values[c * p->rows + w];
}

tt_yield_options tt_yield_options_default(void) {
    tt_yield_options o = {0};
    o.low_wavelength_value = .95;
    o.maximum_temperature = DBL_MAX;
    o.acetone_branches = 2;
    o.maximum_wavelength = DBL_MAX;
    return o;
}
const char *tt_quantum_yield_name(tt_quantum_yield_kind kind) {
    static const char *const names[] = {"base",
                                        "tint",
                                        "NO2 tint",
                                        "O3+hv->O2+O(1D)",
                                        "O3+hv->O2+O(3P)",
                                        "HO2",
                                        "NO3-_(aq)+hv->NO2(aq)+O-",
                                        "CH2O",
                                        "CH3CHO+hv->CH3+HCO",
                                        "C2H5CHO",
                                        "CH2CHCHO+hv->Products",
                                        "MVK+hv->Products",
                                        "CH3COCH3+hv->CH3CO+CH3",
                                        "CH3COCH2CH3+hv->CH3CO+CH2CH3",
                                        "CH3COCHO+hv->CH3CO+HCO",
                                        "ClO+hv->Cl+O(1D)",
                                        "ClO+hv->Cl+O(3P)",
                                        "ClONO2+hv->Cl+NO3",
                                        "ClONO2+hv->ClO+NO2",
                                        "H2SO4 Mills",
                                        "Taylor series"};
    return (unsigned)kind < TT_QY_COUNT ? names[kind] : NULL;
}

static double ozone(double w, double temp) {
    if (w <= 305)
        return .90;
    if (w > 328)
        return w <= 340 ? .08 : 0;
    double kt = .695 * temp, q2 = exp(-825.518 / kt), qfac1 = 1 / (1 + q2), qfac2 = q2 / (1 + q2),
           t300 = temp / 300;
    double fourth = sq(sq((304.225 - w) / 5.576));
    return .0765 + .8036 * qfac1 * exp(-fourth) +
           8.9061 * t300 * t300 * qfac2 * exp(-sq((314.957 - w) / 6.601)) +
           .1192 * pow(t300, 1.5) * exp(-sq((310.737 - w) / 2.187));
}

static double acetone(double w, double temp, double density, const tt_yield_options *o) {
    if (w < 279)
        return o->low_wavelength_value;
    if (w > 327)
        return o->high_wavelength_value;
    double t = clamp(temp, o->minimum_temperature, o->maximum_temperature) / 295;
    double a0 = .350 * pow(t, -1.28), b0 = .068 * pow(t, -2.65), arg = b0 * (w - 248);
    double ca0 = arg > 80 ? 5.e34 : exp(arg) * a0 / (1 - a0), fco = 1 / (1 + ca0), wadj = 1.e7 / w,
           fac;
    if (w < 302) {
        double a1 = 1.600e-19 * pow(t, -2.38), b1 = .55e-3 * pow(t, -3.19);
        double ca1 = a1 * exp(-b1 * (wadj - 33113));
        fac = (1 - fco) / (1 + ca1 * density);
    } else {
        double a2 = 1.62e-17 * pow(t, -10.03), b2 = 1.79e-3 * pow(t, -1.364);
        double ca2 = a2 * exp(-b2 * (wadj - 30488));
        double a3 = 26.29 * pow(t, -6.59), b3 = 5.72e-7 * pow(t, -2.93), c3 = 30006 * pow(t, -.064);
        double ca3 = a3 * exp(-b3 * sq(wadj - c3));
        double a4 = 1.67e-15 * pow(t, -7.25), b4 = 2.08e-3 * pow(t, -1.16);
        double ca4 = a4 * exp(-b4 * (wadj - 30488));
        fac = (1 - fco) * (1 + ca3 + ca4 * density) /
              ((1 + ca3 + ca2 * density) * (1 + ca4 * density));
    }
    double qy = 0;
    if (o->acetone_branches & 1)
        qy = qy + fco;
    if (o->acetone_branches & 2)
        qy = qy + fac;
    return qy;
}

tt_status tt_quantum_yield_validate(tt_quantum_yield_kind kind, size_t nw, size_t count,
                                    const tt_spectral_table *p, const tt_yield_options *o,
                                    int *needs_data_out) {
    if ((unsigned)kind >= TT_QY_COUNT || !nw || !o || (count && !p))
        return TT_INVALID;
    int needs_data = kind == TT_QY_BASE || kind == TT_QY_CH2O || kind == TT_QY_CH3CHO ||
                     kind == TT_QY_C2H5CHO || kind == TT_QY_TINT || kind == TT_QY_NO2 ||
                     kind == TT_QY_H2SO4 || kind == TT_QY_TAYLOR;
    if (needs_data && !count)
        return TT_INVALID;
    for (size_t i = 0; i < count; ++i) {
        if (!p[i].values || p[i].rows != nw || !p[i].columns)
            return TT_INVALID;
        if (kind == TT_QY_TINT || kind == TT_QY_NO2) {
            if (!p[i].temperatures || p[i].temperature_count < 2 ||
                p[i].columns < p[i].temperature_count)
                return TT_INVALID;
            for (size_t j = 0; j < p[i].temperature_count; ++j)
                if (!isfinite(p[i].temperatures[j]) ||
                    (j && !(p[i].temperatures[j] > p[i].temperatures[j - 1])))
                    return TT_INVALID;
        }
    }
    if ((kind == TT_QY_CH2O || kind == TT_QY_CH3CHO) && p[0].columns < 2)
        return TT_INVALID;
    if (kind == TT_QY_ACETONE &&
        (o->acetone_branches < 1 || o->acetone_branches > 3 ||
         o->minimum_temperature > o->maximum_temperature || o->maximum_temperature <= 0))
        return TT_INVALID;
    if (kind == TT_QY_TAYLOR && (!o->coefficient_count || !o->coefficients))
        return TT_INVALID;
    if (kind == TT_QY_H2SO4 &&
        (!o->collision_rates || o->molecular_weight <= 0 || o->molecular_diameter <= 0))
        return TT_INVALID;
    if (needs_data_out)
        *needs_data_out = needs_data;
    return TT_OK;
}

tt_status tt_quantum_yield(tt_quantum_yield_kind kind, const tt_spectral_grid *grid, size_t count,
                           const tt_spectral_table *p, const tt_yield_options *options,
                           double *out) {
    if ((unsigned)kind >= TT_QY_COUNT || !grid || !out || !grid->wavelengths || !grid->levels ||
        !grid->midpoints || !grid->temperature || !grid->air_density ||
        grid->levels > SIZE_MAX / grid->wavelengths || (count && !p))
        return TT_INVALID;
    tt_yield_options defaults = tt_yield_options_default();
    const tt_yield_options *o = options ? options : &defaults;
    size_t nw = grid->wavelengths, nz = grid->levels;
    int needs_data = 0;
    tt_status status = tt_quantum_yield_validate(kind, nw, count, p, o, &needs_data);
    if (status != TT_OK)
        return status;
    /* Validate the profiles a kind uses; as in TUV-x, the others are ignored. */
    int thermal = kind == TT_QY_TINT || kind == TT_QY_NO2 || kind == TT_QY_O3_O1D ||
                  kind == TT_QY_O3_O3P || kind == TT_QY_NO3_AQ || kind == TT_QY_CH2O ||
                  kind == TT_QY_ACETONE || kind == TT_QY_H2SO4;
    int density = kind == TT_QY_CH2O || kind == TT_QY_CH3CHO || kind == TT_QY_C2H5CHO ||
                  kind == TT_QY_CH2CHCHO || kind == TT_QY_MVK || kind == TT_QY_ACETONE ||
                  kind == TT_QY_MEK || kind == TT_QY_CH3COCHO || kind == TT_QY_H2SO4;
    for (size_t z = 0; z < nz; ++z) {
        double temp = grid->temperature[z], m = grid->air_density[z];
        if ((thermal && (!isfinite(temp) || temp <= 0)) || (density && (!isfinite(m) || m < 0)))
            return TT_INVALID;
    }
    for (size_t w = 0; w < nw; ++w) {
        double lambda = grid->midpoints[w];
        if (!isfinite(lambda) || lambda <= 0)
            return TT_INVALID;
        for (size_t z = 0; z < nz; ++z) {
            double temp = grid->temperature[z], m = grid->air_density[z], value = 0,
                   base = needs_data ? parameter(p, w, 0) : 0;
            switch (kind) {
            case TT_QY_BASE:
                value = base;
                break;
            case TT_QY_TINT:
            case TT_QY_NO2:
                for (size_t i = 0; i < count; ++i) {
                    const tt_spectral_table *q = p + i;
                    double t = kind == TT_QY_TINT ? clamp(temp, q->temperatures[0],
                                                          q->temperatures[q->temperature_count - 1])
                                                  : temp;
                    size_t j = 0;
                    while (j + 2 < q->temperature_count && t > q->temperatures[j + 1])
                        ++j;
                    double f =
                        (t - q->temperatures[j]) / (q->temperatures[j + 1] - q->temperatures[j]);
                    value = value + parameter(q, w, j) +
                            f * (parameter(q, w, j + 1) - parameter(q, w, j));
                }
                if (kind == TT_QY_NO2)
                    value = fmax(value, 0);
                break;
            case TT_QY_O3_O1D:
                value = ozone(lambda, temp);
                break;
            case TT_QY_O3_O3P:
                value = 1 - ozone(lambda, temp);
                break;
            case TT_QY_HO2:
                value = lambda >= 248 ? 1 : fmax(0, (1 + 14 * (lambda - 193) / 55) / 15);
                break;
            case TT_QY_NO3_AQ:
                value = exp(-2400 / temp + 3.6);
                break;
            case TT_QY_CH2O: {
                double q1 = base, q2 = parameter(p, w, 1), qt = 1 - q1;
                value = q2;
                if (lambda >= 330 && lambda < 360 && q2 > 0) {
                    double f = (300 - temp) / 80;
                    double q = (1 - (q1 + q2)) / (2.45e19 * q2 * qt);
                    q = q * (1 + .05 * (lambda - 329) * f);
                    value = 1 / (1 / qt + q * m);
                }
                break;
            }
            case TT_QY_CH3CHO: {
                double q1 = parameter(p, w, 1), q2 = 1 - base, q = q1 > 0 ? q2 / q1 - 1 : 0;
                double fac = m / 2.465e19;
                value = clamp(q1 * (1 + q) / (1 + q * fac), 0, 1);
                break;
            }
            case TT_QY_C2H5CHO:
                value = base < 10.0 / 1.e36 ? 0 : fmin(1, 1 / (1 + (1 / base - 1) * (m / 2.45e19)));
                break;
            case TT_QY_CH2CHCHO:
                if (m > 2.6e19)
                    value = .004;
                else {
                    double phi = m <= 8.e17 ? .086 + 1.613e-17 * 8.e17 : .086 + 1.613e-17 * m;
                    value = .004 + 1 / phi;
                }
                break;
            case TT_QY_MVK:
                value = fmin(exp(-.055 * (lambda - 308)) / (5.5 + 9.2e-19 * m), 1);
                break;
            case TT_QY_ACETONE:
                value = acetone(lambda, temp, m, o);
                break;
            case TT_QY_MEK: {
                double ptorr = 760 * m / 2.69e19;
                value = fmin(1 / (.96 + 2.22e-3 * ptorr), 1);
                break;
            }
            case TT_QY_CH3COCHO: {
                double airfac = m * 760 / 2.456e19, phi = clamp(1 - (lambda - 380) / 60, 0, 1),
                       kq = 1.36e8 * exp(-8793 / lambda);
                value = phi > 0 && lambda >= 380 && lambda <= 440 ? phi / (phi + kq * airfac) : phi;
                break;
            }
            case TT_QY_CLO_O1D:
                value = lambda < 263.4 ? 1 : 0;
                break;
            case TT_QY_CLO_O3P:
                value = lambda < 263.4 ? 0 : 1;
                break;
            case TT_QY_CLONO2_CL:
            case TT_QY_CLONO2_CLO:
                value = lambda < 308 ? .6 : lambda <= 364 ? 7.143e-3 * lambda - 1.6 : 1;
                if (kind == TT_QY_CLONO2_CLO)
                    value = 1 - value;
                break;
            case TT_QY_H2SO4:
                value = base;
                if (o->collision_rates[w] >= 0) {
                    double pi = 3.1415926535898, rate = o->collision_rates[w];
                    value =
                        rate /
                        (rate + sqrt(16.0 * 8.31446261815324 * temp / (pi * o->molecular_weight)) *
                                    (pi * sq(o->molecular_diameter) * m * 1.e6));
                    if (z + 1 == nz)
                        value = 1;
                }
                break;
            case TT_QY_TAYLOR:
                value = base;
                if (lambda >= o->minimum_wavelength && lambda <= o->maximum_wavelength) {
                    value = o->coefficients[0];
                    for (size_t j = 1; j < o->coefficient_count; ++j)
                        value = value + o->coefficients[j] * integer_power(lambda, j);
                    value = clamp(value, 0, 1);
                }
                break;
            case TT_QY_COUNT:
                return TT_INVALID;
            }
            out[w * nz + z] = value;
        }
    }
    return TT_OK;
}
