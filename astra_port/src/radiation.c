/* Delta-Eddington equations and boundary conditions from TUV-x 0.17.0.
 * All arrays presented by the API ascend in altitude. Internal layer work
 * descends from the top, following the reference's arithmetic order. */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static double integer_power(double x, size_t n) {
    double y = 1.0;
    while (n) {
        if (n & 1)
            y *= x;
        n >>= 1;
        if (n)
            x *= x;
    }
    return y;
}

tt_status tt_mix_radiators(size_t nr, const tt_radiator *r, size_t nl, size_t nw, size_t nm,
                           double *tau, double *ssa, double *g) {
    if ((nr && !r) || !nl || !nw || !nm || !tau || !ssa || !g || nw > SIZE_MAX / nl ||
        nm > SIZE_MAX / (nl * nw))
        return TT_INVALID;
    size_t count = nl * nw;
    for (size_t j = 0; j < nr; ++j)
        if (!r[j].optical_depth || !r[j].single_scattering_albedo || !r[j].asymmetry)
            return TT_INVALID;
    for (size_t i = 0; i < count; ++i) {
        double scatter = 0.0, absorb = 0.0;
        for (size_t m = 0; m < nm; ++m)
            g[m * count + i] = 0.0;
        for (size_t j = 0; j < nr; ++j) {
            double od = r[j].optical_depth[i], a = r[j].single_scattering_albedo[i],
                   h = r[j].asymmetry[i];
            if (!isfinite(od) || od < 0 || !isfinite(a) || a < 0 || a > 1 || !isfinite(h) ||
                fabs(h) > 1)
                return TT_INVALID;
            double ds = od * a;
            scatter = scatter + ds;
            absorb = absorb + od * (1.0 - a);
            if (!r[j].is_air) {
                for (size_t m = 0; m < nm; ++m)
                    g[m * count + i] = g[m * count + i] + integer_power(h, m + 1) * ds;
            } else if (nm >= 2)
                g[count + i] = g[count + i] + 0.1 * ds;
            else
                g[i] = g[i] + h * ds;
        }
        scatter = fmax(scatter, 1.0 / 1.e36);
        absorb = fmax(absorb, 1.0 / 1.e36);
        tau[i] = scatter + absorb;
        ssa[i] = scatter == 1.0 / 1.e36 ? 1.0 / 1.e36 : scatter / tau[i];
        if (nm >= 2) {
            ssa[i] = fmax(fmin(ssa[i], 1.0 - 1.e-7), 1.e-7);
            tau[i] = fmax(tau[i], 1.e-7);
        }
        for (size_t m = 0; m < nm; ++m) {
            g[m * count + i] = g[m * count + i] / scatter;
        }
    }
    return TT_OK;
}

static double slant_depth(size_t layer, int crossed, const double *path, const double *tau) {
    if (crossed < 0)
        return (double)1.e36f; /* reference literal has default real kind */
    size_t single = (size_t)crossed < layer ? (size_t)crossed : layer;
    double result = 0.0;
    for (size_t j = 0; j < single; ++j)
        result = result + tau[j] * path[j];
    for (size_t j = single; j < (size_t)crossed; ++j)
        result = result + 2.0 * tau[j] * path[j];
    return result;
}

tt_status tt_delta_eddington(const tt_atmosphere *atm, double zenith, const double *surface,
                             tt_radiation *out) {
    if (!atm || !out || !surface || !atm->height || !atm->optical_depth ||
        !atm->single_scattering_albedo || !atm->legendre || !atm->moments || !out->direct_flux ||
        !out->down_flux || !out->up_flux || !out->direct_irradiance || !out->down_irradiance ||
        !out->up_irradiance)
        return TT_INVALID;
    size_t n = atm->layers, nw = atm->wavelengths;
    if (!n || !nw || n > INT32_MAX || n > SIZE_MAX / n ||
        n * n > SIZE_MAX / sizeof(double) - 32 * n - 3 || nw > SIZE_MAX / (n + 1))
        return TT_INVALID;
    double *storage = calloc(n * (n + 1) + 32 * n + 3, sizeof(double));
    int *nid = malloc((n + 1) * sizeof(int));
    if (!storage || !nid) {
        free(storage);
        free(nid);
        return TT_MEMORY;
    }
    double *path = storage, *taun = path + n * (n + 1), *g = taun + n, *om = g + n;
    double *slant = om + n, *tauc = slant + n + 1, *mu2 = tauc + n + 1;
    double *e1 = mu2 + n + 1, *e2 = e1 + n, *e3 = e2 + n, *e4 = e3 + n;
    double *cup = e4 + n, *cdn = cup + n, *cuptn = cdn + n, *cdntn = cuptn + n;
    double *a = cdntn + n, *b = a + 2 * n, *d = b + 2 * n, *rhs = d + 2 * n, *y = rhs + 2 * n,
           *work = y + 2 * n;
    tt_status status = tt_spherical_geometry(n + 1, atm->height, zenith, nid, path);
    if (status)
        goto done;
    double mu = cos(zenith * (3.1415926535898 / 180.0));
    for (size_t w = 0; w < nw; ++w) {
        double rsfc = surface[w];
        if (!isfinite(rsfc) || rsfc < 0 || rsfc > 1) {
            status = TT_INVALID;
            goto done;
        }
        for (size_t i = 0; i < n; ++i) {
            size_t k = w * n + n - 1 - i;
            double gu = atm->legendre[k], omu = atm->single_scattering_albedo[k],
                   tau = atm->optical_depth[k];
            if (!isfinite(tau) || tau < 0 || !isfinite(omu) || omu < 0 || omu > 1 ||
                !isfinite(gu) || fabs(gu) >= 1) {
                status = TT_INVALID;
                goto done;
            }
            double f = gu * gu;
            g[i] = (gu - f) / (1.0 - f);
            om[i] = (1.0 - f) * omu / (1.0 - omu * f);
            taun[i] = (1.0 - omu * f) * tau;
        }
        slant[0] = mu < 0 ? slant_depth(0, nid[0], path, taun) : 0.0;
        tauc[0] = 0.0;
        for (size_t i = 0; i < n; ++i) {
            tauc[i + 1] = tauc[i] + taun[i];
            double gi = copysign(fmin(fabs(g[i]), 1.0 - 1.e-7), g[i]);
            double omi = fmin(om[i], 1.0 - 1.e-7);
            slant[i + 1] = slant_depth(i + 1, nid[i + 1], path + (i + 1) * n, taun);
            double m = 1.0 / sqrt(1.e36);
            if (nid[i + 1] >= 0) {
                if (slant[i + 1] == slant[i])
                    m = sqrt(1.e36);
                else {
                    m = (tauc[i + 1] - tauc[i]) / (slant[i + 1] - slant[i]);
                    m = copysign(fmax(fabs(m), 1.0 / sqrt(1.e36)), m);
                }
            }
            mu2[i + 1] = m;
            double gam1 = (7.0 - omi * (4.0 + 3.0 * gi)) / 4.0;
            double gam2 = -(1.0 - omi * (4.0 - 3.0 * gi)) / 4.0;
            double gam3 = (2.0 - 3.0 * gi * mu) / 4.0, gam4 = 1.0 - gam3;
            double lam = sqrt(gam1 * gam1 - gam2 * gam2);
            double bgam = gam2 != 0.0 ? (gam1 - lam) / gam2 : 0.0;
            double expon = exp(-lam * taun[i]);
            e1[i] = 1.0 + bgam * expon;
            e2[i] = 1.0 - bgam * expon;
            e3[i] = bgam + expon;
            e4[i] = bgam - expon;
            double expon0 = exp(-slant[i]), expon1 = exp(-slant[i + 1]);
            double divisor = lam * lam - 1.0 / (m * m);
            divisor = copysign(fmax(1.e-3, fabs(divisor)), divisor);
            double up = omi * ((gam1 - 1.0 / m) * gam3 + gam4 * gam2) / divisor;
            double dn = omi * ((gam1 + 1.0 / m) * gam4 + gam2 * gam3) / divisor;
            cup[i] = up * expon0;
            cdn[i] = dn * expon0;
            cuptn[i] = up * expon1;
            cdntn[i] = dn * expon1;
        }
        double ssfc = rsfc * mu * exp(-slant[n]);
        a[0] = 0.0;
        b[0] = e1[0];
        d[0] = -e2[0];
        rhs[0] = -cdn[0];
        for (size_t i = 0; i + 1 < n; ++i) {
            size_t row = 2 * i + 2;
            a[row] = e2[i] * e3[i] - e4[i] * e1[i];
            b[row] = e1[i] * e1[i + 1] - e3[i] * e3[i + 1];
            d[row] = e3[i] * e4[i + 1] - e1[i] * e2[i + 1];
            rhs[row] = e3[i] * (cup[i + 1] - cuptn[i]) + e1[i] * (cdntn[i] - cdn[i + 1]);
            row = 2 * i + 1;
            a[row] = e2[i + 1] * e1[i] - e3[i] * e4[i + 1];
            b[row] = e2[i] * e2[i + 1] - e4[i] * e4[i + 1];
            d[row] = e1[i + 1] * e4[i + 1] - e2[i + 1] * e3[i + 1];
            rhs[row] = (cup[i + 1] - cuptn[i]) * e2[i + 1] - (cdn[i + 1] - cdntn[i]) * e4[i + 1];
        }
        a[2 * n - 1] = e1[n - 1] - rsfc * e3[n - 1];
        b[2 * n - 1] = e2[n - 1] - rsfc * e4[n - 1];
        d[2 * n - 1] = 0.0;
        rhs[2 * n - 1] = ssfc - cuptn[n - 1] + rsfc * cdntn[n - 1];
        status = tt_tridiagonal(2 * n, a, b, d, rhs, y, work);
        if (status)
            goto done;
        for (size_t lev = 0; lev <= n; ++lev) {
            size_t k = w * (n + 1) + n - lev;
            double up, down;
            if (!lev) {
                down = 0.0;
                up = y[0] * e3[0] - y[1] * e4[0] + cup[0];
            } else {
                size_t i = lev - 1, row = 2 * i;
                down = y[row] * e3[i] + y[row + 1] * e4[i] + cdntn[i];
                up = y[row] * e1[i] + y[row + 1] * e2[i] + cuptn[i];
            }
            out->direct_flux[k] = exp(-slant[lev]);
            out->direct_irradiance[k] = mu * out->direct_flux[k];
            out->down_irradiance[k] = down;
            out->up_irradiance[k] = up;
            out->down_flux[k] = down / 0.5;
            out->up_flux[k] = up / 0.5;
        }
    }
done:
    free(nid);
    free(storage);
    return status;
}
