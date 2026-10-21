/* Numerical conventions retained from TUV-x (Apache-2.0). */
#include "turbotuv.h"
#include <math.h>

const char *tt_status_string(tt_status status) {
    static const char *const messages[] = {
        "success",           "invalid argument", "outside supported range",
        "allocation failed", "singular system",  "not found"};
    return (unsigned)status < sizeof messages / sizeof *messages ? messages[status]
                                                                 : "unknown status";
}

static int increasing(size_t n, const double *x) {
    if (!x || n < 2 || !isfinite(x[0]))
        return 0;
    for (size_t i = 1; i < n; ++i)
        if (!isfinite(x[i]) || x[i] <= x[i - 1])
            return 0;
    return 1;
}

tt_status tt_interpolate(tt_interpolation method, int fold, size_t nx, const double *x,
                         const double *y, size_t nt, const double *target, double *out) {
    int valid_target =
        method == TT_LINEAR && nt == 1 ? target && isfinite(target[0]) : increasing(nt, target);
    if (!y || !out || !increasing(nx, x) || !valid_target || method < TT_LINEAR ||
        method > TT_FRACTIONAL_TARGET || (fold && method < TT_FRACTIONAL_SOURCE))
        return TT_INVALID;
    if (method == TT_CONSERVING && (x[0] > target[0] || x[nx - 1] < target[nt - 1]))
        return TT_RANGE;
    size_t start = 0, j = 0;
    if (method == TT_LINEAR) {
        for (size_t i = 0; i < nt; ++i) {
            out[i] = 0.0;
            for (j = start; j + 1 < nx; ++j) {
                if (target[i] < x[j] || target[i] >= x[j + 1])
                    continue;
                double slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
                out[i] = y[j] + slope * (target[i] - x[j]);
                start = j;
                break;
            }
        }
        return TT_OK;
    }
    for (size_t i = 0; i + 1 < nt; ++i) {
        double sum = 0.0;
        j = start;
        if (method == TT_CONSERVING) {
            while (j + 1 < nx && x[j + 1] <= target[i])
                ++j;
            start = j;
            for (; j + 1 < nx && x[j] < target[i + 1]; ++j) {
                double a1 = fmax(x[j], target[i]);
                double a2 = fmin(x[j + 1], target[i + 1]);
                double slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
                double b1 = y[j] + slope * (a1 - x[j]);
                double b2 = y[j] + slope * (a2 - x[j]);
                sum = sum + 0.5 * (a2 - a1) * (b2 + b1);
            }
            out[i] = sum / (target[i + 1] - target[i]);
        } else {
            while (j + 1 < nx && x[j + 1] < target[i])
                ++j;
            start = j;
            for (; j + 1 < nx && x[j] <= target[i + 1]; ++j) {
                double a1 = fmax(x[j], target[i]);
                double a2 = fmin(x[j + 1], target[i + 1]);
                if (method == TT_FRACTIONAL_SOURCE)
                    sum = sum + y[j] * (a2 - a1) / (x[j + 1] - x[j]);
                else
                    sum = sum + y[j] * (a2 - a1);
            }
            out[i] = method == TT_FRACTIONAL_SOURCE ? sum : sum / (target[i + 1] - target[i]);
        }
    }
    if (fold) {
        /* Deliberately preserve TUV-x's unusual mixed bin-density/integral
         * convention. Changing it changes existing radiation calculations. */
        if (!j)
            return TT_RANGE;
        --j;
        double a1 = target[nt - 1], a2 = x[j + 1];
        if (a2 > a1 || j + 2 < nx) {
            double tail = y[j] * (a2 - a1) / (x[j + 1] - x[j]);
            for (size_t k = j + 1; k + 1 < nx; ++k)
                tail = tail + y[k] * (x[k + 1] - x[k]);
            out[nt - 2] = out[nt - 2] + tail;
        }
    }
    return TT_OK;
}

tt_status tt_tridiagonal(size_t n, const double *a, const double *b, const double *c,
                         const double *r, double *u, double *cp) {
    if (!n || !a || !b || !c || !r || !u || !cp)
        return TT_INVALID;
    if (b[0] == 0.0)
        return TT_SINGULAR;
    cp[0] = c[0] / b[0];
    u[0] = r[0] / b[0];
    for (size_t i = 1; i < n; ++i) {
        double pivot = b[i] - a[i] * cp[i - 1];
        if (pivot == 0.0)
            return TT_SINGULAR;
        double denom = 1.0 / pivot;
        cp[i] = c[i] * denom;
        u[i] = (r[i] - a[i] * u[i - 1]) * denom;
    }
    for (size_t i = n - 1; i > 0; --i)
        u[i - 1] = u[i - 1] - cp[i - 1] * u[i];
    return TT_OK;
}
