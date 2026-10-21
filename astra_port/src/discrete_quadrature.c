/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "discrete_internal.h"
#include <float.h>
#include <math.h>

tt_status tt_discrete_quadrature(size_t n, double *nodes, double *weights) {
    if (!n || n > 1000 || !nodes || !weights)
        return TT_INVALID;
    if (n == 1) {
        nodes[0] = .5;
        weights[0] = 1;
        return TT_OK;
    }
    double en = (double)n, np1 = (double)(n * (n + 1)), pi = 2 * asin(1),
           cona = (double)(n - 1) / (double)(8 * n * n * n);
    size_t limit = n / 2;
    for (size_t k = 1; k <= limit; ++k) {
        double t = (double)(4 * k - 1) * pi / (double)(4 * n + 2), x = cos(t + cona / tan(t)),
               pm2 = 0, tmp = 0;
        for (size_t iter = 1;; ++iter) {
            pm2 = 1;
            double pm1 = x, p = 0;
            for (size_t j = 2; j <= n; ++j) {
                p = ((double)(2 * j - 1) * x * pm1 - (double)(j - 1) * pm2) / (double)j;
                pm2 = pm1;
                pm1 = p;
            }
            tmp = 1 / (1 - x * x);
            double ppr = en * (pm2 - x * p) * tmp, p2pri = (2 * x * ppr - np1 * p) * tmp;
            double xi = x - (p / ppr) * (1 + (p / ppr) * p2pri / (2 * ppr));
            if (fabs(xi - x) <= 10 * DBL_EPSILON)
                break;
            if (iter > 1000)
                return TT_RANGE;
            x = xi;
        }
        nodes[k - 1] = -x;
        weights[k - 1] = 2 / (tmp * ((en * pm2) * (en * pm2)));
        nodes[n - k] = -nodes[k - 1];
        weights[n - k] = weights[k - 1];
    }
    if (n % 2) {
        nodes[limit] = 0;
        double product = 1;
        for (size_t k = 3; k <= n; k += 2)
            product = product * (double)k / (double)(k - 1);
        weights[limit] = 2 / (product * product);
    }
    for (size_t k = 0; k < n; ++k) {
        nodes[k] = .5 * nodes[k] + .5;
        weights[k] = .5 * weights[k];
    }
    return TT_OK;
}

void tt_discrete_legendre(size_t degree, size_t count, const double *cosines, double *polynomials) {
    for (size_t i = 0; i < count; ++i) {
        double *p = polynomials + i * (degree + 1);
        p[0] = 1;
        if (degree)
            p[1] = cosines[i];
        for (size_t l = 2; l <= degree; ++l)
            p[l] = ((double)(2 * l - 1) * cosines[i] * p[l - 1] - (double)(l - 1) * p[l - 2]) /
                   (double)l;
    }
}
