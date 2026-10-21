#include "../src/discrete_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef TT_HAVE_REFERENCE
void ref_eigen(int, const double *, double *, double *, int *);
void ref_quadrature(int, double *, double *);
void ref_legendre(int, int, const double *, double *);
#endif
static size_t comparisons;
static double worst;
static void check(const char *name, size_t i, double a, double b) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(fabs(a), fabs(b)) + DBL_MIN);
    ++comparisons;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "%s[%zu]: %.17g != %.17g (%.3g eps)\n", name, i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static uint64_t seed = 0x9d2093110025ULL;
static double random_value(void) {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return (seed >> 11) * 0x1p-53;
}
int main(void) {
    const size_t orders[] = {1, 2, 3, 4, 8, 16, 32};
    for (size_t ni = 0; ni < 7; ++ni)
        for (size_t trial = 0; trial < 24; ++trial) {
            size_t n = orders[ni];
            double *work = calloc(4 * n * n + 4 * n, sizeof(double));
            if (!work)
                return 1;
            double *a = work, *v = a + n * n, *refv = v + n * n, *scratch = refv + n * n;
            double *values = scratch + n * n + 2 * n, *refvalues = values + n;
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0; i <= j; ++i) {
                    double value = random_value() - .5;
                    if (i == j)
                        value += (double)n;
                    a[j * n + i] = value;
                    a[i * n + j] = value;
                }
            if (trial % 4 == 1)
                for (size_t j = 0; j < n; ++j)
                    for (size_t i = 0; i < n; ++i)
                        a[j * n + i] *= pow(2, (double)i - j);
            if (trial % 4 == 2 && n > 2)
                for (size_t j = 0; j < n; ++j)
                    for (size_t i = j + 1; i < n; ++i)
                        a[j * n + i] = 0;
            if (trial % 4 == 3 && n > 2) {
                for (size_t j = 1; j < n; ++j)
                    a[j * n] = a[j] = 0;
            }
            if (tt_discrete_eigen(n, a, values, v, scratch)) {
                fprintf(stderr, "eigen failed n=%zu trial=%zu\n", n, trial);
                return 1;
            }
#ifdef TT_HAVE_REFERENCE
            int info;
            ref_eigen((int)n, a, refvalues, refv, &info);
            if (info)
                return 1;
            for (size_t i = 0; i < n; ++i)
                check("eigenvalue", i, values[i], refvalues[i]);
            for (size_t i = 0; i < n * n; ++i)
                check("eigenvector", i, v[i], refv[i]);
#endif
            double norm = 0;
            for (size_t i = 0; i < n; ++i) {
                double row = 0;
                for (size_t j = 0; j < n; ++j)
                    row += fabs(a[j * n + i]);
                norm = fmax(norm, row);
            }
            for (size_t j = 0; j < n; ++j) {
                double vn = 0;
                for (size_t i = 0; i < n; ++i)
                    vn = fmax(vn, fabs(v[j * n + i]));
                if (vn == 0)
                    return 1;
                for (size_t i = 0; i < n; ++i) {
                    double av = 0;
                    for (size_t k = 0; k < n; ++k)
                        av += a[k * n + i] * v[j * n + k];
                    double residual = fabs(av - values[j] * v[j * n + i]);
                    if (!isfinite(residual) ||
                        residual > 64 * n * DBL_EPSILON * (norm + fabs(values[j])) * vn) {
                        fprintf(stderr, "eigen residual n=%zu trial=%zu vector=%zu row=%zu\n", n,
                                trial, j, i);
                        return 1;
                    }
                }
            }
            free(work);
        }
    double nodes[256], weights[256], refnodes[256], refweights[256], poly[65 * 256],
        refpoly[65 * 256];
    for (size_t n = 1; n <= 256; n = n < 32 ? n + 1 : n * 2) {
        if (tt_discrete_quadrature(n, nodes, weights))
            return 1;
#ifdef TT_HAVE_REFERENCE
        ref_quadrature((int)n, refnodes, refweights);
        for (size_t i = 0; i < n; ++i) {
            check("Gauss node", i, nodes[i], refnodes[i]);
            check("Gauss weight", i, weights[i], refweights[i]);
        }
#endif
        for (size_t p = 0; p < 2 * n && p < 16; ++p) {
            double integral = 0;
            for (size_t i = 0; i < n; ++i)
                integral += weights[i] * pow(nodes[i], (double)p);
            double exact = 1 / (double)(p + 1);
            ++comparisons;
            /* The original Newton stopping rule and n-term accumulation
             * contribute rounding beyond a single scalar operation. */
            if (fabs(integral - exact) > 64 * n * DBL_EPSILON * exact) {
                fprintf(stderr, "quadrature integral n=%zu power=%zu\n", n, p);
                return 1;
            }
        }
        tt_discrete_legendre(64, n, nodes, poly);
#ifdef TT_HAVE_REFERENCE
        ref_legendre(64, (int)n, nodes, refpoly);
        for (size_t i = 0; i < 65 * n; ++i)
            check("Legendre polynomial", i, poly[i], refpoly[i]);
#endif
    }
    double zero[4] = {0}, values[2], vectors[4], scratch[8];
    if (tt_discrete_eigen(2, zero, values, vectors, scratch))
        return 1;
    for (size_t i = 0; i < 4; ++i)
        check("zero matrix vectors", i, vectors[i], i % 3 == 0 ? 1 : 0);
    const double complex_matrix[] = {0, 1, -1, 0};
    if (tt_discrete_eigen(2, complex_matrix, values, vectors, scratch) != TT_RANGE)
        return 1;
    printf("Discrete-ordinate helpers: %zu comparisons; worst %.3g relative epsilon-scaled error\n",
           comparisons, worst);
    return 0;
}
