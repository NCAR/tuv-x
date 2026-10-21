#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef TT_HAVE_REFERENCE
void ref_lu(int, int, int, int, int, double *, size_t *, double *, double *, double *, int);
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
static uint64_t seed = 0xe932653c51005aULL;
static double random_value(void) {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return (seed >> 11) * 0x1p-53;
}
int main(void) {
    const size_t orders[] = {1, 2, 3, 7, 16, 31, 81};
    for (size_t ni = 0; ni < 7; ++ni)
        for (int layout = 0; layout < 2; ++layout)
            for (size_t trial = 0; trial < 24; ++trial) {
                size_t n = orders[ni], ml = trial == 23 ? n + 1 : trial % n,
                       mu = trial == 23 ? n + 2 : (trial * 3) % n;
                size_t ld = layout ? 2 * ml + mu + 1 : n + 3;
                double *work = calloc(3 * n * ld + n * n + 6 * n, sizeof(double));
                size_t *pivots = calloc(2 * n, sizeof(size_t));
                if (!work || !pivots)
                    return 1;
                double *original = work, *fact = original + n * ld, *expected = fact + n * ld,
                       *full = expected + n * ld;
                double *z = full + n * n, *refz = z + n, *rhs = refz + n, *b = rhs + n,
                       *refb = b + n;
                for (size_t j = 0; j < n; ++j)
                    for (size_t i = 0; i < n; ++i) {
                        if (layout && ((j > i && j - i > mu) || (i > j && i - j > ml)))
                            continue;
                        double v = 2 * random_value() - 1;
                        if (i == j)
                            v += trial % 3 == 0 ? n : 1.e-3;
                        if (trial % 6 == 1)
                            v *= pow(10, (double)j / n * 8 - 4);
                        if (trial % 6 == 2)
                            v *= 1.e-100;
                        if (trial % 6 == 3)
                            v *= 1.e100;
                        full[j * n + i] = v;
                        original[j * ld + (layout ? ml + mu + i - j : i)] = v;
                    }
                for (size_t i = 0; i < n; ++i)
                    rhs[i] = 2 * random_value() - 1;
                memcpy(fact, original, n * ld * sizeof(double));
                tt_matrix matrix = {(tt_matrix_layout)layout, n, ld, ml, mu, fact};
                double rcond;
                tt_status status = tt_lu_factor(&matrix, pivots, &rcond, z);
                if (status || !isfinite(rcond) || rcond < 0 || rcond > 1 + 32 * DBL_EPSILON) {
                    fprintf(stderr, "factor n=%zu layout=%d trial=%zu: %s\n", n, layout, trial,
                            tt_status_string(status));
                    return 1;
                }
                for (int trans = 0; trans < 2; ++trans) {
                    memcpy(b, rhs, n * sizeof(double));
                    if (tt_lu_solve(&matrix, pivots, trans, b))
                        return 1;
#ifdef TT_HAVE_REFERENCE
                    memcpy(expected, original, n * ld * sizeof(double));
                    memcpy(refb, rhs, n * sizeof(double));
                    double refrcond;
                    ref_lu(layout, (int)n, (int)ld, (int)ml, (int)mu, expected, pivots + n,
                           &refrcond, refz, refb, trans);
                    for (size_t i = 0; i < n * ld; ++i)
                        check("LU factors", i, fact[i], expected[i]);
                    for (size_t i = 0; i < n; ++i) {
                        if (pivots[i] != pivots[n + i])
                            return 1;
                        check("condition vector", i, z[i], refz[i]);
                        check("solution", i, b[i], refb[i]);
                    }
                    check("reciprocal condition", 0, rcond, refrcond);
#endif
                    double matrix_norm = 0, solution_norm = 0, rhs_norm = 0;
                    for (size_t i = 0; i < n; ++i) {
                        double row_sum = 0;
                        for (size_t j = 0; j < n; ++j)
                            row_sum += fabs(full[trans ? i * n + j : j * n + i]);
                        matrix_norm = fmax(matrix_norm, row_sum);
                        solution_norm = fmax(solution_norm, fabs(b[i]));
                        rhs_norm = fmax(rhs_norm, fabs(rhs[i]));
                    }
                    double scale = matrix_norm * solution_norm + rhs_norm;
                    for (size_t i = 0; i < n; ++i) {
                        double sum = 0;
                        for (size_t j = 0; j < n; ++j) {
                            double term = full[trans ? i * n + j : j * n + i] * b[j];
                            sum += term;
                        }
                        if (fabs(sum - rhs[i]) > 32 * n * DBL_EPSILON * scale) {
                            fprintf(stderr,
                                    "large backward error n=%zu row=%zu layout=%d trial=%zu "
                                    "trans=%d got=%.17g expected=%.17g scale=%.17g\n",
                                    n, i, layout, trial, trans, sum, rhs[i], scale);
                            return 1;
                        }
                    }
                }
                free(work);
                free(pivots);
            }
    double singular[] = {1, 2, 2, 4}, condition, work[2], rhs[] = {1, 2};
    size_t pivots[2];
    tt_matrix a = {TT_MATRIX_DENSE, 2, 2, 0, 0, singular};
    if (tt_lu_factor(&a, pivots, &condition, work) != TT_SINGULAR || condition != 0)
        return 1;
    if (tt_lu_solve(&a, pivots, 0, rhs) != TT_SINGULAR)
        return 1;
    a.leading_dimension = 1;
    if (tt_lu_factor(&a, pivots, NULL, NULL) != TT_INVALID)
        return 1;
    printf("Dense/banded LINPACK: %zu comparisons; worst %.3g relative epsilon-scaled error\n",
           comparisons, worst);
    return 0;
}
