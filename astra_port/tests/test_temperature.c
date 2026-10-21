#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef TT_HAVE_REFERENCE
void ref_temperature(const tt_temperature_model *, double, int, const double *, double *);
void ref_temperature_cross(const tt_temperature_model *, int, const double *, const double *, int,
                           int, const double *, const double *, const double *, double *);
#endif
static size_t comparisons;
static double worst;
static void check(size_t i, double a, double b) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(fabs(a), fabs(b)) + DBL_MIN);
    ++comparisons;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "temperature[%zu]: %.17g != %.17g (%.3g eps)\n", i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
int main(void) {
    enum { N = 127, NW = 313, NZ = 29, NT = 5 };
    double wl[N], base[N], out[NW * NZ], expected[NW * NZ], edges[NW + 1], mids[NW], temp[NZ];
    double a[NT * N], b[N], powers[NT] = {0, 1, 2, .5, 3}, sigma[N];
    for (size_t i = 0; i < N; ++i) {
        wl[i] = 180 + 3 * i;
        base[i] = 1.e-20 * (1 + i * .1);
        sigma[i] = base[i] * 2;
    }
    for (size_t i = 0; i <= NW; ++i)
        edges[i] = 150 + 1.5 * i;
    for (size_t i = 0; i < NW; ++i)
        mids[i] = .5 * (edges[i] + edges[i + 1]);
    for (size_t i = 0; i < NZ; ++i)
        temp[i] = 150 + 10 * i;
    tt_spectral_grid grid = {NW, NZ, edges, mids, temp, NULL};
    const tt_temperature_range ranges[] = {{0, 200, 210, 1}, {200, 300, 0, 0}, {300, 500, 290, 1}};
    for (int kind = 0; kind < 4; ++kind)
        for (size_t variant = 0; variant < 8; ++variant) {
            tt_temperature_model p = {0};
            p.kind = (tt_temperature_kind)kind;
            p.first = variant % 2 ? 10 : 0;
            p.count = N - 2 * p.first;
            p.terms = NT;
            p.a = a;
            p.b = b;
            p.powers = powers;
            p.sigma = sigma;
            p.base_temperature =
                kind == TT_TEMPERATURE_BURKHOLDER || kind == TT_TEMPERATURE_HARWOOD ? 0 : 273;
            p.base_wavelength = 100;
            p.base_10 = (variant & 2) != 0;
            p.invert_temperature = (variant & 4) != 0;
            p.burkholder_a = 300;
            p.burkholder_b = -.695;
            p.ranges = ranges;
            p.range_count = variant % 2 ? 3 : 0;
            for (size_t i = 0; i < N; ++i) {
                a[i] = .5 + .01 * i;
                b[i] = .2 + .03 * i;
            }
            if (kind == TT_TEMPERATURE_EXPONENTIAL) {
                const double aa[NT] = {-21, .001, -1.e-6, .002, 1.e-9},
                             bb[NT] = {.001, 1.e-5, 1.e-9, 1.e-6, 1.e-12};
                memcpy(a, aa, sizeof aa);
                memcpy(b, bb, sizeof bb);
            } else if (kind == TT_TEMPERATURE_TAYLOR) {
                for (size_t j = 0; j < NT; ++j)
                    for (size_t i = 0; i < p.count; ++i)
                        a[j * p.count + i] = pow(.001, (double)j + 1) * (1 + i * .01);
            } else if (kind == TT_TEMPERATURE_HARWOOD) {
                for (size_t i = 0; i < N; ++i) {
                    a[i] = -21 + .001 * i;
                    b[i] = -.2 - i;
                }
            }
            for (size_t z = 0; z < NZ; ++z) {
                memcpy(out, base, sizeof base);
                memcpy(expected, base, sizeof base);
                if (tt_temperature_apply(&p, temp[z], N, wl, out))
                    return 1;
#ifdef TT_HAVE_REFERENCE
                ref_temperature(&p, temp[z], N, wl, expected);
                for (size_t i = 0; i < N; ++i)
                    check(i, out[i], expected[i]);
#endif
                for (size_t i = 0; i < N; ++i) {
                    if (!isfinite(out[i]) || out[i] < 0)
                        return 1;
                    if (i < p.first || i >= p.first + p.count)
                        check(i, out[i], base[i]);
                }
            }
            if (tt_temperature_cross_section(&p, N, wl, base, &grid, out))
                return 1;
#ifdef TT_HAVE_REFERENCE
            ref_temperature_cross(&p, N, wl, base, NW, NZ, edges, mids, temp, expected);
            for (size_t i = 0; i < NW * NZ; ++i)
                check(i, out[i], expected[i]);
#endif
        }
    tt_temperature_model p = {0};
    p.kind = TT_TEMPERATURE_EXPONENTIAL;
    p.count = N;
    p.terms = 1;
    p.a = a;
    p.b = b;
    p.powers = powers;
    p.base_10 = 1;
    a[0] = -20;
    b[0] = 0;
    memset(out, 0, N * sizeof(double));
    if (tt_temperature_apply(&p, 300, N, wl, out))
        return 1;
    for (size_t i = 0; i < N; ++i)
        check(i, out[i], 1.e-20);
    p.first = N;
    if (tt_temperature_apply(&p, 300, N, wl, out) != TT_INVALID)
        return 1;
    printf("4 temperature parameterizations: %zu comparisons; worst %.3g relative epsilon-scaled "
           "error\n",
           comparisons, worst);
    return 0;
}
