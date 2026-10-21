#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef TT_HAVE_REFERENCE
void ref_yield(int, int, int, const double *, const double *, const double *, const double *, int,
               const tt_spectral_table *, const tt_yield_options *, double *);
void ref_weight(int, int, const double *, const double *, const double *, double, double, double *);
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
static int order(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
int main(void) {
    enum { NW = 2048, NZ = 30, NP = 3 };
    double edges[NW + 1], mids[NW], data[NP][4 * NW], out[NW * NZ], expected[NW * NZ], rates[NW];
    double temp[NZ], air[NZ];
    const double knots[] = {193, 248, 263.4, 279, 290, 298, 300, 302, 305, 308, 313, 327,
                            328, 330, 340,   360, 364, 366, 380, 390, 400, 440, 700};
    size_t nw = 0;
    for (double w = 140; w < 850; w += .5)
        mids[nw++] = w;
    for (size_t i = 0; i < sizeof knots / sizeof *knots; ++i) {
        mids[nw++] = nextafter(knots[i], -INFINITY);
        mids[nw++] = knots[i];
        mids[nw++] = nextafter(knots[i], INFINITY);
    }
    qsort(mids, nw, sizeof(double), order);
    size_t unique = 1;
    for (size_t i = 1; i < nw; ++i)
        if (mids[i] != mids[unique - 1])
            mids[unique++] = mids[i];
    nw = unique;
    for (size_t w = 0; w < nw; ++w) {
        edges[w] = mids[w];
        rates[w] = w % 3 == 0 ? -1 : 1.e4 * (1 + w % 13);
    }
    edges[nw] = 851;
    double temperatures[NP][4] = {{180, 240, 300, 360}, {190, 250, 310, 370}, {200, 260, 320, 380}};
    const double densities[] = {0,       1.e8,   1.e15,   8.e17, 8.1e17,
                                2.45e19, 2.6e19, 2.61e19, 1.e20, 1.e22};
    for (size_t z = 0; z < NZ; ++z) {
        temp[z] = 150 + 15 * z;
        air[z] = densities[z % 10];
    }
    tt_spectral_grid grid = {nw, NZ, edges, mids, temp, air};
    tt_spectral_table tables[NP];
    tt_yield_options options = tt_yield_options_default();
    options.molecular_diameter = 4.e-10;
    options.molecular_weight = .098;
    options.collision_rates = rates;
    const double coefficients[] = {1.1, -.003, 2.e-6, 4.e-9};
    options.coefficients = coefficients;
    options.coefficient_count = 4;
    options.minimum_wavelength = 279;
    options.maximum_wavelength = 440;
    for (int kind = 0; kind < TT_QY_COUNT; ++kind) {
        for (size_t p = 0; p < NP; ++p) {
            tables[p] = (tt_spectral_table){nw, 4, 4, data[p], temperatures[p]};
            for (size_t c = 0; c < 4; ++c)
                for (size_t w = 0; w < nw; ++w)
                    data[p][c * nw + w] =
                        kind == TT_QY_TAYLOR ? .2 : .05 * (1 + c) + .001 * (w % 17) + .003 * p;
        }
        for (unsigned variant = 1; variant <= (kind == TT_QY_ACETONE ? 6 : 1); ++variant) {
            options.acetone_branches = 1 + (variant - 1) % 3;
            options.minimum_temperature = variant > 3 ? 220 : 0;
            options.maximum_temperature = variant > 3 ? 310 : DBL_MAX;
            tt_status status =
                tt_quantum_yield((tt_quantum_yield_kind)kind, &grid, NP, tables, &options, out);
            if (status) {
                fprintf(stderr, "yield %d: %s\n", kind, tt_status_string(status));
                return 1;
            }
            for (size_t i = 0; i < nw * NZ; ++i)
                if (!isfinite(out[i]))
                    return 1;
#ifdef TT_HAVE_REFERENCE
            ref_yield(kind, (int)nw, NZ, edges, mids, temp, air, NP, tables, &options, expected);
            for (size_t i = 0; i < nw * NZ; ++i)
                check(tt_quantum_yield_name((tt_quantum_yield_kind)kind), i, out[i], expected[i]);
#endif
        }
    }
    for (int kind = 0; kind < TT_WEIGHT_COUNT; ++kind) {
        tt_status status =
            tt_spectral_weight((tt_weight_kind)kind, nw, edges, mids, data[0], 305, 340, out);
        if (status) {
            fprintf(stderr, "weight %d: %s\n", kind, tt_status_string(status));
            return 1;
        }
        for (size_t i = 0; i < nw; ++i)
            if (!isfinite(out[i]))
                return 1;
#ifdef TT_HAVE_REFERENCE
        ref_weight(kind, (int)nw, edges, mids, data[0], 305, 340, expected);
        for (size_t i = 0; i < nw; ++i)
            check("spectral weight", kind * nw + i, out[i], expected[i]);
#endif
        if (kind == TT_WEIGHT_GAUSSIAN) {
            double sum = 0;
            for (size_t i = 0; i < nw; ++i)
                sum += out[i];
            check("gaussian normalization", 0, sum, 1);
        }
        if (kind == TT_WEIGHT_EPPLEY) {
            double sum = 0;
            for (size_t i = 0; i < nw; ++i)
                sum += out[i] * (edges[i + 1] - edges[i]);
            check("Eppley normalization", 0, sum, 90);
        }
    }
    /* Complementary branches sum to unity at every temperature and boundary. */
    for (int branch = 0; branch < 3; ++branch) {
        tt_quantum_yield_kind a = branch == 0   ? TT_QY_O3_O1D
                                  : branch == 1 ? TT_QY_CLO_O1D
                                                : TT_QY_CLONO2_CL;
        if (tt_quantum_yield(a, &grid, 0, NULL, NULL, out) ||
            tt_quantum_yield((tt_quantum_yield_kind)(a + 1), &grid, 0, NULL, NULL, expected))
            return 1;
        for (size_t i = 0; i < nw * NZ; ++i)
            check("complementary branches", i, out[i] + expected[i], 1);
    }
    printf("21 quantum yields and 13 weights: %zu comparisons; worst %.3g relative epsilon-scaled "
           "error\n",
           comparisons, worst);
    return 0;
}
