#include "fixtures.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef TT_HAVE_REFERENCE
void ref_aerosol(int, int, const double *, const double *, int, const double *, double, double,
                 double, double *);
void ref_gas(int, int, const double *, const double *, const double *, const double *, int,
             double *);
void ref_radiator_dataset(double *);
#endif
static size_t checked;
static int scenario;
static double worst;
static void ok(tt_status status) {
    if (status) {
        fprintf(stderr, "case %d: %s\n", scenario, tt_status_string(status));
        exit(1);
    }
}
static void same(const char *name, size_t i, double a, double b) {
    double error = b ? fabs((a - b) / b) / DBL_EPSILON : (a == b ? 0 : INFINITY);
    ++checked;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "case %d %s[%zu]: %.17g != %.17g (%.3g eps)\n", scenario, name, i, a, b,
                error);
        exit(1);
    }
    worst = fmax(worst, error);
}
int main(void) {
    const size_t layers[] = {1, 2, 5, 16, 81, 120}, samples[] = {2, 7, 51, 121};
    enum { NW = 11, NL = 120 };
    double z[NL + 1], w[NW + 1], depth[NL + 1], out[3 * NW * NL], expected[3 * NW * NL];
    double column[NL], sigma[NW], cross[NL * NW];
    for (size_t i = 0; i <= NW; ++i)
        w[i] = 120 + (double)i * 100;
    for (size_t ni = 0; ni < 6; ++ni)
        for (size_t di = 0; di < 4; ++di)
            for (int variant = 0; variant < 6; ++variant)
                for (int top = 0; top < 3; ++top) {
                    ++scenario;
                    size_t n = layers[ni], nd = samples[di], size = n * NW;
                    double ceiling = top == 0 ? 1.3 : top == 1 ? 17.5 : 120;
                    double a = (double)variant / 5, g = (double)variant / 5 * 2 - 1;
                    double tau550 = variant == 0   ? -1
                                    : variant == 1 ? 0
                                    : variant == 2 ? 1.e-20
                                                   : .235;
                    for (size_t i = 0; i <= n; ++i)
                        z[i] = ceiling * (double)i / (double)n;
                    for (size_t i = 0; i < nd; ++i)
                        depth[i] = variant == 4 ? 0 : exp(-(double)i / 8) * .1;
                    tt_grid height = {0}, wave = {0};
                    ok(tt_grid_from_edges(n + 1, z, &height));
                    ok(tt_grid_from_edges(NW + 1, w, &wave));
                    ok(tt_aerosol_optics(&height, &wave, nd, depth, a, g, tau550, 1, out,
                                         out + size, out + 2 * size));
#ifdef TT_HAVE_REFERENCE
                    ref_aerosol((int)n, NW, z, w, (int)nd, depth, a, g, tau550, expected);
                    for (size_t i = 0; i < 3 * size; ++i)
                        same("aerosol", i, out[i], expected[i]);
#endif
                    for (size_t j = 0; j < size; ++j) {
                        if (out[j] < 0 || out[size + j] < 0 ||
                            out[size + j] > 1 + 32 * DBL_EPSILON ||
                            fabs(out[2 * size + j]) > 1 + 32 * DBL_EPSILON)
                            return 1;
                    }
                    for (size_t i = 0; i < n; ++i)
                        column[i] = exp(-(double)i / 8) * 1.e24;
                    for (size_t j = 0; j < NW; ++j) {
                        sigma[j] = (double)j * 3.e-20;
                        for (size_t i = 0; i < n; ++i)
                            cross[j * n + i] = sigma[j];
                    }
                    for (int air = 0; air < 2; ++air) {
                        ok(tt_gas_optics(n, NW, column, cross, air, out, out + size,
                                         out + 2 * size));
#ifdef TT_HAVE_REFERENCE
                        ref_gas((int)n, NW, z, w, column, sigma, air, expected);
                        for (size_t i = 0; i < 3 * size; ++i)
                            same("gas", i, out[i], expected[i]);
#endif
                        for (size_t i = 0; i < size; ++i) {
                            same("gas albedo", i, out[size + i], air);
                            same("gas asymmetry", i, out[2 * size + i], 0);
                        }
                    }
                    tt_grid_destroy(&height);
                    tt_grid_destroy(&wave);
                }
    tt_radiator radiator;
    const tt_dataset *data = tt_test_dataset("test/data/radiator.nc");
    ok(tt_radiator_from_dataset(data, 3, 4, &radiator));
#ifdef TT_HAVE_REFERENCE
    ref_radiator_dataset(expected);
    const double *fields[] = {radiator.optical_depth, radiator.single_scattering_albedo,
                              radiator.asymmetry};
    for (size_t j = 0; j < 3; ++j)
        for (size_t i = 0; i < 12; ++i)
            same("embedded radiator", i, fields[j][i], expected[j * 12 + i]);
#endif
    /* The original test's independent expected values. */
    const double stored[3][12] = {{0, 10, 100, 1, 11, 101, 2, 12, 102, 3, 13, 103},
                                  {.1, .11, .111, .2, .12, .112, .3, .13, .113, .4, .14, .114},
                                  {20, 10, 4, 19, 9, 3, 18, 8, 2, 17, 7, 1}};
    for (size_t i = 0; i < 12; ++i) {
        same("stored depth", i, radiator.optical_depth[i], stored[0][i]);
        same("stored albedo", i, radiator.single_scattering_albedo[i], stored[1][i]);
        same("stored asymmetry", i, radiator.asymmetry[i], stored[2][i]);
    }
    if (tt_radiator_from_dataset(data, 4, 3, &radiator) != TT_INVALID)
        return 1;
    tt_grid height = {0}, wave = {0};
    ok(tt_grid_equal(0, 20, 2, &height));
    ok(tt_grid_equal(500, 600, 100, &wave));
    double scalar = .3;
    ok(tt_aerosol_optics(&height, &wave, 1, &scalar, .9, .7, -1, 1, out, out + 10, out + 20));
    for (size_t i = 0; i < 10; ++i) {
        same("uniform aerosol depth", i, out[i], scalar * (340. / 550.));
        same("uniform aerosol albedo", i, out[10 + i], .9);
        same("uniform aerosol asymmetry", i, out[20 + i], .7);
    }
    tt_grid_destroy(&height);
    tt_grid_destroy(&wave);
    printf(
        "Optical properties: %d cases, %zu comparisons; worst %.3g relative epsilon-scaled error\n",
        scenario, checked, worst);
    return 0;
}
