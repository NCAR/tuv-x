#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef TT_HAVE_REFERENCE
void ref_rates(int, int, int, const double *, const double *, const double *, const double *,
               const double *, const double *, const double *, double, double, double *);
#endif
static size_t comparisons;
static double worst;
static void check(const char *label, size_t i, double a, double b) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(fabs(a), fabs(b)) + DBL_MIN);
    ++comparisons;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "%s[%zu]: %.17g != %.17g (%.3g eps)\n", label, i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static uint64_t seed = 0xc5f05732004bULL;
static double random_value(void) {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return (seed >> 11) * 0x1p-53;
}
int main(void) {
    const size_t wavelengths[] = {1, 3, 17, 151, 1023}, levels[] = {2, 9, 81};
    for (size_t ni = 0; ni < 5; ++ni)
        for (size_t zi = 0; zi < 3; ++zi) {
            size_t nw = wavelengths[ni], nz = levels[zi], n = nw * nz;
            double *work = calloc(10 * n + 4 * nw + 1 + 8 * nz, sizeof(double));
            if (!work)
                return 1;
            double *field = work, *flux = field + 6 * n, *irrad = flux + n, *cross = irrad + n,
                   *yield = cross + n;
            double *edges = yield + n, *mids = edges + nw + 1, *solar = mids + nw,
                   *weights = solar + nw;
            double *temp = weights + nw, *air = temp + nz, *out = air + nz,
                   *expected = out + 3 * nz;
            tt_radiation radiation = {field,         field + n,     field + 2 * n,
                                      field + 3 * n, field + 4 * n, field + 5 * n};
            tt_spectral_grid grid = {nw, nz, edges, mids, temp, air};
            for (size_t z = 0; z < nz; ++z) {
                temp[z] = 170 + 3.2 * z;
                air[z] = 2.45e19 * exp(-.1 * z);
            }
            for (size_t w = 0; w <= nw; ++w)
                edges[w] = 140 + 700 * (double)w / nw;
            for (size_t w = 0; w < nw; ++w) {
                mids[w] = .5 * (edges[w] + edges[w + 1]);
                solar[w] = pow(10, 5 + 10 * random_value());
                weights[w] = random_value();
            }
            for (int trial = 0; trial < 12; ++trial) {
                for (size_t i = 0; i < 6 * n; ++i)
                    field[i] = trial == 0 ? 0 : 3 * random_value() - 1;
                if (tt_actinic_flux(nw, nz, &radiation, solar, flux) ||
                    tt_spectral_irradiance(nw, nz, &radiation, solar, mids, irrad) ||
                    tt_cross_section(TT_CS_RAYLEIGH, &grid, 0, NULL, cross))
                    return 1;
                for (size_t i = 0; i < n; ++i)
                    if (flux[i] < 0 || irrad[i] < 0)
                        return 1;
                for (int kind = 0; kind < 3; ++kind) {
                    tt_quantum_yield_kind qkind = kind == 0   ? TT_QY_NO3_AQ
                                                  : kind == 1 ? TT_QY_HO2
                                                              : TT_QY_O3_O1D;
                    double threshold = 250 + 40 * trial, scale = trial % 2 ? .37 : 2.3;
                    if (tt_quantum_yield(qkind, &grid, 0, NULL, NULL, yield) ||
                        tt_photolysis_rate(nw, nz, flux, cross, yield, scale, out) ||
                        tt_heating_rate(nw, nz, flux, cross, yield, mids, threshold, scale,
                                        out + nz) ||
                        tt_dose_rate(nw, nz, irrad, weights, out + 2 * nz))
                        return 1;
#ifdef TT_HAVE_REFERENCE
                    ref_rates(kind, (int)nw, (int)nz, edges, mids, temp, air, solar, field, weights,
                              threshold, scale, expected);
                    for (size_t i = 0; i < 3 * nz; ++i)
                        check("rate", i, out[i], expected[i]);
#endif
                    if (trial == 0)
                        for (size_t i = 0; i < 3 * nz; ++i)
                            check("darkness", i, out[i], 0);
                }
            }
            free(work);
        }
    /* A single band has an independently known reaction rate. */
    double f = 1.e12, c = 2.e-19, q = .5, result;
    if (tt_photolysis_rate(1, 1, &f, &c, &q, 3, &result))
        return 1;
    check("single band", 0, result, 3.e-7);
    if (tt_heating_rate(1, 1, &f, &c, &q, &f, -1, 1, &result) != TT_INVALID)
        return 1;
    printf("Photolysis, heating, dose: %zu comparisons; worst %.3g relative epsilon-scaled error\n",
           comparisons, worst);
    return 0;
}
