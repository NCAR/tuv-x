#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/generated/oxygen.inc"

#ifdef TT_HAVE_REFERENCE
void ref_oxygen(int, int, const double *, const double *, const double *, const double *, double,
                double *, double *, double *);
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
int main(void) {
    const double edges[] = {100,   121.4, 121.9, 150,   175.4, 177,   178.6, 180.2,
                            181.8, 183.5, 185.2, 186.9, 188.7, 190.5, 192.3, 194.2,
                            196.1, 198,   200,   202,   204.1, 206.2, 220,   350};
    enum { NW = 23, NZ = 127 };
    double temp[NZ], vertical[NZ - 1], slant[NZ], oxygen[NZ], depth[NW * (NZ - 1)], cross[NW * NZ];
    double expected_depth[NW * (NZ - 1)], expected_cross[NW * NZ], coefficients[680];
    tt_oxygen_bands bands;
    if (tt_oxygen_band_grid(NW, edges, &bands) || bands.lyman_alpha != 1 ||
        bands.schumann_runge != 4)
        return 1;
    const size_t levels[] = {2, 3, 11, 81, 127};
    for (size_t li = 0; li < 5; ++li)
        for (int trial = 0; trial < 20; ++trial) {
            size_t nz = levels[li];
            double scale = trial % 2 ? .2095 : 1;
            for (size_t z = 0; z < nz; ++z) {
                temp[z] = 150 + (trial % 5) * 40 + .5 * z;
                double log_column = 60 - 25 * (double)z / (nz - 1);
                if (trial % 5 == 1)
                    log_column = 55 - 15 * (double)z / (nz - 1);
                if (trial % 5 == 2)
                    log_column = 57 + 2 * (1 - (double)z / (nz - 1));
                if (trial % 5 == 3)
                    log_column = 45 - 1.e-8 * z;
                if (trial % 5 == 4)
                    log_column = 56 - 18 * (double)z / (nz - 1);
                slant[z] = exp(log_column) / scale;
                if (trial >= 15 && z == 0)
                    slant[z] = 1.e36;
                oxygen[z] = scale * slant[z];
                if (z + 1 < nz)
                    vertical[z] = slant[z] / (1 + .05 * z + .7 * trial);
            }
            for (size_t i = 0; i < NW * (nz - 1); ++i)
                depth[i] = expected_depth[i] = .125;
            for (size_t i = 0; i < NW * nz; ++i)
                cross[i] = expected_cross[i] = 1.1e-21;
            if (tt_oxygen_optics(&bands, NW, nz, temp, vertical, slant, oxygen, depth, cross))
                return 1;
#ifdef TT_HAVE_REFERENCE
            ref_oxygen(NW, (int)nz, edges, temp, vertical, slant, scale, expected_depth,
                       expected_cross, coefficients);
            for (size_t i = 0; i < NW * (nz - 1); ++i)
                check("oxygen depth", i, depth[i], expected_depth[i]);
            for (size_t i = 0; i < NW * nz; ++i)
                check("oxygen cross section", i, cross[i], expected_cross[i]);
            for (size_t w = 0; w < 17; ++w)
                for (size_t p = 0; p < 20; ++p) {
                    check("Chebyshev A", w * 20 + p, oxygen_a[w][p], coefficients[w * 20 + p]);
                    check("Chebyshev B", w * 20 + p, oxygen_b[w][p],
                          coefficients[340 + w * 20 + p]);
                }
#endif
            for (size_t w = 0; w < NW; ++w) {
                for (size_t z = 0; z < nz; ++z)
                    if (!isfinite(cross[w * nz + z]) || cross[w * nz + z] <= 0)
                        return 1;
                if (w != 1 && (w < 4 || w >= 21)) {
                    for (size_t z = 0; z < nz; ++z)
                        check("unchanged cross section", w * nz + z, cross[w * nz + z], 1.1e-21);
                    for (size_t z = 0; z + 1 < nz; ++z)
                        check("unchanged depth", w * (nz - 1) + z, depth[w * (nz - 1) + z], .125);
                }
            }
        }
    /* The original indexes outside its array when every column is tiny.
     * The C implementation defines this case using its tabulated limits. */
    for (size_t z = 0; z < NZ; ++z) {
        temp[z] = 220;
        slant[z] = oxygen[z] = 1.e10;
        if (z + 1 < NZ)
            vertical[z] = 1.e10;
    }
    if (tt_oxygen_optics(&bands, NW, NZ, temp, vertical, slant, oxygen, depth, cross))
        return 1;
    check("thin SR limit", 0, cross[4 * NZ], 6.2180730e-21);
    for (size_t w = 4; w < 21; ++w)
        for (size_t z = 0; z < NZ; ++z)
            check("thin uniform limit", z, cross[w * NZ + z], cross[w * NZ]);
    const double absent[] = {250, 300, 500}, bad[] = {100, 122, 200, 220};
    if (tt_oxygen_band_grid(2, absent, &bands) || bands.lyman_alpha != SIZE_MAX ||
        bands.schumann_runge != SIZE_MAX)
        return 1;
    if (tt_oxygen_band_grid(3, bad, &bands) != TT_RANGE)
        return 1;
    printf("O2 bands: %zu comparisons; worst %.3g relative epsilon-scaled error\n", comparisons,
           worst);
    return 0;
}
