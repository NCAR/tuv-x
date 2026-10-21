#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void ref_interpolate(int, int, int, const double *, const double *, int, const double *, double *);
void ref_tridiagonal(int, const double *, const double *, const double *, const double *, double *);
void ref_solar(int, int, double, double, double, double *, double *);
void ref_geometry(int, const double *, double, const double *, int *, double *, double *, double *);

static size_t comparisons;
static double worst;
static uint64_t rng = UINT64_C(0x106689d45497fdb5);
static double uniform(void) {
    rng ^= rng >> 12;
    rng ^= rng << 25;
    rng ^= rng >> 27;
    return (double)((rng * UINT64_C(0x2545f4914f6cdd1d)) >> 11) * 0x1p-53;
}
static void check(const char *label, double actual, double reference, double scale) {
    ++comparisons;
    double error =
        fabs(actual - reference) / (DBL_EPSILON * fmax(fabs(reference), scale) + DBL_MIN);
    if (!isfinite(actual) || !isfinite(reference) || error > 32.0) {
        fprintf(stderr, "%s: %.17g != %.17g (%.3g eps)\n", label, actual, reference, error);
        exit(1);
    }
    if (error > worst)
        worst = error;
}
static void ok(tt_status status) {
    if (status) {
        fprintf(stderr, "%s\n", tt_status_string(status));
        exit(1);
    }
}

int main(void) {
    enum { N = 129, M = 34 };
    double x[N], y[N], t[M], expected[N], actual[N], work[N];
    for (int trial = 0; trial < 200; ++trial) {
        x[0] = -100;
        for (int j = 1; j < N; ++j)
            x[j] = x[j - 1] + .1 + uniform();
        for (int j = 0; j < N; ++j)
            y[j] = exp(-.02 * j) * (uniform() + .1);
        for (int j = 0; j < M; ++j)
            t[j] = x[2] + (x[N - 3] - x[2]) * j / (M - 1);
        for (int method = 0; method < 4; ++method)
            for (int fold = 0; fold <= (method >= 2); ++fold) {
                ref_interpolate(method, fold, N, x, y, M, t, expected);
                ok(tt_interpolate((tt_interpolation)method, fold, N, x, y, M, t, actual));
                for (int j = 0; j < M - (method != 0); ++j)
                    check("interpolation", actual[j], expected[j], 0);
            }
    }
    double a[N], b[N], c[N], rhs[N], truth[N];
    for (int trial = 0; trial < 200; ++trial) {
        int n = 1 + trial % N;
        for (int j = 0; j < n; ++j) {
            a[j] = uniform() - .5;
            c[j] = uniform() - .5;
            b[j] = 2 + fabs(a[j]) + fabs(c[j]);
            truth[j] = uniform();
        }
        for (int j = 0; j < n; ++j)
            rhs[j] = b[j] * truth[j] + (j ? a[j] * truth[j - 1] : 0) +
                     (j + 1 < n ? c[j] * truth[j + 1] : 0);
        ref_tridiagonal(n, a, b, c, rhs, expected);
        ok(tt_tridiagonal((size_t)n, a, b, c, rhs, actual, work));
        for (int j = 0; j < n; ++j) {
            check("tridiagonal", actual[j], expected[j], 0);
            check("Fortran tridiagonal residual", expected[j], truth[j], 1);
        }
    }
    for (int year = 1950; year <= 2050; ++year)
        for (int day = 1; day <= 365; day += 13) {
            double hour = 36 * uniform() - 12, lat = 180 * uniform() - 90,
                   lon = 360 * uniform() - 180, dist, elev;
            ref_solar(year, day, hour, lat, lon, &dist, &elev);
            check("solar distance", tt_earth_sun_distance(year, day, hour), dist, 1);
            check("solar elevation", tt_solar_elevation(year, day, hour, lat, lon), elev, 90);
            if (dist < .98 || dist > 1.02 || fabs(elev) > 90)
                return 1;
        }
    enum { Z = 49 };
    double z[Z], air[Z], p[Z * (Z - 1)], pr[Z * (Z - 1)], v[Z], vr[Z], s[Z], sr[Z];
    int nid[Z], nr[Z];
    const double angles[] = {0, 20, 60, 85, 89.9, 90, 90.001, 92, 95, 100, 110, 150, 180};
    for (int trial = 0; trial < 25; ++trial) {
        z[0] = 5 * uniform();
        air[0] = 1e24;
        for (int j = 1; j < Z; ++j) {
            z[j] = z[j - 1] + .1 + 5 * uniform();
            air[j] = 1e24 * exp(-z[j] / 8);
        }
        for (size_t k = 0; k < sizeof angles / sizeof *angles; ++k) {
            ref_geometry(Z, z, angles[k], air, nr, pr, vr, sr);
            ok(tt_spherical_geometry(Z, z, angles[k], nid, p));
            ok(tt_air_mass(Z, air, nid, p, v, s));
            for (int j = 0; j < Z; ++j) {
                if (nid[j] != nr[j])
                    return 1;
                if (j + 1 < Z)
                    check("vertical column", v[j], vr[j], 0);
                check("slant column", s[j], sr[j], 0);
            }
            for (int j = 0; j < Z * (Z - 1); ++j)
                check("spherical geometry", p[j], pr[j], 0);
        }
    }
    printf("%zu C/Fortran and analytic comparisons; worst %.3g epsilon-scaled error\n", comparisons,
           worst);
    return 0;
}
