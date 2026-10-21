#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_host_grid(int, int, const double *, const double *, double *);
void ref_host_profile(int, int, const double *, const double *, const double *, double, double,
                      double *);
void ref_grid(int, int, const double *, double, double, double, const char *, int, double *);
void ref_profile(int, int, const double *, const double *, double, double, double, double *);
void ref_solar_profiles(int, const double *, int, int, int, double, double, double, double *,
                        double *);
void ref_solar_flux(int, const double *, int, const int *, double *);
#endif
static size_t checked;
static double worst;
static int scenario;
static void ok(tt_status status) {
    if (status) {
        fprintf(stderr, "case %d: %s\n", scenario, tt_status_string(status));
        exit(1);
    }
}
static void same(const char *name, size_t i, double a, double b, double scale) {
    double norm = fmax(fabs(b), scale),
           error = norm ? fabs(a - b) / norm / DBL_EPSILON : (a == b ? 0 : INFINITY);
    ++checked;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "case %d %s[%zu]: %.17g != %.17g (%.3g eps)\n", scenario, name, i, a, b,
                error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static const tt_dataset *embedded(const char *path) {
    char name[256];
    snprintf(name, sizeof name, "data/%s", path);
    return tt_dataset_find(name);
}
static void check_grid(const tt_grid *grid, int kind, double first, double last, double delta,
                       const char *path) {
    size_t n = grid->cells;
    ++scenario;
#ifdef TT_HAVE_REFERENCE
    double *expected = calloc(3 * (n + 1), sizeof(double));
    ref_grid(kind, (int)n, grid->edges, first, last, delta, path, (int)strlen(path), expected);
    for (size_t i = 0; i <= n; ++i)
        same("grid edge", i, grid->edges[i], expected[i], 0);
    for (size_t i = 0; i < n; ++i) {
        same("grid midpoint", i, grid->midpoints[i], expected[n + 1 + i], 0);
        same("grid width", i, grid->widths[i], expected[2 * (n + 1) + i], 0);
    }
    free(expected);
#endif
    for (size_t i = 0; i < n; ++i) {
        same("grid width identity", i, grid->edges[i] + grid->widths[i], grid->edges[i + 1], 1);
        if (!(grid->edges[i] < grid->edges[i + 1]))
            exit(1);
    }
}
static void check_profile(const tt_grid *grid, const tt_profile *p, const double *expected,
                          int kind, double angle_scale) {
    size_t n = grid->cells;
    const double *fields[] = {p->edges,   p->midpoints,        p->differences,
                              p->columns, p->extended_columns, p->burden};
#ifdef TT_HAVE_REFERENCE
    for (size_t f = 0; f < 6; ++f) {
        if ((kind == 2 && f == 5) || (kind == 3 && f >= 4) || (kind == 5 && f >= 3))
            continue;
        size_t count = n + ((f == 0 || f == 4) ? 1 : 0);
        for (size_t i = 0; i < count; ++i)
            same("profile", f * (n + 1) + i, fields[f][i], expected[f * (n + 1) + i], angle_scale);
    }
#endif
    for (size_t i = 0; i < n; ++i) {
        if (kind != 6)
            same("profile midpoint", i, p->midpoints[i], .5 * (p->edges[i] + p->edges[i + 1]), 0);
        same("profile difference", i, p->differences[i], p->edges[i + 1] - p->edges[i], 0);
    }
    if (n && kind != 5) {
        same("top column", n, p->columns[n - 1],
             p->extended_columns[n - 1] + p->extended_columns[n], 0);
        double total = 0;
        for (size_t i = n; i-- > 0;)
            total += p->columns[i];
        same("column burden", 0, p->burden[0], total, 0);
    }
}
static void host_checks(void) {
    for (size_t n = 0; n < 245; n += 7)
        for (int mode = 0; mode < 5; ++mode) {
            if (!n && mode)
                continue;
            ++scenario;
            double *z = malloc((n + 1) * sizeof(double)), *v = malloc((n + 1) * sizeof(double));
            double *m = malloc((n + 1) * sizeof(double)), *c = malloc((n + 1) * sizeof(double));
            double *expected = calloc(6 * (n + 1), sizeof(double));
            for (size_t i = 0; i <= n; ++i) {
                z[i] = (double)i * 1.03;
                v[i] = 1.e18 * exp(-(double)i / 8);
                m[i] = 1.e17 / (double)(i + 1);
                c[i] = m[i] * 2.e5;
            }
            tt_grid grid = {0};
            tt_profile p = {0};
            ok(tt_grid_from_edges(n + 1, z, &grid));
            ok(tt_grid_update(&grid, z, mode % 2 ? m : NULL));
#ifdef TT_HAVE_REFERENCE
            ref_host_grid((int)n, mode % 2, z, m, expected);
            for (size_t i = 0; i <= n; ++i)
                same("host grid edge", i, grid.edges[i], expected[i], 0);
            for (size_t i = 0; i < n; ++i) {
                same("host grid midpoint", i, grid.midpoints[i], expected[n + 1 + i], 0);
                same("host grid width", i, grid.widths[i], expected[2 * (n + 1) + i], 0);
            }
#endif
            double scale = 8.01, exo = 3.7e16;
            ok(tt_profile_create(n, &p));
            ok(tt_profile_update(&p, v, mode == 1 ? m : NULL, mode ? c : NULL,
                                 mode == 2 || mode == 4 ? &scale : NULL, mode == 3 ? &exo : NULL));
            if (mode == 4)
                ok(tt_profile_update(&p, v, m, NULL, NULL, NULL));
#ifdef TT_HAVE_REFERENCE
            ref_host_profile((int)n, mode, v, m, c, scale, exo, expected);
#endif
            check_profile(&grid, &p, expected, mode == 1 || mode == 4 ? 6 : 0, 0);
            tt_profile_destroy(&p);
            tt_grid_destroy(&grid);
            free(z);
            free(v);
            free(m);
            free(c);
            free(expected);
        }
}

static void solar_boundaries(const tt_dataset *const *solar) {
    double edges[15] = {0};
    size_t count = 1;
    const double cuts[] = {150.01, 200.07, 1000.99};
    for (size_t i = 0; i < 3; ++i) {
        edges[count++] = nextafter(cuts[i], 0);
        edges[count++] = cuts[i];
        edges[count++] = nextafter(cuts[i], INFINITY);
    }
    edges[count++] = 1500;
    tt_grid grid = {0};
    ok(tt_grid_from_edges(count, edges, &grid));
    double actual[15], expected[15];
    for (int method = 1; method <= 3; ++method) {
        ++scenario;
        tt_interpolation methods[4] = {(tt_interpolation)method, (tt_interpolation)method,
                                       (tt_interpolation)method, TT_FRACTIONAL_TARGET};
        int m[4] = {method, method, method, 3};
        ok(tt_extraterrestrial_flux(&grid, 4, solar, methods, actual));
#ifdef TT_HAVE_REFERENCE
        ref_solar_flux((int)grid.cells, grid.edges, 4, m, expected);
        for (size_t w = 0; w < grid.cells; ++w)
            same("solar splice boundary", w, actual[w], expected[w], 0);
#endif
    }
    tt_grid_destroy(&grid);
}

int main(void) {
    const char *grid_names[] = {"grids/wavelength/combined.grid", "grids/wavelength/cam.csv",
                                "grids/wavelength/kockarts.grid"};
    const char *solar_names[] = {
        "profiles/solar/susim_hi.flx", "profiles/solar/atlas3_1994_317_a.dat",
        "profiles/solar/sao2010.solref.converted", "profiles/solar/neckel.flx"};
    const tt_dataset *solar[4];
    for (size_t i = 0; i < 4; ++i)
        solar[i] = embedded(solar_names[i]);
    for (size_t i = 0; i < 3; ++i) {
        tt_grid grid = {0};
        ok(tt_grid_from_dataset(embedded(grid_names[i]), &grid));
        check_grid(&grid, 2, 0, 0, 0, grid_names[i]);
        double *actual = malloc(grid.cells * sizeof(double)),
               *expected = malloc(grid.cells * sizeof(double));
        for (size_t nf = 1; nf <= 4; ++nf) {
            ++scenario;
            tt_interpolation methods[] = {TT_CONSERVING, TT_CONSERVING, TT_CONSERVING,
                                          TT_FRACTIONAL_TARGET};
            int m[] = {1, 1, 1, 3};
            ok(tt_extraterrestrial_flux(&grid, nf, solar, methods, actual));
#ifdef TT_HAVE_REFERENCE
            ref_solar_flux((int)grid.cells, grid.edges, (int)nf, m, expected);
            for (size_t w = 0; w < grid.cells; ++w)
                same("solar spectrum", w, actual[w], expected[w], 0);
#endif
            for (size_t w = 0; w < grid.cells; ++w)
                if (!isfinite(actual[w]) || actual[w] < 0)
                    exit(1);
        }
        free(actual);
        free(expected);
        tt_grid_destroy(&grid);
    }
    for (size_t i = 0; i < 80; ++i) {
        double first = (double)i * .03125 - 1, delta = .1 + (double)i * .0078125;
        double last = first + delta * (double)(i * 3) + ((i % 3) ? delta * .37 : 0);
        tt_grid grid = {0};
        ok(tt_grid_equal(first, last, delta, &grid));
        check_grid(&grid, 1, first, last, delta, "");
        check_grid(&grid, 0, 0, 0, 0, "");
        tt_grid_destroy(&grid);
    }
    const size_t sizes[] = {0, 1, 2, 5, 16, 81, 120, 240};
    for (size_t si = 0; si < 8; ++si)
        for (int kind = 0; kind < 5; ++kind)
            for (int variant = 0; variant < 6; ++variant) {
                size_t n = sizes[si];
                if (!n && kind != 0)
                    continue;
                ++scenario;
                double *z = malloc((n + 1) * sizeof(double)), *v = malloc((n + 1) * sizeof(double));
                double *expected = calloc(6 * (n + 1), sizeof(double));
                for (size_t i = 0; i <= n; ++i) {
                    z[i] = 2.3 +
                           (n ? ((variant == 5 ? 150. : 120.) - 2.3) * (double)i / (double)n : 0);
                    v[i] = variant % 2 ? 4.5e19 : sin((double)i) * 2.5e17;
                }
                tt_grid grid = {0};
                tt_profile p = {0};
                ok(tt_grid_from_edges(n + 1, z, &grid));
                tt_profile_options options = tt_profile_options_default(
                    kind == 4 ? TT_PROFILE_LINEAR : (tt_profile_kind)kind);
                options.length_conversion = kind == 0 && variant % 2 ? 1 : 1.e5;
                if (kind != 1)
                    options.scale_height = variant == 0 ? 0 : variant == 1 ? .001 : 4.5;
                options.reference_column = variant == 0   ? 1
                                           : variant == 1 ? 0
                                                          : 50 + (double)variant * 100;
                if (kind == 0)
                    ok(tt_profile_from_edges(&grid, v, &options, &p));
                else {
                    const char *path = kind == 3   ? "profiles/atmosphere/ussa.ozone"
                                       : kind == 4 ? "profiles/atmosphere/ussa.temp"
                                                   : "profiles/atmosphere/ussa.dens";
                    ok(tt_profile_from_dataset(&grid, embedded(path), &options, &p));
                }
#ifdef TT_HAVE_REFERENCE
                ref_profile(kind, (int)n, z, v, options.scale_height, options.reference_column,
                            options.length_conversion, expected);
#endif
                check_profile(&grid, &p, expected, kind, 0);
                if (kind == 3 && options.reference_column != 1) {
                    double total = 0;
                    for (size_t i = 0; i < n; ++i)
                        total += p.columns[i];
                    same("O3 Dobson column", 0, total / 2.687e16, options.reference_column, 0);
                }
                tt_grid_destroy(&grid);
                tt_profile_destroy(&p);
                free(z);
                free(v);
                free(expected);
            }
    for (int year = 1950; year <= 2050; year += 10)
        for (int month = 1; month <= 12; ++month) {
            ++scenario;
            tt_grid time = {0};
            tt_profile z = {0}, d = {0};
            ok(tt_grid_equal(-1, 25, .25, &time));
            double *zr = calloc(6 * (time.cells + 1), sizeof(double)),
                   *dr = calloc(6 * (time.cells + 1), sizeof(double));
            double zone = (double)month - 6, lat = (double)(year - 1950) * 1.8 - 90,
                   lon = (double)month * 30 - 180;
            ok(tt_solar_profiles(&time, year, month, 12, zone, lat, lon, &z, &d));
#ifdef TT_HAVE_REFERENCE
            ref_solar_profiles((int)time.cells, time.edges, year, month, 12, zone, lat, lon, zr,
                               dr);
#endif
            check_profile(&time, &z, zr, 5, 90);
            check_profile(&time, &d, dr, 5, 0);
            free(zr);
            free(dr);
            tt_profile_destroy(&z);
            tt_profile_destroy(&d);
            tt_grid_destroy(&time);
        }
    host_checks();
    solar_boundaries(solar);
    {
        /* 29.4 - 13.4 rounds below 16: TUV-x appends 17 one-km points, not 16,
         * so the top edge is interpolated rather than zero. Values from the
         * original profile_o3_t constructor. */
        ++scenario;
        const double edges[] = {0, 10, 20, 29.4}, z[] = {0, 13.4}, v[] = {1.e12, 5.e11};
        const double plain[] = {1.e12, 626865671641.79102, 116041053960.47238, 14282750392.275186};
        const double rescaled[] = {6435590952944.716, 4034251045129.523, 746792757038.1855,
                                   91917939207.69377};
        tt_grid grid = {0};
        ok(tt_grid_from_edges(4, edges, &grid));
        for (int pass = 0; pass < 2; ++pass) {
            tt_profile p = {0};
            tt_profile_options options = tt_profile_options_default(TT_PROFILE_O3);
            if (!pass)
                options.reference_column = 1;
            ok(tt_profile_from_samples(&grid, 2, z, v, &options, &p));
            for (size_t i = 0; i < 4; ++i)
                same("O3 extension edge", i, p.edges[i], (pass ? rescaled : plain)[i], 0);
            tt_profile_destroy(&p);
        }
        tt_grid_destroy(&grid);
    }
    tt_grid bad = {0};
    if (tt_grid_equal(0, 1, 0, &bad) != TT_INVALID ||
        tt_grid_from_edges(0, NULL, &bad) != TT_INVALID)
        return 1;
    tt_grid_destroy(&bad);
    tt_profile_destroy(NULL);
    printf("Grids and profiles: %d cases, %zu comparisons; worst %.3g epsilon-scaled error\n",
           scenario, checked, worst);
    return 0;
}
