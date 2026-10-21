#include "turbotuv_model.h"
#include "fixtures.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_profile_config(int, const char *, int, const char *, int, int *, double *);
#endif
#include "generated/profile_cases.inc"
static size_t cases, compared;
static double worst;
static const tt_dataset *resolve(void *context, const char *path) {
    (void)context;
    return tt_test_dataset(path);
}
static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "profile config %zu: %s\n", cases, message);
        exit(1);
    }
}
static tt_config *parse(const char *text) {
    tt_error error;
    tt_config *c = NULL;
    require(tt_config_parse(text, strlen(text), &c, &error) == TT_OK, error.message);
    return c;
}
static void same(double a, double b, double scale, size_t f, size_t i) {
    double norm = fmax(fmax(fabs(a), fabs(b)), scale);
    double error = norm ? fabs(a - b) / norm / DBL_EPSILON : (a == b ? 0 : INFINITY);
    ++compared;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "profile config %zu field %zu[%zu]: %.17g != %.17g (%.3g eps)\n", cases, f,
                i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static void invalid_profiles(void) {
    const char *const text[] = {
        "{type: unknown}",
        "{type: from csv file, file path: missing, grid: {name: height, units: km}}",
        "{type: from csv file, interpolator: conserving, file path: "
        "data/profiles/atmosphere/ussa.temp, grid: {name: height, units: km}}",
        "{type: from config file, values: [1], grid: {name: height, units: km}}",
        "{type: from config file, grid: {name: height, units: km}}",
        "{type: from config file, uniform value: 1, grid: {name: absent, units: km}}",
        "{type: Earth-Sun distance, year: 2000, month: 2, day: 30}",
        "{type: extraterrestrial flux, file path: []}",
        "{type: extraterrestrial flux, file path: [data/profiles/solar/susim_hi.flx], "
        "interpolator: [bad]}"};
    tt_grid g = {0};
    require(tt_grid_equal(0, 10, 1, &g) == TT_OK, "negative test grid");
    tt_named_grid grids[] = {{"height", "km", &g}, {"time", "hours", &g}, {"wavelength", "nm", &g}};
    for (size_t i = 0; i < sizeof text / sizeof *text; ++i) {
        tt_config *c = parse(text[i]);
        tt_profile p = {0};
        tt_error error;
        tt_status status = tt_profile_from_config(tt_config_root(c), 3, grids, NULL, &p, &error);
        require(status != TT_OK && status == error.status && error.message[0],
                "invalid profile rejected");
        require(!p.edges && !p.midpoints, "failed profile preserves output");
        tt_config_destroy(c);
    }
    tt_grid_destroy(&g);
}
int main(void) {
    tt_data_source source = {resolve, NULL};
    for (size_t c = 0; c < sizeof profile_cases / sizeof *profile_cases; ++c) {
        ++cases;
        const struct profile_case *test = profile_cases + c;
        tt_config *gc = parse(test->grids), *pc = parse(test->config);
        tt_config_node gn = tt_config_root(gc), pn = tt_config_root(pc);
        size_t ng = tt_config_size(gn);
        tt_grid *grid = calloc(ng, sizeof *grid);
        tt_named_grid *named = calloc(ng, sizeof *named);
        tt_error error;
        require(grid && named, "allocation");
        tt_status grid_status = TT_OK;
        for (size_t j = 0; j < ng; ++j) {
            tt_config_node node = tt_config_at(gn, j);
            tt_status gs = tt_grid_from_config(node, &source, grid + j, &error);
            if (gs != TT_OK)
                grid_status = gs;
            named[j] = (tt_named_grid){tt_config_string(tt_config_get(node, "name")),
                                       tt_config_string(tt_config_get(node, "units")), grid + j};
        }
        if (test->valid == -1) {
            require(grid_status == TT_INVALID, "invalid upstream grid type rejected");
            for (size_t j = 0; j < ng; ++j)
                tt_grid_destroy(grid + j);
            free(grid);
            free(named);
            tt_config_destroy(gc);
            tt_config_destroy(pc);
            continue;
        }
        require(grid_status == TT_OK, error.message);
        tt_profile p = {0};
        tt_status status = tt_profile_from_config(pn, ng, named, &source, &p, &error);
        if (!test->valid) {
            /* Dormant upstream inputs omit required coordinates or a grid. */
            require(status == TT_INVALID, "missing required profile settings rejected");
            for (size_t j = 0; j < ng; ++j)
                tt_grid_destroy(grid + j);
            free(grid);
            free(named);
            tt_config_destroy(gc);
            tt_config_destroy(pc);
            continue;
        }
        if (status)
            fprintf(stderr, "%s\n%s\n", test->origin, test->config);
        require(status == TT_OK, error.message);
        size_t n = p.cells;
        double *expected = calloc(6 * (n + 1), sizeof(double));
        int fields[6] = {1, 1, 1, 1, 1, 1};
        require(expected != NULL, "reference allocation");
        const double *v[] = {p.edges,   p.midpoints,        p.differences,
                             p.columns, p.extended_columns, p.burden};
#ifdef TT_HAVE_REFERENCE
        if (getenv("TT_VERBOSE"))
            fprintf(stderr, "profile %zu %s %s\n", cases, test->origin, test->config);
        ref_profile_config((int)strlen(test->grids), test->grids, (int)strlen(test->config),
                           test->config, (int)n, fields, expected);
#else
        for (size_t f = 0; f < 6; ++f)
            memcpy(expected + f * (n + 1), v[f], (n + (f == 0 || f == 4)) * sizeof(double));
#endif
        const char *type = tt_config_string(tt_config_get(pn, "type"));
        for (size_t f = 0; f < 6; ++f)
            if (fields[f])
                for (size_t i = 0; i < n + (f == 0 || f == 4); ++i)
                    same(v[f][i], expected[f * (n + 1) + i],
                         !strcmp(type, "solar zenith angle") ? 1 : 0, f, i);
        for (size_t j = 0; j < ng; ++j)
            tt_grid_destroy(grid + j);
        free(grid);
        free(named);
        free(expected);
        tt_config_destroy(gc);
        tt_config_destroy(pc);
        tt_profile_destroy(&p);
    }
    invalid_profiles();
    printf("Profile configuration: %zu cases, %zu comparisons, worst %.3g relative eps\n", cases,
           compared, worst);
    return 0;
}
