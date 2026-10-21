#include "turbotuv_model.h"
#include "fixtures.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_spectrum(int, int, const char *, int, int, int, const double *, const double *,
                  const double *, const double *, double *);
#endif
#include "generated/spectrum_cases.inc"
static size_t cases, comparisons;
static double worst;
static const tt_dataset *resolve(void *context, const char *path) {
    (void)context;
    return tt_test_dataset(path);
}
static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "spectrum case %zu: %s\n", cases, message);
        exit(1);
    }
}
static void compare(double a, double b, size_t i) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(fabs(a), fabs(b)) + DBL_MIN);
    ++comparisons;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "spectrum %zu [%zu] %.17g != %.17g (%.4g eps)\n", cases, i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static tt_config *parse(const char *text) {
    tt_config *c = NULL;
    tt_error error;
    require(tt_config_parse(text, strlen(text), &c, &error) == TT_OK, error.message);
    return c;
}
static void run(int role, const char *grid_json, const char *json, const char *origin) {
    ++cases;
    tt_data_source source = {resolve, NULL};
    tt_config *config = parse(json), *grid_config = parse(grid_json);
    tt_grid grid = {0};
    tt_error error;
    require(tt_grid_from_config(tt_config_root(grid_config), &source, &grid, &error) == TT_OK,
            error.message);
    tt_spectrum *s = NULL;
    tt_status status = tt_spectrum_create(role, tt_config_root(config), &grid, &source, &s, &error);
    if (status)
        fprintf(stderr, "%s\n%s\n", origin, json);
    require(status == TT_OK, error.message);
    enum { NZ = 23 };
    double temp[NZ] = {150, 173, 180, 194, 200, 203, 210, 218, 233, 235, 240, 250,
                       273, 295, 296, 298, 300, 310, 320, 350, 400, 450, 500};
    if (strstr(origin, "test/data/cross_sections/"))
        for (size_t i = 0; i < NZ; ++i)
            temp[i] = 298.15 + .001 * ((double)i - 11);
    double air[NZ];
    for (size_t i = 0; i < NZ; ++i)
        air[i] = 2.5e19 * exp(-.65 * (double)i);
    size_t n = grid.cells * (role == TT_SPECTRUM_WEIGHT ? 1 : NZ);
    double *out = malloc(n * sizeof(double)), *reference = malloc(n * sizeof(double));
    require(out && reference, "allocation");
    const char *type = tt_config_string(tt_config_get(tt_config_root(config), "type"));
    int repeats = role == TT_SPECTRUM_CROSS_SECTION && strcmp(type, "HNO3+hv->OH+NO2") &&
                          strcmp(type, "H2O2+hv->OH+OH")
                      ? 2
                      : 1;
    for (int midpoint = 0; midpoint < repeats; ++midpoint) {
        status = tt_spectrum_calculate(s, NZ, temp, air, out);
        if (status)
            fprintf(stderr, "%s\n%s\n", origin, json);
        require(status == TT_OK, "calculate");
#ifdef TT_HAVE_REFERENCE
        if (getenv("TT_VERBOSE"))
            fprintf(stderr, "case %zu mid %d %s %s\n", cases, midpoint, origin, type);
        ref_spectrum(role, (int)strlen(json), json, (int)grid.cells, NZ, midpoint, grid.edges,
                     grid.midpoints, temp, air, reference);
#else
        require(tt_spectrum_calculate(s, NZ, temp, air, reference) == TT_OK, "repeat");
#endif
        for (size_t i = 0; i < n; ++i)
            compare(out[i], reference[i], i);
    }
    /* The object owns its configuration and wavelength storage. */
    tt_config_destroy(config);
    tt_config_destroy(grid_config);
    tt_grid_destroy(&grid);
    require(tt_spectrum_calculate(s, NZ, temp, air, reference) == TT_OK, "owned inputs");
    for (size_t i = 0; i < n; ++i)
        compare(out[i], reference[i], i);
    tt_spectrum_destroy(s);
    free(out);
    free(reference);
}
static void invalid_inputs(void) {
    static const char *const configs[] = {
        "{}", "{type: unknown}", "{type: base}",
        "{type: base, netcdf files: [{file path: missing.nc}]}",
        "{type: base, netcdf files: [{file path: data/cross_sections/O2_1.nc, lower extrapolation: "
        "{type: unknown}}]}",
        "{type: base, data: {point values: [{wavelength: 42, value: 1}]}}",
        "{type: base, data: {}, override bands: [{band: unknown, value: 1}]}",
        "{type: base, data: {}, override bands: [{band: range, minimum wavelength: 1000, value: "
        "1}]}",
        "{type: temperature based, parameterization: {type: unknown}}",
        "{type: temperature based, parameterization: {AA: [1], BB: [2,3], lp: [0], base "
        "temperature: 298, base wavelength: 0, logarithm: natural}}",
        /* TUV-x accepts merged files and point data only for base cross sections. */
        "{type: tint, merge data: true, netcdf files: [{file path: "
        "test/data/cross_sections/cross_section.tint.nc}, {file path: "
        "test/data/cross_sections/cross_section.tint.nc}]}",
        "{type: CH3ONO2+hv->CH3O+NO2, netcdf files: [{file path: "
        "data/cross_sections/CH3ONO2_1.nc}], data: {point values: []}}",
        /* A zero Burkholder slope fails every evaluation, so it fails here. */
        "{type: temperature based, parameterization: {type: BURKHOLDER, netcdf file: {file path: "
        "data/cross_sections/HO2NO2_temp_JPL06.nc}, A: 1, B: 0}}"};
    tt_grid grid = {0};
    tt_data_source source = {resolve, NULL};
    require(tt_grid_equal(200, 400, 10, &grid) == TT_OK, "invalid test grid");
    for (size_t i = 0; i < sizeof configs / sizeof *configs; ++i) {
        tt_config *config = parse(configs[i]);
        tt_spectrum *out = NULL;
        tt_error error;
        tt_status status = tt_spectrum_create(TT_SPECTRUM_CROSS_SECTION, tt_config_root(config),
                                              &grid, &source, &out, &error);
        require(status != TT_OK && status == error.status && error.message[0] && error.line,
                "structured error");
        require(out == NULL, "failed constructor preserves output");
        if (i == 3)
            require(status == TT_NOT_FOUND, "embedded data only");
        tt_config_destroy(config);
    }
    tt_config *config = parse("{type: base}");
    tt_spectrum *out = NULL;
    tt_error error;
    require(tt_spectrum_create(TT_SPECTRUM_WEIGHT, tt_config_root(config), &grid, NULL, &out,
                               &error) == TT_INVALID,
            "weight needs data");
    tt_config_destroy(config);
    config = parse("{type: Taylor series, coefficients: []}");
    require(tt_spectrum_create(TT_SPECTRUM_QUANTUM_YIELD, tt_config_root(config), &grid, NULL, &out,
                               &error) == TT_INVALID,
            "Taylor needs coefficients");
    tt_config_destroy(config);
    require(tt_spectrum_table_count(NULL) == 0 && !tt_spectrum_table(NULL, 0), "null views");
    require(tt_spectrum_calculate(NULL, 0, NULL, NULL, NULL) == TT_INVALID, "null evaluator");
    tt_spectrum_destroy(NULL);
    tt_grid_destroy(&grid);
}
/* Datasets supplied by a resolver are validated like compiled ones. */
static const double malformed_wavelength[] = {120, 140, 160, 180};
static const double malformed_values[] = {1.e-20, 2.e-20, 3.e-20, 4.e-20, 5.e-20, 6.e-20,
                                          7.e-20, 8.e-20, 9.e-20, 1.e-19, 2.e-19, 3.e-19};
static const double nan_temperature[] = {200, NAN, 300}, nan_first[] = {NAN};
static const tt_data_variable short_variables[] = {
    {"wavelength", 1, {4}, 4, malformed_wavelength, NULL},
    /* shape {3, 4} needs 12 values; the variable claims 4 */
    {"temperature_parameters", 2, {3, 4}, 4, malformed_values, NULL}};
static const tt_data_variable nan_variables[] = {
    {"wavelength", 1, {4}, 4, malformed_wavelength, NULL},
    {"cross_section_parameters", 2, {3, 4}, 12, malformed_values, NULL},
    {"quantum_yield_parameters", 2, {3, 4}, 12, malformed_values, NULL},
    {"temperature", 1, {3}, 3, nan_temperature, NULL}};
/* OClO reads one temperature per file. */
static const tt_data_variable nan_first_variables[] = {
    {"wavelength", 1, {4}, 4, malformed_wavelength, NULL},
    {"cross_section_parameters", 2, {1, 4}, 4, malformed_values, NULL},
    {"temperature", 1, {1}, 1, nan_first, NULL}};
static const tt_dataset malformed_sets[] = {
    {"short.nc", NULL, 2, short_variables, NULL, 0},
    {"nan.nc", NULL, 4, nan_variables, NULL, 0},
    {"nan_first.nc", NULL, 3, nan_first_variables, NULL, 0}};
static const tt_dataset *malformed(void *context, const char *path) {
    (void)context;
    for (size_t i = 0; i < sizeof malformed_sets / sizeof *malformed_sets; ++i)
        if (!strcmp(path, malformed_sets[i].path))
            return malformed_sets + i;
    return NULL;
}
static void malformed_inputs(void) {
    static const struct {
        int role;
        const char *config;
    } cases[] = {{TT_SPECTRUM_CROSS_SECTION,
                  "{type: temperature based, parameterization: {type: "
                  "TAYLOR_SERIES, netcdf file: {file path: short.nc}, base "
                  "temperature: 298}}"},
                 {TT_SPECTRUM_CROSS_SECTION, "{type: tint, netcdf files: [{file path: nan.nc}]}"},
                 {TT_SPECTRUM_CROSS_SECTION, "{type: OClO+hv->Products, netcdf files: [{file path: "
                                             "nan_first.nc}]}"},
                 {TT_SPECTRUM_QUANTUM_YIELD, "{type: tint, netcdf files: [nan.nc]}"}};
    tt_data_source source = {malformed, NULL};
    tt_grid grid = {0};
    require(tt_grid_equal(130, 170, 10, &grid) == TT_OK, "malformed test grid");
    for (size_t i = 0; i < sizeof cases / sizeof *cases; ++i) {
        tt_config *config = parse(cases[i].config);
        tt_spectrum *out = NULL;
        tt_error error;
        tt_status status =
            tt_spectrum_create(cases[i].role, tt_config_root(config), &grid, &source, &out, &error);
        if (getenv("TT_VERBOSE"))
            fprintf(stderr, "malformed %zu: %s %s\n", i, tt_status_string(status), error.message);
        require(status == TT_INVALID && !out, "malformed dataset rejected");
        tt_config_destroy(config);
    }
    tt_grid_destroy(&grid);
}
int main(void) {
    for (size_t i = 0; i < sizeof spectrum_cases / sizeof *spectrum_cases; ++i) {
        const struct spectrum_case *c = spectrum_cases + i;
        run(c->role, c->grid, c->config, c->origin);
    }
    invalid_inputs();
    malformed_inputs();
    printf("Prepared spectra: %zu cases, %zu comparisons, worst %.3g relative eps\n", cases,
           comparisons, worst);
    return 0;
}
