#include "turbotuv_model.h"
#include "fixtures.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_model_create(int, const char *, int *);
void ref_model_run(int, int, int, int, int, double, double, int, const double *, double *, double *,
                   double *, double *);
void ref_model_label(int, int, int, char *);
void ref_model_destroy(void);
#endif
#include "generated/model_cases.inc"
static size_t cases, runs, compared;
static double worst;
static const char *origin;
static const tt_dataset *resolve(void *context, const char *path) {
    (void)context;
    return tt_test_dataset(path);
}
static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "model %zu %s run %zu: %s\n", cases, origin, runs, message);
        exit(1);
    }
}
static void same(double a, double b, double scale, const char *field, size_t i) {
    double norm = fmax(fmax(fabs(a), fabs(b)), scale);
    double error = norm ? fabs(a - b) / norm / DBL_EPSILON : (a == b ? 0 : INFINITY);
    ++compared;
    if (!isfinite(a) || !isfinite(b) || error > 64) {
        fprintf(stderr, "model %zu %s run %zu %s[%zu]: %.17g != %.17g (%.3g eps)\n", cases, origin,
                runs, field, i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
int main(void) {
    const double angles[] = {0, 30, 60, 85, 89.9, 90, 90.01, 93, 99, 110, 180, 30};
    tt_data_source source = {resolve, NULL};
    for (size_t c = 0; c < sizeof model_cases / sizeof *model_cases; ++c) {
        ++cases;
        const struct model_case *test = model_cases + c;
        origin = test->origin;
        if (getenv("TT_VERBOSE"))
            fprintf(stderr, "Model %zu %s\n", cases, origin);
        tt_config *config = NULL;
        tt_error error;
        require(tt_config_parse(test->config, strlen(test->config), &config, &error) == TT_OK,
                error.message);
        tt_model *model = NULL;
        require(tt_model_create(tt_config_root(config), &source, &model, &error) == TT_OK,
                error.message);
        const char *solver = tt_config_string(tt_config_get(
            tt_config_get(tt_config_get(tt_config_root(config), "radiative transfer"), "solver"),
            "type"));
        int host_solver = !strcmp(solver, "from host");
        tt_config_destroy(config);
        if (!test->runnable) {
            /* This original reader fixture deliberately uses SSA/G above 1. */
            const tt_model_result *result = NULL;
            require(tt_model_run(model, 30, 1, &result, &error) == TT_INVALID && !result,
                    "unphysical compiled radiator rejected by solver");
            tt_model_destroy(model);
            continue;
        }
        tt_model_dimensions s = tt_model_shape(model);
        size_t nz = s.levels, nw = s.wavelengths, n = nz * nw,
               nr = nz * (s.photolysis + s.heating + s.dose);
        size_t counts[] = {s.photolysis, s.heating, s.dose};
#ifdef TT_HAVE_REFERENCE
        int dimensions[5];
        const char *ref_config = test->reference ? test->reference : test->config;
        ref_model_create((int)strlen(ref_config), ref_config, dimensions);
        require((size_t)dimensions[0] == nz && (size_t)dimensions[1] == nw &&
                    (size_t)dimensions[2] == s.photolysis && (size_t)dimensions[3] == s.heating &&
                    (size_t)dimensions[4] == s.dose,
                "reference dimensions");
#endif
        for (int k = 0; k < 3; ++k) {
            for (size_t j = 0; j < counts[k]; ++j) {
                const char *name = tt_model_rate_name(model, (tt_rate_kind)k, j);
                require(name && *name, "rate name");
#ifdef TT_HAVE_REFERENCE
                char label[256];
                ref_model_label(k, (int)j, sizeof label, label);
                if (strcmp(name, label))
                    fprintf(stderr, "kind %d index %zu: C '%s' Fortran '%s'\n", k, j, name, label);
                require(!strcmp(name, label), "reference rate order and name");
#endif
            }
            require(!tt_model_rate_name(model, (tt_rate_kind)k, counts[k]),
                    "out-of-bounds rate name");
        }
        double *host = calloc(6 * n, sizeof(double)),
               *expected = calloc(6 * n + nr + 1, sizeof(double));
        require(host && expected, "test allocation");
        double *rp = expected + 6 * n, *rh = rp + nz * s.photolysis, *rd = rh + nz * s.heating;
        const tt_profile *za = tt_model_profile(model, "solar zenith angle", "degrees"),
                         *ds = tt_model_profile(model, "Earth-Sun distance", "AU");
        size_t configured = za && ds && za->cells == ds->cells ? za->cells + 1 : 0;
        for (size_t a = 0; a < sizeof angles / sizeof *angles + configured; ++a) {
            ++runs;
            double angle = a < sizeof angles / sizeof *angles
                               ? angles[a]
                               : za->edges[a - sizeof angles / sizeof *angles];
            double scale = a < sizeof angles / sizeof *angles
                               ? (a % 2 ? 0.9833 : 1.0167)
                               : ds->edges[a - sizeof angles / sizeof *angles];
            int mode = host_solver && a ? (int)(a % 3) + 1 : 0;
            if (host_solver) {
                for (size_t f = 0; f < 6; ++f)
                    for (size_t w = 0; w < nw; ++w)
                        for (size_t z = 0; z < nz; ++z)
                            host[f * n + w * nz + z] =
                                (0.02 + 0.13 * f) * (1.0 + 0.1 * w) * (0.3 + 0.7 * z / (nz - 1));
                if (!mode)
                    memset(host, 0, 6 * n * sizeof(double));
                if (mode == 3)
                    for (size_t i = 0; i < 6 * n; ++i)
                        if (i % 3 == 0)
                            host[i] *= -2;
                tt_radiation r = {host,
                                  host + n,
                                  host + 2 * n,
                                  mode == 1 ? NULL : host + 3 * n,
                                  mode == 1 ? NULL : host + 4 * n,
                                  mode == 1 ? NULL : host + 5 * n};
                if (mode)
                    require(tt_model_set_radiation(model, &r) == TT_OK, "host radiation input");
            }
            const tt_model_result *result = NULL;
            require(tt_model_run(model, angle, scale, &result, &error) == TT_OK, error.message);
#ifdef TT_HAVE_REFERENCE
            ref_model_run((int)nz, (int)nw, (int)s.photolysis, (int)s.heating, (int)s.dose, angle,
                          scale, mode, host, expected, rp, rh, rd);
#endif
            const double *fields[] = {
                result->radiation.direct_flux,     result->radiation.down_flux,
                result->radiation.up_flux,         result->radiation.direct_irradiance,
                result->radiation.down_irradiance, result->radiation.up_irradiance};
            for (size_t f = 0; f < 6; ++f)
                for (size_t i = 0; i < n; ++i) {
#ifdef TT_HAVE_REFERENCE
                    same(fields[f][i], expected[f * n + i], 1, "radiation", f * n + i);
#else
                    require(isfinite(fields[f][i]), "finite radiation");
#endif
                    if (host_solver)
                        same(fields[f][i], mode == 1 && f >= 3 ? 0 : host[f * n + i] * scale, 1,
                             "host radiation", f * n + i);
                }
            const double *rates[] = {result->photolysis, result->heating, result->dose};
            const double *references[] = {rp, rh, rd};
            for (int k = 0; k < 3; ++k)
                for (size_t i = 0; i < counts[k] * nz; ++i) {
#ifdef TT_HAVE_REFERENCE
                    same(rates[k][i], references[k][i], 0,
                         k == 0   ? "photolysis"
                         : k == 1 ? "heating"
                                  : "dose",
                         i);
#else
                    require(isfinite(rates[k][i]), "finite rate");
#endif
                }
        }
        const tt_model_result *sentinel = NULL;
        require(tt_model_run(model, NAN, 1, &sentinel, &error) != TT_OK && !sentinel,
                "invalid zenith rejected");
        require(tt_model_run(model, 0, NAN, &sentinel, &error) != TT_OK && !sentinel,
                "invalid scale rejected");
        free(host);
        free(expected);
        tt_model_destroy(model);
#ifdef TT_HAVE_REFERENCE
        ref_model_destroy();
#endif
    }
    origin = "invalid configuration";
    for (size_t i = 0; i < sizeof invalid_model_cases / sizeof *invalid_model_cases; ++i) {
        tt_config *config = NULL;
        tt_model *model = NULL;
        tt_error error;
        const char *text = invalid_model_cases[i];
        require(tt_config_parse(text, strlen(text), &config, &error) == TT_OK,
                "negative test parse");
        tt_status status = tt_model_create(tt_config_root(config), &source, &model, &error);
        require(status != TT_OK && error.status == status && error.message[0] && !model,
                "invalid model rejected without changing output");
        tt_config_destroy(config);
    }
    origin = "configuration defaults";
    {
        /* TUV-x's default 550 nm depth of 0 normalizes the aerosol column to zero. */
        static const char *const aerosol =
            "{grids: [{name: height, type: equal interval, units: km, begins at: 0, ends at: 3,"
            " cell delta: 1}, {name: wavelength, type: from config file, units: nm, values: [300,"
            " 400, 500]}], profiles: [{name: temperature, type: from config file, units: K,"
            " grid: {name: height, units: km}, values: [290, 280, 270, 260]}, {name: air, type:"
            " from config file, units: molecule cm-3, grid: {name: height, units: km}, values:"
            " [2.5e19, 2.2e19, 2.0e19, 1.8e19]}, {name: surface albedo, type: from config file,"
            " units: none, grid: {name: wavelength, units: nm}, uniform value: 0.1}],"
            " radiative transfer: {solver: {type: delta eddington}, cross sections: [],"
            " radiators: [{name: aerosols, type: aerosol, optical depths: [0.2, 0.1, 0.05, 0.02],"
            " single scattering albedo: 0.9, asymmetry factor: 0.6%s}]}, O2 absorption: {%s}}";
        const char *variants[][2] = {{"", ""},
                                     {", 550 nm optical depth: 0.235", ""},
                                     {"", "O2 estimate: {scale factor: 0}"}};
        char text[2048];
        for (int v = 0; v < 3; ++v) {
            snprintf(text, sizeof text, aerosol, variants[v][0], variants[v][1]);
            tt_config *config = NULL;
            tt_model *model = NULL;
            tt_error error;
            require(tt_config_parse(text, strlen(text), &config, &error) == TT_OK, error.message);
            if (v == 2) {
                /* No O2 bands on this grid: the estimate is unused but still checked
                 * where bands exist; use the combined grid to reach the check. */
                tt_config_destroy(config);
                const char *bands = strstr(text, "{name: wavelength");
                char modified[2048];
                snprintf(modified, sizeof modified,
                         "%.*s{name: wavelength, type: from csv file,"
                         " units: nm, file path: data/grids/wavelength/combined.grid}%s",
                         (int)(bands - text), text, strstr(bands, "}]") + 1);
                require(tt_config_parse(modified, strlen(modified), &config, &error) == TT_OK,
                        error.message);
                require(tt_model_create(tt_config_root(config), &source, &model, &error) ==
                                TT_RANGE &&
                            !model,
                        "zero O2 scale factor rejected");
                tt_config_destroy(config);
                continue;
            }
            require(tt_model_create(tt_config_root(config), &source, &model, &error) == TT_OK,
                    error.message);
            const tt_radiator *r = tt_model_radiator(model, "aerosols");
            tt_model_dimensions shape = tt_model_shape(model);
            double total = 0;
            for (size_t i = 0; i < (shape.levels - 1) * shape.wavelengths; ++i)
                total += r->optical_depth[i];
            require(v == 0 ? total == 0 : total > 0, "aerosol 550 nm default");
            tt_model_destroy(model);
            tt_config_destroy(config);
        }
    }
    origin = "configuration keys";
    {
        /* As TUV-x's config%validate: unknown keys fail, "__" keys are comments. */
        const char *example = model_cases[0].config; /* examples/tuv_5_4.json */
        static const char *const edits[][2] = {
            {"\"quantum yield\"", "\"scaling factr\":2.5,\"quantum yield\""},
            {"\"treat as air\"", "\"treat as aire\""},
            {"\"type\":\"O3\",", "\"type\":\"O3\",\"scale heigth\":4.5,"},
            {"\"550 nm optical depth\"", "\"550nm optical depth\""},
            {"\"cell delta\"", "\"cell size\""},
            {"\"photolysis\"", "\"__photolysis\""}};
        require(!strcmp(model_cases[0].origin, "examples/tuv_5_4.json"), "example case");
        for (size_t i = 0; i < sizeof edits / sizeof *edits; ++i) {
            const char *at = strstr(example, edits[i][0]);
            require(at != NULL, "key to edit");
            size_t before = (size_t)(at - example), old = strlen(edits[i][0]);
            size_t length = strlen(example) - old + strlen(edits[i][1]);
            char *text = malloc(length + 1);
            require(text != NULL, "allocation");
            memcpy(text, example, before);
            strcpy(text + before, edits[i][1]);
            strcat(text, at + old);
            tt_config *config = NULL;
            tt_model *model = NULL;
            tt_error error;
            require(tt_config_parse(text, length, &config, &error) == TT_OK, error.message);
            tt_status status = tt_model_create(tt_config_root(config), &source, &model, &error);
            /* Renaming the photolysis section to a comment removes the rates. */
            require(i + 1 < sizeof edits / sizeof *edits
                        ? status == TT_INVALID && !model && strstr(error.message, "Unknown")
                        : status == TT_OK && tt_model_shape(model).photolysis == 0,
                    "configuration key validation");
            tt_model_destroy(model);
            tt_config_destroy(config);
            free(text);
        }
    }
    printf("Whole model: %zu configurations, %zu runs, %zu comparisons, worst %.3g eps\n", cases,
           runs, compared, worst);
    return 0;
}
