#include "turbotuv_model.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_model_host_create(int, const char *, int, int, const double *, const double *,
                           const double *, const double *, const double *, double);
void ref_model_host_update(int, int, int, const double *, const double *, const double *,
                           const double *, const double *, const double *, const double *, double,
                           const double *, const double *, const double *);
void ref_model_run(int, int, int, int, int, double, double, int, const double *, double *, double *,
                   double *, double *);
void ref_model_destroy(void);
#endif
#include "generated/host_model_cases.inc"
static size_t scenario, compared;
static double worst;
static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "host model %zu: %s\n", scenario, message);
        exit(1);
    }
}
static void same(double a, double b, double scale, size_t i) {
    double norm = fmax(fmax(fabs(a), fabs(b)), scale);
    double error = norm ? fabs(a - b) / norm / DBL_EPSILON : (a == b ? 0 : INFINITY);
    ++compared;
    if (!isfinite(a) || !isfinite(b) || error > 64) {
        fprintf(stderr, "host model %zu[%zu]: %.17g != %.17g (%.3g eps)\n", scenario, i, a, b,
                error);
        exit(1);
    }
    worst = fmax(worst, error);
}
int main(void) {
    const double angles[] = {0, 30, 60, 85, 89.9, 90, 90.01, 93, 99, 110, 180};
    enum { NZ = 17, NL = NZ - 1 };
    for (size_t c = 0; c < sizeof host_model_cases / sizeof *host_model_cases; ++c) {
        tt_config *config = NULL;
        tt_error error;
        require(tt_config_parse(host_model_cases[c], strlen(host_model_cases[c]), &config,
                                &error) == TT_OK,
                error.message);
        tt_grid height = {0}, wave = {0};
        require(tt_grid_equal(0, 120, 120.0 / NL, &height) == TT_OK, "height grid");
        require(
            tt_grid_from_dataset(tt_dataset_find(c == 1 ? "data/grids/wavelength/cam.csv"
                                                        : "data/grids/wavelength/combined.grid"),
                                 &wave) == TT_OK,
            "wavelength grid");
        size_t nw = wave.cells, n = NZ * nw, layer = NL * nw;
        double z[NZ], mid[NL], air[NZ], am[NL], t[NZ], tm[NL], col[NL], top = 8.01;
        for (size_t i = 0; i < NZ; ++i) {
            z[i] = height.edges[i];
            air[i] = 2.5e19 * exp(-z[i] / 8);
            t[i] = 220 + 60 * exp(-z[i] / 10);
            if (i < NL)
                col[i] = 1e5 * height.widths[i] * 2.5e19 * exp(-height.midpoints[i] / 8);
        }
        tt_profile ap = {0}, tp = {0};
        require(tt_profile_create(NL, &ap) == TT_OK && tt_profile_create(NL, &tp) == TT_OK,
                "host profiles");
        require(tt_profile_update(&ap, air, NULL, col, &top, NULL) == TT_OK, "initial air");
        require(tt_profile_update(&tp, t, NULL, NULL, NULL, NULL) == TT_OK, "initial temperature");
        double *rad = calloc(3 * layer, sizeof(double));
        require(rad != NULL, "host radiator allocation");
        tt_named_grid grids[] = {{"height", "km", &height}, {"wavelength", "nm", &wave}};
        tt_named_profile profiles[] = {{"air", "molecule cm-3", &ap}, {"temperature", "K", &tp}};
        tt_named_radiator cloud = {"host cloud", NL, nw, {rad, rad + layer, rad + 2 * layer, 0}};
        tt_model_inputs inputs = {2, grids, 2, profiles, 1, &cloud};
        tt_model *model = NULL;
        require(tt_model_create_with_inputs(tt_config_root(config), NULL, &inputs, &model,
                                            &error) == TT_OK,
                error.message);
        /* Conflicts and invalid input arrays must not damage a live model. */
        tt_model *bad = NULL;
        tt_model_inputs invalid = inputs;
        invalid.grids = NULL;
        require(tt_model_create_with_inputs(tt_config_root(config), NULL, &invalid, &bad, &error) ==
                        TT_INVALID &&
                    !bad,
                "missing host array");
        tt_named_grid duplicates[] = {grids[0], grids[0]};
        invalid = inputs;
        invalid.grids = duplicates;
        require(tt_model_create_with_inputs(tt_config_root(config), NULL, &invalid, &bad, &error) ==
                        TT_INVALID &&
                    !bad,
                "duplicate host name");
        tt_model_dimensions s = tt_model_shape(model);
        size_t rates = NZ * (s.photolysis + s.heating + s.dose);
        double *ref = calloc(6 * n + rates + 1, sizeof(double));
        double *generic = calloc(8 * n + rates + 1, sizeof(double));
        require(ref != NULL && generic != NULL, "reference allocation");
#ifdef TT_HAVE_REFERENCE
        ref_model_host_create((int)strlen(host_model_cases[c]), host_model_cases[c], NZ, (int)nw, z,
                              wave.edges, air, t, col, top);
#endif
        tt_config_destroy(config);
        tt_grid_destroy(&height);
        tt_grid_destroy(&wave);
        tt_profile_destroy(&ap);
        tt_profile_destroy(&tp);
        const tt_grid *borrowed_grid = tt_model_grid(model, "height", "km");
        const tt_profile *borrowed_profile = tt_model_profile(model, "air", "molecule cm-3");
        const tt_radiator *borrowed_radiator = tt_model_radiator(model, "host cloud");
        for (size_t step = 0; step < 12; ++step) {
            int mode = (int)(step % 3);
            top = mode == 0 ? 6.5 + 0.1 * step : 3e14 * (1 + step);
            for (size_t j = 0; j < NZ; ++j) {
                z[j] = 0.5 + 0.03 * step + (100 + step) * pow((double)j / NL, 1.1);
                air[j] = (2.4e19 + step * 1e17) * exp(-z[j] / 8.5);
                t[j] = 215 + 45 * exp(-z[j] / 10) + step;
            }
            for (size_t j = 0; j < NL; ++j) {
                mid[j] = 0.4 * z[j] + 0.6 * z[j + 1];
                am[j] = 0.45 * air[j] + 0.55 * air[j + 1];
                tm[j] = 0.47 * t[j] + 0.53 * t[j + 1];
                col[j] = 1e5 * (z[j + 1] - z[j]) * sqrt(air[j]) * sqrt(air[j + 1]);
            }
            for (size_t w = 0; w < nw; ++w)
                for (size_t j = 0; j < NL; ++j) {
                    size_t i = w * NL + j;
                    rad[i] = 0.03 * (1 + step) * (1 + 0.2 * cos(w)) * exp(-z[j] / 5);
                    rad[layer + i] = 0.5 + 0.04 * step;
                    rad[2 * layer + i] = 0.2 + 0.01 * step;
                }
            require(tt_model_update_grid(model, "height", "km", NL, z, mode == 1 ? mid : NULL) ==
                        TT_OK,
                    "height update");
            require(tt_model_update_profile(model, "air", "molecule cm-3", NL, air,
                                            mode == 1 ? am : NULL, mode == 2 ? NULL : col,
                                            mode == 0 ? &top : NULL,
                                            mode == 1 ? &top : NULL) == TT_OK,
                    "air update");
            require(tt_model_update_profile(model, "temperature", "K", NL, t, mode == 1 ? tm : NULL,
                                            NULL, NULL, NULL) == TT_OK,
                    "temperature update");
            cloud.radiator.single_scattering_albedo = mode == 1 ? rad + layer : NULL;
            cloud.radiator.asymmetry = mode == 0 ? NULL : rad + 2 * layer;
            require(tt_model_set_radiator(model, &cloud) == TT_OK, "cloud update");
            require(borrowed_grid == tt_model_grid(model, "height", "km") &&
                        borrowed_profile == tt_model_profile(model, "air", "molecule cm-3") &&
                        borrowed_radiator == tt_model_radiator(model, "host cloud"),
                    "stable borrowed views");
#ifdef TT_HAVE_REFERENCE
            ref_model_host_update(NZ, (int)nw, mode, z, mid, air, am, t, tm, col, top, rad,
                                  rad + layer, rad + 2 * layer);
#endif
            for (size_t a = 0; a < sizeof angles / sizeof *angles; ++a) {
                ++scenario;
                const tt_model_result *result = NULL;
                /* The generic sequence must give the prepared results exactly. */
                require(tt_model_set_execution(model, TT_EXECUTION_GENERIC) == TT_OK, "generic");
                require(tt_model_run(model, angles[a], 1.003, &result, &error) == TT_OK,
                        error.message);
                size_t outputs = 8 * n + rates;
                memcpy(generic, result->radiation.direct_flux, 6 * n * sizeof(double));
                memcpy(generic + 6 * n, result->actinic_flux, n * sizeof(double));
                memcpy(generic + 7 * n, result->spectral_irradiance, n * sizeof(double));
                memcpy(generic + 8 * n, result->photolysis, rates * sizeof(double));
                require(tt_model_set_execution(model, TT_EXECUTION_PREPARED) == TT_OK, "prepared");
                require(tt_model_run(model, angles[a], 1.003, &result, &error) == TT_OK,
                        error.message);
                require(
                    !memcmp(generic, result->radiation.direct_flux, 6 * n * sizeof(double)) &&
                        !memcmp(generic + 6 * n, result->actinic_flux, n * sizeof(double)) &&
                        !memcmp(generic + 7 * n, result->spectral_irradiance, n * sizeof(double)) &&
                        !memcmp(generic + 8 * n, result->photolysis, rates * sizeof(double)),
                    "prepared and generic runs differ");
                (void)outputs;
                double *rp = ref + 6 * n, *rh = rp + NZ * s.photolysis, *rd = rh + NZ * s.heating;
#ifdef TT_HAVE_REFERENCE
                ref_model_run(NZ, (int)nw, (int)s.photolysis, (int)s.heating, (int)s.dose,
                              angles[a], 1.003, 0, ref, ref, rp, rh, rd);
#endif
                const double *fields[] = {
                    result->radiation.direct_flux,     result->radiation.down_flux,
                    result->radiation.up_flux,         result->radiation.direct_irradiance,
                    result->radiation.down_irradiance, result->radiation.up_irradiance};
                for (size_t f = 0; f < 6; ++f)
                    for (size_t i = 0; i < n; ++i) {
#ifdef TT_HAVE_REFERENCE
                        same(fields[f][i], ref[f * n + i], 1, f * n + i);
#else
                        require(isfinite(fields[f][i]), "finite radiation");
#endif
                    }
                const double *values[] = {result->photolysis, result->heating, result->dose};
                const double *expected[] = {rp, rh, rd};
                size_t counts[] = {s.photolysis, s.heating, s.dose};
                for (size_t k = 0; k < 3; ++k)
                    for (size_t i = 0; i < NZ * counts[k]; ++i) {
#ifdef TT_HAVE_REFERENCE
                        same(values[k][i], expected[k][i], 0, i);
#else
                        require(isfinite(values[k][i]), "finite rate");
#endif
                    }
            }
        }
        double old = borrowed_grid->edges[0];
        require(tt_model_update_grid(model, "height", "km", NL - 1, z, NULL) == TT_INVALID &&
                    borrowed_grid->edges[0] == old,
                "shape rejection preserves state");
        require(tt_model_update_grid(model, "wavelength", "nm", nw, z, NULL) == TT_INVALID,
                "fixed wavelength grid");
        require(tt_model_update_profile(model, "missing", "K", NL, t, NULL, NULL, NULL, NULL) ==
                    TT_NOT_FOUND,
                "missing profile");
        cloud.layers = NL - 1;
        require(tt_model_set_radiator(model, &cloud) == TT_INVALID, "radiator shape rejection");
        cloud.layers = NL;
        rad[0] = NAN;
        require(tt_model_set_radiator(model, &cloud) == TT_INVALID, "nonfinite radiator rejection");
        tt_model_destroy(model);
#ifdef TT_HAVE_REFERENCE
        ref_model_destroy();
#endif
        free(rad);
        free(ref);
        free(generic);
    }
    printf("Host model: %zu runs, %zu comparisons, worst %.3g eps\n", scenario, compared, worst);
    return 0;
}
