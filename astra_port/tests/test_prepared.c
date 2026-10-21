/* Prepared runs must reproduce the generic runs bit for bit: every spectrum
 * plan against its generic kernel, and every model configuration on varied
 * columns, including twilight, night, and inputs the plans decline. */
#include "model_plan.h"
#include "spectrum_plan.h"
#include "fixtures.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "generated/model_cases.inc"
#include "generated/spectrum_cases.inc"

static size_t compared, spectra, models, columns;
static const char *origin;
static const tt_dataset *resolve(void *context, const char *path) {
    (void)context;
    return tt_test_dataset(path);
}
static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "prepared (%s): %s\n", origin, message);
        exit(1);
    }
}
static uint64_t state = 2026;
static double uniform(double lo, double hi) {
    state = state * 6364136223846793005u + 1442695040888963407u;
    return lo + (hi - lo) * ((double)(state >> 11) / 9007199254740992.0);
}
static void same(const double *a, const double *b, size_t n, const char *what) {
    compared += n;
    if (n && memcmp(a, b, n * sizeof(double))) {
        for (size_t i = 0; i < n; ++i)
            if (memcmp(a + i, b + i, sizeof(double))) {
                fprintf(stderr, "prepared (%s): %s[%zu] %.17g != %.17g\n", origin, what, i, a[i],
                        b[i]);
                break;
            }
        exit(1);
    }
}
static tt_config *parse(const char *text) {
    tt_config *config = NULL;
    tt_error error;
    require(tt_config_parse(text, strlen(text), &config, &error) == TT_OK, error.message);
    return config;
}

/* Temperatures at the kernels' clamp and branch boundaries, and between. */
static const double edges[] = {100,
                               173,
                               190,
                               194,
                               199.999999999999,
                               200,
                               203,
                               209.999999999999,
                               210,
                               218,
                               219.99,
                               220,
                               233,
                               234.999999999999,
                               235,
                               240,
                               250,
                               260,
                               273,
                               279.99,
                               280,
                               295,
                               295.00000000001,
                               296,
                               298,
                               300,
                               300.00000000001,
                               305,
                               320,
                               350,
                               350.00000000001,
                               400,
                               500,
                               1200,
                               9999};
static void column(size_t nz, int pattern, double *t, double *m) {
    for (size_t z = 0; z < nz; ++z) {
        t[z] = pattern % 3 == 0 ? edges[(size_t)uniform(0, sizeof edges / sizeof *edges)]
                                : uniform(150, 350) + (pattern % 3 == 2 ? uniform(0, 2000) : 0);
        m[z] = pattern % 4 == 0 ? (uniform(0, 1) < .2 ? 0 : 8.e17) : uniform(1e10, 3e19);
    }
}

/* Kinds the generated cases do not include. */
static const struct spectrum_case extra_cases[] = {
    {0, "O3 JPL06",
     "{name: wavelength, type: from csv file, units: nm, file path: "
     "data/grids/wavelength/combined.grid}",
     "{type: O3 JPL06, base netcdf file: {file path: data/cross_sections/O3_JPL06_base.nc},"
     " low-temperature netcdf file: {file path: data/cross_sections/O3_JPL06_218K.nc},"
     " high-temperature netcdf file: {file path: data/cross_sections/O3_JPL06_298K.nc}}"}};

/* Every spectrum against its generic kernel at levels of varied state. */
static void spectrum_plans(void) {
    tt_data_source source = {resolve, NULL};
    enum { NZ = 41, TRIALS = 12 };
    double t[NZ], m[NZ];
    size_t generated = sizeof spectrum_cases / sizeof *spectrum_cases;
    for (size_t i = 0; i < generated + sizeof extra_cases / sizeof *extra_cases; ++i) {
        const struct spectrum_case *c =
            i < generated ? spectrum_cases + i : extra_cases + i - generated;
        if (c->role == TT_SPECTRUM_WEIGHT)
            continue;
        origin = c->origin;
        tt_config *config = parse(c->config), *grid_config = parse(c->grid);
        tt_grid grid = {0};
        tt_spectrum *s = NULL;
        tt_spectrum_plan *plan = NULL;
        tt_error error;
        require(tt_grid_from_config(tt_config_root(grid_config), &source, &grid, &error) == TT_OK,
                error.message);
        require(tt_spectrum_create(c->role, tt_config_root(config), &grid, &source, &s, &error) ==
                    TT_OK,
                error.message);
        require(tt_spectrum_plan_create(s, &plan) == TT_OK, "spectrum plan");
        require(!tt_spectrum_plan_generic(plan), "classification guard fell back");
        size_t nw = grid.cells;
        const size_t *varying;
        size_t count = tt_spectrum_plan_varying(plan, &varying);
        const double *constant = tt_spectrum_plan_constant(plan);
        double *generic = malloc(nw * NZ * sizeof(double)),
               *rows = malloc((count ? count : 1) * NZ * sizeof(double)),
               *scratch = malloc((tt_spectrum_plan_scratch(plan, NZ) + 1) * sizeof(double)),
               *value = malloc(nw * NZ * sizeof(double));
        require(generic && rows && scratch && value, "allocation");
        for (int trial = 0; trial < TRIALS; ++trial) {
            column(NZ, trial, t, m);
            tt_status a = tt_spectrum_calculate(s, NZ, t, m, generic);
            tt_status b = tt_spectrum_plan_evaluate(plan, NZ, t, m, scratch, rows);
            require(a == b, "evaluation status");
            if (a)
                continue;
            for (size_t w = 0, k = 0; w < nw; ++w) {
                int varies = k < count && varying[k] == w;
                for (size_t z = 0; z < NZ; ++z)
                    value[w * NZ + z] = varies ? rows[k * NZ + z] : constant[w];
                k += (size_t)varies;
            }
            same(value, generic, nw * NZ, "spectrum");
        }
        ++spectra;
        free(generic);
        free(rows);
        free(scratch);
        free(value);
        tt_spectrum_plan_destroy(plan);
        tt_spectrum_destroy(s);
        tt_grid_destroy(&grid);
        tt_config_destroy(config);
        tt_config_destroy(grid_config);
    }
}

typedef struct {
    const char *name, *units;
    double *edges, *top;
} density;
/* All outputs and radiator views of the last run. */
static size_t collect(tt_model *model, const tt_model_result *r, double *out) {
    tt_model_dimensions s = tt_model_shape(model);
    size_t n = s.levels * s.wavelengths, k = 0, layers = (s.levels - 1) * s.wavelengths;
    const double *fields[] = {r->radiation.direct_flux,
                              r->radiation.down_flux,
                              r->radiation.up_flux,
                              r->radiation.direct_irradiance,
                              r->radiation.down_irradiance,
                              r->radiation.up_irradiance,
                              r->actinic_flux,
                              r->spectral_irradiance};
    for (size_t f = 0; f < 8; ++f, k += n)
        if (out)
            memcpy(out + k, fields[f], n * sizeof(double));
    const double *rates[] = {r->photolysis, r->heating, r->dose};
    size_t counts[] = {s.photolysis, s.heating, s.dose};
    for (size_t i = 0; i < 3; ++i, k += counts[i - 1] * s.levels)
        if (out)
            memcpy(out + k, rates[i], counts[i] * s.levels * sizeof(double));
    const tt_model_profile_entry *unused = NULL;
    (void)unused;
    return k + 3 * layers * 16; /* room for radiator views, filled by views() */
}
static size_t views(tt_model *model, const char *const *names, size_t count, double *out) {
    tt_model_dimensions s = tt_model_shape(model);
    size_t layers = (s.levels - 1) * s.wavelengths, k = 0;
    for (size_t i = 0; i < count; ++i) {
        const tt_radiator *r = tt_model_radiator(model, names[i]);
        const double *v[] = {r->optical_depth, r->single_scattering_albedo, r->asymmetry};
        for (size_t f = 0; f < 3; ++f, k += layers)
            memcpy(out + k, v[f], layers * sizeof(double));
    }
    return k;
}

static void model_runs(void) {
    tt_data_source source = {resolve, NULL};
    const char *densities[][2] = {
        {"air", "molecule cm-3"}, {"O2", "molecule cm-3"}, {"O3", "molecule cm-3"}};
    for (size_t c = 0; c < sizeof model_cases / sizeof *model_cases; ++c) {
        if (!model_cases[c].runnable)
            continue;
        origin = model_cases[c].origin;
        tt_config *config = parse(model_cases[c].config);
        tt_model *model = NULL;
        tt_error error;
        require(tt_model_create(tt_config_root(config), &source, &model, &error) == TT_OK,
                error.message);
        require(tt_model_set_execution(model, TT_EXECUTION_PREPARED) == TT_OK, "no plan");
        require(tt_model_plan_generic_spectra(model->plan) == 0, "classification guard fell back");
        const char *solver = tt_config_string(tt_config_get(
            tt_config_get(tt_config_get(tt_config_root(config), "radiative transfer"), "solver"),
            "type"));
        int host = !strcmp(solver, "from host");
        tt_config_node radiators =
            tt_config_get(tt_config_get(tt_config_root(config), "radiative transfer"), "radiators");
        size_t nr = tt_config_size(radiators);
        const char **names = calloc(nr + 1, sizeof *names);
        require(names != NULL, "allocation");
        for (size_t i = 0; i < nr; ++i)
            names[i] = tt_config_string(tt_config_get(tt_config_at(radiators, i), "name"));
        tt_model_dimensions s = tt_model_shape(model);
        size_t nz = s.levels, nl = nz - 1, n = nz * s.wavelengths;
        size_t size =
            8 * n + nz * (s.photolysis + s.heating + s.dose) + 3 * nl * s.wavelengths * nr;
        double *prepared = malloc(size * sizeof(double)), *generic = malloc(size * sizeof(double));
        double *radiation = calloc(6 * n, sizeof(double));
        const tt_profile *temperature = tt_model_profile(model, "temperature", "K");
        double *t0 = malloc(nz * sizeof(double)), *t = malloc(nz * sizeof(double));
        require(prepared && generic && radiation && temperature && t0 && t, "allocation");
        memcpy(t0, temperature->edges, nz * sizeof(double));
        /* Copies of the density profiles, scaled for each column. */
        tt_profile saved[3] = {{0}};
        int present[3] = {0};
        for (size_t d = 0; d < 3; ++d) {
            const tt_profile *p = tt_model_profile(model, densities[d][0], densities[d][1]);
            if (!p || p->cells != nl)
                continue;
            present[d] = 1;
            require(tt_profile_create(nl, saved + d) == TT_OK, "allocation");
            memcpy(saved[d].edges, p->edges, nz * sizeof(double));
            memcpy(saved[d].extended_columns, p->extended_columns, nz * sizeof(double));
        }
        double *edges = malloc(nz * sizeof(double)), *cols = malloc(nl * sizeof(double));
        require(edges && cols, "allocation");
        enum { COLUMNS = 14 };
        for (int col = 0; col < COLUMNS; ++col, ++columns) {
            /* Daylight, twilight and night; a hot level beyond the plans'
             * limit and an invalid level make the plan decline. */
            double zenith = col < 8 ? uniform(0, 89.9) : col < 11 ? uniform(89.9, 98) : 140;
            double scale = uniform(.983, 1.017), factor = uniform(.9, 1.1);
            for (size_t z = 0; z < nz; ++z)
                t[z] = t0[z] + uniform(-15, 15);
            if (col == 12)
                t[nz / 2] = 2.e4;
            if (col == 13)
                t[nz - 1] = -5;
            require(tt_model_update_profile(model, "temperature", "K", nl, t, NULL, NULL, NULL,
                                            NULL) == TT_OK,
                    "temperature update");
            for (size_t d = 0; d < 3; ++d) {
                if (!present[d])
                    continue;
                double f = d == 2 ? uniform(.7, 1.3) : factor;
                for (size_t z = 0; z < nz; ++z)
                    edges[z] = f * saved[d].edges[z];
                for (size_t z = 0; z < nl; ++z)
                    cols[z] = f * saved[d].extended_columns[z];
                double top = f * saved[d].extended_columns[nl];
                require(tt_model_update_profile(model, densities[d][0], densities[d][1], nl, edges,
                                                NULL, cols, NULL, &top) == TT_OK,
                        "density update");
            }
            if (host) {
                for (size_t i = 0; i < 6 * n; ++i)
                    radiation[i] = uniform(-.1, 1);
                tt_radiation r = {radiation,         radiation + n,     radiation + 2 * n,
                                  radiation + 3 * n, radiation + 4 * n, radiation + 5 * n};
                require(tt_model_set_radiation(model, &r) == TT_OK, "host radiation");
            }
            const tt_model_result *result = NULL;
            tt_error e1, e2;
            require(tt_model_set_execution(model, TT_EXECUTION_PREPARED) == TT_OK, "prepared");
            tt_status a = tt_model_run(model, zenith, scale, &result, &e1);
            size_t k = 0;
            if (!a) {
                k = collect(model, result, prepared) - 3 * nl * s.wavelengths * 16;
                k += views(model, names, nr, prepared + k);
            }
            require(tt_model_set_execution(model, TT_EXECUTION_GENERIC) == TT_OK, "generic");
            tt_status b = tt_model_run(model, zenith, scale, &result, &e2);
            require(a == b && !strcmp(e1.message, e2.message), "run status and message");
            if (a)
                continue;
            size_t j = collect(model, result, generic) - 3 * nl * s.wavelengths * 16;
            j += views(model, names, nr, generic + j);
            require(j == k, "output size");
            same(prepared, generic, k, "outputs and radiator views");
        }
        ++models;
        for (size_t d = 0; d < 3; ++d)
            tt_profile_destroy(saved + d);
        free(edges);
        free(cols);
        free(prepared);
        free(generic);
        free(radiation);
        free(t0);
        free(t);
        free(names);
        tt_model_destroy(model);
        tt_config_destroy(config);
    }
}

/* Blocked Delta-Eddington against mixing plus the generic solver, on the
 * atmospheres of test_radiation.c, with partial and multiple blocks. */
static size_t atmospheres;
static void blocked_transfer(void) {
    enum { NL = 81, NR = 3 };
    const double angles[] = {0, 30, 60, 85, 89.9, 90, 90.01, 93, 99, 110, 180};
    const int sizes[] = {1, 2, 5, 16, 81}, air[] = {0, 1, 0};
    const size_t widths[] = {7, 37};
    origin = "blocked Delta-Eddington";
    for (size_t wi = 0; wi < 2; ++wi) {
        size_t nw = widths[wi], field = (NL + 1) * nw;
        double *tau = malloc(3 * NR * NL * nw * sizeof(double)), *ssa = tau + NR * NL * nw,
               *g = ssa + NR * NL * nw, *mix = malloc(3 * NL * nw * sizeof(double)),
               *expected = malloc(12 * field * sizeof(double)), *actual = expected + 6 * field,
               *work = malloc(tt_block_delta_eddington_work(NL) * sizeof(double)),
               *path = malloc(NL * (NL + 1) * sizeof(double)),
               *surface = malloc(nw * sizeof(double));
        int *nid = malloc((NL + 1) * sizeof(int));
        double z[NL + 1];
        require(tau && mix && expected && work && path && surface && nid, "allocation");
        for (size_t ni = 0; ni < 5; ++ni)
            for (size_t ai = 0; ai < 11; ++ai)
                for (int kind = 0; kind < 5; ++kind, ++atmospheres) {
                    size_t n = (size_t)sizes[ni], count = n * nw, levels = (n + 1) * nw;
                    for (size_t j = 0; j <= n; ++j)
                        z[j] = 2.3 + 120.0 * (double)j / (double)n;
                    for (size_t w = 0; w < nw; ++w)
                        surface[w] = (double)w / (double)(nw - 1);
                    for (size_t r = 0; r < NR; ++r)
                        for (size_t w = 0; w < nw; ++w)
                            for (size_t j = 0; j < n; ++j) {
                                size_t k = r * count + w * n + j;
                                tau[k] =
                                    exp(-(z[j] - 2.3) / 8.0) * (1 + w % 7) * (.1 + r) / (double)n;
                                tau[k] *= kind == 1 ? 1e-10 : kind == 2 ? 1e4 : kind == 4 ? 0 : 1;
                                ssa[k] = kind == 3 ? 1.0 : (double)w / (double)(nw - 1);
                                g[k] = r == 1 ? 0.0 : (-0.2 + 1.1 * (double)w / (double)(nw - 1));
                            }
                    tt_radiator radiators[NR];
                    for (size_t r = 0; r < NR; ++r)
                        radiators[r] =
                            (tt_radiator){tau + r * count, ssa + r * count, g + r * count, air[r]};
                    double scale = kind % 2 ? 0.987 : 1;
                    require(tt_mix_radiators(NR, radiators, n, nw, 1, mix, mix + count,
                                             mix + 2 * count) == TT_OK,
                            "mixing");
                    tt_atmosphere atm = {n, nw, 1, z, mix, mix + count, mix + 2 * count};
                    tt_radiation r1 = {expected,
                                       expected + levels,
                                       expected + 2 * levels,
                                       expected + 3 * levels,
                                       expected + 4 * levels,
                                       expected + 5 * levels};
                    tt_radiation r2 = {actual,
                                       actual + levels,
                                       actual + 2 * levels,
                                       actual + 3 * levels,
                                       actual + 4 * levels,
                                       actual + 5 * levels};
                    require(tt_delta_eddington(&atm, angles[ai], surface, &r1) == TT_OK, "generic");
                    for (size_t i = 0; i < 6 * levels; ++i)
                        expected[i] = expected[i] * scale;
                    require(tt_spherical_geometry(n + 1, z, angles[ai], nid, path) == TT_OK,
                            "geometry");
                    require(tt_block_delta_eddington(NR, radiators, n, nw, nid, path, angles[ai],
                                                     surface, scale, work, &r2) == TT_OK,
                            "blocked");
                    same(actual, expected, 6 * levels, "blocked radiation");
                }
        free(tau);
        free(mix);
        free(expected);
        free(work);
        free(path);
        free(surface);
        free(nid);
    }
}

int main(void) {
    spectrum_plans();
    blocked_transfer();
    model_runs();
    printf("Prepared execution: %zu spectra, %zu atmospheres, %zu models, %zu columns, %zu values "
           "identical\n",
           spectra, atmospheres, models, columns, compared);
    return 0;
}
