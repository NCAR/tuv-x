/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Column benchmark: one model evaluates many columns, as a host model does.
 * Every column changes the temperature, air, O2 and O3 profiles, the solar
 * zenith angle and the radiation scale through the public host-update API,
 * so no column-dependent work can be reused between columns. The program
 * prints the elapsed time and an FNV-1a hash of every output value of every
 * column; equal hashes from two builds mean bit-identical results. Build it
 * without floating-point contraction (-ffp-contract=off, as CMake does), so
 * that the generated column inputs are the same with every compiler.
 *
 * --sweep instead holds the profiles fixed and runs column c at solar zenith
 * angle 89 c / N degrees with unit scale, timing only the runs, and hashes the
 * outputs the original TUV-x core also returns: the six radiation fields and
 * the rates. bench/tuvx_bench.F90 runs the same sweep through TUV-x itself. */
#define _POSIX_C_SOURCE 200809L
#include "turbotuv_model.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}
/* Deterministic and platform-independent column states. */
static uint64_t next(uint64_t *state) {
    *state = *state * 6364136223846793005u + 1442695040888963407u;
    return *state >> 11;
}
static double uniform(uint64_t *state, double lo, double hi) {
    return lo + (hi - lo) * ((double)next(state) / 9007199254740992.0);
}
static uint64_t hash(uint64_t h, size_t n, const double *v) {
    const unsigned char *p = (const unsigned char *)v;
    for (size_t i = 0; i < n * sizeof(double); ++i)
        h = (h ^ p[i]) * 1099511628211u;
    return h;
}
typedef struct {
    const char *name, *units;
    tt_profile original;
    double *edges, *columns;
} profile_state;
static int save(tt_model *m, const char *name, const char *units, profile_state *s) {
    const tt_profile *p = tt_model_profile(m, name, units);
    s->name = name;
    s->units = units;
    if (!p)
        return 0;
    size_t n = p->cells;
    if (tt_profile_create(n, &s->original))
        return -1;
    memcpy(s->original.edges, p->edges, (n + 1) * sizeof(double));
    memcpy(s->original.columns, p->columns, n * sizeof(double));
    memcpy(s->original.extended_columns, p->extended_columns, (n + 1) * sizeof(double));
    s->edges = malloc((2 * n + 1) * sizeof(double));
    s->columns = s->edges ? s->edges + n + 1 : NULL;
    return s->edges ? 1 : -1;
}
/* Scale a density profile, keeping its extra column above the model top. */
static tt_status scale_profile(tt_model *m, profile_state *s, double factor) {
    size_t n = s->original.cells;
    for (size_t i = 0; i <= n; ++i)
        s->edges[i] = factor * s->original.edges[i];
    double top = factor * s->original.extended_columns[n];
    for (size_t i = 0; i < n; ++i)
        s->columns[i] = factor * s->original.extended_columns[i];
    return tt_model_update_profile(m, s->name, s->units, n, s->edges, NULL, s->columns, NULL, &top);
}
static int usage(void) {
    fputs("Usage: turbotuv_bench CONFIG [--columns N] [--repeat R] [--seed S] [--output FILE]\n"
          "                      [--generic] [--sweep]\n",
          stderr);
    return 2;
}
int main(int argc, char **argv) {
    const char *path = NULL, *output = NULL;
    long columns = 1000, repeat = 1;
    unsigned long long seed = 2026;
    int generic = 0, sweep = 0;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--columns") && i + 1 < argc)
            columns = strtol(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--repeat") && i + 1 < argc)
            repeat = strtol(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc)
            seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--output") && i + 1 < argc)
            output = argv[++i];
        else if (!strcmp(argv[i], "--generic"))
            generic = 1;
        else if (!strcmp(argv[i], "--sweep"))
            sweep = 1;
        else if (argv[i][0] == '-' || path)
            return usage();
        else
            path = argv[i];
    }
    if (!path || columns < 1 || repeat < 1)
        return usage();
    tt_config *config = NULL;
    tt_model *model = NULL;
    tt_error error = {0};
    if (tt_config_read_file(path, &config, &error) ||
        tt_model_create(tt_config_root(config), NULL, &model, &error)) {
        fprintf(stderr, "turbotuv_bench: %s\n", error.message);
        return 1;
    }
    tt_config_destroy(config);
    tt_status execution =
        tt_model_set_execution(model, generic ? TT_EXECUTION_GENERIC : TT_EXECUTION_PREPARED);
    tt_model_dimensions shape = tt_model_shape(model);
    size_t nz = shape.levels, nl = nz - 1, n = nz * shape.wavelengths;
    profile_state air, o2, o3;
    const tt_profile *temperature = tt_model_profile(model, "temperature", "K");
    double *t0 = malloc(nz * sizeof(double)), *t = malloc(nz * sizeof(double));
    if (save(model, "air", "molecule cm-3", &air) < 0 ||
        save(model, "O2", "molecule cm-3", &o2) < 0 ||
        save(model, "O3", "molecule cm-3", &o3) < 0 || !temperature || !t0 || !t) {
        fputs("turbotuv_bench: the configuration needs temperature and air profiles\n", stderr);
        return 1;
    }
    memcpy(t0, temperature->edges, nz * sizeof(double));
    FILE *file = output ? fopen(output, "wb") : NULL;
    if (output && !file) {
        perror(output);
        return 1;
    }
    double (*column)[7] = malloc((size_t)columns * sizeof *column);
    if (!column) {
        fputs("turbotuv_bench: cannot allocate column states\n", stderr);
        return 1;
    }
    uint64_t state = seed;
    for (long c = 0; c < columns; ++c) {
        static const double range[7][2] = {{0, 95},   {0.983, 1.017},         {-10, 10},
                                           {0, 8},    {0, 6.283185307179586}, {0.95, 1.05},
                                           {0.7, 1.3}};
        for (size_t k = 0; k < 7; ++k)
            column[c][k] = uniform(&state, range[k][0], range[k][1]);
        if (sweep) {
            column[c][0] = 89.0 * (double)c / (double)columns;
            column[c][1] = 1;
        }
    }
    uint64_t digest = 14695981039346656037u;
    double elapsed = 0;
    /* Column 0 runs once untimed first, so the first timed column is not cold. */
    for (long r = 0; r < repeat; ++r)
        for (long c = -(r == 0); c < columns; ++c) {
            const double *v = column[c < 0 ? 0 : c];
            double start = seconds();
            const tt_model_result *result = NULL;
            tt_status status = TT_OK;
            if (!sweep) {
                for (size_t i = 0; i < nz; ++i)
                    t[i] = t0[i] + v[2] + v[3] * sin(v[4] + 0.1 * (double)i);
                status = tt_model_update_profile(model, "temperature", "K", nl, t, NULL, NULL, NULL,
                                                 NULL);
                if (!status && air.edges)
                    status = scale_profile(model, &air, v[5]);
                if (!status && o2.edges)
                    status = scale_profile(model, &o2, v[5]);
                if (!status && o3.edges)
                    status = scale_profile(model, &o3, v[6]);
            }
            if (!status)
                status = tt_model_run(model, v[0], v[1], &result, &error);
            if (status) {
                fprintf(stderr, "turbotuv_bench: column %ld: %s\n", c, tt_status_string(status));
                return 1;
            }
            if (c >= 0)
                elapsed += seconds() - start;
            if (c < 0 || r)
                continue;
            const double *fields[] = {result->radiation.direct_flux,
                                      result->radiation.down_flux,
                                      result->radiation.up_flux,
                                      result->radiation.direct_irradiance,
                                      result->radiation.down_irradiance,
                                      result->radiation.up_irradiance,
                                      result->actinic_flux,
                                      result->spectral_irradiance,
                                      result->photolysis,
                                      result->heating,
                                      result->dose};
            /* TUV-x returns no actinic flux or irradiance to compare in a sweep. */
            size_t counts[] = {n,
                               n,
                               n,
                               n,
                               n,
                               n,
                               sweep ? 0 : n,
                               sweep ? 0 : n,
                               nz * shape.photolysis,
                               nz * shape.heating,
                               nz * shape.dose};
            for (size_t k = 0; k < 11; ++k) {
                digest = hash(digest, counts[k], fields[k]);
                if (file && fwrite(fields[k], sizeof(double), counts[k], file) != counts[k]) {
                    perror(output);
                    return 1;
                }
            }
        }
    if (file && fclose(file)) {
        perror(output);
        return 1;
    }
    long total = columns * repeat;
    printf("config=%s columns=%ld levels=%zu wavelengths=%zu photolysis=%zu heating=%zu dose=%zu "
           "execution=%s\n",
           path, total, nz, shape.wavelengths, shape.photolysis, shape.heating, shape.dose,
           generic     ? "generic"
           : execution ? "generic (not prepared)"
                       : "prepared");
    printf("elapsed_ms=%.3f per_column_us=%.2f hash=%016llx\n", 1e3 * elapsed,
           1e6 * elapsed / (double)total, (unsigned long long)digest);
    tt_profile_destroy(&air.original);
    tt_profile_destroy(&o2.original);
    tt_profile_destroy(&o3.original);
    free(air.edges);
    free(o2.edges);
    free(o3.edges);
    free(t0);
    free(t);
    free(column);
    tt_model_destroy(model);
    return 0;
}
