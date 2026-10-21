/* Independent model instances share immutable configuration and data only. */
#include "turbotuv_model.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum { WORKERS = 4, RUNS = 4 };
typedef struct {
    tt_config_node config;
    size_t worker, count;
    double *expected;
    int success;
} worker;
static const double angles[WORKERS] = {0, 30, 89.9, 93};
static size_t flatten(const tt_model_result *r, tt_model_dimensions s, double *target,
                      const double *expected) {
    const double *arrays[] = {r->radiation.direct_flux,
                              r->radiation.down_flux,
                              r->radiation.up_flux,
                              r->radiation.direct_irradiance,
                              r->radiation.down_irradiance,
                              r->radiation.up_irradiance,
                              r->actinic_flux,
                              r->spectral_irradiance,
                              r->photolysis,
                              r->heating,
                              r->dose};
    size_t n = s.levels * s.wavelengths, counts[] = {n,
                                                     n,
                                                     n,
                                                     n,
                                                     n,
                                                     n,
                                                     n,
                                                     n,
                                                     s.photolysis * s.levels,
                                                     s.heating * s.levels,
                                                     s.dose * s.levels},
           offset = 0;
    for (size_t k = 0; k < 11; ++k) {
        if (target)
            memcpy(target + offset, arrays[k], counts[k] * sizeof(double));
        if (expected && memcmp(expected + offset, arrays[k], counts[k] * sizeof(double)))
            return 0;
        offset += counts[k];
    }
    return offset;
}
static void *calculate(void *argument) {
    worker *w = argument;
    tt_error error;
    tt_model *m = NULL;
    w->success = tt_model_create(w->config, NULL, &m, &error) == TT_OK;
    for (size_t i = 0; i < RUNS && w->success; ++i) {
        const tt_model_result *result = NULL;
        w->success = tt_model_run(m, angles[w->worker] + 0.01 * i, 0.99 + 0.001 * i, &result,
                                  &error) == TT_OK;
        if (w->success)
            w->success =
                flatten(result, tt_model_shape(m), NULL, w->expected + i * w->count) == w->count;
    }
    tt_model_destroy(m);
    return NULL;
}
int main(int argc, char **argv) {
    size_t compared = 0;
    for (int c = 1; c < argc; ++c) {
        tt_config *config = NULL;
        tt_error error;
        if (tt_config_read_file(argv[c], &config, &error)) {
            fprintf(stderr, "%s\n", error.message);
            return 1;
        }
        tt_model *serial = NULL;
        if (tt_model_create(tt_config_root(config), NULL, &serial, &error)) {
            fprintf(stderr, "%s\n", error.message);
            return 1;
        }
        tt_model_dimensions s = tt_model_shape(serial);
        size_t count = s.levels * (8 * s.wavelengths + s.photolysis + s.heating + s.dose);
        worker workers[WORKERS];
        pthread_t threads[WORKERS];
        for (size_t w = 0; w < WORKERS; ++w) {
            workers[w] =
                (worker){tt_config_root(config), w, count, calloc(RUNS * count, sizeof(double)), 0};
            if (!workers[w].expected)
                return 1;
            for (size_t i = 0; i < RUNS; ++i) {
                const tt_model_result *result = NULL;
                if (tt_model_run(serial, angles[w] + 0.01 * i, 0.99 + 0.001 * i, &result, &error))
                    return 1;
                if (flatten(result, s, workers[w].expected + i * count, NULL) != count)
                    return 1;
            }
        }
        tt_model_destroy(serial);
        for (size_t w = 0; w < WORKERS; ++w)
            if (pthread_create(threads + w, NULL, calculate, workers + w))
                return 1;
        for (size_t w = 0; w < WORKERS; ++w) {
            if (pthread_join(threads[w], NULL) || !workers[w].success) {
                fprintf(stderr, "Concurrent worker %zu differs from serial result\n", w);
                return 1;
            }
            compared += count * RUNS;
            free(workers[w].expected);
        }
        tt_config_destroy(config);
    }
    printf("Concurrent models: %zu bit-identical values across %d workers\n", compared, WORKERS);
    return 0;
}
