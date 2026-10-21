/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#define _POSIX_C_SOURCE 200809L
#include "output.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static void usage(FILE *f) {
    fputs("Usage: turbotuv CONFIG [--output FILE.json | --netcdf FILE.nc]\n"
          "                      [--zenith DEGREES] [--scale FACTOR]\n"
          "\nRuns the configured time/solar profiles, or one supplied zenith angle.\n"
          "JSON goes to stdout unless an output path is supplied. Scientific data\n"
          "are compiled in; CONFIG may be an original TUV-x JSON or YAML file.\n",
          f);
}
static int number(const char *s, double *out) {
    char *end = NULL;
    errno = 0;
    *out = strtod(s, &end);
    return s != end && !*end && errno != ERANGE && isfinite(*out);
}
int main(int argc, char **argv) {
    const char *config_path = NULL, *output_path = NULL;
    double zenith = NAN, scale = NAN;
    int netcdf = 0;
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (!strcmp(arg, "--help") || !strcmp(arg, "-h")) {
            usage(stdout);
            return 0;
        }
        if (!strcmp(arg, "--version")) {
            puts("TurboTUV " TT_VERSION);
            return 0;
        }
        if (!strcmp(arg, "--output") || !strcmp(arg, "-o") || !strcmp(arg, "--netcdf")) {
            if (output_path || ++i == argc) {
                usage(stderr);
                return 2;
            }
            netcdf = !strcmp(arg, "--netcdf");
            output_path = argv[i];
        } else if (!strcmp(arg, "--zenith") || !strcmp(arg, "--scale")) {
            double *value = !strcmp(arg, "--zenith") ? &zenith : &scale;
            if (++i == argc || !number(argv[i], value)) {
                fputs("Invalid numeric argument\n", stderr);
                return 2;
            }
        } else if (arg[0] == '-' || config_path) {
            usage(stderr);
            return 2;
        } else
            config_path = arg;
    }
    if (!config_path || (!isnan(zenith) && (zenith < 0 || zenith > 180))) {
        usage(stderr);
        return 2;
    }
    tt_config *config = NULL;
    tt_model *model = NULL;
    tt_error error = {0};
    tt_output *writer = NULL;
    FILE *file = stdout;
    char *temporary = NULL;
    int success = 0;
    if (tt_config_read_file(config_path, &config, &error) ||
        tt_model_create(tt_config_root(config), NULL, &model, &error))
        goto done;
    tt_config_destroy(config);
    config = NULL;
    const tt_grid *time = tt_model_grid(model, "time", "hours");
    const tt_profile *angles = tt_model_profile(model, "solar zenith angle", "degrees"),
                     *distances = tt_model_profile(model, "Earth-Sun distance", "AU");
    size_t steps = 1;
    if (isnan(zenith)) {
        if (!time || !angles || !distances || time->cells != angles->cells ||
            time->cells != distances->cells) {
            snprintf(error.message, sizeof error.message,
                     "Configure matching time, solar zenith angle and Earth-Sun distance profiles, "
                     "or supply --zenith");
            goto done;
        }
        steps = time->cells + 1;
    }
    if (output_path) {
        temporary = malloc(strlen(output_path) + 16);
        if (!temporary) {
            snprintf(error.message, sizeof error.message, "Cannot allocate output path");
            goto done;
        }
        sprintf(temporary, "%s.tmp.XXXXXX", output_path);
        int fd = mkstemp(temporary);
        if (fd < 0) {
            snprintf(error.message, sizeof error.message, "Cannot create output: %s",
                     strerror(errno));
            goto done;
        }
        /* mkstemp creates mode 0600; publish with the replaced file's mode,
         * or the mode a new file would have. */
        struct stat existing;
        mode_t mask = umask(0);
        umask(mask);
        mode_t mode = stat(output_path, &existing) == 0 && S_ISREG(existing.st_mode)
                          ? existing.st_mode & 07777
                          : 0666 & ~mask;
        if (fchmod(fd, mode)) {
            close(fd);
            snprintf(error.message, sizeof error.message, "Cannot set output mode: %s",
                     strerror(errno));
            goto done;
        }
        if (netcdf) {
            close(fd);
            file = NULL;
        } else {
            file = fdopen(fd, "wb");
            if (!file) {
                close(fd);
                snprintf(error.message, sizeof error.message, "Cannot open output stream");
                goto done;
            }
        }
    }
    if (!tt_output_open(file, netcdf ? temporary : NULL, model, &writer, &error))
        goto done;
    for (size_t i = 0; i < steps; ++i) {
        double angle = isnan(zenith) ? angles->edges[i] : zenith,
               factor = isnan(scale) ? (isnan(zenith) ? distances->edges[i] : 1) : scale,
               hour = isnan(zenith) ? time->edges[i] : 0;
        const tt_model_result *result = NULL;
        if (tt_model_run(model, angle, factor, &result, &error) ||
            !tt_output_step(writer, hour, angle, factor, result, &error))
            goto done;
    }
    success = tt_output_close(writer, &error);
    writer = NULL;
    if (file && file != stdout) {
        if (fclose(file)) {
            success = 0;
            snprintf(error.message, sizeof error.message, "Cannot close output");
        }
        file = NULL;
    }
    if (success && temporary && rename(temporary, output_path)) {
        success = 0;
        snprintf(error.message, sizeof error.message, "Cannot publish output: %s", strerror(errno));
    }
done:
    if (writer) /* leave the JSON unterminated so it cannot pass as complete */
        tt_output_abandon(writer);
    if (file && file != stdout)
        fclose(file);
    if (temporary) {
        if (!success)
            unlink(temporary);
        free(temporary);
    }
    tt_model_destroy(model);
    tt_config_destroy(config);
    if (!success) {
        fprintf(stderr, "TurboTUV: %s", error.message);
        if (error.line)
            fprintf(stderr, " at line %zu, column %zu", error.line, error.column);
        fputc('\n', stderr);
    }
    return success ? 0 : 1;
}
