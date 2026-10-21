/* Configuration and data parsing must not depend on the process locale. A
 * host that selects a locale with a ',' decimal point gets the same model. */
#include "turbotuv_model.h"
#include <locale.h>
#include <stdio.h>
#include <string.h>
static int run(const char *path, double *out, size_t capacity, size_t *count) {
    tt_config *config = NULL;
    tt_model *model = NULL;
    tt_error error;
    const tt_model_result *r = NULL;
    if (tt_config_read_file(path, &config, &error) ||
        tt_model_create(tt_config_root(config), NULL, &model, &error) ||
        tt_model_run(model, 35, 1, &r, &error)) {
        fprintf(stderr, "locale %s: %s\n", setlocale(LC_NUMERIC, NULL), error.message);
        tt_model_destroy(model);
        tt_config_destroy(config);
        return 0;
    }
    tt_model_dimensions s = tt_model_shape(model);
    size_t n = s.levels * s.photolysis, m = s.levels * s.wavelengths;
    *count = n + m;
    if (*count <= capacity) {
        memcpy(out, r->photolysis, n * sizeof(double));
        memcpy(out + n, r->actinic_flux, m * sizeof(double));
    }
    tt_model_destroy(model);
    tt_config_destroy(config);
    return *count <= capacity;
}
int main(int argc, char **argv) {
    static double reference[200000], localized[200000];
    size_t n = 0, m = 0;
    if (argc != 2 || !run(argv[1], reference, 200000, &n))
        return 1;
    const char *names[] = {"de_DE.UTF-8", "de_DE", "fr_FR.UTF-8", "fr_FR", "ru_RU.UTF-8"};
    for (size_t i = 0; i < sizeof names / sizeof *names; ++i) {
        if (!setlocale(LC_NUMERIC, names[i]) || strcmp(localeconv()->decimal_point, ","))
            continue;
        if (!run(argv[1], localized, 200000, &m) || m != n ||
            memcmp(reference, localized, n * sizeof(double))) {
            fprintf(stderr, "locale %s changed the model\n", names[i]);
            return 1;
        }
        printf("Locale %s: %zu values identical to the C locale\n", names[i], n);
        return 0;
    }
    printf("No locale with a comma decimal point is installed; skipped\n");
    return 0;
}
