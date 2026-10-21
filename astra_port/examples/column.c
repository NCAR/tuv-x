/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include <turbotuv_model.h>
#include <stdio.h>
int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: column CONFIG\n");
        return 2;
    }
    tt_config *config = NULL;
    tt_model *model = NULL;
    const tt_model_result *result = NULL;
    tt_error error;
    tt_status status = tt_config_read_file(argv[1], &config, &error);
    if (!status)
        status = tt_model_create(tt_config_root(config), NULL, &model, &error);
    tt_config_destroy(config); /* The model owns everything needed to run. */
    if (!status)
        status = tt_model_run(model, 30.0, 1.0, &result, &error);
    if (status)
        fprintf(stderr, "%s\n", error.message);
    else {
        tt_model_dimensions shape = tt_model_shape(model);
        for (size_t i = 0; i < shape.photolysis; ++i)
            printf("%s: %.17g s-1 at the surface\n",
                   tt_model_rate_name(model, TT_RATE_PHOTOLYSIS, i),
                   result->photolysis[i * shape.levels]);
    }
    tt_model_destroy(model);
    return status ? 1 : 0;
}
