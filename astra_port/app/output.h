#ifndef TURBOTUV_OUTPUT_H
#define TURBOTUV_OUTPUT_H
#include "turbotuv_model.h"
#include <stdio.h>
typedef struct tt_output tt_output;
int tt_output_open(FILE *json, const char *netcdf_path, const tt_model *model, tt_output **out,
                   tt_error *error);
int tt_output_step(tt_output *out, double time, double zenith, double scale,
                   const tt_model_result *result, tt_error *error);
int tt_output_close(tt_output *out, tt_error *error);
/* Release a writer after a failure without completing its document. */
void tt_output_abandon(tt_output *out);
#endif
