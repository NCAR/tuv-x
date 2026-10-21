/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Prepared model execution: the generic run's steps and arithmetic, with
 * column-independent spectral values computed once, varying values from
 * spectrum plans, fused rate integration, and Delta-Eddington transfer in
 * blocks of wavelengths. Results are bit-identical to the generic run. */
#ifndef TT_MODEL_PLAN_H
#define TT_MODEL_PLAN_H
#include "model_internal.h"

typedef struct tt_model_plan tt_model_plan;
/* Build after the model is complete. On failure the model runs generically. */
tt_status tt_model_plan_create(const tt_model *model, tt_model_plan **out);
void tt_model_plan_destroy(tt_model_plan *plan);
/* Run one column. Returns 1 with the model's results written, or 0 when the
 * column needs the generic run: invalid or unusual inputs, which the generic
 * run reports or evaluates. Model state is unchanged except result and work
 * storage, which the generic run overwrites. */
int tt_model_plan_run(tt_model_plan *plan, tt_model *model, double zenith, double scale);
/* Spectrum plans whose classification guard fell back to generic evaluation
 * (expected zero; tests require it). */
size_t tt_model_plan_generic_spectra(const tt_model_plan *plan);
/* Delta-Eddington with the radiator mixing fused in, for all wavelengths in
 * blocks: per wavelength the arithmetic of tt_mix_radiators followed by
 * tt_delta_eddington. Outputs are multiplied by scale, as the model does.
 * Returns TT_INVALID or TT_SINGULAR where those would. */
size_t tt_block_delta_eddington_work(size_t layers);
tt_status tt_block_delta_eddington(size_t radiator_count, const tt_radiator *radiators,
                                   size_t layers, size_t wavelengths, const int *nid,
                                   const double *dsdh, double zenith, const double *surface,
                                   double scale, double *work, tt_radiation *out);
#endif
