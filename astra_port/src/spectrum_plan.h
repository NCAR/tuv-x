/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Prepared evaluation of one spectrum on a model's fixed wavelength grid.
 * A plan separates the wavelengths whose values cannot change with the
 * column's temperature or air density, which are evaluated once with the
 * generic kernel, from those that can. Varying wavelengths are evaluated with
 * the generic kernel's expressions, rearranged only by hoisting subexpressions
 * that depend on wavelength or level alone. Every value is therefore the
 * generic kernel's value bit for bit. */
#ifndef TT_SPECTRUM_PLAN_H
#define TT_SPECTRUM_PLAN_H
#include "spectrum_internal.h"

typedef struct tt_spectrum_plan tt_spectrum_plan;

/* The plan borrows the spectrum, which must outlive it. */
tt_status tt_spectrum_plan_create(const tt_spectrum *spectrum, tt_spectrum_plan **out);
void tt_spectrum_plan_destroy(tt_spectrum_plan *plan);
/* Values at wavelengths that never vary; entries at varying wavelengths are
 * the reference column's values and must not be used. */
const double *tt_spectrum_plan_constant(const tt_spectrum_plan *plan);
/* Ascending indices of the wavelengths whose values vary between columns. */
size_t tt_spectrum_plan_varying(const tt_spectrum_plan *plan, const size_t **indices);
/* The constant classification assumes temperatures at most this value. */
double tt_spectrum_plan_maximum_temperature(const tt_spectrum_plan *plan);
/* Scratch doubles needed to evaluate levels levels. */
size_t tt_spectrum_plan_scratch(const tt_spectrum_plan *plan, size_t levels);
/* Evaluate the varying wavelengths: rows[k * levels + z] for the k-th varying
 * wavelength. Temperatures must be finite, positive and at most the maximum
 * above; air densities must be finite and nonnegative. Parameterized cross
 * sections may return TT_RANGE where the generic kernel does. */
tt_status tt_spectrum_plan_evaluate(const tt_spectrum_plan *plan, size_t levels,
                                    const double *temperature, const double *air, double *scratch,
                                    double *rows);
/* True if the classification guard rejected this spectrum's classification;
 * every wavelength then uses the generic kernel. Tests require false. */
int tt_spectrum_plan_generic(const tt_spectrum_plan *plan);
/* Plans with equal prepared data produce equal values; models evaluate them once. */
int tt_spectrum_plan_equal(const tt_spectrum_plan *a, const tt_spectrum_plan *b);
#endif
