/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef TT_DISCRETE_INTERNAL_H
#define TT_DISCRETE_INTERNAL_H
#include "turbotuv.h"

/* Internal double-precision DISORT helpers. Matrices and eigenvectors are
 * column-major. Eigen work has n*n+2*n entries. Eigenvalues are real, unordered;
 * eigenvector normalization follows the original solver. */
tt_status tt_discrete_eigen(size_t n, const double *matrix, double *values, double *vectors,
                            double *work);
tt_status tt_discrete_quadrature(size_t n, double *nodes, double *weights);
void tt_discrete_legendre(size_t degree, size_t count, const double *cosines, double *polynomials);
#endif
