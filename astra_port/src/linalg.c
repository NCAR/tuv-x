/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv.h"
#include <math.h>
#include <stdint.h>

static size_t smaller(size_t a, size_t b) {
    return a < b ? a : b;
}
static size_t larger(size_t a, size_t b) {
    return a > b ? a : b;
}
static double *entry(const tt_matrix *a, size_t row, size_t col) {
    size_t index = a->layout == TT_MATRIX_DENSE ? row : a->lower + a->upper + row - col;
    return a->values + col * a->leading_dimension + index;
}
static size_t below(const tt_matrix *a, size_t col) {
    size_t count = a->order - col - 1;
    return a->layout == TT_MATRIX_DENSE ? count : smaller(count, a->lower);
}
static size_t above(const tt_matrix *a, size_t col) {
    return a->layout == TT_MATRIX_DENSE ? col : smaller(col, a->lower + a->upper);
}
static int valid(const tt_matrix *a) {
    if (!a || !a->order || !a->values || (unsigned)a->layout > TT_MATRIX_BAND ||
        !a->leading_dimension || a->order > SIZE_MAX / sizeof(double) / a->leading_dimension)
        return 0;
    if (a->layout == TT_MATRIX_DENSE)
        return a->leading_dimension >= a->order;
    return a->upper < SIZE_MAX && a->lower <= (SIZE_MAX - a->upper - 1) / 2 &&
           a->leading_dimension >= 2 * a->lower + a->upper + 1;
}
static void swap(double *a, double *b) {
    double t = *a;
    *a = *b;
    *b = t;
}
static void scale(size_t n, double factor, double *x) {
    for (size_t i = 0; i < n; ++i)
        x[i] = factor * x[i];
}
static void axpy(size_t n, double factor, const double *x, double *y) {
    if (factor != 0)
        for (size_t i = 0; i < n; ++i)
            y[i] = y[i] + factor * x[i];
}
static double dot(size_t n, const double *x, const double *y) {
    double sum = 0;
    for (size_t i = 0; i < n; ++i)
        sum = sum + x[i] * y[i];
    return sum;
}
static double absolute_sum(size_t n, const double *x) {
    double sum = 0;
    for (size_t i = 0; i < n; ++i)
        sum = sum + fabs(x[i]);
    return sum;
}
static double normalize(size_t n, double *x) {
    double s = 1 / absolute_sum(n, x);
    scale(n, s, x);
    return s;
}

/* Dense and packed matrices share the same elimination and condition
 * estimator. Only the extent and location of each contiguous column differ. */
static tt_status factor(tt_matrix *a, size_t *pivots) {
    size_t n = a->order, ju = 0, fill = 0;
    tt_status status = TT_OK;
    if (a->layout == TT_MATRIX_BAND) {
        size_t m = a->lower + a->upper;
        for (size_t j = a->upper + 1; j < smaller(n, m + 1) - 1; ++j)
            for (size_t i = m - j; i < a->lower; ++i)
                a->values[j * a->leading_dimension + i] = 0;
        fill = smaller(n, m + 1) - 1;
    }
    for (size_t k = 0; k + 1 < n; ++k) {
        if (a->layout == TT_MATRIX_BAND && fill < n) {
            for (size_t i = 0; i < a->lower; ++i)
                a->values[fill * a->leading_dimension + i] = 0;
            ++fill;
        }
        size_t count = below(a, k), pivot = 0;
        double *column = entry(a, k, k);
        double maximum = 0;
        for (size_t i = 0; i <= count; ++i)
            if (fabs(column[i]) > maximum) {
                maximum = fabs(column[i]);
                pivot = i;
            }
        pivots[k] = k + pivot;
        if (column[pivot] == 0) {
            status = TT_SINGULAR;
            continue;
        }
        swap(column + pivot, column);
        scale(count, -1 / column[0], column + 1);
        ju =
            a->layout == TT_MATRIX_DENSE ? n - 1 : smaller(larger(ju, a->upper + pivots[k]), n - 1);
        for (size_t j = k + 1; j <= ju; ++j) {
            double *target = entry(a, k, j), *selected = entry(a, pivots[k], j);
            double t = *selected;
            swap(target, selected);
            axpy(count, t, column + 1, target + 1);
        }
    }
    pivots[n - 1] = n - 1;
    if (*entry(a, n - 1, n - 1) == 0)
        status = TT_SINGULAR;
    return status;
}

static double condition(const tt_matrix *a, const size_t *pivots, double norm, double *z) {
    size_t n = a->order, ju = 0;
    double ek = 1;
    for (size_t i = 0; i < n; ++i)
        z[i] = 0;
    /* Solve U^T w = e, choosing e to maximize local growth. */
    for (size_t k = 0; k < n; ++k) {
        double diagonal = *entry(a, k, k);
        if (z[k] != 0)
            ek = copysign(fabs(ek), -z[k]);
        if (fabs(ek - z[k]) > fabs(diagonal)) {
            double s = fabs(diagonal) / fabs(ek - z[k]);
            scale(n, s, z);
            ek = s * ek;
        }
        double wk = ek - z[k], wkm = -ek - z[k], s = fabs(wk), sm = fabs(wkm);
        if (diagonal != 0) {
            wk = wk / diagonal;
            wkm = wkm / diagonal;
        } else {
            wk = 1;
            wkm = 1;
        }
        ju =
            a->layout == TT_MATRIX_DENSE ? n - 1 : smaller(larger(ju, a->upper + pivots[k]), n - 1);
        for (size_t j = k + 1; j <= ju; ++j) {
            double u = *entry(a, k, j);
            sm = sm + fabs(z[j] + wkm * u);
            z[j] = z[j] + wk * u;
            s = s + fabs(z[j]);
        }
        if (s < sm) {
            double t = wkm - wk;
            wk = wkm;
            for (size_t j = k + 1; j <= ju; ++j)
                z[j] = z[j] + t * (*entry(a, k, j));
        }
        z[k] = wk;
    }
    normalize(n, z);
    /* Solve L^T y = w, then L v = y. */
    for (size_t k = n; k-- > 0;) {
        size_t count = below(a, k);
        if (k + 1 < n)
            z[k] = z[k] + dot(count, entry(a, k, k) + 1, z + k + 1);
        if (fabs(z[k]) > 1)
            scale(n, 1 / fabs(z[k]), z);
        swap(z + pivots[k], z + k);
    }
    normalize(n, z);
    double ynorm = 1;
    for (size_t k = 0; k < n; ++k) {
        double t = z[pivots[k]];
        swap(z + pivots[k], z + k);
        if (k + 1 < n)
            axpy(below(a, k), t, entry(a, k, k) + 1, z + k + 1);
        if (fabs(z[k]) > 1) {
            double s = 1 / fabs(z[k]);
            scale(n, s, z);
            ynorm = s * ynorm;
        }
    }
    double s = normalize(n, z);
    ynorm = s * ynorm;
    /* Solve U z = v, scaling to avoid overflow. */
    for (size_t k = n; k-- > 0;) {
        double diagonal = *entry(a, k, k);
        if (fabs(z[k]) > fabs(diagonal)) {
            s = fabs(diagonal) / fabs(z[k]);
            scale(n, s, z);
            ynorm = s * ynorm;
        }
        z[k] = diagonal != 0 ? z[k] / diagonal : 1;
        size_t count = above(a, k);
        axpy(count, -z[k], entry(a, k - count, k), z + k - count);
    }
    s = normalize(n, z);
    ynorm = s * ynorm;
    return norm != 0 ? ynorm / norm : 0;
}

tt_status tt_lu_factor(tt_matrix *a, size_t *pivots, double *reciprocal_condition, double *work) {
    if (!valid(a) || !pivots || (reciprocal_condition && !work))
        return TT_INVALID;
    double norm = 0;
    if (reciprocal_condition) {
        size_t length = a->lower + 1, first = a->lower + a->upper;
        for (size_t j = 0; j < a->order; ++j) {
            if (a->layout == TT_MATRIX_DENSE)
                norm = fmax(norm, absolute_sum(a->order, entry(a, 0, j)));
            else {
                norm =
                    fmax(norm, absolute_sum(length, a->values + j * a->leading_dimension + first));
                if (first > a->lower)
                    --first;
                if (j < a->upper)
                    ++length;
                if (j + a->lower + 1 >= a->order)
                    --length;
            }
        }
    }
    tt_status status = factor(a, pivots);
    if (reciprocal_condition)
        *reciprocal_condition = condition(a, pivots, norm, work);
    return status;
}

tt_status tt_lu_solve(const tt_matrix *a, const size_t *pivots, int transpose, double *rhs) {
    if (!valid(a) || !pivots || !rhs)
        return TT_INVALID;
    size_t n = a->order;
    for (size_t k = 0; k < n; ++k) {
        if (pivots[k] < k || pivots[k] >= n || pivots[k] - k > below(a, k))
            return TT_INVALID;
        if (*entry(a, k, k) == 0)
            return TT_SINGULAR;
    }
    if (!transpose) {
        for (size_t k = 0; k + 1 < n; ++k) {
            double t = rhs[pivots[k]];
            swap(rhs + pivots[k], rhs + k);
            axpy(below(a, k), t, entry(a, k, k) + 1, rhs + k + 1);
        }
        for (size_t k = n; k-- > 0;) {
            rhs[k] = rhs[k] / (*entry(a, k, k));
            size_t count = above(a, k);
            axpy(count, -rhs[k], entry(a, k - count, k), rhs + k - count);
        }
    } else {
        for (size_t k = 0; k < n; ++k) {
            size_t count = above(a, k);
            double t = dot(count, entry(a, k - count, k), rhs + k - count);
            rhs[k] = (rhs[k] - t) / (*entry(a, k, k));
        }
        for (size_t k = n - 1; k-- > 0;) {
            rhs[k] = rhs[k] + dot(below(a, k), entry(a, k, k) + 1, rhs + k + 1);
            swap(rhs + pivots[k], rhs + k);
        }
    }
    return TT_OK;
}
