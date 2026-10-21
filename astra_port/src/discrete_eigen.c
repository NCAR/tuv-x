/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Real-eigenvalue specialization of TUV-x ASYMTX, adapted there from EISPACK:
 * Parlett-Reinsch balancing and Martin-Wilkinson Hessenberg/QR reduction. */
#include "discrete_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

typedef ptrdiff_t index_t;
/* Local matrix notation keeps the row/column convention visible in equations. */
#define A(i, j) a[(size_t)(j) * order + (size_t)(i)]
#define V(i, j) v[(size_t)(j) * order + (size_t)(i)]
static void exchange(double *x, double *y) {
    double t = *x;
    *x = *y;
    *y = t;
}

static void balance(size_t order, double *a, double *work, index_t *low, index_t *high) {
    index_t m = (index_t)order, l = 0, k = m - 1;
    for (;;) {
        int found = 0;
        for (index_t j = k; j >= 0; --j) {
            double row = 0;
            for (index_t i = 0; i <= k; ++i)
                if (i != j)
                    row = row + fabs(A(j, i));
            if (row == 0) {
                work[k] = (double)j;
                if (j != k) {
                    for (index_t i = 0; i <= k; ++i)
                        exchange(&A(i, j), &A(i, k));
                    for (index_t i = l; i < m; ++i)
                        exchange(&A(j, i), &A(k, i));
                }
                --k;
                found = 1;
                break;
            }
        }
        if (!found)
            break;
    }
    for (;;) {
        int found = 0;
        for (index_t j = l; j <= k; ++j) {
            double col = 0;
            for (index_t i = l; i <= k; ++i)
                if (i != j)
                    col = col + fabs(A(i, j));
            if (col == 0) {
                work[l] = (double)j;
                if (j != l) {
                    for (index_t i = 0; i <= k; ++i)
                        exchange(&A(i, j), &A(i, l));
                    for (index_t i = l; i < m; ++i)
                        exchange(&A(j, i), &A(l, i));
                }
                ++l;
                found = 1;
                break;
            }
        }
        if (!found)
            break;
    }
    for (index_t i = l; i <= k; ++i)
        work[i] = 1;
    for (;;) {
        int converged = 1;
        for (index_t i = l; i <= k; ++i) {
            double col = 0, row = 0;
            for (index_t j = l; j <= k; ++j)
                if (j != i) {
                    col = col + fabs(A(j, i));
                    row = row + fabs(A(i, j));
                }
            double f = 1, g = row / 16, h = col + row;
            while (col < g) {
                f = f * 16;
                col = col * 256;
            }
            g = row * 16;
            while (col >= g) {
                f = f / 16;
                col = col / 256;
            }
            if ((col + row) / f < .95 * h) {
                work[i] = work[i] * f;
                converged = 0;
                for (index_t j = l; j < m; ++j)
                    A(i, j) = A(i, j) / f;
                for (index_t j = 0; j <= k; ++j)
                    A(j, i) = A(j, i) * f;
            }
        }
        if (converged)
            break;
    }
    *low = l;
    *high = k;
}

static void hessenberg(size_t order, double *a, double *v, double *work, index_t l, index_t k) {
    index_t m = (index_t)order;
    for (index_t n = l + 1; n <= k - 1; ++n) {
        double h = 0, scale = 0;
        work[n + m] = 0;
        for (index_t i = n; i <= k; ++i)
            scale = scale + fabs(A(i, n - 1));
        if (scale == 0)
            continue;
        for (index_t i = k; i >= n; --i) {
            work[i + m] = A(i, n - 1) / scale;
            h = h + work[i + m] * work[i + m];
        }
        double g = -copysign(sqrt(h), work[n + m]);
        h = h - work[n + m] * g;
        work[n + m] = work[n + m] - g;
        for (index_t j = n; j < m; ++j) {
            double f = 0;
            for (index_t i = k; i >= n; --i)
                f = f + work[i + m] * A(i, j);
            for (index_t i = n; i <= k; ++i)
                A(i, j) = A(i, j) - work[i + m] * f / h;
        }
        for (index_t i = 0; i <= k; ++i) {
            double f = 0;
            for (index_t j = k; j >= n; --j)
                f = f + work[j + m] * A(i, j);
            for (index_t j = n; j <= k; ++j)
                A(i, j) = A(i, j) - work[j + m] * f / h;
        }
        work[n + m] = scale * work[n + m];
        A(n, n - 1) = scale * g;
    }
    for (index_t n = k - 2; n >= l; --n) {
        double f = A(n + 1, n);
        if (f == 0)
            continue;
        f = f * work[n + 1 + m];
        for (index_t i = n + 2; i <= k; ++i)
            work[i + m] = A(i, n);
        for (index_t j = 0; j < m; ++j) {
            double g = 0;
            for (index_t i = n + 1; i <= k; ++i)
                g = g + work[i + m] * V(i, j);
            g = g / f;
            for (index_t i = n + 1; i <= k; ++i)
                V(i, j) = V(i, j) + g * work[i + m];
        }
    }
}

static tt_status qr(size_t order, double *a, double *v, double *values, index_t l, index_t k,
                    double norm) {
    index_t m = (index_t)order, n = k;
    double shift = 0;
    while (n >= l) {
        unsigned iterations = 0;
        for (;;) {
            index_t lb = n;
            for (; lb > l; --lb) {
                double s = fabs(A(lb - 1, lb - 1)) + fabs(A(lb, lb));
                if (s == 0)
                    s = norm;
                if (fabs(A(lb, lb - 1)) <= DBL_EPSILON * s)
                    break;
            }
            double x = A(n, n);
            if (lb == n) {
                A(n, n) = x + shift;
                values[n] = A(n, n);
                --n;
                break;
            }
            index_t n1 = n - 1, n2 = n - 2;
            double y = A(n1, n1), w = A(n, n1) * A(n1, n);
            if (lb == n1) {
                double p = (y - x) * .5, q = p * p + w, z = sqrt(fabs(q));
                if (q < 0)
                    return TT_RANGE;
                A(n, n) = x + shift;
                x = A(n, n);
                A(n1, n1) = y + shift;
                z = p + copysign(z, p);
                values[n1] = x + z;
                values[n] = values[n1];
                if (z != 0)
                    values[n] = x - w / z;
                x = A(n, n1);
                double r = sqrt(x * x + z * z);
                if (r == 0)
                    return TT_RANGE;
                p = x / r;
                q = z / r;
                for (index_t j = n1; j < m; ++j) {
                    z = A(n1, j);
                    A(n1, j) = q * z + p * A(n, j);
                    A(n, j) = q * A(n, j) - p * z;
                }
                for (index_t i = 0; i <= n; ++i) {
                    z = A(i, n1);
                    A(i, n1) = q * z + p * A(i, n);
                    A(i, n) = q * A(i, n) - p * z;
                }
                for (index_t i = l; i <= k; ++i) {
                    z = V(i, n1);
                    V(i, n1) = q * z + p * V(i, n);
                    V(i, n) = q * V(i, n) - p * z;
                }
                n = n2;
                break;
            }
            if (iterations == 30)
                return TT_RANGE;
            if (iterations == 10 || iterations == 20) {
                shift = shift + x;
                for (index_t i = l; i <= n; ++i)
                    A(i, i) = A(i, i) - x;
                double s = fabs(A(n, n1)) + fabs(A(n1, n2));
                x = .75 * s;
                y = x;
                w = -.4375 * (s * s);
            }
            ++iterations;
            index_t start = n2;
            double p = 0, q = 0, r = 0, s = 0, z = 0;
            for (; start >= lb; --start) {
                z = A(start, start);
                r = x - z;
                s = y - z;
                p = (r * s - w) / A(start + 1, start) + A(start, start + 1);
                q = A(start + 1, start + 1) - z - r - s;
                r = A(start + 2, start + 1);
                s = fabs(p) + fabs(q) + fabs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if (start == lb)
                    break;
                double u = fabs(A(start, start - 1)) * (fabs(q) + fabs(r));
                double vtest = fabs(p) * (fabs(A(start - 1, start - 1)) + fabs(z) +
                                          fabs(A(start + 1, start + 1)));
                if (u <= DBL_EPSILON * vtest)
                    break;
            }
            A(start + 2, start) = 0;
            for (index_t j = start + 3; j <= n; ++j) {
                A(j, j - 2) = 0;
                A(j, j - 3) = 0;
            }
            for (index_t ka = start; ka <= n1; ++ka) {
                int not_last = ka != n1;
                if (ka == start) {
                    s = copysign(sqrt(p * p + q * q + r * r), p);
                    if (lb != start)
                        A(ka, ka - 1) = -A(ka, ka - 1);
                } else {
                    p = A(ka, ka - 1);
                    q = A(ka + 1, ka - 1);
                    r = not_last ? A(ka + 2, ka - 1) : 0;
                    x = fabs(p) + fabs(q) + fabs(r);
                    if (x == 0)
                        continue;
                    p = p / x;
                    q = q / x;
                    r = r / x;
                    s = copysign(sqrt(p * p + q * q + r * r), p);
                    A(ka, ka - 1) = -s * x;
                }
                p = p + s;
                x = p / s;
                y = q / s;
                z = r / s;
                q = q / p;
                r = r / p;
                for (index_t j = ka; j < m; ++j) {
                    p = A(ka, j) + q * A(ka + 1, j);
                    if (not_last) {
                        p = p + r * A(ka + 2, j);
                        A(ka + 2, j) = A(ka + 2, j) - p * z;
                    }
                    A(ka + 1, j) = A(ka + 1, j) - p * y;
                    A(ka, j) = A(ka, j) - p * x;
                }
                for (index_t i = 0; i <= n && i <= ka + 3; ++i) {
                    p = x * A(i, ka) + y * A(i, ka + 1);
                    if (not_last) {
                        p = p + z * A(i, ka + 2);
                        A(i, ka + 2) = A(i, ka + 2) - p * r;
                    }
                    A(i, ka + 1) = A(i, ka + 1) - p * q;
                    A(i, ka) = A(i, ka) - p;
                }
                for (index_t i = l; i <= k; ++i) {
                    p = x * V(i, ka) + y * V(i, ka + 1);
                    if (not_last) {
                        p = p + z * V(i, ka + 2);
                        V(i, ka + 2) = V(i, ka + 2) - p * r;
                    }
                    V(i, ka + 1) = V(i, ka + 1) - p * q;
                    V(i, ka) = V(i, ka) - p;
                }
            }
        }
    }
    return TT_OK;
}

static void back_substitute(size_t order, double *a, double *v, const double *values,
                            const double *work, index_t l, index_t k, double norm) {
    index_t m = (index_t)order;
    if (norm != 0) {
        for (index_t n = m - 1; n >= 0; --n) {
            index_t n2 = n;
            A(n, n) = 1;
            for (index_t i = n - 1; i >= 0; --i) {
                double w = A(i, i) - values[n];
                if (w == 0)
                    w = DBL_EPSILON * norm;
                double r = A(i, n);
                for (index_t j = n2; j < n; ++j)
                    r = r + A(i, j) * A(j, n);
                A(i, n) = -r / w;
                n2 = i;
            }
        }
        for (index_t i = 0; i < m; ++i)
            if (i < l || i > k)
                for (index_t j = i; j < m; ++j)
                    V(i, j) = A(i, j);
        if (k >= 0)
            for (index_t j = m - 1; j >= l; --j)
                for (index_t i = l; i <= k; ++i) {
                    double z = 0;
                    for (index_t n = l; n <= j && n <= k; ++n)
                        z = z + V(i, n) * A(n, j);
                    V(i, j) = z;
                }
    }
    for (index_t i = l; i <= k; ++i)
        for (index_t j = 0; j < m; ++j)
            V(i, j) = V(i, j) * work[i];
    for (index_t i = l - 1; i >= 0; --i) {
        index_t j = (index_t)work[i];
        if (i != j)
            for (index_t n = 0; n < m; ++n)
                exchange(&V(i, n), &V(j, n));
    }
    for (index_t i = k + 1; i < m; ++i) {
        index_t j = (index_t)work[i];
        if (i != j)
            for (index_t n = 0; n < m; ++n)
                exchange(&V(i, n), &V(j, n));
    }
}

tt_status tt_discrete_eigen(size_t order, const double *input, double *values, double *v,
                            double *work) {
    if (!order || !input || !values || !v || !work || order > PTRDIFF_MAX ||
        order > SIZE_MAX / sizeof(double) / (order + 2))
        return TT_INVALID;
    if (order == 1) {
        values[0] = input[0];
        v[0] = 1;
        return TT_OK;
    }
    if (order == 2) {
        double difference = input[0] - input[3];
        double discriminant = difference * difference + 4 * input[2] * input[1];
        if (discriminant < 0)
            return TT_RANGE;
        double sign = input[0] < input[3] ? -1 : 1;
        values[0] = .5 * (input[0] + input[3] + sign * sqrt(discriminant));
        values[1] = .5 * (input[0] + input[3] - sign * sqrt(discriminant));
        v[0] = v[3] = 1;
        if (input[0] == input[3] && (input[1] == 0 || input[2] == 0)) {
            double norm = fabs(input[0]) + fabs(input[2]) + fabs(input[1]) + fabs(input[3]);
            double w = DBL_EPSILON * norm;
            v[1] = w == 0 ? 0 : input[1] / w;
            v[2] = w == 0 ? 0 : -input[2] / w;
        } else {
            v[1] = input[1] / (values[0] - input[3]);
            v[2] = input[2] / (values[1] - input[0]);
        }
        return TT_OK;
    }
    double *a = work, *scratch = a + order * order;
    memcpy(a, input, order * order * sizeof(double));
    memset(v, 0, order * order * sizeof(double));
    memset(values, 0, order * sizeof(double));
    for (size_t i = 0; i < order; ++i)
        V(i, i) = 1;
    index_t l, k;
    balance(order, a, scratch, &l, &k);
    if (k - 1 >= l + 1)
        hessenberg(order, a, v, scratch, l, k);
    double norm = 0;
    size_t start = 0;
    for (size_t i = 0; i < order; ++i) {
        for (size_t j = start; j < order; ++j)
            norm = norm + fabs(A(i, j));
        start = i;
        if ((index_t)i < l || (index_t)i > k)
            values[i] = A(i, i);
    }
    tt_status status = qr(order, a, v, values, l, k, norm);
    if (!status)
        back_substitute(order, a, v, values, scratch, l, k, norm);
    return status;
}
