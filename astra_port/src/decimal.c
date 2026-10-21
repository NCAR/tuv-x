/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Configuration and data files always write '.' as the decimal point, but
 * strtod reads the current locale's, so a host that sets a locale with ','
 * would otherwise fail to read them. The number is rewritten as its digits
 * and a decimal exponent (123.456e7 as 123456e4), a form without a decimal
 * point that strtod reads identically in every locale. Both forms denote
 * the same exact value, so the correctly rounded result is unchanged. */
#include "decimal_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int tt_decimal_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}
static int digit(char c) {
    return c >= '0' && c <= '9';
}

double tt_decimal(const char *text, const char **end) {
    const char *p = text;
    *end = text;
    while (tt_decimal_space(*p))
        ++p;
    int negative = *p == '-';
    if (*p == '+' || *p == '-')
        ++p;
    const char *integer = p;
    size_t before = 0, after = 0;
    while (digit(p[before]))
        ++before;
    p += before;
    const char *fraction = p;
    if (*p == '.') {
        fraction = ++p;
        while (digit(p[after]))
            ++after;
        p += after;
    }
    if (!before && !after)
        return 0;
    /* Saturate huge exponents; they overflow or underflow either way. */
    const long long limit = 1000000000000000LL;
    long long exponent = 0;
    if (*p == 'e' || *p == 'E') {
        const char *q = p + 1;
        int minus = *q == '-';
        if (*q == '+' || *q == '-')
            ++q;
        if (digit(*q)) {
            for (; digit(*q); ++q)
                if (exponent < limit)
                    exponent = exponent * 10 + (*q - '0');
            exponent = minus ? -exponent : exponent;
            p = q;
        }
    }
    exponent -= after < (size_t)limit ? (long long)after : limit;
    char small[256], *buffer = small;
    size_t length = before + after + 32;
    if (length > sizeof small && !(buffer = malloc(length)))
        return 0;
    size_t n = 0;
    buffer[n++] = negative ? '-' : '+';
    memcpy(buffer + n, integer, before);
    n += before;
    memcpy(buffer + n, fraction, after);
    n += after;
    snprintf(buffer + n, length - n, "e%lld", exponent);
    double value = strtod(buffer, NULL);
    if (buffer != small)
        free(buffer);
    *end = p;
    return value;
}
