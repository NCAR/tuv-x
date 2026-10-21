/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef TT_DECIMAL_INTERNAL_H
#define TT_DECIMAL_INTERNAL_H
/* Convert a decimal number independently of the process locale: optional
 * leading whitespace and sign, digits with an optional '.', and an optional
 * e/E exponent. *end is set past the number, or to text if there is none.
 * The value is the correctly rounded double, as strtod gives in the "C" locale. */
double tt_decimal(const char *text, const char **end);
/* Whitespace as classified in the "C" locale. */
int tt_decimal_space(char c);
#endif
