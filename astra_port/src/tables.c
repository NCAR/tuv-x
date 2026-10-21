/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "table_internal.h"
#include "decimal_internal.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const char *line_end(const char *p) {
    const char *end = strchr(p, '\n');
    return end ? end : p + strlen(p);
}

tt_status tt_text_columns(const tt_dataset *data, size_t columns, int skip_header, size_t *rows,
                          double **values) {
    if (!data || !data->text || !columns || !rows || !values)
        return TT_INVALID;
    const char *start = data->text;
    if (skip_header) {
        start = line_end(start);
        if (*start)
            ++start;
    }
    size_t count = 0;
    for (const char *p = start; *p;) {
        const char *end = line_end(p);
        while (p < end && tt_decimal_space(*p))
            ++p;
        if (p < end && !strchr("#!$%*", *p))
            ++count;
        p = *end ? end + 1 : end;
    }
    if (!count || count > SIZE_MAX / sizeof(double) / columns)
        return TT_INVALID;
    double *result = malloc(count * columns * sizeof(double));
    if (!result)
        return TT_MEMORY;
    size_t row = 0;
    for (const char *p = start; *p;) {
        const char *end = line_end(p);
        while (p < end && tt_decimal_space(*p))
            ++p;
        if (p < end && !strchr("#!$%*", *p)) {
            for (size_t col = 0; col < columns; ++col) {
                while (p < end && (tt_decimal_space(*p) || *p == ','))
                    ++p;
                const char *last = p;
                while (last < end && !tt_decimal_space(*last) && *last != ',')
                    ++last;
                size_t length = (size_t)(last - p);
                char token[128];
                const char *after;
                if (!length || length >= sizeof token) {
                    free(result);
                    return TT_INVALID;
                }
                memcpy(token, p, length);
                token[length] = 0;
                for (size_t i = 0; i < length; ++i)
                    if (token[i] == 'D' || token[i] == 'd')
                        token[i] = 'e';
                double value = tt_decimal(token, &after);
                if (*after || !isfinite(value)) {
                    free(result);
                    return TT_INVALID;
                }
                result[col * count + row] = value;
                p = last;
            }
            ++row;
        }
        p = *end ? end + 1 : end;
    }
    *rows = count;
    *values = result;
    return TT_OK;
}
