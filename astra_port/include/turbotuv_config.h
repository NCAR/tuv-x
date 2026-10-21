/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef TURBOTUV_CONFIG_H
#define TURBOTUV_CONFIG_H
#include "turbotuv.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct tt_config tt_config;
typedef struct {
    const tt_config *document;
    size_t index;
} tt_config_node;
typedef enum {
    TT_CONFIG_MISSING,
    TT_CONFIG_NULL,
    TT_CONFIG_SCALAR,
    TT_CONFIG_SEQUENCE,
    TT_CONFIG_MAPPING
} tt_config_kind;
typedef struct {
    tt_status status;
    size_t line, column; /* one-based; zero when unavailable */
    char message[256];
} tt_error;

/* Parse one JSON or YAML document. The document owns all nodes/strings;
 * numbers are read with '.' decimal points whatever the process locale.
 * borrowed node handles remain valid until destroy. Output is unchanged on
 * failure. Cycles, duplicate/non-scalar mapping keys and depth >128 are
 * rejected; non-cyclic YAML aliases retain shared immutable nodes. */
tt_status tt_config_parse(const char *text, size_t length, tt_config **out, tt_error *error);
tt_status tt_config_read_file(const char *path, tt_config **out, tt_error *error);
void tt_config_destroy(tt_config *config);
tt_config_node tt_config_root(const tt_config *config);
tt_config_kind tt_config_type(tt_config_node node);
size_t tt_config_size(tt_config_node node);
tt_config_node tt_config_get(tt_config_node node, const char *key);
/* Sequence item or mapping value in source order; missing if out of range. */
tt_config_node tt_config_at(tt_config_node node, size_t index);
const char *tt_config_key(tt_config_node node, size_t index);
const char *tt_config_string(tt_config_node node);
tt_status tt_config_number(tt_config_node node, double *out);
tt_status tt_config_integer(tt_config_node node, int *out);
tt_status tt_config_boolean(tt_config_node node, int *out);
void tt_config_location(tt_config_node node, size_t *line, size_t *column);

#ifdef __cplusplus
}
#endif
#endif
