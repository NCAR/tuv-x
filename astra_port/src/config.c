/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "turbotuv_config.h"
#include "decimal_internal.h"
#include "namespace.h"
#include "yaml.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct tt_config {
    yaml_document_t yaml;
    unsigned char *null_nodes;
};

static int null_word(const char *s) {
    return !*s || !strcmp(s, "~") || !strcmp(s, "null") || !strcmp(s, "Null") || !strcmp(s, "NULL");
}

/* LibYAML's composed nodes omit whether a scalar tag was explicit. For
 * documents containing null-like scalars, retain that distinction from the
 * event stream, so !!str null and untagged null remain different values. */
static int classify_nulls(tt_config *config, const char *text, size_t length, size_t count) {
    int needed = 0;
    for (size_t i = 0; i < count; ++i) {
        const yaml_node_t *node = config->yaml.nodes.start + i;
        if (node->type == YAML_SCALAR_NODE && null_word((const char *)node->data.scalar.value))
            needed = 1;
    }
    if (!needed)
        return 1;
    yaml_parser_t parser;
    if (!yaml_parser_initialize(&parser))
        return 0;
    yaml_parser_set_input_string(&parser, (const unsigned char *)text, length);
    size_t index = 0;
    int success = 1, done = 0;
    while (!done) {
        yaml_event_t event;
        if (!yaml_parser_parse(&parser, &event)) {
            success = 0;
            break;
        }
        if (event.type == YAML_SCALAR_EVENT || event.type == YAML_SEQUENCE_START_EVENT ||
            event.type == YAML_MAPPING_START_EVENT)
            ++index;
        if (event.type == YAML_SCALAR_EVENT && index <= count) {
            const char *tag = (const char *)event.data.scalar.tag;
            int untyped = event.data.scalar.plain_implicit || (tag && !strcmp(tag, "!"));
            config->null_nodes[index] =
                (unsigned char)(null_word((const char *)event.data.scalar.value) &&
                                (untyped || (tag && !strcmp(tag, YAML_NULL_TAG))));
        }
        done = event.type == YAML_STREAM_END_EVENT;
        yaml_event_delete(&event);
    }
    yaml_parser_delete(&parser);
    return success;
}

static tt_status fail(tt_error *error, tt_status status, size_t line, size_t column,
                      const char *message) {
    if (error) {
        error->status = status;
        error->line = line;
        error->column = column;
        snprintf(error->message, sizeof error->message, "%s", message);
    }
    return status;
}

static const yaml_node_t *node_at(tt_config_node node) {
    if (!node.document || !node.index)
        return NULL;
    const yaml_document_t *doc = &node.document->yaml;
    size_t count = (size_t)(doc->nodes.top - doc->nodes.start);
    return node.index <= count ? doc->nodes.start + node.index - 1 : NULL;
}

typedef struct {
    const char *key;
    size_t index;
} key_entry;
static int compare_keys(const void *a, const void *b) {
    const key_entry *x = a, *y = b;
    int c = strcmp(x->key, y->key);
    return c ? c : (x->index > y->index) - (x->index < y->index);
}
/* Index of the first key, in source order, that repeats an earlier key of
 * this mapping (SIZE_MAX if none). Only the keys before the first invalid
 * key are considered, because validation stops there. Large mappings are
 * sorted rather than compared pairwise. Returns nonzero if out of memory. */
static int first_duplicate(const tt_config *doc, const yaml_node_t *node, size_t *out) {
    size_t n = 0, total = (size_t)(node->data.mapping.pairs.top - node->data.mapping.pairs.start);
    const char **keys = NULL;
    key_entry *entries = NULL;
    for (; n < total; ++n) {
        const yaml_node_t *key =
            node_at((tt_config_node){doc, (size_t)node->data.mapping.pairs.start[n].key});
        if (!key || key->type != YAML_SCALAR_NODE ||
            memchr(key->data.scalar.value, 0, key->data.scalar.length))
            break;
    }
    *out = SIZE_MAX;
    if (n < 2)
        return 0;
    keys = malloc(n * sizeof *keys);
    if (!keys)
        return 1;
    for (size_t i = 0; i < n; ++i)
        keys[i] = (const char *)node_at(
                      (tt_config_node){doc, (size_t)node->data.mapping.pairs.start[i].key})
                      ->data.scalar.value;
    if (n <= 16) {
        for (size_t i = 1; i < n && *out == SIZE_MAX; ++i)
            for (size_t j = 0; j < i; ++j)
                if (!strcmp(keys[i], keys[j])) {
                    *out = i;
                    break;
                }
    } else {
        entries = malloc(n * sizeof *entries);
        if (!entries) {
            free(keys);
            return 1;
        }
        for (size_t i = 0; i < n; ++i)
            entries[i] = (key_entry){keys[i], i};
        qsort(entries, n, sizeof *entries, compare_keys);
        /* In each run of equal keys, the second entry is the first repeat. */
        for (size_t i = 1; i < n; ++i)
            if (!strcmp(entries[i].key, entries[i - 1].key) &&
                (i < 2 || strcmp(entries[i - 1].key, entries[i - 2].key)) &&
                entries[i].index < *out)
                *out = entries[i].index;
    }
    free(entries);
    free(keys);
    return 0;
}

static tt_status validate_node(tt_config *doc, size_t index, size_t depth, unsigned char *visited,
                               tt_error *error) {
    const yaml_node_t *node = node_at((tt_config_node){doc, index});
    if (!node)
        return fail(error, TT_INVALID, 0, 0, "Invalid configuration node");
    if (depth > 128 || visited[index] == 1)
        return fail(error, TT_INVALID, node->start_mark.line + 1, node->start_mark.column + 1,
                    depth > 128 ? "Configuration nesting exceeds 128"
                                : "Cyclic configuration alias");
    if (visited[index] >= 2)
        return depth + visited[index] - 2 <= 128
                   ? TT_OK
                   : fail(error, TT_INVALID, node->start_mark.line + 1, node->start_mark.column + 1,
                          "Configuration nesting exceeds 128");
    visited[index] = 1;
    unsigned char height = 0;
    if (node->type == YAML_MAPPING_NODE) {
        size_t duplicate = SIZE_MAX;
        if (first_duplicate(doc, node, &duplicate))
            return fail(error, TT_MEMORY, 0, 0, "Cannot allocate configuration validation");
        for (const yaml_node_pair_t *pair = node->data.mapping.pairs.start;
             pair < node->data.mapping.pairs.top; ++pair) {
            const yaml_node_t *key = node_at((tt_config_node){doc, (size_t)pair->key});
            if (!key || key->type != YAML_SCALAR_NODE)
                return fail(error, TT_INVALID, node->start_mark.line + 1,
                            node->start_mark.column + 1,
                            "Configuration keys must be scalar strings");
            if (memchr(key->data.scalar.value, 0, key->data.scalar.length))
                return fail(error, TT_INVALID, key->start_mark.line + 1, key->start_mark.column + 1,
                            "NUL characters are not supported in configuration keys");
            if ((size_t)(pair - node->data.mapping.pairs.start) == duplicate)
                return fail(error, TT_INVALID, key->start_mark.line + 1, key->start_mark.column + 1,
                            "Duplicate configuration key");
            tt_status status = validate_node(doc, (size_t)pair->value, depth + 1, visited, error);
            if (status)
                return status;
            unsigned char child_height = (unsigned char)(visited[pair->value] - 1);
            if (child_height > height)
                height = child_height;
        }
    } else if (node->type == YAML_SEQUENCE_NODE) {
        for (const yaml_node_item_t *item = node->data.sequence.items.start;
             item < node->data.sequence.items.top; ++item) {
            tt_status status = validate_node(doc, (size_t)*item, depth + 1, visited, error);
            if (status)
                return status;
            unsigned char child_height = (unsigned char)(visited[*item] - 1);
            if (child_height > height)
                height = child_height;
        }
    } else if (node->type == YAML_SCALAR_NODE) {
        if (memchr(node->data.scalar.value, 0, node->data.scalar.length))
            return fail(error, TT_INVALID, node->start_mark.line + 1, node->start_mark.column + 1,
                        "NUL characters are not supported in configuration strings");
    }
    visited[index] = (unsigned char)(height + 2);
    return TT_OK;
}

void tt_config_destroy(tt_config *config) {
    if (config) {
        yaml_document_delete(&config->yaml);
        free(config->null_nodes);
        free(config);
    }
}

tt_status tt_config_parse(const char *text, size_t length, tt_config **out, tt_error *error) {
    if (error)
        *error = (tt_error){0};
    if (!text || !out)
        return fail(error, TT_INVALID, 0, 0, "Missing configuration text or output");
    tt_config *config = calloc(1, sizeof(*config));
    yaml_parser_t parser;
    if (!config || !yaml_parser_initialize(&parser)) {
        free(config);
        return fail(error, TT_MEMORY, 0, 0, "Cannot allocate configuration parser");
    }
    yaml_parser_set_input_string(&parser, (const unsigned char *)text, length);
    tt_status status = TT_OK;
    if (!yaml_parser_load(&parser, &config->yaml)) {
        status = fail(error, parser.error == YAML_MEMORY_ERROR ? TT_MEMORY : TT_INVALID,
                      parser.problem_mark.line + 1, parser.problem_mark.column + 1,
                      parser.problem ? parser.problem : "Cannot parse configuration");
        /* LibYAML releases the partial document on a parse error. */
        free(config);
        config = NULL;
    } else {
        yaml_document_t extra;
        if (!yaml_parser_load(&parser, &extra))
            status = fail(error, parser.error == YAML_MEMORY_ERROR ? TT_MEMORY : TT_INVALID,
                          parser.problem_mark.line + 1, parser.problem_mark.column + 1,
                          parser.problem ? parser.problem : "Trailing configuration content");
        else {
            if (yaml_document_get_root_node(&extra))
                status = fail(error, TT_INVALID, extra.start_mark.line + 1,
                              extra.start_mark.column + 1, "Expected one configuration document");
            yaml_document_delete(&extra);
        }
    }
    yaml_parser_delete(&parser);
    if (!status) {
        size_t count = config->yaml.nodes.start
                           ? (size_t)(config->yaml.nodes.top - config->yaml.nodes.start)
                           : 0;
        if (!count)
            status = fail(error, TT_INVALID, 0, 0, "Empty configuration document");
        else {
            unsigned char *visited = calloc(count + 1, 1);
            if (!visited)
                status = fail(error, TT_MEMORY, 0, 0, "Cannot allocate configuration validation");
            else {
                status = validate_node(config, 1, 0, visited, error);
                config->null_nodes = visited;
                memset(visited, 0, count + 1);
                if (!status && !classify_nulls(config, text, length, count))
                    status = fail(error, TT_MEMORY, 0, 0, "Cannot classify configuration scalars");
            }
        }
    }
    if (status)
        tt_config_destroy(config);
    else
        *out = config;
    return status;
}

tt_status tt_config_read_file(const char *path, tt_config **out, tt_error *error) {
    if (error)
        *error = (tt_error){0};
    if (!path || !out)
        return fail(error, TT_INVALID, 0, 0, "Missing configuration path or output");
    FILE *file = fopen(path, "rb");
    if (!file)
        return fail(error, TT_NOT_FOUND, 0, 0, "Cannot open configuration file");
    size_t size = 0, capacity = 4096;
    char *text = malloc(capacity);
    tt_status status = text ? TT_OK : TT_MEMORY;
    while (!status) {
        size_t got = fread(text + size, 1, capacity - size, file);
        size += got;
        if (size < capacity) {
            if (ferror(file))
                status = TT_INVALID;
            break;
        }
        if (capacity > SIZE_MAX / 2) {
            status = TT_RANGE;
            break;
        }
        capacity *= 2;
        char *new_text = realloc(text, capacity);
        if (!new_text) {
            status = TT_MEMORY;
            break;
        }
        text = new_text;
    }
    fclose(file);
    if (status)
        fail(error, status, 0, 0, "Cannot read configuration file");
    else
        status = tt_config_parse(text, size, out, error);
    free(text);
    return status;
}

tt_config_node tt_config_root(const tt_config *config) {
    return (tt_config_node){config, config ? 1 : 0};
}

tt_config_kind tt_config_type(tt_config_node value) {
    const yaml_node_t *node = node_at(value);
    if (!node)
        return TT_CONFIG_MISSING;
    if (node->type == YAML_MAPPING_NODE)
        return TT_CONFIG_MAPPING;
    if (node->type == YAML_SEQUENCE_NODE)
        return TT_CONFIG_SEQUENCE;
    if (node->type != YAML_SCALAR_NODE)
        return TT_CONFIG_NULL;
    if (value.document->null_nodes[value.index])
        return TT_CONFIG_NULL;
    return TT_CONFIG_SCALAR;
}

size_t tt_config_size(tt_config_node value) {
    const yaml_node_t *node = node_at(value);
    if (!node)
        return 0;
    if (node->type == YAML_MAPPING_NODE)
        return (size_t)(node->data.mapping.pairs.top - node->data.mapping.pairs.start);
    if (node->type == YAML_SEQUENCE_NODE)
        return (size_t)(node->data.sequence.items.top - node->data.sequence.items.start);
    return 0;
}

const char *tt_config_string(tt_config_node value) {
    const yaml_node_t *node = node_at(value);
    return node && tt_config_type(value) == TT_CONFIG_SCALAR ? (const char *)node->data.scalar.value
                                                             : NULL;
}

tt_config_node tt_config_get(tt_config_node value, const char *key) {
    const yaml_node_t *node = node_at(value);
    if (node && key && node->type == YAML_MAPPING_NODE)
        for (const yaml_node_pair_t *p = node->data.mapping.pairs.start;
             p < node->data.mapping.pairs.top; ++p) {
            const yaml_node_t *k = node_at((tt_config_node){value.document, (size_t)p->key});
            if (!strcmp(key, (const char *)k->data.scalar.value))
                return (tt_config_node){value.document, (size_t)p->value};
        }
    return (tt_config_node){0};
}

tt_config_node tt_config_at(tt_config_node value, size_t index) {
    const yaml_node_t *node = node_at(value);
    if (!node || index >= tt_config_size(value))
        return (tt_config_node){0};
    if (node->type == YAML_SEQUENCE_NODE)
        return (tt_config_node){value.document, (size_t)node->data.sequence.items.start[index]};
    return (tt_config_node){value.document, (size_t)node->data.mapping.pairs.start[index].value};
}

const char *tt_config_key(tt_config_node value, size_t index) {
    const yaml_node_t *node = node_at(value);
    if (!node || node->type != YAML_MAPPING_NODE || index >= tt_config_size(value))
        return NULL;
    const yaml_node_t *key = node_at(
        (tt_config_node){value.document, (size_t)node->data.mapping.pairs.start[index].key});
    return (const char *)key->data.scalar.value;
}

tt_status tt_config_number(tt_config_node node, double *out) {
    const char *text = tt_config_string(node);
    if (!text || !*text || !out)
        return TT_INVALID;
    const char *end;
    double value = tt_decimal(text, &end);
    if (end == text)
        return TT_INVALID;
    while (tt_decimal_space(*end))
        ++end;
    if (end == text || *end || !isfinite(value))
        return TT_INVALID;
    *out = value;
    return TT_OK;
}

tt_status tt_config_integer(tt_config_node node, int *out) {
    const char *text = tt_config_string(node);
    if (!text || !*text || !out)
        return TT_INVALID;
    errno = 0;
    char *end;
    long value = strtol(text, &end, 0);
    if (end == text)
        return TT_INVALID;
    while (tt_decimal_space(*end))
        ++end;
    if (errno == ERANGE || end == text || *end || value < INT_MIN || value > INT_MAX)
        return TT_INVALID;
    *out = (int)value;
    return TT_OK;
}

tt_status tt_config_boolean(tt_config_node node, int *out) {
    const char *text = tt_config_string(node);
    if (!text || !out || strlen(text) > 5)
        return TT_INVALID;
    char lower[6];
    int all_lower = 1, all_upper = 1, capitalized = isupper((unsigned char)text[0]);
    for (size_t j = 0; text[j]; ++j) {
        all_lower = all_lower && islower((unsigned char)text[j]);
        all_upper = all_upper && isupper((unsigned char)text[j]);
        if (j)
            capitalized = capitalized && islower((unsigned char)text[j]);
    }
    if (!all_lower && !all_upper && !capitalized)
        return TT_INVALID;
    size_t i = 0;
    do {
        lower[i] = (char)tolower((unsigned char)text[i]);
    } while (text[i++]);
    if (!strcmp(lower, "true") || !strcmp(lower, "yes") || !strcmp(lower, "y") ||
        !strcmp(lower, "on"))
        *out = 1;
    else if (!strcmp(lower, "false") || !strcmp(lower, "no") || !strcmp(lower, "n") ||
             !strcmp(lower, "off"))
        *out = 0;
    else
        return TT_INVALID;
    return TT_OK;
}

void tt_config_location(tt_config_node value, size_t *line, size_t *column) {
    const yaml_node_t *node = node_at(value);
    if (line)
        *line = node ? node->start_mark.line + 1 : 0;
    if (column)
        *column = node ? node->start_mark.column + 1 : 0;
}
