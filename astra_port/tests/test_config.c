#include "turbotuv_config.h"
#include "fixtures.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_HAVE_REFERENCE
void ref_configuration(int, const char *, int, double *, int *, int *, double *, char *);
#endif
static size_t checked, documents;
static void require(int test, const char *message) {
    ++checked;
    if (!test) {
        fprintf(stderr, "configuration check: %s\n", message);
        exit(1);
    }
}
static tt_config *parse(const char *text, size_t size) {
    tt_config *config = NULL;
    tt_error error;
    tt_status status = tt_config_parse(text, size, &config, &error);
    if (status) {
        fprintf(stderr, "configuration %zu (%zu:%zu): %s\n", documents, error.line, error.column,
                error.message);
        exit(1);
    }
    ++documents;
    return config;
}
static size_t visit(tt_config_node node) {
    tt_config_kind kind = tt_config_type(node);
    require(kind != TT_CONFIG_MISSING, "present node");
    size_t result = 1, count = tt_config_size(node);
    for (size_t i = 0; i < count; ++i) {
        tt_config_node child = tt_config_at(node, i);
        if (kind == TT_CONFIG_MAPPING) {
            const char *key = tt_config_key(node, i);
            require(key != NULL, "mapping key");
            require(tt_config_get(node, key).index == child.index, "mapping lookup");
        }
        result += visit(child);
    }
    require(tt_config_type(tt_config_at(node, count)) == TT_CONFIG_MISSING, "range checked access");
    return result;
}
int main(void) {
    for (size_t i = 0; i < tt_test_dataset_count() + tt_dataset_count(); ++i) {
        const tt_dataset *d = i < tt_test_dataset_count()
                                  ? tt_test_dataset_at(i)
                                  : tt_dataset_at(i - tt_test_dataset_count());
        if (!d->text || (!strstr(d->path, ".json") && !strstr(d->path, ".yml")))
            continue;
        const char *text = d->text;
        char *wrapped = NULL;
        size_t length = d->text_length;
        if (!strncmp(d->path, "data/", 5)) {
            /* These two upstream JSON files are fragments inserted by the
             * original regression scripts into an enclosing object. */
            wrapped = malloc(length + 3);
            wrapped[0] = '{';
            memcpy(wrapped + 1, text, length);
            wrapped[length + 1] = '}';
            wrapped[length + 2] = 0;
            text = wrapped;
            length += 2;
        }
        tt_config *config = parse(text, length);
        free(wrapped);
        require(visit(tt_config_root(config)) > 0, d->path);
        tt_config_destroy(config);
    }
    for (int k = 0; k < 200; ++k) {
        char input[1024];
        double wanted = exp((double)(k - 100) * 5.5);
        int integer = k * 12345 - 98572, flag = k % 2;
        snprintf(input, sizeof input,
                 "number: %.17g\ninteger: %d\nflag: %s\nvalues: [1e-300, -7.425e28, %.17g, "
                 "-0.0]\nstring: \"a\\nb\\t\\u03bb\"\n",
                 wanted, integer, flag ? "True" : "OFF", wanted);
        tt_config *config = parse(input, strlen(input));
        tt_config_node root = tt_config_root(config);
        double number, values[4];
        int iv, bv;
        require(tt_config_number(tt_config_get(root, "number"), &number) == TT_OK,
                "number conversion");
        require(tt_config_integer(tt_config_get(root, "integer"), &iv) == TT_OK,
                "integer conversion");
        require(tt_config_boolean(tt_config_get(root, "flag"), &bv) == TT_OK, "boolean conversion");
        for (size_t i = 0; i < 4; ++i)
            require(tt_config_number(tt_config_at(tt_config_get(root, "values"), i), values + i) ==
                        TT_OK,
                    "array value");
        const char *string = tt_config_string(tt_config_get(root, "string"));
        require(number == wanted && iv == integer && bv == flag, "typed values");
        require(!strcmp(string, "a\nb\tλ"), "Unicode/escape string");
#ifdef TT_HAVE_REFERENCE
        double nr, vr[4];
        int ir, br;
        char sr[128];
        ref_configuration((int)strlen(input), input, 4, &nr, &ir, &br, vr, sr);
        require(number == nr && iv == ir && bv == br && !strcmp(string, sr),
                "Fortran scalar parity");
        for (size_t i = 0; i < 4; ++i)
            require(values[i] == vr[i] && signbit(values[i]) == signbit(vr[i]),
                    "Fortran array parity");
#endif
        tt_config_destroy(config);
    }
    const char *equivalent[] = {"{\"x\": [2.5, {\"unicode\": \"\\u03bb\"}], \"nil\": null}",
                                "x: [2.5, {unicode: λ}]\nnil: null\n"};
    for (size_t i = 0; i < 2; ++i) {
        tt_config *c = parse(equivalent[i], strlen(equivalent[i]));
        tt_config_node root = tt_config_root(c);
        double number;
        require(tt_config_number(tt_config_at(tt_config_get(root, "x"), 0), &number) == TT_OK &&
                    number == 2.5,
                "flow sequence");
        require(!strcmp(tt_config_string(
                            tt_config_get(tt_config_at(tt_config_get(root, "x"), 1), "unicode")),
                        "λ"),
                "nested Unicode");
        require(tt_config_type(tt_config_get(root, "nil")) == TT_CONFIG_NULL, "null node");
        tt_config_destroy(c);
    }
    const char *bad[] = {"",
                         "[1,2",
                         "{a: 1, a: 2}",
                         "a: &a [*a]",
                         "[a]\n---\n[b]",
                         "? [a,b]\n: c",
                         "x: \"\\u0000\"",
                         "\"\\u0000\": value"};
    for (size_t i = 0; i < sizeof bad / sizeof *bad; ++i) {
        tt_config *c = NULL;
        tt_error e;
        require(tt_config_parse(bad[i], strlen(bad[i]), &c, &e) == TT_INVALID && !c && e.message[0],
                "parse error");
    }
    const char *aliases = "base: &base [1, 2, 3]\ncopy: *base\n";
    tt_config *c = parse(aliases, strlen(aliases));
    tt_config_node root = tt_config_root(c);
    require(tt_config_get(root, "base").index == tt_config_get(root, "copy").index,
            "shared immutable aliases");
    tt_config_destroy(c);
    const char *invalid = "numeric: '  '\nboolean: tRuE\ninteger: 2147483648\n";
    c = parse(invalid, strlen(invalid));
    root = tt_config_root(c);
    double d;
    int iv;
    require(tt_config_number(tt_config_get(root, "numeric"), &d) == TT_INVALID,
            "empty numeric scalar");
    require(tt_config_boolean(tt_config_get(root, "boolean"), &iv) == TT_INVALID,
            "invalid boolean case");
    require(tt_config_integer(tt_config_get(root, "integer"), &iv) == TT_INVALID,
            "integer overflow");
    tt_config_destroy(c);
    const char *tags = "literal: !!str null\nplain: null\nquoted: 'null'\nempty: ''\n";
    c = parse(tags, strlen(tags));
    root = tt_config_root(c);
    require(!strcmp(tt_config_string(tt_config_get(root, "literal")), "null"),
            "explicit string tag");
    require(tt_config_type(tt_config_get(root, "plain")) == TT_CONFIG_NULL, "implicit null");
    require(!strcmp(tt_config_string(tt_config_get(root, "quoted")), "null"), "quoted null string");
    require(!strcmp(tt_config_string(tt_config_get(root, "empty")), ""), "empty string");
    tt_config_destroy(c);
    char nested[600];
    memset(nested, '[', 150);
    nested[150] = '0';
    memset(nested + 151, ']', 150);
    nested[301] = 0;
    c = NULL;
    tt_error error;
    require(tt_config_parse(nested, 301, &c, &error) == TT_INVALID && !c, "depth bound");
#ifdef TT_HAVE_REFERENCE
    const char *file_name = "test-config-reference.tmp";
#else
    const char *file_name = "test-config-invariants.tmp";
#endif
    FILE *file = fopen(file_name, "wb");
    require(file != NULL, "temporary config file");
    require(fwrite(tags, 1, strlen(tags), file) == strlen(tags), "write temporary config");
    require(fclose(file) == 0, "close temporary config");
    require(tt_config_read_file(file_name, &c, &error) == TT_OK, "read configuration file");
    tt_config_destroy(c);
    require(remove(file_name) == 0, "remove temporary config");
    c = NULL;
    require(tt_config_read_file(file_name, &c, &error) == TT_NOT_FOUND && !c,
            "missing configuration file");
    /* Deliberately exercise truncated and malformed parser inputs under sanitizers. */
    const char *fragment = "{x: [1, 2, {s: \"a\\nb\\t\\u03bb\"}], y: false}";
    for (size_t i = 0; i < strlen(fragment); ++i) {
        c = NULL;
        tt_error e;
        tt_status status = tt_config_parse(fragment, i, &c, &e);
        require(status == TT_OK || status == TT_INVALID, "truncated input result");
        tt_config_destroy(c);
    }
    printf("Configuration: %zu documents, %zu checks\n", documents, checked);
    return 0;
}
