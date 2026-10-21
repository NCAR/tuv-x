#ifndef TT_BUILD_INTERNAL_H
#define TT_BUILD_INTERNAL_H
#include "turbotuv_model.h"
#include <stddef.h>

typedef struct {
    tt_data_source source;
    tt_error *error;
    tt_status status;
} tt_builder;
void tt_build_fail(tt_builder *b, tt_config_node node, tt_status status, const char *format, ...);
const char *tt_build_string(tt_builder *b, tt_config_node node, const char *key,
                            const char *fallback);
double tt_build_number(tt_builder *b, tt_config_node node, const char *key, double fallback);
int tt_build_integer(tt_builder *b, tt_config_node node, const char *key);
int tt_build_boolean(tt_builder *b, tt_config_node node, const char *key, int fallback);
int tt_build_has(tt_config_node node, const char *key);
void *tt_build_alloc(tt_builder *b, size_t count, size_t size);
double *tt_build_array(tt_builder *b, tt_config_node node, const char *key, size_t *count);
const tt_dataset *tt_build_dataset(tt_builder *b, tt_config_node node, const char *path);
tt_interpolation tt_build_interpolator(tt_builder *b, tt_config_node node, int *fold);
int tt_build_equal(double a, double b);
/* Check a mapping as TUV-x's config%validate does: every required key is
 * present, and every other key is optional or begins with "__". Lists are
 * NULL-terminated. A missing node is accepted; its readers report it. */
void tt_build_keys(tt_builder *b, tt_config_node node, const char *const *required,
                   const char *const *optional);
#endif
