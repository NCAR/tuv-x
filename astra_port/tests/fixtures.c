#include "fixtures.h"
#include <string.h>
#include "generated/data.inc"
size_t tt_test_dataset_count(void) {
    return sizeof tt_catalogue / sizeof *tt_catalogue;
}
const tt_dataset *tt_test_dataset_at(size_t index) {
    return index < tt_test_dataset_count() ? &tt_catalogue[index] : NULL;
}
const tt_dataset *tt_test_dataset(const char *path) {
    for (size_t i = 0; i < sizeof tt_catalogue / sizeof *tt_catalogue; ++i)
        if (!strcmp(path, tt_catalogue[i].path))
            return &tt_catalogue[i];
    return tt_dataset_find(path);
}
