#include "turbotuv.h"
#include <string.h>

#include "generated/data.inc"

size_t tt_dataset_count(void) {
    return sizeof tt_catalogue / sizeof *tt_catalogue;
}

const tt_dataset *tt_dataset_at(size_t index) {
    return index < tt_dataset_count() ? &tt_catalogue[index] : NULL;
}

const tt_dataset *tt_dataset_find(const char *path) {
    if (!path)
        return NULL;
    while (path[0] == '.' && path[1] == '/')
        path += 2;
    size_t lo = 0, hi = tt_dataset_count();
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp = strcmp(path, tt_catalogue[mid].path);
        if (!cmp)
            return &tt_catalogue[mid];
        if (cmp < 0)
            hi = mid;
        else
            lo = mid + 1;
    }
    return NULL;
}

const tt_data_variable *tt_dataset_variable(const tt_dataset *data, const char *name) {
    if (!data || !name)
        return NULL;
    for (size_t i = 0; i < data->variable_count; ++i)
        if (!strcmp(data->variables[i].name, name))
            return &data->variables[i];
    return NULL;
}
