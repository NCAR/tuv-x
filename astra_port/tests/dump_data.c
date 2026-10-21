#include "turbotuv.h"
#include <stdio.h>
#ifdef TT_TEST_FIXTURES
#include "fixtures.h"
#define tt_dataset_count tt_test_dataset_count
#define tt_dataset_at tt_test_dataset_at
#endif

/* Binary output for the independent extraction audit in verify_data.py. */
int main(void) {
    for (size_t i = 0; i < tt_dataset_count(); ++i) {
        const tt_dataset *d = tt_dataset_at(i);
        if (d->text && fwrite(d->text, 1, d->text_length, stdout) != d->text_length)
            return 1;
        for (size_t j = 0; j < d->variable_count; ++j) {
            const tt_data_variable *v = &d->variables[j];
            if (fwrite(v->values, sizeof(double), v->count, stdout) != v->count)
                return 1;
        }
    }
    return ferror(stdout) != 0;
}
