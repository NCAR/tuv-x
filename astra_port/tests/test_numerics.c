#include "turbotuv.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x);                                \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)

int main(void) {
    double x[] = {0, 1, 2, 3, 4}, y[] = {1, 3, 5, 7, 9}, target[] = {.5, 1.5, 2.5, 3.5};
    double out[5];
    CHECK(tt_interpolate(TT_CONSERVING, 0, 5, x, y, 4, target, out) == TT_OK);
    CHECK(out[0] == 3 && out[1] == 5 && out[2] == 7);
    CHECK(tt_interpolate(TT_LINEAR, 0, 5, x, y, 5, x, out) == TT_OK);
    CHECK(out[0] == 1 && out[3] == 7 && out[4] == 0);
    double bad[] = {0, 1, 1};
    CHECK(tt_interpolate(TT_CONSERVING, 0, 3, bad, y, 4, target, out) == TT_INVALID);
    CHECK(tt_day_of_year(2000, 2, 29) == 60);
    CHECK(tt_day_of_year(2001, 2, 29) == 0);
    CHECK(tt_day_of_year(2024, 12, 31) == 366);
    CHECK(tt_day_of_year(2024, 1, 0) == 0);
    CHECK(isnan(tt_earth_sun_distance(1900, 1, 0)));
    double a[] = {0, -1, -1}, b[] = {2, 2, 2}, c[] = {-1, -1, 0}, r[] = {1, 0, 1}, work[3];
    CHECK(tt_tridiagonal(3, a, b, c, r, out, work) == TT_OK);
    CHECK(fabs(out[0] - 1) < 1e-15 && fabs(out[2] - 1) < 1e-15);
    b[0] = 0;
    CHECK(tt_tridiagonal(3, a, b, c, r, out, work) == TT_SINGULAR);
    double z[] = {0, 1, 10}, paths[6];
    int nid[3];
    CHECK(tt_spherical_geometry(3, z, 0, nid, paths) == TT_OK);
    CHECK(nid[0] == 0 && nid[2] == 2 && paths[4] == 1 && paths[5] == 1);
    CHECK(tt_spherical_geometry(3, z, 180, nid, paths) == TT_OK);
    CHECK(nid[0] == -1 && nid[2] == -1);
    CHECK(tt_dataset_count() == 190);
    for (size_t i = 0; i < tt_dataset_count(); ++i) {
        const tt_dataset *d = tt_dataset_at(i);
        CHECK(d && tt_dataset_find(d->path) == d && strlen(d->sha256) == 64);
        if (d->text)
            CHECK(d->text[d->text_length] == 0);
        for (size_t j = 0; j < d->variable_count; ++j) {
            const tt_data_variable *v = &d->variables[j];
            size_t count = 1;
            for (size_t k = 0; k < v->rank; ++k)
                count *= v->shape[k];
            CHECK(count == v->count);
            for (size_t k = 0; k < count; ++k)
                CHECK(isfinite(v->values[k]));
        }
    }
    CHECK(!tt_dataset_find("missing") && !tt_dataset_at(tt_dataset_count()));
    puts("numerical invariants, error handling, and embedded catalogue passed");
    return 0;
}
