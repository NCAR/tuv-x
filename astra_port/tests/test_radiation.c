#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef TT_HAVE_REFERENCE
void ref_delta(int, int, int, const double *, double, const double *, const double *,
               const double *, const double *, const int *, double *, double *);
#endif
static size_t checked;
static double worst;
static int scenario;
static void same(const char *name, size_t i, double a, double b) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(1.0, fabs(b)));
    ++checked;
    if (!isfinite(a) || !isfinite(b) || error > 64) {
        fprintf(stderr, "case %d %s[%zu]: %.17g != %.17g (%.3g eps)\n", scenario, name, i, a, b,
                error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static void ok(tt_status s) {
    if (s) {
        fprintf(stderr, "case %d: %s\n", scenario, tt_status_string(s));
        exit(1);
    }
}
int main(void) {
    enum { NL = 81, NW = 7, NR = 3 };
    double z[NL + 1], surface[NW], tau[NR * NL * NW], ssa[NR * NL * NW], g[NR * NL * NW];
    double mix[3 * NL * NW], mix_ref[3 * NL * NW], out[6 * (NL + 1) * NW],
        expected[6 * (NL + 1) * NW];
    const double angles[] = {0, 30, 60, 85, 89.9, 90, 90.01, 93, 99, 110, 180};
    const int sizes[] = {1, 2, 5, 16, 81}, air[] = {0, 1, 0};
    for (size_t ni = 0; ni < 5; ++ni)
        for (size_t ai = 0; ai < 11; ++ai)
            for (int kind = 0; kind < 5; ++kind) {
                ++scenario;
                size_t n = (size_t)sizes[ni], count = n * NW, field = (n + 1) * NW;
                for (size_t j = 0; j <= n; ++j)
                    z[j] = 2.3 + 120.0 * (double)j / (double)n;
                for (size_t w = 0; w < NW; ++w)
                    surface[w] = (double)w / (NW - 1);
                for (size_t r = 0; r < NR; ++r)
                    for (size_t w = 0; w < NW; ++w)
                        for (size_t j = 0; j < n; ++j) {
                            size_t k = r * count + w * n + j;
                            tau[k] = exp(-(z[j] - 2.3) / 8.0) * (1 + w) * (.1 + r) / (double)n;
                            if (kind == 1)
                                tau[k] *= 1e-10;
                            if (kind == 2)
                                tau[k] *= 1e4;
                            if (kind == 4)
                                tau[k] = 0;
                            ssa[k] = kind == 3 ? 1.0 : (double)w / (NW - 1);
                            g[k] = r == 1 ? 0.0 : (-0.2 + 1.1 * (double)w / (NW - 1));
                        }
                tt_radiator radiators[NR];
                for (size_t r = 0; r < NR; ++r)
                    radiators[r] =
                        (tt_radiator){tau + r * count, ssa + r * count, g + r * count, air[r]};
#ifdef TT_HAVE_REFERENCE
                ref_delta((int)n, NW, NR, z, angles[ai], surface, tau, ssa, g, air, mix_ref,
                          expected);
#endif
                ok(tt_mix_radiators(NR, radiators, n, NW, 1, mix, mix + count, mix + 2 * count));
#ifdef TT_HAVE_REFERENCE
                for (size_t j = 0; j < 3 * count; ++j)
                    same("mixed optics", j, mix[j], mix_ref[j]);
#endif
                tt_atmosphere atm = {n, NW, 1, z, mix, mix + count, mix + 2 * count};
                tt_radiation radiation = {out,
                                          out + field,
                                          out + 2 * field,
                                          out + 3 * field,
                                          out + 4 * field,
                                          out + 5 * field};
                ok(tt_delta_eddington(&atm, angles[ai], surface, &radiation));
                for (size_t j = 0; j < 6 * field; ++j)
                    if (!isfinite(out[j]))
                        return 1;
#ifdef TT_HAVE_REFERENCE
                for (size_t j = 0; j < 6 * field; ++j)
                    same("radiation", j, out[j], expected[j]);
#endif
                if (angles[ai] < 90) {
                    for (size_t w = 0; w < NW; ++w) {
                        same("top direct beam", w, radiation.direct_flux[w * (n + 1) + n], 1.0);
                        same("surface reflection", w, radiation.up_irradiance[w * (n + 1)],
                             surface[w] * (radiation.direct_irradiance[w * (n + 1)] +
                                           radiation.down_irradiance[w * (n + 1)]));
                    }
                }
            }
    printf("%d atmospheres, %zu comparisons; worst %.3g epsilon-scaled error\n", scenario, checked,
           worst);
    return 0;
}
