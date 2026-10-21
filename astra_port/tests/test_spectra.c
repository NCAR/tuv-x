#include "turbotuv.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef TT_HAVE_REFERENCE
void ref_cross(int, int, int, const double *, const double *, const double *, int,
               const tt_spectral_table *, double *);
#endif
static size_t comparisons;
static double worst;
static void check(const char *name, size_t i, double a, double b) {
    double error = fabs(a - b) / (DBL_EPSILON * fmax(fabs(a), fabs(b)) + DBL_MIN);
    ++comparisons;
    if (!isfinite(a) || !isfinite(b) || error > 32) {
        fprintf(stderr, "%s[%zu]: %.17g != %.17g (%.3g eps)\n", name, i, a, b, error);
        exit(1);
    }
    worst = fmax(worst, error);
}
static int order(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
int main(void) {
    enum { NW = 4096, NZ = 25, NP = 4 };
    double edges[NW + 1], mids[NW], data[NP][4 * NW], out[NZ * NW], expected[NZ * NW];
    const double temp[NZ] = {150, 173, 180, 194, 200, 203, 210, 218, 233, 235, 240, 250, 273,
                             295, 296, 298, 300, 310, 320, 350, 400, 450, 500, 550, 600};
    const double knots[] = {173, 184.9, 185, 190, 194, 195, 196.078, 240, 250, 260, 270,
                            284, 290,   306, 330, 335, 340, 342.5,   345, 350, 550};
    size_t nw = 0;
    for (double w = 140; w < 850; w += .5)
        mids[nw++] = w;
    for (size_t i = 0; i < sizeof knots / sizeof *knots; ++i) {
        mids[nw++] = nextafter(knots[i], -INFINITY);
        mids[nw++] = knots[i];
        mids[nw++] = nextafter(knots[i], INFINITY);
    }
    qsort(mids, nw, sizeof(double), order);
    size_t unique = 1;
    for (size_t i = 1; i < nw; ++i)
        if (mids[i] != mids[unique - 1])
            mids[unique++] = mids[i];
    nw = unique;
    for (size_t i = 0; i < nw; ++i)
        edges[i] = mids[i];
    edges[nw] = 851;
    tt_spectral_grid grid = {nw, NZ, edges, mids, temp, NULL};
    double temperatures[NP][4] = {
        {180, 240, 300, 360}, {190, 250, 310, 370}, {200, 260, 320, 380}, {210, 270, 330, 390}};
    tt_spectral_table tables[NP];
    for (int kind = 0; kind < TT_CS_COUNT; ++kind) {
        for (size_t p = 0; p < NP; ++p) {
            tables[p] = (tt_spectral_table){nw, 4, 4, data[p], temperatures[p]};
            for (size_t c = 0; c < 4; ++c)
                for (size_t w = 0; w < nw; ++w) {
                    double value = (1 + w * .001 + p * .1) * (c + 1) * 1.e-20;
                    if (kind == TT_CS_ACETONE && c)
                        value = pow(1.e-3, (double)c) * .2;
                    if (kind == TT_CS_CLONO2 && c)
                        value = pow(1.e-3, (double)c) * .2;
                    if (kind == TT_CS_CH2O && c == 1)
                        value *= 1.e-4;
                    if ((kind == TT_CS_CH3ONO2 || kind == TT_CS_HNO3 || kind == TT_CS_RONO2) &&
                        c == 1)
                        value = .001 * (w % 9);
                    if (kind == TT_CS_N2O5 && p == 1)
                        value = .01 * (1 + w % 7);
                    data[p][c * nw + w] = value;
                }
        }
        double hcfc[15] = {-46,    -.1,    -.001, 1.e-5, -1.e-7, .001,   -1.e-4, 1.e-6,
                           -1.e-8, 1.e-10, 1.e-6, 1.e-8, 1.e-10, 1.e-12, 1.e-14};
        double hcfc_temp[] = {295, 250};
        if (kind == TT_CS_HCFC)
            tables[0] = (tt_spectral_table){5, 3, 2, hcfc, hcfc_temp};
        size_t count = kind == TT_CS_HCFC ? 1 : NP;
        tt_status status = tt_cross_section((tt_cross_section_kind)kind, &grid, count, tables, out);
        if (status) {
            fprintf(stderr, "%s: %s\n", tt_cross_section_name((tt_cross_section_kind)kind),
                    tt_status_string(status));
            return 1;
        }
        for (size_t i = 0; i < nw * NZ; ++i)
            if (!isfinite(out[i]))
                return 1;
#ifdef TT_HAVE_REFERENCE
        ref_cross(kind, (int)nw, NZ, edges, mids, temp, (int)count, tables, expected);
        for (size_t i = 0; i < nw * NZ; ++i)
            check(tt_cross_section_name((tt_cross_section_kind)kind), i, out[i], expected[i]);
#endif
    }
    /* Temperatures are validated only by kernels that use them, as in TUV-x. */
    {
        double cold[NZ];
        memcpy(cold, temp, sizeof cold);
        cold[3] = 0;
        tt_spectral_grid frozen = {nw, NZ, edges, mids, cold, NULL};
        for (size_t p = 0; p < NP; ++p)
            tables[p] = (tt_spectral_table){nw, 4, 4, data[p], temperatures[p]};
        for (int kind = 0; kind < TT_CS_COUNT; ++kind) {
            int independent = kind == TT_CS_BASE || kind == TT_CS_HOBR ||
                              kind == TT_CS_NITROXY_ACETONE || kind == TT_CS_NITROXY_ETHANOL ||
                              kind == TT_CS_RAYLEIGH || kind == TT_CS_T_BUTYL_NITRATE;
            if (kind == TT_CS_HCFC)
                continue;
            tt_status status =
                tt_cross_section((tt_cross_section_kind)kind, &frozen, NP, tables, out);
            if (independent ? status != TT_OK : status != TT_INVALID) {
                fprintf(stderr, "%s with 0 K: %s\n",
                        tt_cross_section_name((tt_cross_section_kind)kind),
                        tt_status_string(status));
                return 1;
            }
        }
    }
    /* Regrid a real embedded spectrum: its integral must be conserved when
     * target bins are the original knots (with zero extension outside). */
    const tt_dataset *d = tt_dataset_find("data/cross_sections/CCl4_1.nc");
    const tt_data_variable *wl = tt_dataset_variable(d, "wavelength");
    const tt_data_variable *values = tt_dataset_variable(d, "cross_section_parameters");
    if (!d || !wl || !values || wl->count > NW)
        return 1;
    if (tt_regrid_column(d, "cross_section_parameters", 0, wl->count - 1, wl->values,
                         TT_EXTRAPOLATE_ZERO, 0, TT_EXTRAPOLATE_ZERO, 0, out))
        return 1;
    for (size_t i = 0; i + 1 < wl->count; ++i)
        check("embedded regrid", i, out[i], .5 * (values->values[i] + values->values[i + 1]));
    printf("26 cross-section kernels: %zu comparisons; worst %.3g relative epsilon-scaled error\n",
           comparisons, worst);
    return 0;
}
