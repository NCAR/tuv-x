#include "turbotuv.h"
#include <math.h>

static double erythema(double w) {
    if (w <= 298)
        return 1;
    if (w <= 328)
        return pow(10, .094 * (298 - w));
    if (w <= 400)
        return pow(10, .015 * (140 - w));
    return 1.e-36;
}
static double scup(double w) {
    const double x[] = {270, 302, 334, 367, 400}, a[] = {-10.91, -.86, -8.60, -9.36, -13.15};
    double p = 0;
    for (size_t i = 0; i < 5; ++i) {
        double numerator = 1, denominator = 1;
        for (size_t j = 0; j < 5; ++j)
            if (j != i) {
                numerator *= w - x[j];
                denominator *= x[i] - x[j];
            }
        p = p + (a[i] / denominator) * numerator;
    }
    return exp(p);
}

tt_status tt_spectral_weight(tt_weight_kind kind, size_t nw, const double *edges,
                             const double *mids, const double *base, double a, double b,
                             double *out) {
    if ((unsigned)kind >= TT_WEIGHT_COUNT || !nw || !edges || !mids || !out ||
        ((kind == TT_WEIGHT_BASE || kind == TT_WEIGHT_EPPLEY) && !base))
        return TT_INVALID;
    double accum = 0, scup_factor = kind == TT_WEIGHT_SCUP ? 1 / scup(300) : 1;
    for (size_t w = 0; w < nw; ++w) {
        double l = mids[w], value = 0;
        if (!isfinite(l) || l <= 0 || !isfinite(edges[w]) || !isfinite(edges[w + 1]) ||
            edges[w + 1] <= edges[w])
            return TT_INVALID;
        switch (kind) {
        case TT_WEIGHT_BASE:
        case TT_WEIGHT_EPPLEY:
            value = base[w];
            break;
        case TT_WEIGHT_NOTCH:
            value = a < l && l < b ? 1 : 0;
            break;
        case TT_WEIGHT_GAUSSIAN:
            value = exp(-(log(2) * .04 * ((l - a) * (l - a))));
            break;
        case TT_WEIGHT_PAR:
            value = 400 < l && l < 700 ? 8.36e-3 * l : 0;
            break;
        case TT_WEIGHT_EXPONENTIAL:
            value = pow(10, (300 - l) / 14);
            break;
        case TT_WEIGHT_SCUP:
            value = scup(l) * scup_factor;
            break;
        case TT_WEIGHT_ERYTHEMA:
            value = erythema(l);
            break;
        case TT_WEIGHT_UV_INDEX:
            value = 40 * erythema(l);
            break;
        case TT_WEIGHT_PHYTOPLANKTON:
            if (l > 290 && l < 400)
                value = fmax(0, -3.17e-6 + exp(112.5 + l * (-.6223 + l * 7.67e-4)));
            break;
        case TT_WEIGHT_PLANT:
            value = 570.25 + l * (-4.70144 + l * (.01274 + l * (-1.13118e-5)));
            if (value < 0 || l > 313)
                value = 0;
            break;
        case TT_WEIGHT_PLANT_FLINT:
        case TT_WEIGHT_PLANT_FLINT_EXT:
            value = exp(4.688272 * exp(-exp(.1703411 * (l - 307.867) / 1.15)) +
                        ((390 - l) / 121.7557 - 4.183832));
            value = value * l / 300;
            if (value < 0 || l > (kind == TT_WEIGHT_PLANT_FLINT ? 366 : 390))
                value = 0;
            break;
        case TT_WEIGHT_COUNT:
            return TT_INVALID;
        }
        out[w] = value;
        if (kind == TT_WEIGHT_GAUSSIAN)
            accum = accum + value;
        if (kind == TT_WEIGHT_EPPLEY)
            accum = accum + value * (edges[w + 1] - edges[w]);
    }
    if (kind == TT_WEIGHT_GAUSSIAN || kind == TT_WEIGHT_EPPLEY) {
        if (accum == 0 || !isfinite(accum))
            return TT_RANGE;
        for (size_t w = 0; w < nw; ++w)
            out[w] = kind == TT_WEIGHT_GAUSSIAN ? out[w] / accum : 90 * out[w] / accum;
    }
    return TT_OK;
}
