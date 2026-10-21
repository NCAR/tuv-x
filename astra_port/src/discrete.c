/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Pseudo-spherical, delta-M discrete ordinates with a Lambertian surface.
 * TUV-x fixes thermal emission off and computes the azimuth-averaged fluxes;
 * this implementation specializes those same choices. */
#include "discrete_internal.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double transfer_pi = 3.1415926535898;
static const double solver_pi = 3.14159265358979323846;

typedef struct {
    size_t layers, streams, half, cut, band_rows;
    int layer_cut;
    double solar_cosine;
    double *paths, *tau, *scaled_tau, *delta, *ssa, *omega, *fraction, *gl;
    double *slant, *original_slant, *beam, *level_tau, *beam_cosine;
    double *cosines, *weights, *polynomials, *solar_polynomials;
    double *eigenvalues, *eigenvectors, *particular, *constants;
    double *band, *rhs, *condition;
    size_t *level_layer, *pivots;
    int *nid;
} discrete_state;

/* One scratch object is reused for every layer and wavelength. The reference
 * supports at most 32 streams, hence at most 16 reduced eigenvectors. */
typedef struct {
    double minus[16 * 16], plus[16 * 16], product[16 * 16], vectors[16 * 16], values[16];
    double eigen_work[16 * 16 + 2 * 16], scattering[32 * 32], matrix[32 * 32], condition[32];
    double rhs[32], saved[32], downward[32];
    size_t pivots[32];
} layer_work;

#define MAT(a, r, c, n) ((a)[(c) * (n) + (r)])

static tt_status layer_eigen(const discrete_state *s, const double *gl, double *kk, double *gc,
                             layer_work *q) {
    size_t ns = s->streams, h = s->half;
    double *cc = q->scattering, *minus = q->minus, *plus = q->plus;
    for (size_t i = 0; i < h; ++i) {
        for (size_t j = 0; j < ns; ++j) {
            double sum = 0;
            for (size_t l = 0; l < ns; ++l)
                sum = sum + gl[l] * s->polynomials[i * ns + l] * s->polynomials[j * ns + l];
            MAT(cc, i, j, ns) = .5 * sum * s->weights[j];
        }
        for (size_t j = 0; j < h; ++j) {
            MAT(cc, i + h, j, ns) = MAT(cc, i, j + h, ns);
            MAT(cc, i + h, j + h, ns) = MAT(cc, i, j, ns);
            double alpha = MAT(cc, i, j, ns) / s->cosines[i],
                   beta = MAT(cc, i, j + h, ns) / s->cosines[i];
            MAT(minus, i, j, h) = alpha - beta;
            MAT(plus, i, j, h) = alpha + beta;
        }
        MAT(minus, i, i, h) = MAT(minus, i, i, h) - 1 / s->cosines[i];
        MAT(plus, i, i, h) = MAT(plus, i, i, h) - 1 / s->cosines[i];
    }
    for (size_t i = 0; i < h; ++i)
        for (size_t j = 0; j < h; ++j) {
            double sum = 0;
            for (size_t k = 0; k < h; ++k)
                sum = sum + MAT(plus, i, k, h) * MAT(minus, k, j, h);
            MAT(q->product, i, j, h) = sum;
        }
    tt_status status = tt_discrete_eigen(h, q->product, q->values, q->vectors, q->eigen_work);
    if (status)
        return status;
    for (size_t i = 0; i < h; ++i) {
        q->values[i] = sqrt(fabs(q->values[i]));
        kk[i + h] = q->values[i];
        kk[h - 1 - i] = -q->values[i];
        if (q->values[i] == 0)
            return TT_SINGULAR;
    }
    for (size_t j = 0; j < h; ++j)
        for (size_t i = 0; i < h; ++i) {
            double sum = 0;
            for (size_t k = 0; k < h; ++k)
                sum = sum + MAT(minus, i, k, h) * MAT(q->vectors, k, j, h);
            MAT(plus, i, j, h) = sum / q->values[j];
        }
    for (size_t j = 0; j < h; ++j)
        for (size_t i = 0; i < h; ++i) {
            double p = MAT(plus, i, j, h), m = MAT(q->vectors, i, j, h);
            MAT(gc, i + h, j + h, ns) = .5 * (p + m);
            MAT(gc, h - 1 - i, j + h, ns) = .5 * (p - m);
            p = -p;
            MAT(gc, i + h, h - 1 - j, ns) = .5 * (p + m);
            MAT(gc, h - 1 - i, h - 1 - j, ns) = .5 * (p - m);
        }
    return TT_OK;
}

static tt_status beam_solution(const discrete_state *s, const double *gl, double mu, double *zz,
                               int *unstable, layer_work *q) {
    size_t ns = s->streams, h = s->half;
    for (size_t i = 0; i < ns; ++i) {
        for (size_t j = 0; j < ns; ++j)
            MAT(q->matrix, i, j, ns) = -MAT(q->scattering, i, j, ns);
        MAT(q->matrix, i, i, ns) = 1 + s->cosines[i] / mu + MAT(q->matrix, i, i, ns);
        double sum = 0;
        for (size_t k = 0; k < ns; ++k)
            sum = sum + gl[k] * s->polynomials[i * ns + k] * s->solar_polynomials[k];
        q->rhs[i] = sum / (4 * solver_pi);
    }
    tt_matrix matrix = {TT_MATRIX_DENSE, ns, ns, 0, 0, q->matrix};
    double rcond;
    tt_status status = tt_lu_factor(&matrix, q->pivots, &rcond, q->condition);
    if (status && status != TT_SINGULAR)
        return status;
    *unstable = fabs(rcond) < 1.e-7;
    if (*unstable)
        return TT_OK;
    status = tt_lu_solve(&matrix, q->pivots, 0, q->rhs);
    if (status)
        return status;
    for (size_t i = 0; i < h; ++i) {
        zz[i + h] = q->rhs[i];
        zz[h - 1 - i] = q->rhs[i + h];
    }
    return TT_OK;
}

static tt_status solve_layer(const discrete_state *s, size_t layer, layer_work *q) {
    size_t ns = s->streams;
    double *gl = s->gl + layer * ns, *kk = s->eigenvalues + layer * ns;
    double *gc = s->eigenvectors + layer * ns * ns, *zz = s->particular + layer * ns;
    double original = s->omega[layer], ssa = original, down_ssa = 0, add = -1.e-7;
    int upward = 0, downward = 0, done = 0, no_up = 0, no_down = 0, small = 0, large = 0;
    memcpy(q->saved, gl, ns * sizeof(double));
    for (size_t attempt = 0; attempt < 20000002; ++attempt) {
        if (downward > 0 || upward > 0) {
            double factor = (ssa + add) / ssa;
            for (size_t k = 0; k < ns; ++k)
                gl[k] = gl[k] * factor;
            ssa = ssa + add;
            if (ssa < 100 * DBL_EPSILON) {
                no_down = 1;
                downward = -1;
                small = 1;
            } else if (ssa > 1 - 100 * DBL_EPSILON) {
                no_up = 1;
                upward = -1;
                large = 1;
            }
        }
        if (!small && !large) {
            tt_status status = layer_eigen(s, gl, kk, gc, q);
            if (status)
                return status;
            int unstable;
            status = beam_solution(s, gl, s->beam_cosine[layer], zz, &unstable, q);
            if (status)
                return status;
            if (!unstable && upward == 0 && downward == 0)
                return TT_OK;
            if (unstable) {
                if (upward == 0)
                    ++downward;
                else if (upward > 0)
                    ++upward;
                continue;
            }
        }
        if (!large) {
            small = 0;
            if (downward != 0 && upward == 0) {
                memcpy(q->downward, gl, ns * sizeof(double));
                memcpy(gl, q->saved, ns * sizeof(double));
                down_ssa = ssa;
                ssa = original;
                add = 1.e-7;
                ++upward;
                continue;
            }
        } else
            large = 0;
        if (done) {
            memcpy(gl, q->saved, ns * sizeof(double));
            return TT_OK;
        }
        if (no_up && no_down) {
            memcpy(gl, q->saved, ns * sizeof(double));
            ssa = original;
            done = 1;
            continue;
        }
        if (no_up) {
            memcpy(gl, q->downward, ns * sizeof(double));
            ssa = down_ssa;
            done = 1;
            continue;
        }
        if (no_down) {
            memcpy(gl, q->saved, ns * sizeof(double));
            return TT_OK;
        }
        if (upward > downward) {
            memcpy(gl, q->downward, ns * sizeof(double));
            ssa = down_ssa;
            done = 1;
            continue;
        }
        memcpy(gl, q->saved, ns * sizeof(double));
        return TT_OK;
    }
    return TT_RANGE;
}

static tt_status setup_layers(discrete_state *s, const tt_atmosphere *atmosphere,
                              size_t wavelength) {
    size_t n = s->layers, ns = s->streams, count = n * atmosphere->wavelengths;
    double absorption = 0;
    s->tau[0] = s->scaled_tau[0] = 0;
    s->cut = n;
    for (size_t l = 0; l < n; ++l) {
        size_t i = wavelength * n + n - 1 - l;
        double tau = atmosphere->optical_depth[i], ssa = atmosphere->single_scattering_albedo[i];
        if (!isfinite(tau) || tau < 0 || !isfinite(ssa) || ssa < 0 || ssa > 1)
            return TT_INVALID;
        if (ssa == 1)
            ssa = 1 - 100 * DBL_EPSILON;
        s->ssa[l] = ssa;
        s->tau[l + 1] = s->tau[l] + tau;
        if (absorption < 10000)
            s->cut = l + 1;
        absorption = absorption + (1 - ssa) * tau;
        double f = atmosphere->legendre[(ns - 1) * count + i];
        if (!isfinite(f) || f < -1 || f >= 1)
            return TT_INVALID;
        s->fraction[l] = f;
        s->omega[l] = ssa * (1 - f) / (1 - f * ssa);
        s->delta[l] = (1 - f * ssa) * tau;
        s->scaled_tau[l + 1] = s->scaled_tau[l] + s->delta[l];
        for (size_t k = 0; k < ns; ++k) {
            double moment = k ? atmosphere->legendre[(k - 1) * count + i] : 1;
            if (!isfinite(moment) || moment < -1 || moment > 1)
                return TT_INVALID;
            s->gl[l * ns + k] = (double)(2 * k + 1) * s->omega[l] * (moment - f) / (1 - f);
        }
    }
    s->layer_cut = absorption >= 10000 && n > 1;
    if (!s->layer_cut)
        s->cut = n;
    s->slant[0] = s->original_slant[0] = 0;
    if (s->solar_cosine < 0) {
        if (s->nid[0] < 0)
            s->slant[0] = s->original_slant[0] = 1.e36;
        else
            for (int j = 0; j < s->nid[0]; ++j) {
                s->slant[0] = s->slant[0] + 2 * s->delta[j] * s->paths[j];
                size_t i = wavelength * n + n - 1 - (size_t)j;
                s->original_slant[0] =
                    s->original_slant[0] + 2 * atmosphere->optical_depth[i] * s->paths[j];
            }
    }
    s->beam[0] = exp(-s->slant[0]);
    for (size_t l = 1; l <= n; ++l) {
        s->beam_cosine[l - 1] = 1 / 1.e36;
        if (s->nid[l] < 0)
            s->slant[l] = s->original_slant[l] = 1.e36;
        else {
            double sum = 0, original = 0;
            size_t end = (size_t)s->nid[l], split = end < l ? end : l;
            for (size_t j = 0; j < split; ++j) {
                sum = sum + s->delta[j] * s->paths[l * n + j];
                original = original + atmosphere->optical_depth[wavelength * n + n - 1 - j] *
                                          s->paths[l * n + j];
            }
            for (size_t j = split; j < end; ++j) {
                sum = sum + 2 * s->delta[j] * s->paths[l * n + j];
                original = original + 2 * atmosphere->optical_depth[wavelength * n + n - 1 - j] *
                                          s->paths[l * n + j];
            }
            s->slant[l] = sum;
            s->original_slant[l] = original;
            if (s->slant[l] == s->slant[l - 1])
                s->beam_cosine[l - 1] = 1.e36;
            else {
                double mu =
                    (s->scaled_tau[l] - s->scaled_tau[l - 1]) / (s->slant[l] - s->slant[l - 1]);
                s->beam_cosine[l - 1] = copysign(fmax(fabs(mu), 1 / 1.e36), mu);
            }
        }
        s->beam[l] = exp(-s->slant[l]);
    }
    for (size_t z = 0; z <= n; ++z) {
        size_t layer = 0;
        while (layer + 1 < n && !(s->tau[z] >= s->tau[layer] && s->tau[z] <= s->tau[layer + 1]))
            ++layer;
        s->level_layer[z] = layer;
        s->level_tau[z] = s->scaled_tau[layer] +
                          (1 - s->ssa[layer] * s->fraction[layer]) * (s->tau[z] - s->tau[layer]);
    }
    for (size_t i = 0; i < s->half; ++i)
        for (size_t l = 0; l < n; ++l)
            if (fabs(s->beam_cosine[l]) < 1.e5 &&
                fabs(1 - fabs(s->beam_cosine[l]) / s->cosines[i]) < (double).05f)
                s->beam_cosine[l] = s->beam_cosine[l] * (double).999f;
    return TT_OK;
}

static tt_status boundary_solution(discrete_state *s, double albedo) {
    size_t ns = s->streams, h = s->half, n = s->cut, ld = s->band_rows, ncol = ns * n;
    double attenuation[16];
    memset(s->band, 0, ld * ncol * sizeof(double));
    for (size_t l = 0; l < n; ++l) {
        const double *kk = s->eigenvalues + l * ns, *gc = s->eigenvectors + l * ns * ns;
        for (size_t i = 0; i < h; ++i)
            attenuation[i] = exp(kk[i] * s->delta[l]);
        for (size_t i = 0; i < ns; ++i)
            for (size_t j = 0; j < ns; ++j) {
                size_t row = ld - 2 * ns - i + j, col = l * ns + i;
                double g = MAT(gc, j, i, ns);
                MAT(s->band, row + ns, col, ld) = i < h ? g : g * attenuation[ns - 1 - i];
                MAT(s->band, row, col, ld) = i < h ? -g * attenuation[i] : -g;
            }
    }
    for (size_t i = 0; i < ns; ++i) {
        double expa = i < h ? exp(s->eigenvalues[i] * s->scaled_tau[1]) : 1;
        size_t row = ld - 2 * ns - i + h;
        for (size_t j = h; j-- > 0;)
            MAT(s->band, row++, i, ld) = MAT(s->eigenvectors, j, i, ns) * expa;
    }
    const double *gc = s->eigenvectors + (n - 1) * ns * ns, *zz = s->particular + (n - 1) * ns;
    for (size_t i = 0; i < ns; ++i) {
        double expa = i < h ? 1 : attenuation[ns - 1 - i];
        size_t row = ld - ns - i;
        for (size_t j = h; j < ns; ++j) {
            double value = MAT(gc, j, i, ns);
            if (!s->layer_cut) {
                double sum = 0;
                for (size_t k = 0; k < h; ++k)
                    sum = sum + s->weights[k] * s->cosines[k] * albedo * MAT(gc, h - 1 - k, i, ns);
                value = value - 2 * sum;
            }
            MAT(s->band, row++, (n - 1) * ns + i, ld) = value * expa;
        }
    }
    for (size_t i = 0; i < h; ++i) {
        s->rhs[i] = -s->particular[h - 1 - i] * s->beam[0];
        double value = -zz[i + h] * s->beam[n];
        if (!s->layer_cut) {
            double sum = 0;
            for (size_t j = 0; j < h; ++j)
                sum = sum + s->weights[j] * s->cosines[j] * albedo * (zz[h - 1 - j] * s->beam[n]);
            value = 2 * sum + (albedo * s->solar_cosine / solver_pi - zz[i + h]) * s->beam[n];
        }
        s->rhs[ncol - h + i] = value;
    }
    for (size_t l = 0; l + 1 < n; ++l)
        for (size_t i = 0; i < ns; ++i)
            s->rhs[h + l * ns + i] =
                (s->particular[(l + 1) * ns + i] - s->particular[l * ns + i]) * s->beam[l + 1];
    tt_matrix matrix = {TT_MATRIX_BAND, ncol, ld, 3 * h - 1, 3 * h - 1, s->band};
    double condition;
    tt_status status = tt_lu_factor(&matrix, s->pivots, &condition, s->condition);
    if (status)
        return status;
    status = tt_lu_solve(&matrix, s->pivots, 0, s->rhs);
    if (!status)
        memcpy(s->constants, s->rhs, ncol * sizeof(double));
    return status;
}

static void radiation_output(const discrete_state *s, size_t wavelength, tt_radiation *out) {
    size_t n = s->layers, ns = s->streams, h = s->half;
    for (size_t z = 0; z <= n; ++z) {
        size_t l = s->level_layer[z], i = wavelength * (n + 1) + n - z;
        double fd = 0, fu = 0, ed = 0, eu = 0, beam = 0, direct = 0;
        if (!s->layer_cut || l < s->cut) {
            const double *kk = s->eigenvalues + l * ns, *gc = s->eigenvectors + l * ns * ns;
            const double *ll = s->constants + l * ns, *zz = s->particular + l * ns;
            double factor = exp(-s->slant[z]);
            direct = s->solar_cosine * exp(-s->original_slant[z]);
            beam = factor / (4 * solver_pi);
            for (size_t iq = 0; iq < ns; ++iq) {
                double intensity = 0;
                for (size_t j = 0; j < h; ++j)
                    intensity =
                        intensity + MAT(gc, iq, j, ns) * ll[j] *
                                        exp(-kk[j] * (s->level_tau[z] - s->scaled_tau[l + 1]));
                for (size_t j = h; j < ns; ++j)
                    intensity = intensity + MAT(gc, iq, j, ns) * ll[j] *
                                                exp(-kk[j] * (s->level_tau[z] - s->scaled_tau[l]));
                intensity = intensity + zz[iq] * factor;
                size_t q = iq < h ? h - 1 - iq : iq - h;
                if (iq < h) {
                    fd = fd + s->weights[q] * intensity;
                    ed = ed + s->weights[q] * s->cosines[q] * intensity;
                } else {
                    fu = fu + s->weights[q] * intensity;
                    eu = eu + s->weights[q] * s->cosines[q] * intensity;
                }
            }
            eu = 2 * solver_pi * eu;
            ed = 2 * solver_pi * ed;
            ed = (ed + s->solar_cosine * factor) - direct;
            fu = (2 * solver_pi * fu) / (4 * solver_pi);
            fd = (2 * solver_pi * fd) / (4 * solver_pi);
        }
        out->direct_flux[i] = beam * (4 * transfer_pi);
        out->down_flux[i] = fd * (4 * transfer_pi);
        out->up_flux[i] = fu * (4 * transfer_pi);
        out->direct_irradiance[i] = direct;
        out->down_irradiance[i] = ed;
        out->up_irradiance[i] = eu;
    }
}

tt_status tt_discrete_ordinates(const tt_atmosphere *atmosphere, size_t streams, double zenith,
                                const double *surface, tt_radiation *out) {
    if (!atmosphere || !surface || !out || streams < 2 || streams > 32 || streams % 2 ||
        !isfinite(zenith) || zenith < 0 || zenith > 180 || atmosphere->moments < streams ||
        !atmosphere->layers || atmosphere->layers > INT_MAX - 1 || !atmosphere->wavelengths ||
        !atmosphere->height || !atmosphere->optical_depth ||
        !atmosphere->single_scattering_albedo || !atmosphere->legendre || !out->direct_flux ||
        !out->down_flux || !out->up_flux || !out->direct_irradiance || !out->down_irradiance ||
        !out->up_irradiance)
        return TT_INVALID;
    size_t n = atmosphere->layers, ns = streams, ld = 9 * (ns / 2) - 2,
           max = SIZE_MAX / sizeof(double);
    size_t fixed = ns * ns + 3 * ns + 32, per_layer = ns * ns + ns * (ld + 8) + 16;
    if (n > max / (n + 1) || n * (n + 1) > max - fixed ||
        n > (max - fixed - n * (n + 1)) / per_layer ||
        atmosphere->wavelengths > SIZE_MAX / (n + 1) / streams)
        return TT_INVALID;
    double *memory = calloc(n * (n + 1) + n * per_layer + fixed, sizeof(double));
    size_t *indices = malloc((n + 1 + ns * n) * sizeof(size_t));
    int *nid = malloc((n + 1) * sizeof(int));
    layer_work *scratch = malloc(sizeof(*scratch));
    if (!memory || !indices || !nid || !scratch) {
        free(memory);
        free(indices);
        free(nid);
        free(scratch);
        return TT_MEMORY;
    }
    discrete_state s = {0};
    s.layers = n;
    s.streams = ns;
    s.half = ns / 2;
    s.band_rows = ld;
    s.solar_cosine = cos(zenith * (transfer_pi / 180));
    s.nid = nid;
    s.level_layer = indices;
    s.pivots = indices + n + 1;
    double *next = memory;
#define TAKE(field, count)                                                                         \
    s.field = next;                                                                                \
    next += (count)
    TAKE(paths, n * (n + 1));
    TAKE(tau, n + 1);
    TAKE(scaled_tau, n + 1);
    TAKE(delta, n);
    TAKE(ssa, n);
    TAKE(omega, n);
    TAKE(fraction, n);
    TAKE(gl, n * ns);
    TAKE(slant, n + 1);
    TAKE(original_slant, n + 1);
    TAKE(beam, n + 1);
    TAKE(level_tau, n + 1);
    TAKE(beam_cosine, n);
    TAKE(cosines, ns);
    TAKE(weights, ns);
    TAKE(polynomials, ns * ns);
    TAKE(solar_polynomials, ns);
    TAKE(eigenvalues, n * ns);
    TAKE(eigenvectors, n * ns * ns);
    TAKE(particular, n * ns);
    TAKE(constants, n * ns);
    TAKE(band, ld * ns * n);
    TAKE(rhs, ns * n);
    TAKE(condition, ns * n);
#undef TAKE
    tt_status status = tt_spherical_geometry(n + 1, atmosphere->height, zenith, nid, s.paths);
    if (!status)
        status = tt_discrete_quadrature(s.half, s.cosines, s.weights);
    if (!status) {
        for (size_t i = 0; i < s.half; ++i) {
            s.cosines[i + s.half] = -s.cosines[i];
            s.weights[i + s.half] = s.weights[i];
        }
        tt_discrete_legendre(ns - 1, s.half, s.cosines, s.polynomials);
        double sign = -1;
        for (size_t l = 0; l < ns; ++l) {
            sign = -sign;
            for (size_t i = s.half; i < ns; ++i)
                s.polynomials[i * ns + l] = sign * s.polynomials[(i - s.half) * ns + l];
        }
        double beam_cosine = -s.solar_cosine;
        tt_discrete_legendre(ns - 1, 1, &beam_cosine, s.solar_polynomials);
        for (size_t w = 0; w < atmosphere->wavelengths; ++w) {
            if (!isfinite(surface[w]) || surface[w] < 0 || surface[w] > 1) {
                status = TT_INVALID;
                break;
            }
            status = setup_layers(&s, atmosphere, w);
            if (status)
                break;
            for (size_t l = 0; l < s.cut; ++l) {
                status = solve_layer(&s, l, scratch);
                if (status)
                    break;
            }
            if (status)
                break;
            status = boundary_solution(&s, surface[w]);
            if (status)
                break;
            radiation_output(&s, w, out);
        }
    }
    free(memory);
    free(indices);
    free(nid);
    free(scratch);
    return status;
}
