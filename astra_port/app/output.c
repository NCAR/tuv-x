/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "output.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef TT_NETCDF_OUTPUT
#include <netcdf.h>
#endif
struct tt_output {
    FILE *json;
    const tt_model *model;
    tt_model_dimensions shape;
    size_t step;
#ifdef TT_NETCDF_OUTPUT
    int file, time, zenith, scale, temperature, radiation[8], groups[3];
    int *rates;
#endif
};
static const char *const field_names[] = {
    "direct_flux",     "down_flux",     "up_flux",      "direct_irradiance",
    "down_irradiance", "up_irradiance", "actinic_flux", "spectral_irradiance"};
static const char *const rate_names[] = {"photolysis", "heating", "dose"};
static int fail(tt_error *error, const char *message) {
    if (error) {
        *error = (tt_error){.status = TT_INVALID};
        snprintf(error->message, sizeof error->message, "%s", message);
    }
    return 0;
}
static void string(FILE *f, const char *s) {
    fputc('"', f);
    for (const unsigned char *p = (const unsigned char *)s; *p; ++p) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', f);
            fputc(*p, f);
        } else if (*p < 32)
            fprintf(f, "\\u%04x", *p);
        else
            fputc(*p, f);
    }
    fputc('"', f);
}
static int array(FILE *f, size_t n, const double *v) {
    fputc('[', f);
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(v[i]))
            return 0;
        fprintf(f, "%s%.17g", i ? "," : "", v[i]);
    }
    fputc(']', f);
    return 1;
}
static int matrix(FILE *f, size_t rows, size_t columns, const double *v) {
    fputc('[', f);
    for (size_t i = 0; i < rows; ++i) {
        if (i)
            fputc(',', f);
        if (!array(f, columns, v + i * columns))
            return 0;
    }
    fputc(']', f);
    return 1;
}
#ifdef TT_NETCDF_OUTPUT
static int nc_error(int status, tt_error *error) {
    return status == NC_NOERR ? 1 : fail(error, nc_strerror(status));
}
#define NC(call)                                                                                   \
    do {                                                                                           \
        if (!nc_error((call), error))                                                              \
            return 0;                                                                              \
    } while (0)
static int variable(int file, const char *name, const char *units, int rank, const int *dims,
                    int *id, tt_error *error) {
    NC(nc_def_var(file, name, NC_DOUBLE, rank, dims, id));
    NC(nc_put_att_text(file, *id, "units", strlen(units), units));
    return 1;
}
static int nc_open_output(tt_output *o, const char *path, tt_error *error) {
    NC(nc_create(path, NC_NETCDF4 | NC_CLOBBER, &o->file));
    const char *schema = "TurboTUV-1";
    NC(nc_put_att_text(o->file, NC_GLOBAL, "format", strlen(schema), schema));
    int td, zd, wd, ed, altitude, wavelength, edges;
    NC(nc_def_dim(o->file, "time", NC_UNLIMITED, &td));
    NC(nc_def_dim(o->file, "vertical_level", o->shape.levels, &zd));
    NC(nc_def_dim(o->file, "wavelength", o->shape.wavelengths, &wd));
    NC(nc_def_dim(o->file, "wavelength_edge", o->shape.wavelengths + 1, &ed));
    int spectral[] = {td, wd, zd}, vertical[] = {td, zd};
    if (!variable(o->file, "time", "hours", 1, &td, &o->time, error) ||
        !variable(o->file, "solar zenith angle", "degrees", 1, &td, &o->zenith, error) ||
        !variable(o->file, "Earth-Sun distance", "AU", 1, &td, &o->scale, error) ||
        !variable(o->file, "altitude", "km", 1, &zd, &altitude, error) ||
        !variable(o->file, "wavelength", "nm", 1, &wd, &wavelength, error) ||
        !variable(o->file, "wavelength edges", "nm", 1, &ed, &edges, error) ||
        !variable(o->file, "temperature", "K", 2, vertical, &o->temperature, error))
        return 0;
    const char *nc_names[] = {"direct radiation",  "downward radiation",  "upward radiation",
                              "direct irradiance", "downward irradiance", "upward irradiance",
                              "actinic flux",      "spectral irradiance"};
    for (size_t i = 0; i < 8; ++i)
        if (!variable(o->file, nc_names[i],
                      i < 6    ? "1"
                      : i == 6 ? "photon cm-2 s-1"
                               : "W m-2",
                      3, spectral, &o->radiation[i], error))
            return 0;
    size_t counts[] = {o->shape.photolysis, o->shape.heating, o->shape.dose}, offset = 0;
    const char *units[] = {"s-1", "J s-1", "W m-2"};
    for (int k = 0; k < 3; ++k) {
        NC(nc_def_grp(o->file, rate_names[k], &o->groups[k]));
        for (size_t i = 0; i < counts[k]; ++i) {
            const char *label = tt_model_rate_name(o->model, (tt_rate_kind)k, i);
            char name[NC_MAX_NAME + 1];
            if (strlen(label) > NC_MAX_NAME)
                return fail(error, "Rate label exceeds NetCDF name limit");
            strcpy(name, label);
            for (char *p = name; *p; ++p)
                if (*p == '/')
                    *p = '_';
            if (!variable(o->groups[k], name, units[k], 2, vertical, o->rates + offset + i, error))
                return 0;
            NC(nc_put_att_text(o->groups[k], o->rates[offset + i], "label", strlen(label), label));
        }
        offset += counts[k];
    }
    NC(nc_enddef(o->file));
    NC(nc_put_var_double(o->file, altitude, tt_model_grid(o->model, "height", "km")->edges));
    NC(nc_put_var_double(o->file, wavelength,
                         tt_model_grid(o->model, "wavelength", "nm")->midpoints));
    NC(nc_put_var_double(o->file, edges, tt_model_grid(o->model, "wavelength", "nm")->edges));
    return 1;
}
static int nc_step(tt_output *o, double time, double zenith, double scale,
                   const double *const *fields, const double *const *rates, tt_error *error) {
    size_t start[] = {o->step, 0, 0}, one[] = {1}, vertical[] = {1, o->shape.levels},
           spectral[] = {1, o->shape.wavelengths, o->shape.levels};
    NC(nc_put_vara_double(o->file, o->time, start, one, &time));
    NC(nc_put_vara_double(o->file, o->zenith, start, one, &zenith));
    NC(nc_put_vara_double(o->file, o->scale, start, one, &scale));
    NC(nc_put_vara_double(o->file, o->temperature, start, vertical,
                          tt_model_profile(o->model, "temperature", "K")->edges));
    for (size_t i = 0; i < 8; ++i)
        NC(nc_put_vara_double(o->file, o->radiation[i], start, spectral, fields[i]));
    size_t counts[] = {o->shape.photolysis, o->shape.heating, o->shape.dose}, offset = 0;
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < counts[k]; ++i)
            NC(nc_put_vara_double(o->groups[k], o->rates[offset + i], start, vertical,
                                  rates[k] + i * o->shape.levels));
        offset += counts[k];
    }
    return 1;
}
#endif
int tt_output_open(FILE *json, const char *path, const tt_model *model, tt_output **out,
                   tt_error *error) {
    tt_output *o = calloc(1, sizeof *o);
    if (!o)
        return fail(error, "Cannot allocate output writer");
    o->json = json;
    o->model = model;
    o->shape = tt_model_shape(model);
    if (path) {
#ifdef TT_NETCDF_OUTPUT
        o->file = -1;
        size_t count = o->shape.photolysis + o->shape.heating + o->shape.dose;
        o->rates = calloc(count ? count : 1, sizeof(int));
        if (!o->rates || !nc_open_output(o, path, error)) {
            if (!o->rates)
                fail(error, "Cannot allocate NetCDF variables");
            if (o->file >= 0)
                nc_close(o->file);
            free(o->rates);
            free(o);
            return 0;
        }
#else
        free(o);
        return fail(error, "NetCDF output is disabled; configure TURBOTUV_NETCDF_OUTPUT=ON");
#endif
    } else {
        fputs("{\"format\":\"TurboTUV-1\",\"altitude_km\":", json);
        array(json, o->shape.levels, tt_model_grid(model, "height", "km")->edges);
        fputs(",\"wavelength_nm\":", json);
        array(json, o->shape.wavelengths, tt_model_grid(model, "wavelength", "nm")->midpoints);
        fputs(",\"wavelength_edges_nm\":", json);
        array(json, o->shape.wavelengths + 1, tt_model_grid(model, "wavelength", "nm")->edges);
        size_t counts[] = {o->shape.photolysis, o->shape.heating, o->shape.dose};
        for (int k = 0; k < 3; ++k) {
            fprintf(json, ",\"%s_labels\":[", rate_names[k]);
            for (size_t i = 0; i < counts[k]; ++i) {
                if (i)
                    fputc(',', json);
                string(json, tt_model_rate_name(model, (tt_rate_kind)k, i));
            }
            fputc(']', json);
        }
        fputs(",\"steps\":[\n", json);
    }
    *out = o;
    return 1;
}
int tt_output_step(tt_output *o, double time, double zenith, double scale, const tt_model_result *r,
                   tt_error *error) {
    const double *fields[] = {r->radiation.direct_flux,
                              r->radiation.down_flux,
                              r->radiation.up_flux,
                              r->radiation.direct_irradiance,
                              r->radiation.down_irradiance,
                              r->radiation.up_irradiance,
                              r->actinic_flux,
                              r->spectral_irradiance};
    const double *rates[] = {r->photolysis, r->heating, r->dose};
    size_t counts[] = {o->shape.photolysis, o->shape.heating, o->shape.dose};
    if (!o->json) {
#ifdef TT_NETCDF_OUTPUT
        if (!nc_step(o, time, zenith, scale, fields, rates, error))
            return 0;
#endif
    } else {
        FILE *f = o->json;
        fprintf(f,
                "%s{\"time_hours\":%.17g,\"solar_zenith_angle\":%.17g,\"radiation_scale\":%.17g,"
                "\"temperature_K\":",
                o->step ? ",\n" : "", time, zenith, scale);
        if (!array(f, o->shape.levels, tt_model_profile(o->model, "temperature", "K")->edges))
            return fail(error, "Nonfinite temperature");
        for (size_t i = 0; i < 8; ++i) {
            fprintf(f, ",\"%s\":", field_names[i]);
            if (!matrix(f, o->shape.wavelengths, o->shape.levels, fields[i]))
                return fail(error, "Nonfinite radiation");
        }
        for (size_t i = 0; i < 3; ++i) {
            fprintf(f, ",\"%s\":", rate_names[i]);
            if (!matrix(f, counts[i], o->shape.levels, rates[i]))
                return fail(error, "Nonfinite rate");
        }
        fputc('}', f);
        if (ferror(f))
            return fail(error, "Cannot write output");
    }
    ++o->step;
    return 1;
}
void tt_output_abandon(tt_output *o) {
    if (!o)
        return;
    if (o->json)
        fflush(o->json);
#ifdef TT_NETCDF_OUTPUT
    else
        nc_close(o->file);
    free(o->rates);
#endif
    free(o);
}
int tt_output_close(tt_output *o, tt_error *error) {
    if (!o)
        return 1;
    int success = 1;
    if (o->json) {
        fputs("\n]}\n", o->json);
        if (fflush(o->json) || ferror(o->json))
            success = fail(error, "Cannot finish JSON output");
    }
#ifdef TT_NETCDF_OUTPUT
    else
        success = nc_error(nc_close(o->file), error);
    free(o->rates);
#endif
    free(o);
    return success;
}
