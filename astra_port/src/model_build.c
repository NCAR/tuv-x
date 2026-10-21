/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "model_internal.h"
#include "model_plan.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static char *string_copy(tt_builder *b, const char *value) {
    if (!value) {
        tt_build_fail(b, (tt_config_node){0}, TT_INVALID, "Missing host input name or units");
        return NULL;
    }
    char *out = tt_build_alloc(b, strlen(value) + 1, 1);
    if (out)
        strcpy(out, value);
    return out;
}
static char *name_copy(tt_builder *b, tt_config_node node, const char *key) {
    return string_copy(b, tt_build_string(b, node, key, NULL));
}
/* Keys TUV-x accepts: src/core.F90, src/radiative_transfer/,
 * src/photolysis_rates.F90, src/heating_rates.F90 and src/dose_rates.F90. */
static const char *const root_required[] = {"radiative transfer", "grids", "profiles",
                                            "O2 absorption", NULL};
static const char *const root_optional[] = {"photolysis", "dose rates", "enable diagnostics", NULL};
static const char *const transfer_required[] = {"cross sections", "radiators", NULL};
static const char *const solver_key[] = {"solver", NULL};
static const char *const type_key[] = {"type", NULL};
static const char *const streams_keys[] = {"type", "number of streams", NULL};
static const char *const radiator_required[] = {
    "name", "type", "cross section", "vertical profile", "vertical profile units", NULL};
static const char *const radiator_optional[] = {"enable diagnostics", "treat as air", NULL};
static const char *const aerosol_required[] = {
    "type", "optical depths", "single scattering albedo", "asymmetry factor", "name", NULL};
static const char *const aerosol_optional[] = {"550 nm optical depth", "enable diagnostics", NULL};
static const char *const netcdf_radiator[] = {"type", "name", "file path", NULL};
static const char *const reactions_key[] = {"reactions", NULL};
static const char *const rates_key[] = {"rates", NULL};
static const char *const diagnostics_key[] = {"enable diagnostics", NULL};
static const char *const reaction_required[] = {"name", "cross section", "quantum yield", NULL};
static const char *const reaction_optional[] = {"scaling factor", "heating", NULL};
static const char *const energy_key[] = {"energy term", NULL};
static const char *const dose_required[] = {"name", "weights", NULL};

static size_t combined_count(tt_builder *b, size_t configured, size_t supplied) {
    if (supplied > SIZE_MAX - configured) {
        tt_build_fail(b, (tt_config_node){0}, TT_RANGE, "Too many host inputs");
        return 0;
    }
    return configured + supplied;
}
size_t tt_model_grid_index(const tt_model *m, const char *name, const char *units) {
    for (size_t i = 0; i < m->grid_count; ++i)
        if (m->grids[i].name && m->grids[i].units && !strcmp(name, m->grids[i].name) &&
            !strcmp(units, m->grids[i].units))
            return i;
    return SIZE_MAX;
}
size_t tt_model_profile_index(const tt_model *m, const char *name, const char *units) {
    for (size_t i = 0; i < m->profile_count; ++i)
        if (m->profiles[i].name && m->profiles[i].units && !strcmp(name, m->profiles[i].name) &&
            !strcmp(units, m->profiles[i].units))
            return i;
    return SIZE_MAX;
}
static size_t require_profile(tt_builder *b, tt_config_node node, const tt_model *m,
                              const char *name, const char *units, size_t cells) {
    size_t i = tt_model_profile_index(m, name, units);
    if (i == SIZE_MAX)
        tt_build_fail(b, node, TT_NOT_FOUND, "Missing profile '%s' with units '%s'", name, units);
    else if (m->profiles[i].value.cells != cells)
        tt_build_fail(b, node, TT_INVALID, "Profile '%s' has the wrong grid size", name);
    return i;
}
static void build_grids(tt_builder *b, tt_config_node config, const tt_model_inputs *inputs,
                        tt_model *m) {
    tt_config_node list = tt_config_get(config, "grids");
    size_t configured = tt_config_size(list);
    m->grid_count = combined_count(b, configured, inputs->grid_count);
    m->grids = tt_build_alloc(b, m->grid_count, sizeof *m->grids);
    for (size_t i = 0; i < m->grid_count && b->status == TT_OK; ++i) {
        tt_config_node node = tt_config_at(list, i);
        tt_error error;
        if (i < configured) {
            m->grids[i].name = name_copy(b, node, "name");
            m->grids[i].units = name_copy(b, node, "units");
            tt_status status = tt_grid_from_config(node, &b->source, &m->grids[i].value, &error);
            if (status)
                tt_build_fail(b, node, status, "%s", error.message);
        } else {
            const tt_named_grid *g = inputs->grids + i - configured;
            m->grids[i].name = string_copy(b, g->name);
            m->grids[i].units = string_copy(b, g->units);
            tt_status status =
                g->grid && g->grid->cells < SIZE_MAX
                    ? tt_grid_from_edges(g->grid->cells + 1, g->grid->edges, &m->grids[i].value)
                    : TT_INVALID;
            if (!status)
                status = tt_grid_update(&m->grids[i].value, g->grid->edges, g->grid->midpoints);
            if (status)
                tt_build_fail(b, node, status, "Invalid host grid");
        }
        if (b->status == TT_OK && tt_model_grid_index(m, m->grids[i].name, m->grids[i].units) != i)
            tt_build_fail(b, node, TT_INVALID, "Duplicate grid '%s'", m->grids[i].name);
    }
    if (b->status != TT_OK)
        return;
    m->height = tt_model_grid_index(m, "height", "km");
    m->wavelength = tt_model_grid_index(m, "wavelength", "nm");
    if (m->height == SIZE_MAX || m->wavelength == SIZE_MAX) {
        tt_build_fail(b, config, TT_NOT_FOUND,
                      "Model requires height (km) and wavelength (nm) grids");
        return;
    }
    size_t nl = m->grids[m->height].value.cells, nw = m->grids[m->wavelength].value.cells;
    if (!nl || !nw || nl == SIZE_MAX || nw > SIZE_MAX / (nl + 1) / sizeof(double) / 32 ||
        nl > SIZE_MAX / (nl + 1) / sizeof(double)) {
        tt_build_fail(b, config, TT_RANGE, "Invalid or oversized model grid");
        return;
    }
    m->shape.levels = nl + 1;
    m->shape.wavelengths = nw;
}
static void build_profiles(tt_builder *b, tt_config_node config, const tt_model_inputs *inputs,
                           tt_model *m) {
    tt_config_node list = tt_config_get(config, "profiles");
    size_t configured = tt_config_size(list);
    m->profile_count = combined_count(b, configured, inputs->profile_count);
    m->profiles = tt_build_alloc(b, m->profile_count, sizeof *m->profiles);
    tt_named_grid *grids = tt_build_alloc(b, m->grid_count, sizeof *grids);
    if (!grids)
        return;
    for (size_t i = 0; i < m->grid_count; ++i)
        grids[i] = (tt_named_grid){m->grids[i].name, m->grids[i].units, &m->grids[i].value};
    for (size_t i = 0; i < m->profile_count && b->status == TT_OK; ++i) {
        tt_config_node node = tt_config_at(list, i);
        tt_error error;
        if (i < configured) {
            m->profiles[i].name = name_copy(b, node, "name");
            m->profiles[i].units = name_copy(b, node, "units");
            tt_status status = tt_profile_from_config(node, m->grid_count, grids, &b->source,
                                                      &m->profiles[i].value, &error);
            if (status)
                tt_build_fail(b, node, status, "%s", error.message);
        } else {
            const tt_named_profile *p = inputs->profiles + i - configured;
            m->profiles[i].name = string_copy(b, p->name);
            m->profiles[i].units = string_copy(b, p->units);
            tt_status status = tt_model_copy_profile(p->profile, &m->profiles[i].value);
            if (status)
                tt_build_fail(b, node, status, "Invalid host profile");
        }
        if (b->status == TT_OK &&
            tt_model_profile_index(m, m->profiles[i].name, m->profiles[i].units) != i)
            tt_build_fail(b, node, TT_INVALID, "Duplicate profile '%s'", m->profiles[i].name);
    }
    free(grids);
    if (b->status != TT_OK)
        return;
    m->temperature = require_profile(b, config, m, "temperature", "K", m->shape.levels - 1);
    m->air = require_profile(b, config, m, "air", "molecule cm-3", m->shape.levels - 1);
    m->solar = tt_model_profile_index(m, "extraterrestrial flux", "photon cm-2 s-1");
    if (m->solar != SIZE_MAX && m->profiles[m->solar].value.cells != m->shape.wavelengths)
        tt_build_fail(b, config, TT_INVALID, "Extraterrestrial flux has the wrong grid size");
}
static tt_spectrum *spectrum(tt_builder *b, tt_config_node node, tt_spectrum_role role,
                             const tt_model *m) {
    tt_spectrum *out = NULL;
    tt_error error;
    if (b->status != TT_OK)
        return NULL;
    tt_status status =
        tt_spectrum_create(role, node, &m->grids[m->wavelength].value, &b->source, &out, &error);
    if (status)
        tt_build_fail(b, node, status, "%s", error.message);
    return out;
}
static void build_transfer(tt_builder *b, tt_config_node config, const tt_model_inputs *inputs,
                           tt_model *m) {
    tt_config_node node = tt_config_get(config, "radiative transfer");
    tt_build_keys(b, node, transfer_required, solver_key);
    tt_config_node solver = tt_config_get(node, "solver");
    const char *type = tt_build_string(b, solver, "type", NULL);
    if (!strcmp(type, "delta eddington") || !strcmp(type, "from host"))
        tt_build_keys(b, solver, type_key, NULL);
    else if (!strcmp(type, "discrete ordinate"))
        tt_build_keys(b, solver, streams_keys, NULL);
    if (!strcmp(type, "delta eddington"))
        m->streams = 1;
    else if (!strcmp(type, "discrete ordinate")) {
        int streams = tt_build_integer(b, solver, "number of streams");
        if (streams < 2 || streams > 32 || streams % 2)
            tt_build_fail(b, solver, TT_RANGE,
                          "Discrete ordinate streams must be even, 2 through 32");
        else
            m->streams = (size_t)streams;
    } else if (!strcmp(type, "from host")) {
        m->host_solver = 1;
        m->streams = 1;
    } else
        tt_build_fail(b, solver, TT_INVALID, "Unknown solver '%s'", type);
    if (!m->host_solver)
        m->surface = require_profile(b, node, m, "surface albedo", "none", m->shape.wavelengths);
    tt_config_node list = tt_config_get(node, "cross sections");
    m->spectrum_count = tt_config_size(list);
    m->spectra = tt_build_alloc(b, m->spectrum_count, sizeof *m->spectra);
    for (size_t i = 0; i < m->spectrum_count && b->status == TT_OK; ++i) {
        tt_config_node cs = tt_config_at(list, i);
        m->spectra[i].name = name_copy(b, cs, "name");
        m->spectra[i].spectrum = spectrum(b, cs, TT_SPECTRUM_CROSS_SECTION, m);
    }
    list = tt_config_get(node, "radiators");
    size_t configured = tt_config_size(list);
    m->radiator_count = combined_count(b, configured, inputs->radiator_count);
    m->radiators = tt_build_alloc(b, m->radiator_count, sizeof *m->radiators);
    m->radiator_views = tt_build_alloc(b, m->radiator_count, sizeof *m->radiator_views);
    const tt_grid *height = &m->grids[m->height].value,
                  *wavelength = &m->grids[m->wavelength].value;
    size_t nl = height->cells, nw = wavelength->cells, n = nl * nw;
    for (size_t i = 0; i < m->radiator_count && b->status == TT_OK; ++i) {
        tt_config_node rad = tt_config_at(list, i);
        tt_model_radiator_entry *r = m->radiators + i;
        r->name = i < configured ? name_copy(b, rad, "name")
                                 : string_copy(b, inputs->radiators[i - configured].name);
        type = i < configured ? tt_build_string(b, rad, "type", NULL) : "from host";
        r->data = tt_build_alloc(b, n, 3 * sizeof(double));
        if (!r->data)
            break;
        r->view = (tt_radiator){r->data, r->data + n, r->data + 2 * n, 0};
        if (i >= configured) {
            tt_status status =
                tt_model_copy_radiator(inputs->radiators + i - configured, r, nl, nw);
            if (status)
                tt_build_fail(b, rad, status, "Invalid host radiator");
        } else if (!strcmp(type, "base")) {
            tt_build_keys(b, rad, radiator_required, radiator_optional);
            r->gas = 1;
            const char *pn = tt_build_string(b, rad, "vertical profile", NULL),
                       *pu = tt_build_string(b, rad, "vertical profile units", NULL);
            r->profile = require_profile(b, rad, m, pn, pu, nl);
            const char *cs = tt_build_string(b, rad, "cross section", NULL);
            r->spectrum = SIZE_MAX;
            for (size_t j = 0; j < m->spectrum_count; ++j)
                if (!strcmp(cs, m->spectra[j].name)) {
                    r->spectrum = j;
                    break;
                }
            if (r->spectrum == SIZE_MAX)
                tt_build_fail(b, rad, TT_NOT_FOUND, "Missing radiator cross section '%s'", cs);
            r->view.is_air = tt_build_boolean(b, rad, "treat as air", 0);
        } else if (!strcmp(type, "aerosol")) {
            tt_build_keys(b, rad, aerosol_required, aerosol_optional);
            size_t samples = 0;
            double *depths = tt_build_array(b, rad, "optical depths", &samples);
            double albedo = tt_build_number(b, rad, "single scattering albedo", NAN);
            double asymmetry = tt_build_number(b, rad, "asymmetry factor", NAN);
            /* TUV-x defaults to 0, which normalizes the aerosol column to zero. */
            double tau550 = tt_build_number(b, rad, "550 nm optical depth", 0);
            double alpha = tt_build_number(b, rad, "Alpha", 1);
            tt_status status =
                b->status == TT_OK
                    ? tt_aerosol_optics(height, wavelength, samples, depths, albedo, asymmetry,
                                        tau550, alpha, r->data, r->data + n, r->data + 2 * n)
                    : b->status;
            if (status)
                tt_build_fail(b, rad, status, "Cannot construct aerosol '%s'", r->name);
            free(depths);
        } else if (!strcmp(type, "from netcdf file")) {
            tt_build_keys(b, rad, netcdf_radiator, NULL);
            const char *path = tt_build_string(b, rad, "file path", NULL);
            const tt_dataset *data = tt_build_dataset(b, rad, path);
            tt_radiator view;
            tt_status status =
                b->status == TT_OK ? tt_radiator_from_dataset(data, nl, nw, &view) : b->status;
            if (status)
                tt_build_fail(b, rad, status, "Invalid compiled radiator '%s'", path);
            else {
                memcpy(r->data, view.optical_depth, n * sizeof(double));
                memcpy(r->data + n, view.single_scattering_albedo, n * sizeof(double));
                memcpy(r->data + 2 * n, view.asymmetry, n * sizeof(double));
            }
        } else
            tt_build_fail(b, rad, TT_INVALID, "Unknown radiator type '%s'", type);
        m->radiator_views[i] = r->view;
        if (r->name && !strcmp(r->name, "O2"))
            m->oxygen_radiator = i;
        for (size_t j = 0; j < i && b->status == TT_OK; ++j)
            if (!strcmp(r->name, m->radiators[j].name))
                tt_build_fail(b, rad, TT_INVALID, "Duplicate radiator '%s'", r->name);
    }
}
static void build_rates(tt_builder *b, tt_config_node config, tt_model *m) {
    tt_build_keys(b, tt_config_get(config, "photolysis"), reactions_key, diagnostics_key);
    tt_build_keys(b, tt_config_get(config, "dose rates"), rates_key, diagnostics_key);
    tt_config_node list = tt_config_get(tt_config_get(config, "photolysis"), "reactions");
    m->shape.photolysis = tt_config_size(list);
    m->reactions = tt_build_alloc(b, m->shape.photolysis, sizeof *m->reactions);
    for (size_t i = 0; i < m->shape.photolysis && b->status == TT_OK; ++i) {
        tt_config_node node = tt_config_at(list, i), cs = tt_config_get(node, "cross section");
        tt_build_keys(b, node, reaction_required, reaction_optional);
        tt_build_keys(b, tt_config_get(node, "heating"), energy_key, NULL);
        tt_model_reaction *r = m->reactions + i;
        r->name = name_copy(b, node, "name");
        r->cross_section = spectrum(b, cs, TT_SPECTRUM_CROSS_SECTION, m);
        r->yield = spectrum(b, tt_config_get(node, "quantum yield"), TT_SPECTRUM_QUANTUM_YIELD, m);
        r->oxygen = tt_build_boolean(b, cs, "apply O2 bands", 0);
        r->scale = tt_build_number(b, node, "scaling factor", 1);
        r->heating_index = SIZE_MAX;
        if (tt_build_has(node, "heating")) {
            r->heating_index = m->shape.heating++;
            r->energy_wavelength =
                tt_build_number(b, tt_config_get(node, "heating"), "energy term", NAN);
            if (r->energy_wavelength <= 0)
                tt_build_fail(b, node, TT_RANGE, "Heating energy wavelength must be positive");
        }
    }
    list = tt_config_get(tt_config_get(config, "dose rates"), "rates");
    m->shape.dose = tt_config_size(list);
    m->doses = tt_build_alloc(b, m->shape.dose, sizeof *m->doses);
    for (size_t i = 0; i < m->shape.dose && b->status == TT_OK; ++i) {
        tt_config_node node = tt_config_at(list, i);
        tt_build_keys(b, node, dose_required, NULL);
        m->doses[i].name = name_copy(b, node, "name");
        tt_spectrum *weights = spectrum(b, tt_config_get(node, "weights"), TT_SPECTRUM_WEIGHT, m);
        m->doses[i].weights = tt_build_alloc(b, m->shape.wavelengths, sizeof(double));
        tt_status status = b->status == TT_OK
                               ? tt_spectrum_calculate(weights, 0, NULL, NULL, m->doses[i].weights)
                               : b->status;
        if (status)
            tt_build_fail(b, node, status, "Cannot construct dose weights");
        tt_spectrum_destroy(weights);
    }
    if (m->shape.photolysis || m->shape.dose)
        m->solar = require_profile(b, config, m, "extraterrestrial flux", "photon cm-2 s-1",
                                   m->shape.wavelengths);
}
static void build_oxygen(tt_builder *b, tt_config_node config, tt_model *m) {
    const tt_grid *w = &m->grids[m->wavelength].value;
    tt_status status = tt_oxygen_band_grid(w->cells, w->edges, &m->bands);
    if (status) {
        tt_build_fail(b, config, status, "Wavelength grid does not match the O2 bands");
        return;
    }
    if (m->bands.lyman_alpha == SIZE_MAX && m->bands.schumann_runge == SIZE_MAX)
        return;
    tt_config_node node = tt_config_get(config, "O2 absorption");
    m->scaled_oxygen = tt_build_has(node, "O2 estimate");
    if (m->scaled_oxygen) {
        m->oxygen_scale =
            tt_build_number(b, tt_config_get(node, "O2 estimate"), "scale factor", NAN);
        if (!(m->oxygen_scale > 0))
            tt_build_fail(b, node, TT_RANGE, "O2 scale factor must be positive");
    } else
        m->oxygen = require_profile(b, node, m, "O2", "molecule cm-3", m->shape.levels - 1);
    /* The coefficient values are compiled. Validate a configured path so a
     * typo or unsupported replacement cannot silently select different data. */
    if (m->bands.schumann_runge != SIZE_MAX) {
        const char *path = tt_build_string(b, node, "cross section parameters file",
                                           "data/cross_sections/O2_parameters.txt");
        const tt_dataset *data = tt_build_dataset(b, node, path);
        const tt_dataset *expected = tt_dataset_find("data/cross_sections/O2_parameters.txt");
        if (data && (!data->sha256 || strcmp(data->sha256, expected->sha256)))
            tt_build_fail(b, node, TT_INVALID,
                          "O2 band coefficients differ from the compiled reference");
    }
}
static void workspace(tt_builder *b, tt_model *m) {
    size_t nz = m->shape.levels, nw = m->shape.wavelengths, nl = nz - 1, n = nz * nw,
           layer = nw * nl;
    m->nid = tt_build_alloc(b, nz, sizeof(int));
    m->path = tt_build_alloc(b, nz, nl * sizeof(double));
    m->vertical = tt_build_alloc(b, nz, 4 * sizeof(double));
    m->radiation_storage = tt_build_alloc(b, n, 6 * sizeof(double));
    m->work = tt_build_alloc(b, n, 4 * sizeof(double));
    m->optical_storage = tt_build_alloc(b, layer, (2 + m->streams) * sizeof(double));
    size_t rates = m->shape.photolysis + m->shape.heating + m->shape.dose;
    m->rate_storage = tt_build_alloc(b, rates, nz * sizeof(double));
    if (m->host_solver)
        m->host_radiation = tt_build_alloc(b, n, 6 * sizeof(double));
    if (b->status != TT_OK)
        return;
    m->slant_air = m->vertical + nz;
    m->slant_oxygen = m->slant_air + nz;
    m->oxygen_vertical = m->slant_oxygen + nz;
    double *r = m->radiation_storage;
    m->result.radiation = (tt_radiation){r, r + n, r + 2 * n, r + 3 * n, r + 4 * n, r + 5 * n};
    m->result.actinic_flux = m->work;
    m->result.spectral_irradiance = m->work + n;
    m->cross_section = m->work + 2 * n;
    m->yield = m->work + 3 * n;
    m->result.photolysis = m->rate_storage;
    m->result.heating = m->rate_storage + nz * m->shape.photolysis;
    m->result.dose = m->rate_storage + nz * (m->shape.photolysis + m->shape.heating);
    m->atmosphere = (tt_atmosphere){nl,
                                    nw,
                                    m->streams,
                                    m->grids[m->height].value.edges,
                                    m->optical_storage,
                                    m->optical_storage + layer,
                                    m->optical_storage + 2 * layer};
}
tt_status tt_model_create(tt_config_node config, const tt_data_source *source, tt_model **out,
                          tt_error *error) {
    return tt_model_create_with_inputs(config, source, NULL, out, error);
}
tt_status tt_model_create_with_inputs(tt_config_node config, const tt_data_source *source,
                                      const tt_model_inputs *inputs, tt_model **out,
                                      tt_error *error) {
    if (error)
        memset(error, 0, sizeof *error);
    tt_builder b = {source ? *source : (tt_data_source){0}, error, TT_OK};
    const tt_model_inputs empty = {0};
    if (!inputs)
        inputs = &empty;
    if ((inputs->grid_count && !inputs->grids) || (inputs->profile_count && !inputs->profiles) ||
        (inputs->radiator_count && !inputs->radiators)) {
        tt_build_fail(&b, config, TT_INVALID, "Missing host input array");
        return b.status;
    }
    if (!out || tt_config_type(config) != TT_CONFIG_MAPPING) {
        tt_build_fail(&b, config, TT_INVALID, "Expected model configuration mapping and output");
        return b.status;
    }
    tt_build_keys(&b, config, root_required, root_optional);
    if (b.status != TT_OK)
        return b.status;
    tt_model *m = tt_build_alloc(&b, 1, sizeof *m);
    if (!m)
        return b.status;
    m->height = m->wavelength = m->temperature = m->air = m->oxygen = m->surface = m->solar =
        m->oxygen_radiator = SIZE_MAX;
    build_grids(&b, config, inputs, m);
    if (b.status == TT_OK)
        build_profiles(&b, config, inputs, m);
    if (b.status == TT_OK)
        build_transfer(&b, config, inputs, m);
    if (b.status == TT_OK)
        build_rates(&b, config, m);
    if (b.status == TT_OK)
        build_oxygen(&b, config, m);
    if (b.status == TT_OK)
        workspace(&b, m);
    if (b.status != TT_OK) {
        tt_model_destroy(m);
        return b.status;
    }
    /* Prepared runs are an optimization: without a plan, runs are generic. */
    if (tt_model_plan_create(m, &m->plan))
        m->plan = NULL;
    *out = m;
    return TT_OK;
}
void tt_model_destroy(tt_model *m) {
    if (!m)
        return;
    tt_model_plan_destroy(m->plan);
    for (size_t i = 0; m->grids && i < m->grid_count; ++i) {
        free(m->grids[i].name);
        free(m->grids[i].units);
        tt_grid_destroy(&m->grids[i].value);
    }
    for (size_t i = 0; m->profiles && i < m->profile_count; ++i) {
        free(m->profiles[i].name);
        free(m->profiles[i].units);
        tt_profile_destroy(&m->profiles[i].value);
    }
    for (size_t i = 0; m->spectra && i < m->spectrum_count; ++i) {
        free(m->spectra[i].name);
        tt_spectrum_destroy(m->spectra[i].spectrum);
    }
    for (size_t i = 0; m->radiators && i < m->radiator_count; ++i) {
        free(m->radiators[i].name);
        free(m->radiators[i].data);
    }
    for (size_t i = 0; m->reactions && i < m->shape.photolysis; ++i) {
        free(m->reactions[i].name);
        tt_spectrum_destroy(m->reactions[i].cross_section);
        tt_spectrum_destroy(m->reactions[i].yield);
    }
    for (size_t i = 0; m->doses && i < m->shape.dose; ++i) {
        free(m->doses[i].name);
        free(m->doses[i].weights);
    }
    free(m->grids);
    free(m->profiles);
    free(m->spectra);
    free(m->radiators);
    free(m->radiator_views);
    free(m->reactions);
    free(m->doses);
    free(m->nid);
    free(m->path);
    free(m->vertical);
    free(m->work);
    free(m->radiation_storage);
    free(m->host_radiation);
    free(m->optical_storage);
    free(m->rate_storage);
    free(m);
}
tt_model_dimensions tt_model_shape(const tt_model *m) {
    return m ? m->shape : (tt_model_dimensions){0};
}
const tt_grid *tt_model_grid(const tt_model *m, const char *name, const char *units) {
    if (!m || !name || !units)
        return NULL;
    size_t i = tt_model_grid_index(m, name, units);
    return i == SIZE_MAX ? NULL : &m->grids[i].value;
}
const tt_profile *tt_model_profile(const tt_model *m, const char *name, const char *units) {
    if (!m || !name || !units)
        return NULL;
    size_t i = tt_model_profile_index(m, name, units);
    return i == SIZE_MAX ? NULL : &m->profiles[i].value;
}
const tt_radiator *tt_model_radiator(const tt_model *m, const char *name) {
    if (m && name)
        for (size_t i = 0; i < m->radiator_count; ++i)
            if (!strcmp(name, m->radiators[i].name))
                return &m->radiators[i].view;
    return NULL;
}
const char *tt_model_rate_name(const tt_model *m, tt_rate_kind kind, size_t i) {
    if (!m)
        return NULL;
    if (kind == TT_RATE_PHOTOLYSIS)
        return i < m->shape.photolysis ? m->reactions[i].name : NULL;
    if (kind == TT_RATE_DOSE)
        return i < m->shape.dose ? m->doses[i].name : NULL;
    if (kind == TT_RATE_HEATING && i < m->shape.heating)
        for (size_t j = 0; j < m->shape.photolysis; ++j)
            if (m->reactions[j].heating_index == i)
                return m->reactions[j].name;
    return NULL;
}
