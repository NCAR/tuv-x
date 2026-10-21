/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "model_internal.h"
#include "model_plan.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static tt_status failure(tt_error *error, tt_status status, const char *step) {
    if (error) {
        error->status = status;
        snprintf(error->message, sizeof error->message, "%s: %s", step, tt_status_string(status));
    }
    return status;
}
tt_status tt_model_set_radiation(tt_model *m, const tt_radiation *r) {
    if (!m || !m->host_solver || !r || !r->direct_flux || !r->down_flux || !r->up_flux)
        return TT_INVALID;
    size_t n = m->shape.levels * m->shape.wavelengths;
    const double *fields[] = {r->direct_flux,       r->down_flux,       r->up_flux,
                              r->direct_irradiance, r->down_irradiance, r->up_irradiance};
    for (size_t f = 0; f < 6; ++f)
        if (fields[f])
            for (size_t i = 0; i < n; ++i)
                if (!isfinite(fields[f][i]))
                    return TT_INVALID;
    for (size_t f = 0; f < 6; ++f) {
        if (fields[f])
            memcpy(m->host_radiation + f * n, fields[f], n * sizeof(double));
        else
            memset(m->host_radiation + f * n, 0, n * sizeof(double));
    }
    return TT_OK;
}
tt_status tt_model_set_execution(tt_model *m, tt_execution execution) {
    if (!m || (execution != TT_EXECUTION_PREPARED && execution != TT_EXECUTION_GENERIC))
        return TT_INVALID;
    m->generic = execution == TT_EXECUTION_GENERIC;
    return m->generic || m->plan ? TT_OK : TT_NOT_FOUND;
}
tt_status tt_model_run(tt_model *m, double zenith, double scale, const tt_model_result **out,
                       tt_error *error) {
    if (error)
        memset(error, 0, sizeof *error);
    if (!m || !out || !isfinite(scale))
        return failure(error, TT_INVALID, "Invalid model run arguments");
    /* A prepared run gives the generic results; for inputs it does not
     * cover it declines, and the generic sequence decides. */
    if (m->plan && !m->generic && tt_model_plan_run(m->plan, m, zenith, scale)) {
        *out = &m->result;
        return TT_OK;
    }
    const tt_grid *height = &m->grids[m->height].value,
                  *wavelength = &m->grids[m->wavelength].value;
    const tt_profile *temperature = &m->profiles[m->temperature].value,
                     *air = &m->profiles[m->air].value;
    size_t nz = m->shape.levels, nw = m->shape.wavelengths, nl = nz - 1, n = nz * nw,
           layer = nl * nw;
    tt_status status = tt_spherical_geometry(nz, height->edges, zenith, m->nid, m->path);
    if (status)
        return failure(error, status, "Spherical geometry");
    int oxygen = m->bands.lyman_alpha != SIZE_MAX || m->bands.schumann_runge != SIZE_MAX;
    if (oxygen) {
        status = tt_air_mass(nz, air->extended_columns, m->nid, m->path, m->vertical, m->slant_air);
        if (status)
            return failure(error, status, "Air columns");
        if (m->scaled_oxygen)
            for (size_t i = 0; i < nz; ++i)
                m->slant_oxygen[i] = m->oxygen_scale * m->slant_air[i];
        else {
            const tt_profile *o2 = &m->profiles[m->oxygen].value;
            status = tt_air_mass(nz, o2->extended_columns, m->nid, m->path, m->oxygen_vertical,
                                 m->slant_oxygen);
            if (status)
                return failure(error, status, "O2 columns");
        }
    }
    for (size_t i = 0; i < m->radiator_count; ++i) {
        tt_model_radiator_entry *r = m->radiators + i;
        if (r->gas) {
            status =
                tt_spectrum_calculate(m->spectra[r->spectrum].spectrum, nl, temperature->midpoints,
                                      air->midpoints, m->cross_section);
            if (status)
                return failure(error, status, "Radiator cross section");
            status = tt_gas_optics(nl, nw, m->profiles[r->profile].value.columns, m->cross_section,
                                   r->view.is_air, r->data, r->data + layer, r->data + 2 * layer);
            if (status)
                return failure(error, status, "Gas optical properties");
        }
    }
    if (oxygen && m->oxygen_radiator != SIZE_MAX) {
        status = tt_oxygen_optics(&m->bands, nw, nz, temperature->edges, m->vertical, m->slant_air,
                                  m->slant_oxygen, m->radiators[m->oxygen_radiator].data, NULL);
        if (status)
            return failure(error, status, "O2 optical depth");
    }
    if (m->host_solver)
        memcpy(m->radiation_storage, m->host_radiation, 6 * n * sizeof(double));
    else {
        status = tt_mix_radiators(m->radiator_count, m->radiator_views, nl, nw, m->streams,
                                  m->optical_storage, m->optical_storage + layer,
                                  m->optical_storage + 2 * layer);
        if (status)
            return failure(error, status, "Radiator mixing");
        const double *albedo = m->profiles[m->surface].value.midpoints;
        status = m->streams == 1
                     ? tt_delta_eddington(&m->atmosphere, zenith, albedo, &m->result.radiation)
                     : tt_discrete_ordinates(&m->atmosphere, m->streams, zenith, albedo,
                                             &m->result.radiation);
        if (status)
            return failure(error, status, "Radiative transfer");
    }
    /* TUV-x apply_scale_factor multiplies by the supplied distance directly. */
    for (size_t i = 0; i < 6 * n; ++i)
        m->radiation_storage[i] = m->radiation_storage[i] * scale;
    if (m->solar != SIZE_MAX) {
        const double *solar = m->profiles[m->solar].value.midpoints;
        status = tt_actinic_flux(nw, nz, &m->result.radiation, solar, m->work);
        if (status)
            return failure(error, status, "Actinic flux");
        status = tt_spectral_irradiance(nw, nz, &m->result.radiation, solar, wavelength->midpoints,
                                        m->work + n);
        if (status)
            return failure(error, status, "Spectral irradiance");
    }
    for (size_t i = 0; i < m->shape.photolysis; ++i) {
        tt_model_reaction *r = m->reactions + i;
        status = tt_spectrum_calculate(r->cross_section, nz, temperature->edges, air->edges,
                                       m->cross_section);
        if (status)
            return failure(error, status, "Reaction cross section");
        status = tt_spectrum_calculate(r->yield, nz, temperature->edges, air->edges, m->yield);
        if (status)
            return failure(error, status, "Reaction quantum yield");
        if (oxygen && r->oxygen) {
            status = tt_oxygen_optics(&m->bands, nw, nz, temperature->edges, m->vertical,
                                      m->slant_air, m->slant_oxygen, NULL, m->cross_section);
            if (status)
                return failure(error, status, "O2 reaction cross section");
        }
        status = tt_photolysis_rate(nw, nz, m->work, m->cross_section, m->yield, r->scale,
                                    m->rate_storage + i * nz);
        if (status)
            return failure(error, status, "Photolysis rate");
        if (r->heating_index != SIZE_MAX) {
            status =
                tt_heating_rate(nw, nz, m->work, m->cross_section, m->yield, wavelength->midpoints,
                                r->energy_wavelength, r->scale,
                                m->rate_storage + (m->shape.photolysis + r->heating_index) * nz);
            if (status)
                return failure(error, status, "Heating rate");
        }
    }
    for (size_t i = 0; i < m->shape.dose; ++i) {
        status = tt_dose_rate(nw, nz, m->work + n, m->doses[i].weights,
                              m->rate_storage + (m->shape.photolysis + m->shape.heating + i) * nz);
        if (status)
            return failure(error, status, "Dose rate");
    }
    *out = &m->result;
    return TT_OK;
}
