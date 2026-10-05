.. _host-radiation-field:

###################################
Host-Supplied Radiation Field
###################################

Overview
========

TUV-x can calculate photolysis rate constants, dose rates, and heating rates
from a radiation field that the host application supplies. The host application
selects this mode with the ``from host`` radiative transfer solver. TUV-x then
performs no radiative transfer calculation of its own.

Use this mode when the host model already solves radiative transfer across the
TUV-x wavelength range, and you need the photolysis rate constants to stay
consistent with that solution. Read
:ref:`host-radiation-field-resolution` first. A host model that resolves only
the visible and the near infrared cannot drive the ultraviolet reactions.

The other way to couple
-----------------------

TUV-x accepts host data at two different points, and the ``from host`` solver
is the more demanding of the two.

The ``from host`` solver takes the finished radiation field. The host model
must therefore do the ultraviolet radiative transfer itself. It needs
ultraviolet gas absorption cross sections, ultraviolet Rayleigh scattering, and
a wavelength grid fine enough for photolysis.

The ``from host`` radiators take one step earlier in the calculation. The host
model supplies the optical properties of its clouds and its aerosols, which are
the optical depth, the single scattering albedo, and the asymmetry factor for
each vertical layer and each wavelength bin. TUV-x adds its own gas absorbers
and runs its own solver on its own wavelength grid.

The second route asks much less of the host model. Cloud and aerosol optical
properties change smoothly with wavelength, so a host model can usually report
them at ultraviolet wavelengths even when it does not solve radiative transfer
there. Ultraviolet gas absorption does not change smoothly, and that is the
part the ``from host`` solver forces the host model to get right.

Both routes make the photolysis calculation see the same clouds and the same
aerosols that the host model uses for its own radiation calculation. That is
the consistency that usually matters, because the clouds and the aerosols are
what differ between the two models. Choose the radiators unless the host model
truly resolves the ultraviolet. See
:ref:`configuration-radiators-from-host`.

.. _host-radiation-field-contract:

What the host application supplies
==================================

TUV-x forms a photolysis rate constant as:

.. math::

   J_i(z) = \sum_\lambda
            \left( f_{dr} + f_{up} + f_{dn} \right)_{z,\lambda} \,
            F^\infty_\lambda \, \sigma_{i,z,\lambda} \, \phi_{i,z,\lambda}

The host application supplies the three actinic flux components
:math:`f_{dr}`, :math:`f_{up}`, and :math:`f_{dn}`. These are the same
quantities that the internal solvers produce. TUV-x uses their sum for the
photolysis rate constants and for the heating rates.

Dose rates use the irradiance
-----------------------------

Dose rates are the one exception. TUV-x forms a dose rate from the three
**irradiance** components :math:`e_{dr}`, :math:`e_{up}`, and :math:`e_{dn}`,
because a dose rate is an energy flux on a horizontal surface:

.. math::

   D_i(z) = \sum_\lambda
            \left( e_{dr} + e_{up} + e_{dn} \right)_{z,\lambda} \,
            F^\infty_\lambda \, \frac{hc}{\lambda} \, w_{i,\lambda}

The three irradiance components are optional arguments of ``update``, and
TUV-x sets an omitted component to zero. A host application that requests dose
rates must therefore supply all three. If it supplies only the actinic flux
components, the dose rates are zero and TUV-x reports no error.

Actinic flux, not irradiance
----------------------------

The actinic flux integrates the radiance over all directions with equal weight.
The irradiance weights it by the cosine of the angle from the vertical. A
molecule responds to the actinic flux, so TUV-x uses the actinic flux for every
rate calculation.

The internal solvers form one from the other. For the direct beam, the factor is
the cosine of the solar zenith angle. For the two diffuse components, it is the
mean cosine of the diffuse stream, which the two-stream approximation sets. The
Eddington approximation in ``src/radiative_transfer/solvers/delta_eddington.F90``
uses 0.5:

.. code-block:: Fortran

   edr = mu  * fdr      ! mu  = cos( solar zenith angle )
   fdn = edn / mu1      ! mu1 = 0.5 for the Eddington approximation
   fup = eup / mu1

A host application that computes irradiances must divide by the same factors
that its own solver used. If you supply an irradiance where TUV-x expects an
actinic flux, the rate constants are too low, by a factor that varies with the
solar zenith angle.

Units and shape
---------------

- The values are **dimensionless**. TUV-x normalizes the radiation field by the
  extraterrestrial flux, and the internal solvers set the incident flux to
  unity.
- :math:`F^\infty_\lambda` comes from the ``extraterrestrial flux``
  :ref:`profile <configuration-profiles>`, in units of
  photon cm\ :sup:`-2` s\ :sup:`-1` per wavelength bin. Do **not** fold it into
  the values that you supply.
- The values are **per wavelength bin**, not per nm.
- Do **not** fold the Earth-Sun distance factor into the values that you
  supply. TUV-x applies the ``earth_sun_distance`` argument of
  :f:func:`~tuvx_core/core_t%run` to the whole radiation field.
- The array shape is (vertical interface, wavelength bin). The number of
  vertical interfaces is the number of ``height`` grid cells plus one. The
  number of wavelength bins is the number of ``wavelength`` grid cells.
- The vertical index follows the ``height`` grid. Index 1 is the first edge of
  that grid, which is the **lowest altitude**, and the last index is the top of
  the atmosphere. Both internal solvers reverse their own top-down working
  arrays to return the radiation field in this order, so a host application
  that fills the arrays from the top down inverts its whole profile. TUV-x
  cannot detect that, because an inverted profile is a valid set of numbers.

The same units and the same shape apply to the three spectral irradiance
components :math:`e_{dr}`, :math:`e_{up}`, and :math:`e_{dn}`. TUV-x also
reports all six components through
:f:func:`~tuvx_core/core_t%get_radiation_field`.

To convert an absolute photon flux :math:`F` in
photon cm\ :sup:`-2` s\ :sup:`-1` into the value that TUV-x expects, divide by
the extraterrestrial flux in the same bin. The host application reads that
profile with:

.. code-block:: Fortran

   etfl => core%get_profile( "extraterrestrial flux", "photon cm-2 s-1" )
   ! etfl%mid_val_( i_bin ) holds the flux in each wavelength bin
   deallocate( etfl )   ! get_profile returns a copy that the caller owns

.. _configuration-solvers-from-host:

Configuration
=============

Set the radiative transfer solver type to ``from host``:

.. code-block:: JSON

   "radiative transfer": {
     "cross sections": [ ],
     "radiators": [ ],
     "solver": {
       "type": "from host"
     }
   }

There are no other configuration options for this solver. The ``cross
sections`` and ``radiators`` arrays are still required, but they may be empty.
The solver ignores the radiators, the solar zenith angle, and the spherical
geometry, because the host application accounts for all of them when it
computes the radiation field.

Keep the ``O2 absorption`` configuration in place when your wavelength grid
starts at or below 175.4 nm. Keep the ``apply O2 bands`` flag on the
O\ :sub:`2` photolysis cross section as well. TUV-x then applies the
Lyman-alpha and the Schumann-Runge parameterizations to the O\ :sub:`2` **cross
section**. That correction is independent of the solver. See
:ref:`host-radiation-field-la-srb`.

The ``O2`` radiator has no effect in this mode. The ``from host`` solver ignores
every radiator, so TUV-x computes the O\ :sub:`2` optical depth and then
discards it.

Host application call sequence
==============================

.. code-block:: Fortran

   use tuvx_core,             only : core_t
   use tuvx_solver_from_host, only : radiation_field_updater_t

   class(core_t), pointer :: core
   type(radiation_field_updater_t) :: updater
   real(dk), allocatable :: fdr(:,:), fup(:,:), fdn(:,:)
   real(dk), allocatable :: rates(:,:)

   core => core_t( config_path )
   updater = core%get_radiation_field_updater( )

   ! allocate (vertical interface, wavelength bin) arrays and fill them
   ! with the normalized actinic flux components

   call updater%update( direct_actinic_flux   = fdr,                          &
                        upward_actinic_flux   = fup,                          &
                        downward_actinic_flux = fdn )
   call core%run( solar_zenith_angle = sza,                                   &
                  earth_sun_distance = esd,                                   &
                  photolysis_rate_constants = rates )

The host application must call ``update`` before every call to ``run``. The
solver holds the last radiation field that it received, and it starts with a
field of zero. TUV-x reports no error when the host application skips an
update, so a missed update gives either stale rates or rates of zero.

``get_radiation_field_updater`` takes an optional ``found`` flag. TUV-x stops
with an error when the flag is absent and the configured solver is of another
type. When the flag is present, TUV-x sets it to false and returns an updater
that the host application must not use.

The wavelength grid
===================

The host application does not have to use the standard TUV-x wavelength grid.
Supply a ``wavelength`` grid through the
:ref:`grids <configuration-grids>` that a host application provides, and TUV-x
interpolates the cross sections and the quantum yields onto it when it builds
the core. The O\ :sub:`2` parameterization grid below is the one hard
constraint.

.. _host-radiation-field-la-srb:

Lyman-alpha and Schumann-Runge indices
======================================

TUV-x parameterizes O\ :sub:`2` absorption in two spectral regions. Both
parameterizations require the exact wavelength grid edges listed below. The
indices refer to the standard wavelength grid in
``data/grids/wavelength/combined.grid``, which has 157 edges and 156 bins.

.. list-table::
   :header-rows: 1
   :widths: 25 15 20 20 20

   * - Region
     - Bins
     - Bin indices
     - Wavelength range
     - Source
   * - Lyman-alpha
     - 1
     - 2
     - 121.4 - 121.9 nm
     - ``src/la_sr_bands.F90`` ``wlla``
   * - Schumann-Runge
     - 17
     - 21 - 37
     - 175.4 - 206.2 nm
     - ``src/la_sr_bands.F90`` ``wlsrb``

TUV-x checks these edges at construction time. It reports
``Lyman alpha grid mis-match`` when the grid does not match. For wavelengths
below 205.8 nm, only the pre-specified grid is permitted.

TUV-x applies each check only when the host grid spans the region. It enables
the Lyman-alpha parameterization when the grid starts at or below 121.4 nm and
ends at or above 121.9 nm. It enables the Schumann-Runge parameterization when
the grid starts at or below 175.4 nm and ends at or above 206.2 nm. A grid that
starts above 175.4 nm therefore avoids both checks.

TUV-x then computes O\ :sub:`2` photolysis from the plain cross section on the
grid that remains. O\ :sub:`2` photolysis extends to about 242 nm, so the
reaction is still present, but the part of it that the two parameterizations
describe is gone.

.. _host-radiation-field-resolution:

Spectral resolution
===================

The radiation field must resolve the ultraviolet. This is the one requirement in
this document that TUV-x cannot check for you.

Between 300 nm and 340 nm, the product of the cross section and the quantum
yield for O\ :sub:`3` + h\ |nu| |rarr| O\ :sub:`2` + O(\ :sup:`1`\ D) falls by
more than three orders of magnitude, while ozone absorption makes the
transmission rise steeply across the same interval. A rate constant is the
integral of the product of the two. A flux averaged over an interval that wide
records nothing about where in the interval the photons sit, so the integral
that it produces is wrong by orders of magnitude, not by a few percent.

.. |nu| unicode:: U+03BD
.. |rarr| unicode:: U+2192
