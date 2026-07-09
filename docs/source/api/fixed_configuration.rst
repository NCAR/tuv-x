Example Fixed Configuration
===========================

``tuvx/fixed_configuration.hpp`` is a copyable example showing how specific
species are assembled from the general, species-agnostic transform forms. It is
the one place where particular parameter values, species identities, and
literature citations live -- the role a host model or MUSICA's configuration
layer would otherwise fill. The core transform library itself contains no
species names or hard-coded constants.

Cross Sections
--------------

.. doxygenfunction:: tuvx::fixed_configuration::rayleigh
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::hobr
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::t_butyl_nitrate
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::nitroxy_acetone
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::nitroxy_ethanol
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::n2o
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::cl2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::h2o2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::chbr3
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::hno3
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::rono2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::ch3ono2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::ch2o
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::cfc11
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::no2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::ccl4
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::chcl3
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::clono2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::acetone
   :project: tuvx

Spectral Weights
----------------

Dose-rate action spectra: unitless, wavelength-only weights applied to spectral
irradiance during dose-rate calculation.

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::standard_human_erythema
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::uv_index
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::scup_mice
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::exp_decay
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::par
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::gaussian
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::plant_damage
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::plant_damage_flint_caldwell
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::plant_damage_flint_caldwell_ext
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::spectral_weights::phytoplankton_boucher
   :project: tuvx

Quantum Yields
--------------

Unitless photodissociation branching fractions (typically 0-1), applied
alongside the cross-section in photolysis-rate calculation.

.. doxygenfunction:: tuvx::fixed_configuration::quantum_yields::clo_cl_o1d
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::quantum_yields::clo_cl_o3p
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::quantum_yields::clono2_cl_no3
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::quantum_yields::clono2_clo_no2
   :project: tuvx

.. doxygenfunction:: tuvx::fixed_configuration::quantum_yields::ho2_oh_o
   :project: tuvx
