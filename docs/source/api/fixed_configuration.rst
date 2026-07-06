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
