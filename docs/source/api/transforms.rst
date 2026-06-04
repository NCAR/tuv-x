Transforms
==========

The transform system provides a composable API for calculating per-wavelength,
per-height weight arrays used in photolysis rate calculation.  A
``TransformFunc`` is any callable that fills a weight array
``[wavelength x height x column]`` given the current atmospheric state.

Type System
-----------

.. doxygentypedef:: tuvx::TransformFunc
   :project: tuvx

.. doxygenstruct:: tuvx::PiecewiseRegion
   :project: tuvx
   :members:

.. doxygenstruct:: tuvx::SpectralBand
   :project: tuvx
   :members:

Factory Functions
-----------------

.. doxygenfunction:: tuvx::constant
   :project: tuvx

.. doxygenfunction:: tuvx::wrap_analytic
   :project: tuvx

.. doxygenfunction:: tuvx::from_data
   :project: tuvx

.. doxygenfunction:: tuvx::temperature_interpolation
   :project: tuvx

.. doxygenfunction:: tuvx::polynomial_scaling
   :project: tuvx

.. doxygenfunction:: tuvx::exponential_scaling
   :project: tuvx

.. doxygenfunction:: tuvx::linear_correction
   :project: tuvx

.. doxygenfunction:: tuvx::stern_volmer
   :project: tuvx

.. doxygenfunction:: tuvx::parameterized
   :project: tuvx

Analytic Forms
--------------

General, configurable wavelength-only shapes that recur across many
cross-sections, quantum yields, and spectral weights. All parameters are
caller-supplied; no species or constant is hard-coded.

.. doxygenstruct:: tuvx::LogNormalBand
   :project: tuvx
   :members:

.. doxygenstruct:: tuvx::LogNormalBandsParameters
   :project: tuvx
   :members:

.. doxygenfunction:: tuvx::log_normal_bands
   :project: tuvx

.. doxygenstruct:: tuvx::ExpPolynomialParameters
   :project: tuvx
   :members:

.. doxygenfunction:: tuvx::exp_polynomial
   :project: tuvx

Combinators
-----------

.. doxygenfunction:: tuvx::in_region
   :project: tuvx

.. doxygenfunction:: tuvx::multiply
   :project: tuvx

.. doxygenfunction:: tuvx::add
   :project: tuvx

.. doxygenfunction:: tuvx::piecewise
   :project: tuvx

.. doxygenfunction:: tuvx::clamp
   :project: tuvx

.. doxygenfunction:: tuvx::override_band
   :project: tuvx

.. doxygenfunction:: tuvx::scale
   :project: tuvx

.. doxygenfunction:: tuvx::named_band
   :project: tuvx
