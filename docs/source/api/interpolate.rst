Interpolation
=============

Grid interpolation primitives (``tuvx/interpolate.hpp``).  These resample
tabulated data from a source axis onto a target axis -- for example, matching a
molecular cross section, quantum yield, or action spectrum onto the model
wavelength grid.  They are pure numerics (arrays in, arrays out) and do not
extrapolate; use :cpp:func:`tuvx::add_point` to pad the source so it spans the
target.

``interpolate_linear`` treats both axes as discrete points and returns one value
per target point.  The three area/overlap methods treat the axes as bin edges
and return one value per target bin (length ``x_target.Size() - 1``).

.. doxygenfunction:: tuvx::interpolate_linear
   :project: tuvx

.. doxygenfunction:: tuvx::interpolate_conserving
   :project: tuvx

.. doxygenfunction:: tuvx::interpolate_fractional_source
   :project: tuvx

.. doxygenfunction:: tuvx::interpolate_fractional_target
   :project: tuvx

.. doxygenfunction:: tuvx::add_point
   :project: tuvx
