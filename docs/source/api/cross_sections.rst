Cross Sections
==============

Library of pre-built cross-section ``TransformFunc`` factories, ported from the
Fortran ``src/cross_sections/`` algorithms.  Each factory returns a transform
that calculates absorption cross-section weights in m\ :sup:`2` (SI).

Analytic Cross Sections
-----------------------

Temperature- and height-independent cross-sections defined by a closed-form
analytic function of wavelength alone.  The result is broadcast over all height
levels and columns.

.. doxygenfunction:: tuvx::cross_sections::rayleigh
   :project: tuvx

.. doxygenfunction:: tuvx::cross_sections::hobr
   :project: tuvx

.. doxygenfunction:: tuvx::cross_sections::t_butyl_nitrate
   :project: tuvx

.. doxygenfunction:: tuvx::cross_sections::nitroxy_acetone
   :project: tuvx

.. doxygenfunction:: tuvx::cross_sections::nitroxy_ethanol
   :project: tuvx
