TUV-x
=====

TUV-x is a photolysis rate constant calculator being rewritten from Fortran to C++.
Documentation will be expanded as the API stabilizes through each implementation phase.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   api/index

Current C++ Components
----------------------

**Data Structures**

- ``Grid`` — wavelength and altitude grids with midpoint/edge access
- ``Profile`` — 1D atmospheric profiles (temperature, density, etc.)
- ``Array1D``, ``Array2D``, ``Array3D`` — typed multi-dimensional arrays

**Linear Algebra**

- ``TridiagonalMatrix`` — tridiagonal system storage and Thomas algorithm solver

**Radiative Transfer**

- ``ConstituentState`` — optical properties (optical depth, single-scattering albedo, asymmetry parameter)
- ``RadiationField`` — actinic flux, irradiance, and heating rate storage
- ``DeltaEddington`` — two-stream radiative transfer solver (Toon et al., 1989)
