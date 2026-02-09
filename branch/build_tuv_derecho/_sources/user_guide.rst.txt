.. Usage information for TUV-x

###################################
User Guide
###################################

General
=======

Stand-alone TUV-x can be run from the command-line as:

.. code-block:: bash

   ./tuv-x configuration_file.json

The ``configuration_file.json`` contains the TUV-x configuration information described in
:ref:`configuration <configuration>`.

If photolysis rate constants are included in the configuration,
TUV-x will output a file named ``photolysis_rate_constants.nc`` in the working directory. This
file will contain the photolysis rate constants for each reaction at each vertical level
for every time specified in the configuration file.

If dose rates are included in the configuration,
TUV-x will output a file named ``dose_rates.nc`` in the working directory.
This file will contain the calculated dose rates at
each vertical level for every time specified in the configuration.

Some example configurations are available in the |examples|_ folder.
These will be copied to the build folder during the CMake build step.
To run an example, from the build folder run:

.. |examples| replace:: ``examples`` 
.. _examples: https://github.com/NCAR/tuv-x/tree/main/examples

.. code-block:: bash

   ./tuv-x examples/full_config.json

.. _configuration:

Configuration
=============

The TUV-x configuration
`JSON <https://www.json.org/json-en.html>`_
file is used to set up TUV-x to calculate
a set of user-selected photolysis rate constants and/or dose rates.
It also allows you to specify details about how these
are calculated (solver options, which absorbing species should be
considered in determining the radiation field, etc.) and the
conditions they should be calculated for (temperature, solar flux,
solar zenith angle, etc.).

Configuration data for specific components of TUV-x may
include paths to needed data files. These paths can be absolute
or relative to the folder in which you run the ``tuvx`` executable.
(There are no hard-coded paths to data files in the TUV-x source
code; they are all specified in the JSON configuration file.)

TUV-x JSON data include a number of required and optional
elements.
Users can additionally decorate the JSON file with
custom elements by using a "``__``" prefix for the custom key.
TUV-x will ignore all key-value pairs whose key includes this
prefix.
These can be useful for including notes, descriptions,
or references.
They can also be useful if TUV-x is embedded
in an application that needs to attach specific information
to TUV-x objects.

.. code-block:: JSON

   {
     "tuv object": {
       "a TUV-x required key": "foo",
       "a TUV-x optional key": "sometimes bar",
       "__my note": "TUV-x will ignore this",
       "__user provided data" : {
         "will TUV-x ignore this object?": true
       }
     }
   }


The TUV-x configuration JSON file has six objects:

.. code-block:: JSON
   :force:

   {
     "O2 absorption": { ... },
     "grids": [ ... ],
     "profiles": [ ... ],
     "radiative transfer": { ... },
     "photolysis": { ... },
     "dose rates": { ... }
     "enable diagnostics" : false,
   }


We use elipses (``...``) here and throughout this section for
JSON data that are described in a separate sub-section of this
page.

The first four objects
(:ref:`O2 absorption <configuration-o2-absorption>`,
:ref:`grids <configuration-grids>`,
:ref:`profiles <configuration-profiles>`, and
:ref:`radiative transfer <configuration-radiation>`)
are required for every TUV-x run.
The last two objects
(:ref:`photolysis rate constants <configuration-photolysis>` and
:ref:`dose rates <configuration-dose-rates>`)
are optional and allow the user to
calculate photolysis rates or dose rates or both.
Finally, the ``enable diagnostics`` field is used to output diagnostics from
the core of tuv-x. If set to true, a folder called output will be created. 
This flag is optional and defaults to false.

The following sections describe each of these six JSON
object.

.. _configuration-o2-absorption:

O2 Absorption
-------------

This object describes how TUV-x calculates the special
Lyman--Alpha and Shumann--Runge absorption bands for O2.
It has a single required data member that provides the
path to a text file containing parameters needed in
these calculations:

.. code-block:: JSON

   "O2 absorption": {
     "cross section parameters file": "data/cross_sections/O2_parameters.txt"
   }


This file is a hold-over from the original TUV and has
not yet been converted into a NetCDF file. It can be found
in the location shown above relative to the ``tuv-x/`` root
directory.


.. _configuration-grids:

Grids
-----

Grids define axes along which TUV-x data are distributed or calculations
are performed.
Depending on your use case, certain grids may or may not
need to be included.
At minimum, you will need to include a ``height``
grid (units: ``km``) that defines the vertical grid TUV-x operates on,
and a ``wavelength`` grid (units: ``nm``), that defines the wavelength
bins for which optical property data is provided.
For stand-alone TUV-x, you will also need to provide a ``time`` grid
(units: ``hours``) that provides the time of day for which to calculate
the photolysis and/or dose rates.

There are three key-value pairs that are required for every grid, ``type``
, ``units``, and ``name``.
The ``type`` is a string that describes how the grid data is specified,
and must be one of ``equal interval``,
``from csv file``, and ``from config file``.
These types are described below.
The ``units`` are the units for the grid values and are used to ensure
consistency throughout configuration inputs and
with TUV conventions.

Equal Interval
^^^^^^^^^^^^^^

This grid type is defined by start and end points, and a cell width:


.. code-block:: JSON

   {
     "name": "my grid"
     "type": "equal interval",
     "units": "foos",
     "begins at": 12.0,
     "ends at": 112.0,
     "cell delta": 10.0
   }


From CSV File
^^^^^^^^^^^^^

This grid type is defined by data from a text file.
The first line of the file is ignored and can be used for header info.
The remaining lines should each include a single real number
that will be used to define the grid.
Values should be monotonic and increasing.


.. code-block:: JSON

   {
     "name": "my grid"
     "type": "from csv file",
     "units": "bars",
     "file path": "path/to/file"
   }


From Config File
^^^^^^^^^^^^^^^^

The values for each cell of this grid type are included directly in
the configuration file.
Values should be monotonically increasing.


.. code-block:: JSON

   {
     "name": "my grid"
     "type": "from config file",
     "units": "foos",
     "values": [ 1.0, 7.0, 12.5 ]
   }


.. _configuration-profiles:

Profiles
--------

Profiles define parameters on a :ref:`grid <configuration-grids>`.
Similar to grids, depending on your use case, certain profiles may
or may not need to be included.
At minimum, you will need to include the key-value pairs listed
in the following table.

=========================  ==============
profile                    grid
=========================  ==============
``temperature``            ``height``
``solar zenith angle``     ``time``
``Earth-Sun distance``     ``time``
``surface albedo``         ``wavelength``
``extraterrestrial flux``  ``wavelength``
``air``                    ``height``
=========================  ==============

The air profile defines the number density (molecules cm\ :sup:`-3`\ ) of air
as a function of height.

In addition to the required profiles, profiles of atmospheric
constituents that affect the radiation field by absorbing or
scattering light should be included when these exist as
:ref:`radiators <configuration-radiators>` in the
:ref:`radiative transfer <configuration-radiation>` object.

As with grids, there are two key-value pairs required for every
profile, ``type``, ``units``, and ``name``.
The ``units`` are the units for the profile data and are used to ensure
consistency throughout the configuration data and
with expectations in the code.
The ``type`` is a string that describes how the profile data is specified.
There are two general-use types that can be used for any profile,
``from csv file`` and ``from config file``.
In addition, there are several profile types that are useful for describing
specific types of profiles, ``O2``, ``O3``, ``air``, ``solar zenith angle``,
``Earth-Sun distance``, and ``extraterrestrial flux``.

The specific profile types are described below.

.. _configuration-profiles-from-csv:

From CSV file
^^^^^^^^^^^^^

This profile type loads profile data from a text file where the
data is expected to be space-separated with the first column being
the grid-point value and the second column being the value of the
profile at that grid-point.
Any number of header lines can be included at the top of the file
by prefixing them with any of ``#!$%*``.

=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``file path``              required
``name``                   required
``interpolator``           optional
``scale height``           optional
=========================  ==============

The format for the profile is:

.. code-block:: JSON

   {
     "name": "my profile"
     "type": "from csv file",
     "units": "foos",
     "file path": "path/to/file",
     "grid": {
       "name": "bar",
       "units": "bazes"
     }
   }

An optional ``interpolator`` key can be used to specify the
interpolation strategy to apply to the incoming data if needed.
The possible values are:

- ``linear``
  Standard linear interpolation scheme. This is the default
  interpolator used when one is not specified in the configuration
  data.

- ``conserving``
  Linear interpolation that conserves the area under the curve.

- ``fractional source``
  Interpolation scheme for profiles that are integrated values.
  The interpolation is based on the fractional overlap between
  source and target grid sections relative to the source
  grid width.

- ``fractional target``
  Interpolation scheme for profiles that are integrated values.
  The interpolation is based on the fractional overlap between
  source and target grid sections relative to the target
  grid width.


.. _configuration-profiles-from-config:

From Config File
^^^^^^^^^^^^^^^^

This profile type extracts profile data directly from the configuration
data file. There are two options for this profile format. Each option will 
require that ``type`` and ``units`` be present in the configuration file. 

The ``grid`` option is required. In the case of ``values``, the grid ensures
that the number of elements in ``values`` matches the length of the grid.
When using the ``uniform value`` option, the ``edge_val``, or the edges of the 
cells, are determined from the grid.

=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``grid``                   required
``name``                   required
``values``                 optional
``uniform value``          optional
=========================  ==============

The first option specifies a single uniform value at every grid point:

.. code-block:: JSON

   {
     "name": "my profile"
     "type": "from config file",
     "units": "bars",
     "uniform value": 12.3,
     "grid": {
       "name": "foo",
       "units": "bazes"
     }
   }


The second option specifies values at each grid point in an array:


.. code-block:: JSON

     {
       "name": "my profile"
       "type": "from config file",
       "units": "bars",
       "values": [ 12.3, 32.4, 103.2 ],
       "grid": {
         "name": "foo",
         "units": "bazes"
       }
   }

.. _configuration-profiles-from-host:

From Host
^^^^^^^^^

.. todo:: fill this out

Solar Zenith Angle
^^^^^^^^^^^^^^^^^^

This profile is specifically for calculating the solar zenith angle
as a function of time. 

=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``year``                   required
``month``                  required
``day``                    required
``longitude``              required
``latitude``               required
``name``                   required
``time zone``              optional
=========================  ==============

Its configuration takes the form:

.. code-block:: JSON

   {
     "name": "solar zenith angle",
     "type": "solar zenith angle",
     "units": "degrees",
     "year" : 2002,
     "month": 3,
     "day": 21,
     "longitude": 0.0,
     "latitude": 0.0
   }


The latitude and longitude are in degrees. There is an optional
argument ``time zone`` that defaults to 0, and is an offset
in hours to adjust for a specific time zone relative to GMT.


Earth-Sun Distance
^^^^^^^^^^^^^^^^^^

This profile is specifically for calculating the Earth-Sun distance
as a function of time. This profile requires that a grid named ``time`` be
defined in the :ref:`configuration-grids` section.

=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``year``                   required
``month``                  required
``day``                    required
``name``                   required
``time zone``              optional
=========================  ==============

Its configuration takes the form:

.. code-block:: JSON

     {
       "name": "Earth-Sun distance",
       "type": "Earth-Sun distance",
       "units": "AU",
       "year" : 2002,
       "month": 3,
       "day": 21
     }

There is an optional argument ``time zone`` that defaults to 0,
and is an offset in hours to adjust for a specific time zone
relative to GMT.


Other Profiles
^^^^^^^^^^^^^^

The remaining profile types support legacy code from the original
implementation of TUV. These will be gradually removed as configurations
employing the general-use profile types are developed.

These can be used as follows (note that file paths are relative
to the root ``tuv-x/`` folder):


.. code-block:: JSON

      {
         "name": "O3",
         "type": "O3",
         "units": "molecule cm-3",
         "file path": "data/profiles/atmosphere/ussa.ozone"
      },
      {
         "name": "air",
         "type": "air",
         "units": "molecule cm-3",
         "file path": "data/profiles/atmosphere/ussa.dens"
      },
      {
         "name": "O2",
         "type": "O2",
         "units": "molecule cm-3",
         "file path": "data/profiles/atmosphere/ussa.dens"
      },
      {
         "name": "extraterrestrial flux",
         "type": "extraterrestrial flux",
         "units": "photon cm-2 s-1",
         "file path": ["data/profiles/solar/susim_hi.flx",
                      "data/profiles/solar/atlas3_1994_317_a.dat",
                      "data/profiles/solar/sao2010.solref.converted",
                      "data/profiles/solar/neckel.flx"],
         "interpolator": ["","","","fractional target"]
      }

Air Keys
~~~~~~~~
=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``file path``              required
``name``                   required
=========================  ==============

Extraterrestrial Flux Keys
~~~~~~~~~~~~~~~~~~~~~~~~~~
=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``file path``              required
``interpolator``           required
``name``                   required
``enable diagnostics``     optional
=========================  ==============

The regressoin tests compare the new version of TUV-x to the old version. One
way is by directly comparing output. The `enable diagnostics` allows for this
ouptut to be disabled. If this is enabled, a folder named `output` will be 
created in the same directory TUV-x is run from.

O2 Keys
~~~~~~~
=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``file path``              required
``name``                   required
``interpolator``           optional
``scale height``           optional
=========================  ==============

O3 Keys
~~~~~~~
=========================  ==============
keys                       Required/Optional
=========================  ==============
``type``                   required
``units``                  required 
``file path``              required
``name``                   required
``interpolator``           optional
``scale height``           optional
``reference column``       optional
=========================  ==============

.. _configuration-radiation:

Radiative Transfer
------------------

Solvers
^^^^^^^

Radiative transfer specifies how the radiation field is calculated.
The general format for radiative transfer is:

.. code-block:: JSON

   "radiative transfer": {
     "cross sections": [
       {
         "name": "foo",
         "type": "base",
         "netcdf files": [ "my/data/file.nc" ]
       }
     ],
     "radiators": [
       {
         "name": "foo",
         "type": "base",
         "cross section": "foo",
         "vertical profile": "foo",
         "vertical profile units": "molecule cm-3"
       }
     ],
     "solver": {
       "type": "delta eddington"
     }
   }

The ``cross sections`` and ``radiators`` are required arrays
(even if they are of zero length).
The ``cross sections`` define the absorption cross sections
for the radiators.
Cross section configuration formats are described
:ref:`below <configuration-cross-sections>`.

The ``radiators`` define the atmospheric constituents that
should be considered in the calculation of the radiation
field.
Radiator configuration formats are described
:ref:`below <configuration-radiators>`.

The ``solver`` key is required and specifies the solver to
use for radiative transfer.
There are currently two options, described below.

Delta Eddington
~~~~~~~~~~~~~~~

This is a fast 2-stream solver based on:
Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.
The configuration format for the delta-Eddington solver is:


.. code-block:: JSON

   "type" : "delta eddington"


There are no configuration options for this solver.

Discrete Ordinate
~~~~~~~~~~~~~~~~~

This is an ``n``-stream solver, where ``n`` is an
even number between 2 and 32.
The configuration format for the discrete ordinate solver is:


.. code-block:: JSON

   "type": "discrete ordinate",
   "number of streams": 4


The ``number of streams`` is required and specifies the number
of streams to solve for.
The value must be an even number between 2 and 32.

.. _configuration-radiators:

Radiators
^^^^^^^^^

Radiators represent atmospheric constituents that attenuate
solar radiation and will be considered in calculations of the
radiation field.

Base
~~~~

The generic configuration data format for standard radiators is as
follows:

.. code-block:: JSON

   {
     "type": "base",
     "name": "foo",
     "treat as air": true,
     "cross section": "foo",
     "vertical profile": "foo",
     "vertical profile units": "molecule cm-3",
     "enable diagnostics": false
   }
  

===========================  ==============
keys                         Required/Optional
===========================  ==============
``name``                     optional
``type``                     required
``treat as air``             optional
``cross section``            required 
``vertical profile``         required
``vertical profile units``   required
``enable diagnostics``       optional
===========================  ==============

The regression tests compare the new version of TUV-x to the old version. One
way is by directly comparing output. The ``enable diagnostics`` allows for this
ouptut to be disabled. If this is enabled, a folder named ``output`` will be 
created in the same directory TUV-x is run from.

The ``treat as air`` flag can be used to indicate that the radiator
should be treated in a unique way specific to air in the calculation
of optical properties.
The default value is false.

The ``cross section`` must be the name of a cross
section in the list of ``cross sections`` in the
:ref:`configuration-radiation` object.
The vertical profile must be the name of a profile with
the provided units that is present in the list of
:ref:`configuration-profiles`.
These profiles should describe the concentration of
the constituent on the ``height`` grid.

Aerosols
~~~~~~~~
A special radiator type exists for aerosols, which
provides fixed optical depths at each wavelength.
An example of the aerosol configuration is provided
below.


.. code-block:: JSON

   {
     "name": "aerosols",
     "type": "aerosol",
     "optical depths": [2.40e-01, 1.06e-01, 4.56e-02, 1.91e-02, 1.01e-02, 7.63e-03,
                        5.38e-03, 5.00e-03, 5.15e-03, 4.94e-03, 4.82e-03, 4.51e-03,
                        4.74e-03, 4.37e-03, 4.28e-03, 4.03e-03, 3.83e-03, 3.78e-03,
                        3.88e-03, 3.08e-03, 2.26e-03, 1.64e-03, 1.23e-03, 9.45e-04,
                        7.49e-04, 6.30e-04, 5.50e-04, 4.21e-04, 3.22e-04, 2.48e-04,
                        1.90e-04, 1.45e-04, 1.11e-04, 8.51e-05, 6.52e-05, 5.00e-05,
                        3.83e-05, 2.93e-05, 2.25e-05, 1.72e-05, 1.32e-05, 1.01e-05,
                        7.72e-06, 5.91e-06, 4.53e-06, 3.46e-06, 2.66e-06, 2.04e-06,
                        1.56e-06, 1.19e-06, 9.14e-07],
     "single scattering albedo": 0.99,
     "asymmetry factor": 0.61,
     "550 nm optical depth": 0.235,
     "enable diagnostics": false
   }


============================    ==============
keys                            Required/Optional
============================    ==============
``name``                        optional
``type``                        required
``optical depths``              required 
``single scattering albdeo``    required
``asymmetry factor``            required
``550 nm optical depth``        optional
``enable diagnostics``          optional
============================    ==============

The regressoin tests compare the new version of TUV-x to the old version. One
way is by directly comparing output. The `enable diagnostics` allows for this
ouptut to be disabled. If this is enabled, a folder named `output` will be 
created in the same directory TUV-x is run from.

The optical depths are expected to be on the ``wavelength``
grid.

.. _configuration-radiators-from-host:

From Host
~~~~~~~~~

.. todo:: fill this out

.. _configuration-photolysis:

Photolysis Reactions
--------------------

The configuration for photolysis reactions takes the following form:


.. code-block:: JSON
   :force:

   "photolysis": {
     "enable diagnostics": false,
     "reactions": [
      {
        "name": "my first reaction",
        "cross section": { ... },
        "quantum yield": { ... },
        "scaling factor": 1.3
      },
      {
        "name": "my second reaction",
        "cross section": { ... },
        "quantum yield": { ... },
      }
     ]
   }


Each member of ``reactions`` describes a photolysis
reaction that TUV-x will calculate a rate constant for at runtime.
The name for each photolysis reaction is user-defined, 
associated with the calculated rate constant,
and can be used for mapping to a chemistry solver or other
package.
Each reaction must have a ``cross section``, whose configuration
format is described :ref:`here <configuration-cross-sections>`.
Each reaction must also have a ``quantum yield``, whose
configuration is described :ref:`here <configuration-quantum-yields>`.
The ``scaling factor`` is a optional scaling factor that will be
applied to the calculated rate constant.

Diagnostic output can be enabled by setting ``enable diagnostics`` to ``true``.
This keyword is not required and is ``false`` by default. 
When enabled,  a folder named
`output` will be created with some diagnostic output for the cross sections
and quantum yields. This is only used for regression tests and will be removed
in the future.

The file ``data/photolysis_rate_constants.json`` contains
configuration data for every photolysis rate constant that
can be calculated from data available in the ``data/``
folder.

.. _configuration-dose-rates:

Dose Rates
----------

The configuration for dose rates takes the following form:


.. code-block:: JSON
   :force:

   "dose rates": {
     "my first dose rate": {
       "weights": { ... }
     },
     "my second dose rate": {
       "weights": { ... }
     }
   }


Each member of ``dose rates`` describes a dose rate that
TUV-x will calculate at runtime.
The key for each dose rate is a user-defined name that
will be associated with the calculated rate.
Each dose rate must have a ``weights`` object that
defines a spectral weight, whose configuration format
is described :ref:`here <configuration-spectral-weights>`.

Additionally, diagnostic output can be enabled by adding ``enable diagnostics`` to the
json configuration like in the sample below. In this case, a folder named
`output` will be created with some diagnostic output for the cross sections
and quantum yields. This is only used for regression tests and will be removed
in the future.


.. code-block:: JSON
   :force:

   "dose rates": {
     "enable diagnostics" : true,
     "my first dose rate": {
       "weights": { ... }
     },
     "my second dose rate": {
       "weights": { ... }
     }
   }

Additional Objects
------------------

These configuration JSON objects are used in one or more
of the six high-level JSON objects in the TUV-x data file.


.. _configuration-cross-sections:

Cross Sections
^^^^^^^^^^^^^^

The configuration for a standard cross section is as
follows:

.. code-block:: JSON

   {
     "type": "base",
     "netcdf files": [
       {
         "file path": "path/to/my/netcdf/file.nc",
         "lower extrapolation": {
           "type": "boundary"
         },
         "upper extrapotion": {
           "type": "constant",
           "value": 0.0
         },
         "zero below": 215.4,
         "zero above": 768.4
       }
     ],
     "override bands": [
       {
         "band": "lyman-alpha",
         "value": 0.32
       }
     ],
     "apply O2 bands": true
   }


The NetCDF file should be structured as follows::

   dimensions:
	    bins = NUMBER_OF_WAVELENGTH_BINS ;
	    parameters = 1 ;
   variables:
	    double wavelength(bins) ;
		     wavelength:units = "nm" ;
	    double cross_section_parameters(parameters, bins) ;


Here, ``NUMBER_OF_WAVELENGTH_BINS`` is the number of wavelength
bins the cross section values are provided on.

The ``cross_section_parameters`` array should hold the value
of the cross section at each wavelength.
TUV-x will interpolate the cross section
data onto the TUV-x wavelength grid, as specified by the
"wavelength" grid.

For objects in the ``netcdf files`` array, the ``file path``
is required. The remaining keys are optional.

The ``lower extrapolation`` and ``upper extrapolation``
keys will append data beyond the lower and upper limits
of the input data set, respectively. When the ``type``
is ``boundary``, the last data point in the data set is
used for all points beyond the last point. When the
``type`` is ``constant`` a user-specified ``value`` is
used for all points beyond the last point.

The ``zero above`` and ``zero below`` keys are used to
zero all points from the file above or below the
specified wavelength values AFTER INTERPOLATION.

An ``override bands`` array can be included to specify
fixed values over specific bands. These values will be
applied immediately before returning the cross section
and thus override any calculated or interpolated values.
Valid ``band`` values are ``lyman-alpha``,
``schumann-runge``, and ``schumann-runge continuum``.
The wavelength grid must cover the Lyman-alpha and
Schumann-Runge bands if an override is specified,
otherwise a configuration error is returned.

The ``apply O2 bands`` flag can be added to photolysis
rates that dissociate ``O2``. If this flag is set to
``true``, cross section values in the Lyman-alpha and
Schumann-Runge wavelength bands will be overwritten with
TUV-x's custom cross section values for O2 in these bands.

A number of custom cross section types have been developed
when more complex algorithms are needed to calculate
cross sections.
These generally apply to a specific photolysis reaction.
Their configuration data formats are demonstrated
in ``data/photolysis_rate_constants.json``.

.. _configuration-quantum-yields:

Quantum Yields
^^^^^^^^^^^^^^

The configuration for a standard quantum yield can be of two forms.
The first describes a quantum yield that is a constant value at
all wavelengths:


.. code-block:: JSON

   {
     "type": "base",
     "constant value": 1.0
   }


The second describes a quantum yield that is read from a NetCDF
file:


.. code-block:: JSON

   {
     "type": "base",
     "netcdf files": [ "path/to/my/netcdf/file.nc" ]
   }


The NetCDF file should be structured as follows::

   dimensions:
	   bins = NUMBER_OF_WAVELENGTH_BINS ;
	   parameters = 1 ;
   variables:
	   double wavelength(bins) ;
		   wavelength:units = "nm" ;
	   double quantum_yield_parameters(parameters, bins) ;
		   quantum_yield_parameters:units = "fraction" ;


Here, ``NUMBER_OF_WAVELENGTH_BINS`` is the number of wavelength
bins the quantum yield values are provided on.

The ``quantum_yield_parameters`` array should hold the value
of the quantum yield at each wavelength.
TUV-x will perform interpolation of the quantum yield data
to the native wavelength grid.

For either type of quantum yield, an ``override bands``
array can be included:


.. code-block:: JSON

   {
     "type": "base",
     "constant value": 1.0,
     "override bands": [
       {
         "band": "lyman-alpha",
         "value": 0.32
       }
     ]
   }


An ``override bands`` array instructs TUV-x to apply
fixed values over specific bands. These values will be
applied immediately before returning the quantum yield
and thus override any calculated or interpolated values.
Valid ``band`` values are ``lyman-alpha``,
``schumann-runge``, and ``schumann-runge continuum``.
The wavelength grid must cover the Lyman-alpha and
Schumann-Runge bands if an override is specified,
otherwise a configuration error is returned.

A number of custom quantum yield types have been developed
when more complex algorithms are needed to calculate
quantum yields.
These generally apply to a specific photolysis reaction.
Their configuration data formats are demonstrated
in ``data/photolysis_rate_constants.json``.


.. _configuration-spectral-weights:

Spectral Weights
^^^^^^^^^^^^^^^^

Standard spectral weight types load spectral weights from a NetCDF
file.
Their configuration format is as follows:


.. code-block:: JSON

   {
     "type": "base",
     "netcdf files": [ "path/to/my/netcdf/file.nc" ]
   }

The NetCDF file should be structured as follows::

   dimensions:
	   bins = NUMBER_OF_WAVELENGTH_BINS ;
	   parameters = 1 ;
   variables:
	   double wavelength(bins) ;
		   wavelength:units = "nm" ;
	   double spectral_weight_parameters(parameters, bins) ;
		   spectral_weight_parameters:hdr = "" ;


Here, ``NUMBER_OF_WAVELENGTH_BINS`` is the number of wavelength
bins the spectral weight values are provided on.

The ``spectral_weight_parameters`` array should hold the value
of the spectral weight at each wavelength.
TUV-x will perform interpolation of the spectral weight  data
to the native wavelength grid.

The ``type`` and ``netcdf files`` keys are required. There are
also two optional key-value pairs: ``lower extrapolation`` and
``upper extrapolation``.
These both have the same structure.
An example of the ``lower extrapolation`` follows.

.. code-block:: JSON

   "lower extrapolation": {
     "type": "boundary"
   }


The value of ``type`` can be ``boundary`` or ``constant``.
If ``boundary`` is selected, the value at the lower (or upper)
boundary of the input data will be extended to the
lower (or upper) extent of the TUV-x wavelength grid.
If ``constant`` is selected, a key-value pair ``value``
must also be present:

.. code-block:: JSON

   "lower extrapolation": {
     "type": "constant",
     "value": 12.2
   }


This value will be used between the lower (or upper)
boundaries of the input data and TUV-x wavelength grids.
If no lower or upper extrapolation is specified, the
values between the input data and TUV-x wavelength grids
will be 0.


A number of custom spectral weight types have been developed
when more complex algorithms are needed to calculate
spectral weights.
Their configuration data formats are demonstrated
in ``data/dose_rates.json``.
