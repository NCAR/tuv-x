.. Installation instructions for TUV-x

###################################
Getting Started
###################################

Installation
============

Docker
------

The quickest way to get started with TUV-x is with Docker.
The only requirement is that you have `Docker Desktop <https://www.docker.com/get-started>`_
installed and running.
With Docker Desktop running, open a terminal window.

To get the latest release of TUV-x, run the following command to start the TUV-x container:

.. code-block:: bash

   docker run -it ghcr.io/NCAR/tuv-x:release bash


To get the most recent, pre-release version of TUV-x instead run:

.. code-block:: bash

   docker run -it ghcr.io/NCAR/tuv-x:main bash


Inside the container, you can run the TUV-x tests from the ``/build/`` folder:

.. code-block:: bash

   cd build/
   make test


.. _install-local:

Local Installation
------------------

To build TUV-x locally, you will need to have:

- `NetCDF-Fortran <https://github.com/Unidata/netcdf-fortran>`_
- `JSON-Fortran <https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz>`_

If you plan to run the TUV-x tests, you will also need Python 3.x with numpy and scipy.

First, clone the TUV-x repo:

.. code-block:: bash

   git clone https://github.com/NCAR/tuv-x.git
   cd tuv-x


Alternatively, you can download the TUV-x source code for a particular release
`here <https://github.com/NCAR/tuv-x/releases>`_.

To build TUV-x, from the root TUV-x folder run:

.. code-block:: bash

   mkdir build
   cd build/
   export JSON_FORTRAN_HOME=/path/to/json-fortran/installation/
   cmake ..
   make -j 8
   make test

The ``JSON_FORTRAN_HOME`` environment variable should point to the root JSON-Fortran
installation folder, such that the jsonfortran library is located at
``JSON_FORTRAN_HOME/lib/libjsonfortran.so``.

Depending on where you have the NetCDF library installed, you may also
need to set the CMake variables for paths to their include folders and to the library files.

After verifying that the tests pass, you can instal tuv-x with ``make install``.
This will install the tuv-x standalone binary, the static library,
and some data files. The install location will be the default
for your system but you can change the location by modifying
``CMAKE_INSTALL_PREFIX`` during the cmake configure phase.

Assuming your prefix is ``/usr/local`` TUV-x will install into your file system
with these paths:

.. code-block:: bash

   :: tree -L 3 --filelimit=5 /usr/local/ncar 
   /usr/local/ncar
   └── tuvx-0.3.0
      ├── bin
      │   └── tuv-x
      ├── data  [7 entries exceeds filelimit, not opening dir]
      ├── examples
      │   └── full_config.json
      ├── include  [104 entries exceeds filelimit, not opening dir]
      └── lib
         └── libtuvx.a

To include TUV-x in your cmake project after installing TUV-x, you can
find TUV-x with ``find_package`` and add it to your target like this

.. code-block:: CMake

   find_package(tuvx REQUIRED)

   target_link_libraries(your_lib
      PUBLIC 
         musica::tuvx
   )


.. _install-mpi:

MPI Support
-----------

To build TUV-x with MPI support, first clone or download the TUV-x source code as
described in :ref:`install-local`.
Then, from the root TUV-x folder, you can build a Docker image with MPI support and
run it in a container:

.. code-block:: bash

   docker build -t tuv-x . -f Dockerfile.mpi
   docker run -it tuv-x bash
   make test


Alternatively, you can follow the instructions in :ref:`install-local`, replacing
the call to cmake with:

.. code-block:: bash

   cmake -D CMAKE_Fortran_COMPILER=/path/to/mpif90 \
         -D ENABLE_MPI:BOOL=TRUE \
         ..


You should replace ``path/to/mpif90`` with the path to your local Fortran MPI compiler.
