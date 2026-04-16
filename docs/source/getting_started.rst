Getting Started
===============

Requirements
------------

- CMake 3.21 or newer
- A C++20 compiler (GCC 12+, Clang 17+, Apple Clang, MSVC 2022+)

Building
--------

.. code-block:: bash

   git clone https://github.com/NCAR/tuv-x.git
   cd tuv-x
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   ctest --test-dir build --output-on-failure

Using TUV-x as a Library
------------------------

TUV-x is exposed via the ``musica::tuvx`` CMake target:

.. code-block:: cmake

   find_package(tuvx REQUIRED)
   target_link_libraries(my_app PRIVATE musica::tuvx)

.. code-block:: cpp

   #include <tuvx/grid.hpp>
   #include <tuvx/profile.hpp>
   #include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>
