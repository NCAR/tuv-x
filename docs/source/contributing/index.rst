.. _Contributing:

Contributing
============

For all proposed changes (bug fixes, new features, documentation updates, etc.)
please file an `issue <https://github.com/NCAR/tuv-x/issues/new/choose>`_
detailing what you intend to do.

The NSF NCAR software developers will work with you on that issue to help
recommend implementations, work through questions, or provide other background
knowledge.

Testing
-------

Any code you contribute must include associated unit tests and a regression or
integration test. Our tests run across various platforms automatically. New
features must include at least a brief example in the documentation.

Style guide
-----------

TUV-x uses `clang-format` and `clang-tidy` to enforce consistent C++ style.
All warnings are treated as errors. After each PR is merged, an automated
GitHub action runs `clang-format` over the codebase.

Building the documentation
--------------------------

All docs live in the ``docs`` directory. Install Python dependencies first.

Virtualenv
^^^^^^^^^^

.. code-block:: bash

   pip install virtualenv
   virtualenv tuvx_env
   source tuvx_env/bin/activate
   cd docs
   pip install -r requirements.txt
   cd .. && mkdir -p build && cd build
   cmake -DTUVX_BUILD_DOCS=ON ..
   cmake --build . --target docs
   open docs/sphinx/index.html

Venv
^^^^

.. code-block:: bash

   python -m venv tuvx_env
   source tuvx_env/bin/activate
   cd docs
   pip install -r requirements.txt
   cd .. && mkdir -p build && cd build
   cmake -DTUVX_BUILD_DOCS=ON ..
   cmake --build . --target docs
   open docs/sphinx/index.html
