.. TUV-x documentation master file, created by
   sphinx-quickstart on Fri May 20 16:12:10 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. tuv documentation HTML titles
..
.. # (over and under) for module headings
.. = for sections
.. - for subsections
.. ^ for subsubsections
.. ~ for subsubsubsections
.. " for paragraphs

###################################
Welcome to the TUV-x documentation!
###################################

.. todolist::

**TUV-x** (tropospheric ultraviolet-extended) is a Fortran library for 
calculating photolysis rate constants and dose rates. 
TUV-x provides a suite of algorithms for calculating photolysis
rate constants and dose rates that users can apply to their systems of
interest through simple configuration options.

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: Getting started
        :img-top: _static/index_getting_started.svg
        :link: getting_started
        :link-type: doc

        Check out the getting started guide to build and install TUV-x.

    .. grid-item-card::  User guide
        :img-top: _static/index_user_guide.svg
        :link: user_guide
        :link-type: doc

        Using TUV-x is quite striatforward. Check it out here!

    .. grid-item-card::  API reference
        :img-top: _static/index_api.svg
        :link: api/index
        :link-type: doc

        The source code for TUV-x is heavily documented.

    .. grid-item-card::  Contributors guide
        :img-top: _static/index_contribute.svg
        :link: contributing/index
        :link-type: doc

        If you'd like to contribute some new science code or update the docs,
        checkout the contributors guide!

.. note::

   This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   user_guide
   api/index
   contributing/index
   citation

