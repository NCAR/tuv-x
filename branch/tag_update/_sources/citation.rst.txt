.. how to cite TUV-x

############
Citing TUV-x
############

TUV-x is a software library, but the software library is based off of
scientific work. Therefore, we ask that you cite both the scientific publication
detailing the science behing TUV-x as well as the software.

Citing the Science
==================

The original publication that describes the science should be cited. You can read
the work `here <https://doi.org/10.1007/978-3-540-69044-3_1>`_.

.. code-block:: Tex
    :caption: Recommended citation: Madronich, Sasha, and Siri Flocke (1999), The role of solar radiation in atmospheric chemistry, in Handbook of Environmental Chemistry, edited by P. Boule, pp. 1-26, Springer-Verlag, Heidelberg

      @incollection{madronich_role_1999,
        address = {Berlin, Heidelberg},
        series = {The {Handbook} of {Environmental} {Chemistry}},
        title = {The {Role} of {Solar} {Radiation} in {Atmospheric} {Chemistry}},
        isbn = {978-3-540-69044-3},
        url = {https://doi.org/10.1007/978-3-540-69044-3_1},
        language = {en},
        booktitle = {Environmental {Photochemistry}},
        publisher = {Springer},
        author = {Madronich, Sasha and Flocke, Siri},
        editor = {Boule, Pierre},
        year = {1999},
        doi = {10.1007/978-3-540-69044-3_1},
        keywords = {Earth-Sun geometry., photolysis rate coefficients, radiative transfer, solar radiation, spectral actinic flux},
        pages = {1--26},
      }

Citing the Software
===================

The software can be cited in two ways, either the project as a whole, or 
specific versions. We use `Zenodo <https://doi.org/10.5281/zenodo.7126039>`_ to associate a DOI 
with each release of TUV-x.

Cite all versions
^^^^^^^^^^^^^^^^^
The DOI that represents all versions of TUV-x, and, when clicked, will resolve
to the most recent version of TUV-x, should be cited with this DOI.

.. code-block:: Tex
    :caption: This will always point to the most recent version.

        @software{ncar.tuvx,
          author       = {Stacy Walters and
                          Matt Dawson and
                          Kyle Shores},
          title        = {NCAR/tuv-x},
          month        = sep,
          year         = 2022,
          publisher    = {Zenodo},
          doi          = {10.5281/zenodo.7126039},
          url          = {https://doi.org/10.5281/zenodo.7126039}
        }

Version 0.2.0
^^^^^^^^^^^^^

.. code-block:: Tex
    :caption: This will always point to version 0.2.0.

        @software{ncar.tuvx,
          author       = {Stacy Walters and
                          Matt Dawson and
                          Kyle Shores},
          title        = {NCAR/tuv-x: Version 0.2.0},
          month        = sep,
          year         = 2022,
          publisher    = {Zenodo},
          version      = {v0.2.0},
          doi          = {10.5281/zenodo.7126040},
          url          = {https://doi.org/10.5281/zenodo.7126040}
        }