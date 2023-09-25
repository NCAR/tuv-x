

Build the docs
==============

Our documentation is written in `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_.
Read their `primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ 
to get up to speed on writing in the 
`reStructuredText format <https://docutils.sourceforge.io/rst.html>`_.

Install Dependencies
--------------------

All of the dependencies should be installed with `pip` or `conda`, except
for the package which documents TUV-x's source code,  
`sphinx-fortran <https://sphinx-fortran.readthedocs.io/en/latest/>`_. This is
because published `sphinx-fortran` do not work with python3. 

Install each of the requirements in 
`requirments.txt <https://github.com/NCAR/tuv-x/blob/main/docs/requirements.txt>`_.

Build and view the output
-------------------------
From the root directory of tuv-x, run these commands

.. code-block:: console
  
  cd docs
  make html

  #macOS only
  open _build/html/index.html


Navigate to `tuv-x/docs/_build/html/index.html` in your browser and you should
see the docs.

Preview branch changes
^^^^^^^^^^^^^^^^^^^^^^
If you are added as a contributor to the github project, anytime you push a commit
on your branch, a copy of your documentation is hosted on github pages at
``https://ncar.github.io/tuv-x/branch/<your-branch-name>``. You can preview
any changes you make there before submitting a PR to ensure that your changes
work on our github pages.

Style guide
-----------

.. stolen from https://doc.akka.io/docs/akka/2.3/dev/documentation.html#Sections

We use the following identifiers to add logical structure in the documentation.

- ``#`` (over and under) for module headings
- ``=`` for sections
- ``-`` for subsections
- ``^`` for subsubsections
- ``~`` for subsubsubsections
- ``"`` for paragraphs