ETFL
====
|Documentation Status| |Build Status| |Codecov| |Codacy branch grade| |license| |doi|

ETFL: A formulation for flux balance models accounting for expression, thermodynamics, and resource allocation constraints

Paper_: Salvy, P., Hatzimanikatis, V. The ETFL formulation allows multi-omics integration in thermodynamics-compliant metabolism and expression models. Nat Commun 11, 30 (2020) doi:10.1038/s41467-019-13818-7


See `ecETFL <https://github.com/EPFL-LCSB/ecetfl/>`_ for an E. coli model.

See also yETFL_ for a yeast model!

This code is an early release. You will need pyTFA_ to run it.
We recommend using commercial solvers such as CPLEX or Gurobi to run these problems.

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/etfl.git /path/to/etfl
    cd /path/to/etfl
    git lfs install
    git lfs pull

**This module was developed in Python 3.5, and it is recommended to run Python 3.5 
to run commercial solvers such as Gurobi and CPLEX.**
Other Python versions (2.7, 3.4) might also work but are not officially supported (see the `CI builds <https://travis-ci.org/EPFL-LCSB/etfl>`_)


This module requires
`pyTFA <https://github.com/EPFL-LCSB/pytfa/>`_, as well as
`COBRApy <https://github.com/opencobra/cobrapy/>`_, and
`optlang <https://github.com/biosustain/optlang>`_ to work
properly. The installer should take care of that for you. You might also
want to install a dedicated solver. GLPK, CPLEX and Gurobi are
supported.

Installation
------------

The module can be installed like any Python package:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/etfl.git /path/to/etfl
    cd /path/to/etfl
    python3 setup.py install
    
The installation process should not exceed a minute if the requirements are installed. If they are not, it might take longer as the installer installs them first.

Running the code
----------------

You can run the examples in etfl/tutorials:

.. code:: bash

   cd etfl/tutorials
   python test_small.py

You can also run them inside IPython to experiment and play with the
objects:

.. code:: bash

   ipython
   run test_small.py
   m.print_info()

Docker
------

We recommend the use of Docker to set up a container that will have the proper environment and package version to make ETFL work.

Right now, the ETFL Docker is built on top of the pyTFA Docker. 
If you want to use Docker-based install, you will need a working pytfa docker image, with either CPLEX or Gurobi on it. 
You can install them by following the instructions in pyTFA's Documentation_.

More details are available in the `Docker folder <https://github.com/EPFL-LCSB/etfl/tree/master/docker>`_

License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/etfl/blob/master/LICENSE.txt>`_ file for more details.

.. _Paper: https://www.nature.com/articles/s41467-019-13818-7

.. _Preprint: https://www.biorxiv.org/content/10.1101/590992v1

.. _yETFL: https://github.com/EPFL-LCSB/yetfl

.. _Documentation: https://etfl.readthedocs.io/en/latest/solver.html

.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/etfl/blob/master/LICENSE.txt
.. |Documentation Status| image:: https://readthedocs.org/projects/etfl/badge/?version=latest
   :target: http://etfl.readthedocs.io/en/latest/?badge=latest
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/etfl.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/etfl
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/etfl.svg
   :target: https://codecov.io/gh/EPFL-LCSB/etfl
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/57efd28bef86473a8075fde96e132c28
   :target: https://www.codacy.com/app/realLCSB/etfl
.. |doi| image:: https://zenodo.org/badge/DOI/10.1038/s41467-019-13818-7.svg
    :target: https://doi.org/10.1038/s41467-019-13818-7
