ETFL
====
|Build Status| |Codecov| |Codacy branch grade| |license|

ETFL: A formulation for flux balance models accounting for expression, thermodynamics, and resource allocation constraints

Preprint_ on BioRxiv

This code is an early release. You will need pyTFA_ to run it.
We recommend using commercial solvers such as CPLEX or Gurobi to run these problems.

Docker
------

We recommend the use of Docker to set up a container that will have the proper environment and package version to make ETFL work.

Right now, the ETFL Docker is built on top of the pyTFA Docker. 
If you want to use Docker-based install, you will need a working pytfa docker image, with either CPLEX or Gurobi on it. 
You can install them by following the instructions in pyTFA's Documentation_.

.. _Preprint: https://www.biorxiv.org/content/10.1101/590992v1
.. _pyTFA: https://github.com/EPFL-LCSB/pytfa
.. _Documentation: https://pytfa.readthedocs.io/en/latest/solver.html
.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/etfl/blob/master/LICENSE.txt
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/etfl.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/etfl
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/etfl.svg
   :target: https://codecov.io/gh/EPFL-LCSB/etfl
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/46bab484396946a8be07a82276f3e9dc/master.svg
   :target: https://www.codacy.com/app/realLCSB/etfl
