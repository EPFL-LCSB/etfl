ETFL
====

.. image:: https://api.codacy.com/project/badge/Grade/db2121b54ecc4284bf13c7eb430d5809
   :alt: Codacy Badge
   :target: https://app.codacy.com/app/EPFL-LCSB/etfl?utm_source=github.com&utm_medium=referral&utm_content=EPFL-LCSB/etfl&utm_campaign=Badge_Grade_Dashboard

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
