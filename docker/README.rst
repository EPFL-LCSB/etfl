ETFL Docker
===========

This Docker offers a suitable environment to run ETFL.

Requirements
------------

Make sure `docker`_ is installed. 
You will need the pyTFA container already built, as it is where the solver hooks are set-up.
Follow the instructions at `pytfa/docker`_

Running the Docker
------------------
Since you will need to build the pytfa container, we recommend the following architecture:

.. code::

    └───src_folder
        └───pytfa
        └───etfl

The ``run`` script will build pyTFA with the local installation at ``../pytfa``.
if pyTFA is located somewhere else, you can edit the arguments of the run script. 
If you want to download pyTFA from PIP, simply change the last line in the Dockerfile to remove the local installation.

First, build the container with ``build.bat`` or ``. build``. Then start
the container with ``run.bat`` or ``. run``.

.. code:: bash

   . build
   . run

You can run the examples in /etfl/tutorials:

.. code:: bash

   cd /etfl/tutorials
   python test_small.py

You can also run them inside IPython to experiment and play with the
objects:

.. code:: bash

   ipython
   run test_small.py
   m.print_info()

Additional information
----------------------

If you are running your Docker container in a Unix-based environment,
you might get permission errors on the ``.sh`` scripts. This is because
permissions are inherited from the host environment. You can fix this by
running in the ``docker`` folder:

.. code:: bash

   chmod +x utils/*.sh

.. _docker: https://www.docker.com/
.. _solver/instructions.txt: https://github.com/EPFL-LCSB/pytfa/blob/master/docker/solvers/instructions.txt
.. _pytfa/docker: https://github.com/EPFL-LCSB/pytfa/blob/master/docker/README.md

