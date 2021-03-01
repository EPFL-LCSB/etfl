yETFL------------------------------------------------To construct ETFL model for the yeast, Saccharomyces cerevisiae. Paper: Omid Oftadeh, Pierre Salvy, Maria Masid, Maxime Curvat, Ljubisa Miskovic and Vassily Hatzimanikatis. A genome-scale metabolic model of Saccharomyces cerevisiae that integrates expression constraints and reaction thermodynamicsRequirements------------------------------------------------
For this package the following packages are required:
- cobra (https://opencobra.github.io/cobrapy/)
- pytfa (https://github.com/EPFL-LCSB/pytfa)

The code can be run on different types of OS including UNIX and Windows.

You will need Git LFS to properly download some binary files:

git lfs install
git lfs pull

We recommend the use of Docker to set up a container instead of running the codes directly. To solve the genome-scale models, a licensed commercial solver, either Gurobi or CPLEX, should be installed. The free academic versions can be obtained through Gurobi for Acedemics and Researchers (https://www.gurobi.com/academia/academic-program-and-licenses/) or IBM Academic initiative (https://developer.ibm.com/academic/). Alternatively, to test the general performance of the code, the scripts in tests/ subdirectory can be run with the free solver GLPK. For further details, see README.rst.Installing------------------------------------------------Run the following command in a bash terminal:git clone -b dev_yetfl https://github.com/EPFL-LCSB/etfl.git

To install the package the setup.py script should be run. The installation process will be fast (It takes <10 seconds).Generating yETFL models------------------------------------------------
To create the yETFL models, run the script tutorials/helper_gen_models_yeast.py. The output will be saved automatically in json format.Analysing the models------------------------------------------------To benchmark the growth rate at different uptake rates, run the script tutorials/benchmark_yeast_models.py. The output will be automatically saved for different uptake rates in a csv file.To simulate the Crabtree effect, run the script tutorials/crabtree_cheby_parsim.py. The output will be automatically saved for different uptake rates in a csv file.To do gene essentiality, run the script tutorials/timing_gene_essentiality.py. The output will be automatically saved for different genes in a csv file. The run time for the simulations might be different based on the system performance. However, the run time is typically in the order of ~hours. 