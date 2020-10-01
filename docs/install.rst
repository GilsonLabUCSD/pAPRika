************
Installation
************

We recommend installing *pAPRika* in a fresh ``conda`` environment if possible. There are three ways to install this package:

Installing from Conda
---------------------

To install the latest release of *pAPRika* from conda, run::

    conda install -c conda-forge paprika

In order to use all features of *pAPRika*, you must either have ``AmberTools`` (http://ambermd.org/AmberTools.php) in your `$PATH` or separately install ``AmberTools`` with::

    conda install -c conda-forge ambertools=20

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^
You can install ``OpenMM`` from the ``omnia`` channel if you want to run APR simulations with ``OpenMM``::

    conda install -c omnia openmm

If you want to run simulations with `Plumed <https://www.plumed.org/>`_-based restraints (needed for running APR in ``GROMACS``) you can compile Plumed from source or install through conda::

    conda install -c conda-forge plumed

Although ``GROMACS`` is available in conda, a version that is patched with Plumed is currently not available. Therefore, if you want to run ``GROMACS`` simulations in *pAPRika* you will need to compile from source manually (patched with Plumed).

Installing a stable version from source
---------------------------------------

To install the a stable version of *pAPRika* download the source code from `Github <https://github.com/slochower/pAPRika/releases>`_.
Then, unzip the files and change to the ``paprika`` directory::

    tar -xvzf paprika-version.tgz
    cd pAPRika-version/

Change the ``name`` field in ``devtools/conda-envs/test_env.yaml`` to be ``paprika`` and create the environment::

    conda env create -f devtools/conda-envs/test_env.yaml

Activate the environment::

    conda activate paprika

and install *pAPRika* in the environment::

    pip install .



Installing latest from source
-----------------------------

To install *pAPRika* with the latest features, clone the repository from the ``master`` branch on `Github <https://github.com/slochower/pAPRika>`_::

    git clone https://github.com/slochower/pAPRika.git

Change directory to the ``paprika`` folder, change the ``name`` field in ``devtools/conda-envs/test_env.yaml`` to ``paprika`` and create the conda environment::

    conda env create -f devtools/conda-envs/test_env.yaml

Activate the environment::

    conda activate paprika

and install *pAPRika* in the environment::

    pip install .

