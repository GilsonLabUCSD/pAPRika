Install pAPRika
==========================

We recommend installing pAPRika in a fresh ``conda`` environment if possible. There are three ways to install this package:

1. The latest release on ``conda-forge``:
    a. ``conda install -c conda-forge paprika``
    b. To use all features of pAPRika, you must either have _AmberTools: http://ambermd.org/AmberTools.php in your `$PATH` or separately install AmberTools with ``conda install -c http://ambermd.org/downloads/ambertools/conda/ ambertools=19``.
    c. To use OpenMM features: ``conda install -c omnia openmm``.

2. The master branch on GitHub:
    a. Clone this ``git`` repository, then inside the ``paprika`` directory:
    b. Change the ``name`` field in ``devtools/conda-envs/test_env.yaml`` to be ``paprika``.
    c. Create the environment: ``conda env create -f devtools/conda-envs/test_env.yaml``.
    d. Activate the environment: ``conda activate paprika``
    e. Install ``paprika`` in the environment: ``pip install .``
    
3. The latest release on GitHub:
    a. Download _the latest release: https://github.com/slochower/pAPRika/releases, extract it, and change to the ``paprika`` directory:
    b. Change the ``name`` field in ``devtools/conda-envs/test_env.yaml`` to be ``paprika``.
    c. Create the environment: ``conda env create -f devtools/conda-envs/test_env.yaml``.
    d. Activate the environment: ``conda activate paprika``
    e. Install ``paprika`` in the environment: ``pip install .``