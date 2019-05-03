Install pAPRika
==========================

We recommend installing pAPRika in a fresh `conda` environment if possible. There are three ways to install this package:

1. The latest release on `conda-forge`:
    1. `conda install -c conda-forge paprika`
    2. To use all features of pAPRika, you must either have [AmberTools](http://ambermd.org/AmberTools.php) in your `$PATH` or separately install AmberTools with `conda install -c http://ambermd.org/downloads/ambertools/conda/ ambertools=18`.
    3. To use OpenMM features: `conda install -c omnia openmm`.

2. The master branch on GitHub:
    1. Clone this `git` repository, then inside the `paprika` directory:
    2. Change the `name` field in `devtools/conda-envs/test_env.yaml` to be `paprika`.
    3. Create the environment: `conda env create -f devtools/conda-envs/test_env.yaml`.
    4. Activate the environment: `conda activate paprika`
    5. Install `paprika` in the environment: `pip install .`
    
3. The latest release on GitHub:
    1. Download [the latest release](https://github.com/slochower/pAPRika/releases), extract it, and change to the `paprika` directory:
    2. Change the `name` field in `devtools/conda-envs/test_env.yaml` to be `paprika`.
    3. Create the environment: `conda env create -f devtools/conda-envs/test_env.yaml`.
    4. Activate the environment: `conda activate paprika`
    5. Install `paprika` in the environment: `pip install .`