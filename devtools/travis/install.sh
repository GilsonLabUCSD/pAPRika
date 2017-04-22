wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh

bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
conda update conda -y
conda install --yes conda-build pip
conda config --add channels omnia

conda create -y -n myenv python=$PYTHON_VERSION \
       numpy scipy pandas nose openmm coverage nose-timer \
       python-coveralls netCDF4

source activate myenv

wget https://github.com/ParmEd/ParmEd/archive/2.7.3.tar.gz -O parmed.tar.gz
tar xvfz parmed.tar.gz
cd ParmEd-2.7.3
python setup.py install
cd ../../
