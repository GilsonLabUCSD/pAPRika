wget http://repo.continuum.io/miniconda/Miniconda3-4.3.11-Linux-x86_64.sh -O miniconda.sh

bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
conda update conda -y
conda install --yes conda-build pip numpy scipy pandas coverage
conda config --add channels omnia

conda create -y -n myenv python=3
source activate myenv

wget https://github.com/ParmEd/ParmEd/archive/2.7.3.tar.gz -O parmed.tar.gz
tar xvfz parmed.tar.gz
cd ParmEd-2.7.3
python setup.py install
cd ../../
