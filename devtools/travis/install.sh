rm -rf "$HOME/miniconda"
if [ -d "$HOME/miniconda/bin/python" ]; then
      export PATH="$HOME/miniconda/bin:$PATH"
      conda list
      pip list
else
      if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
            wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
      else
            wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
      fi    
      chmod +x ./miniconda.sh
      ./miniconda.sh -b -u -p $HOME/miniconda
      export PATH="$HOME/miniconda/bin:$PATH"
      
      conda config --add channels omnia --add channels conda-forge
      conda create -y -n myenv python=$PYTHON_VERSION
      conda install -y -n myenv \
            openmm numpy scipy pandas pytest pytest-cov codecov mdtraj pymbar
      conda install -y -n myenv \
            ambertools=17.0 -c http://ambermd.org/downloads/ambertools/conda/
      pip install codecov
fi

source activate myenv
python --version
