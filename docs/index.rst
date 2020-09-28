.. paprika documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pAPRika
==================

.. role:: green-font
.. role:: red-font
.. role:: ignore-width

.. |tick|    replace:: :green-font:`✓`
.. |cross|   replace:: :red-font:`✕`
.. |delta|   unicode:: U+0394
.. |ast|     replace:: :ignore-width:`*`
.. |dast|    replace:: :ignore-width:`**`
.. |br| raw:: html

   <br />

*pAPRika is a python toolkit for setting up, running, and analyzing free-energy molecular dynamics simulations.*

===================

Theory
------
todo: explain APR theory (add images/diagram and equations)



===================

Supported Molecular Dynamics Engine
-----------------------------------
Currently, `pAPRika` can be used to setup and analyze APR simulations with `AMBER <https://ambermd.org/>`_,
`OpenMM <http://openmm.org/>`_ and `GROMACS <http://www.gromacs.org/>`_. These MD engines provides their own interface
(e.g., AMBER uses NMR-style restraints). `pAPRika` also provides modules to generate `Plumed <https://www.plumed.org/>`_-based
restraints restraints, which is supported in a number of MD engines.

Future releases of `pAPRika` will include support for other MD engines like ,
`NAMD <https://www.ks.uiuc.edu/Research/namd/>`_ and `LAMMPS <https://lammps.sandia.gov/>`_. Plumed is supported by all
of these but for NAMD only NPT simulations can be performed. Hence, a module converting `pAPRika` restraints to a
`Colvars <https://github.com/Colvars/colvars>`_ module, which is supported in NAMD and LAMMPS, is also in the works.

.. table:: Molecular Dynamics engine that is supported in *pAPRika* (or planned) and the respective restraint modules it supports.
   :widths: auto
   :align: center
   :class: clean-table property-table

   +-----------------+----------------+----------------+-----------------+-----------------+
   || MD Engine      || NMR           || XML           || Plumed         || Colvars\ |ast| |
   +=================+================+================+=================+=================+
   || AMBER          || |tick|        || |cross|       || |tick|         || |cross|        |
   +-----------------+----------------+----------------+-----------------+-----------------+
   || OpenMM         || |cross|       || |tick|        || |tick|\ |dast| || |cross|        |
   +-----------------+----------------+----------------+-----------------+-----------------+
   || GROMACS        || |cross|       || |cross|       || |tick|         || |cross|        |
   +-----------------+----------------+----------------+-----------------+-----------------+
   || NAMD\ |ast|    || |cross|       || |cross|       || |tick|         || |tick|         |
   +-----------------+----------------+----------------+-----------------+-----------------+
   || LAMMPS\ |ast|  || |cross|       || |cross|       || |tick|         || |tick|         |
   +-----------------+----------------+----------------+-----------------+-----------------+

*\* Currently not supported but will be available in future releases.* |br|
*\*\* Supported but have not yet been tested.*

===================

pAPRika Workflow
----------------
todo: workflow-diagram (explain the caveats when using different MD engines).

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   Overview <self>
   install


.. toctree::
   :maxdepth: 10
   :hidden:
   :caption: Tutorials

   01 - Implicit Solvent <tutorials/01-tutorial-cb6-but.ipynb>
   02 - Explicit Solvent <tutorials/02-tutorial-cb6-but-pbc.ipynb>
   03 - K-Cl dissociation <tutorials/03-tutorial-k-cl.ipynb>
   04 - OpenMM with XML <tutorials/04-tutorial-cb6-but-openmm.ipynb>
   05 - AMBER with Plumed <tutorials/05-tutorial-cb6-but-plumed.ipynb>
   06 - GROMACS Simulation <tutorials/06-tutorial-cb6-but-gromacs.ipynb>
..   07 - NAMD Simulation <tutorials/07-tutorial-cb6-but-namd.ipynb>
..   08 - LAMMPS Simulation <tutorials/07-tutorial-cb6-but-lammps.ipynb>

.. toctree::
  :maxdepth: 10
  :hidden:
  :caption: Developer Documentation

  modules/align
  modules/analysis
  modules/dummy
  modules/evaluator
  modules/restraints
  modules/simulation
  modules/tleap
  modules/utils
  releasehistory
