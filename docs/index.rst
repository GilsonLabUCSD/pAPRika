.. paprika documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================
pAPRika
==================

.. role:: green-font
.. role:: red-font
.. role:: ignore-width

.. |tick|    replace:: :green-font:`✓`
.. |cross|   replace:: :red-font:`✕`
.. |delta|   unicode:: U+0394
.. |ast|     replace:: :ignore-width:`*`


*pAPRika is a python toolkit for setting up, running, and analyzing free-energy molecular dynamics simulations.*

==================

Theory
------

==================

Supported Molecular Dynamics Engine
-----------------------------------
Currently, pAPRika can be used to setup and analyze APR simulations with `Amber <https://ambermd.org/>`_ and
`OpenMM <http://openmm.org/>`_. Both of these MD engines provides their own interface for restraints namely NMR and XML
for Amber and OpenMM, respectively. pAPRika also provides modules to generate `Plumed <https://www.plumed.org/>`_-based
restraints restraints, which is supported with both engines.

Future releases of pAPRika will include support for other MD engines like `Gromacs <http://www.gromacs.org/>`_,
`NAMD <https://www.ks.uiuc.edu/Research/namd/>`_ and `LAMMPS <https://lammps.sandia.gov/>`_. Plumed is supported by all
of these but for NAMD only NPT simulations can be performed. Hence, a module converting pAPRika restraints to a
`Colvars <https://github.com/Colvars/colvars>`_ module, which is supported in NAMD and LAMMPS, is also in the works.

.. table:: Molecular Dynamics engine that is supported in pAPRika (planned) and the respective restraint modules it supports.
   :widths: auto
   :align: center
   :class: clean-table property-table

   +-----------------+----------------+----------------+----------------+----------------+
   || MD Engine      || NMR           || XML           || Plumed        || Colvars       |
   +=================+================+================+================+================+
   || Amber          || |tick|        || |cross|       || |tick|        || |cross|       |
   +-----------------+----------------+----------------+----------------+----------------+
   || OpenMM         || |cross|       || |tick|        || |tick|        || |cross|       |
   +-----------------+----------------+----------------+----------------+----------------+
   || Gromacs\ |ast| || |cross|       || |cross|       || |tick|        || |cross|       |
   +-----------------+----------------+----------------+----------------+----------------+
   || NAMD\ |ast|    || |cross|       || |cross|       || |tick|        || |tick|        |
   +-----------------+----------------+----------------+----------------+----------------+
   || LAMMPS\ |ast|  || |cross|       || |cross|       || |tick|        || |tick|        |
   +-----------------+----------------+----------------+----------------+----------------+

*\* Entries marked with an asterisk are currently not supported but will be available in future releases.*

Tutorials
---------

.. Setup the side-pane table of contents.

========


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
   05 - Amber with Plumed <tutorials/05-tutorial-cb6-but-plumed.ipynb>
..   06 - Gromacs Simulation <tutorials/06-tutorial-cb6-but-gromacs.ipynb>
..   07 - NAMD Simulation <tutorials/07-tutorial-cb6-but-namd.ipynb>
..   08 - LAMMPS Simulation <tutorials/07-tutorial-cb6-but-lammps.ipynb>

.. toctree::
  :maxdepth: 10
  :hidden:
  :caption: Developer Documentation

  align
  analysis
  dummy
  evaluator
  restraints
  simulation
  tleap
  utils
  releasehistory
