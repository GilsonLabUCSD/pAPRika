.. paprika documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



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

.. image :: _static/images/paprika.png
   :width: 350px
   :align: left

|
|
|
|
|
|

*pAPRika is a python toolkit for setting up, running, and analyzing free-energy molecular dynamics simulations based on the attach-pull-release (APR) method.*

===================

******
Theory
******


Binding free-energy
===================

todo: explain APR theory (add images/diagram and equations)

.. math ::

   \begin{eqnarray}
      \Delta G^{\circ}_{b} = -(W_{attach} + W_{pull} + W_{release-conf} + W_{release-std})
   \end{eqnarray}

.. image :: _static/images/pmf.png
   :width: 600px
   :align: center


.. math ::

   \begin{eqnarray}
      W_{release-std} & = & RT\,\text{ln}\left( \frac{C^{\circ}}{8\pi^2} \right) \\
                      & + & RT\,\text{ln}\left(\int^{\infty}_{0}\int^{\pi}_{0}\int^{2\pi}_{0} e^{-\beta U(r,\theta,\phi)} r^2 \text{sin}(\theta) \, d\phi d\theta dr \right) \\
                      & + & RT\,\text{ln}\left(\int^{2\pi}_{0}\int^{\pi}_{0}\int^{2\pi}_{0} e^{-\beta U(\alpha,\beta,\gamma)} \text{sin}(\theta) \, d\alpha d\beta d\gamma \right)
   \end{eqnarray}


Binding enthalpy
================

todo: Write about how to calculate binding enthalpy

.. math ::

   \begin{eqnarray}
   \Delta H_{bind} = \left<U_{bound} \right> - \left<U_{unbound} \right>
   \end{eqnarray}



===================



***********************************
Supported Molecular Dynamics Engine
***********************************
Currently, `pAPRika` can be used to setup and analyze APR simulations with `AMBER <https://ambermd.org/>`_,
`OpenMM <http://openmm.org/>`_ and `GROMACS <http://www.gromacs.org/>`_. These MD engines provides their own interface
(e.g., AMBER uses NMR-style restraints). `pAPRika` also provides modules to generate `Plumed <https://www.plumed.org/>`_-based
restraints restraints, which is supported in a number of MD engines.

Future releases of `pAPRika` will include support for other MD engines like ,
`NAMD <https://www.ks.uiuc.edu/Research/namd/>`_ and `LAMMPS <https://lammps.sandia.gov/>`_. Plumed is supported by all
of these but for NAMD only NPT simulations can be performed. Hence, a module converting `pAPRika` restraints to a
`Colvars <https://github.com/Colvars/colvars>`_ module, which is supported in NAMD and LAMMPS, is also in the works.

.. table:: MD engines that are supported in `pAPRika` (or planned) and the respective restraint modules it supports.
   :widths: auto
   :align: center
   :class: clean-table property-table

   +-----------------+---------------------------------------------------------------------+
   ||                || Restraints Module                                                  |
   || MD Engine      +----------------+----------------+-----------------+-----------------+
   ||                || NMR           || XML           || Plumed         || Colvars\ |ast| |
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

****************
pAPRika Workflow
****************
todo: workflow-diagram (explain the caveats when using different MD engines).


===================


**********
References
**********
1. Velez-Vega, C. & Gilson, M. K. Overcoming Dissipation in the Calculation of Standard Binding Free Energies by Ligand Extraction. J. Comput. Chem. 34, 2360–2371 (2013).
2. Fenley, A. T., Henriksen, N. M., Muddana, H. S. & Gilson, M. K. Bridging calorimetry and simulation through precise calculations of cucurbituril-guest binding enthalpies. J. Chem. Theory Comput. 10, 4069–4078 (2014).
3. Henriksen, N. M., Fenley, A. T. & Gilson, M. K. Computational Calorimetry: High-Precision Calculation of Host−Guest Binding Thermodynamics. J. Comput. Theory Nanosci. 11, 4377–4394 (2015).




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

  source/align
  source/analysis
  source/dummy
  source/evaluator
  source/restraints
  source/simulation
  source/tleap
  source/utils
  releasehistory
