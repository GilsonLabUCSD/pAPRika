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

*pAPRika is a python toolkit for setting up, running, and analyzing free energy molecular dynamics simulations.*


===================


Theory
******


Binding free energy
===================

In the attach-pull-release (APR) method, a guest molecule is physically pulled out of the host molecule along a defined
path by applying harmonic restraints. The binding free energy is then computed in terms of the sum of the work required
for (1) `attaching` restraints to the bound complex, (2) `pulling` the guest molecule out of the host molecule and (3)
`releasing` the restraints of the unbound complex:

.. math ::
   \Delta G^\circ_\text{bind} = -(W_\text{attach} + W_\text{pull} + W_\text{release-conf} + W_\text{release-std})
   :label: eq:bind


.. image :: _static/images/pmf.png
   :width: 600px
   :align: center

In the equation above, the work of releasing the restraints is split into two: :math:`W_\text{release-conf}` and
:math:`W_\text{release-std}`. These represents the work for releasing conformational restraints (applied to host or
guest molecules) and the guest's translational and rotational restraints to the standard concentration, respectively.
The translational and rotational degrees of freedom of the guest molecule is defined in spherical coordinates
:math:`(r,\theta,\phi)` and Euler angles :math:`(\alpha,\beta,\gamma)`, respectively. Unlike the work for releasing
conformational restraints, the work for releasing the guest's restraints to the standard concentration can be evaluated
analytically:

.. math ::
   \begin{eqnarray}
      W_\text{release-std} & = & RT\,\text{ln}\left( \frac{C^{\circ}}{8\pi^2} \right) \\
         & + & RT\,\text{ln}\left(\int^{\infty}_{0}\int^{\pi}_{0}\int^{2\pi}_{0} e^{-\beta U(r,\theta,\phi)} r^2 \text{sin}(\theta) \, d\phi \, d\theta \, dr \right) \\
         & + & RT\,\text{ln}\left(\int^{2\pi}_{0}\int^{\pi}_{0}\int^{2\pi}_{0} e^{-\beta U(\alpha,\beta,\gamma)} \text{sin}(\theta) \, d\alpha \, d\beta \, d\gamma \right)
   \end{eqnarray}
   :label: eq:release

For other terms in equation :eq:`eq:bind`, the calculation is broken up into a series of independent windows. During the
attach and release phase, the force constant of the restraints increases and decreases, respectively. In the pull phase,
the distance between the guest and host is discretized from point A (bound) to point B (unbound) by changing the
equilibrium position of the distance restraint. The work for each phase can be estimated using either `thermodynamic integration`
(TI) or the `multistate-Bennett-Acceptance-Ratio` (MBAR).

.. math ::
   \begin{eqnarray}
      W_\text{attach,release} & = & \int_{0}^{1} \left<F\right>_{\lambda} \, d\lambda \\
      W_\text{pull} & = & \int_{A}^{B} \left<F\right>_{r} \, dr
   \end{eqnarray}
   :label: eq:TI

`pAPRika` estimates the standard error of the mean (SEM) with either `block data` analysis or `autocorrelation`.

Binding enthalpy
================

The binding enthalpy is the difference in the partial molar enthalpies of the bound complex and the separated molecules.
Setting up binding enthalpy calculations is more straightforward than the binding free energy if we use the `direct`
method. In the `direct` method, the binding enthalpy is obtained from the difference in the mean potential energy of the
bound and unbound complex.

.. math ::
   \Delta H_\text{bind} = \left<U_\text{bound} \right> - \left<U_\text{unbound} \right>
   :label: eq:enthalpy

In the context of APR simulations, the mean potential energy is estimated from the first and last window. These windows,
however, will need to be run for much longer until the fluctuation of the mean potential energy decreases below a desired
value. Another way to estimate the binding enthalpy is with the `van't Hoff` method:

.. math ::
   \text{ln}\, K_\text{eq} = -\frac{\Delta H}{RT} + \frac{\Delta S}{R}
   :label: eq:vant

Here, the binding enthalpy is obtained by a linear regression of :math:`\text{ln}\,K_\text{eq}` vs :math:`1/T`. The
equilibrium constant will need to be evaluated at different temperatures in order the estimate the binding enthalpy
with the `van't Hoff` method.

===================




Supported Molecular Dynamics Engine
***********************************
Currently, `pAPRika` can be used to setup and analyze APR simulations with `AMBER <https://ambermd.org/>`_,
`OpenMM <http://openmm.org/>`_ and `GROMACS <http://www.gromacs.org/>`_. These MD engines provides their own interface
(e.g., AMBER uses NMR-style restraints). `pAPRika` also provides modules to generate `Plumed <https://www.plumed.org/>`_-based
restraints, which is supported in a number of MD engines.

Future releases of `pAPRika` will include support for other MD engines like ,
`NAMD <https://www.ks.uiuc.edu/Research/namd/>`_ and `LAMMPS <https://lammps.sandia.gov/>`_. Plumed is supported by all
of these but for `NAMD` only NPT simulations can be performed. Hence, a module converting `pAPRika` restraints to a
`Colvars <https://github.com/Colvars/colvars>`_ module, which is supported in `NAMD` and `LAMMPS`, is also in the works.

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




License
*******
`pAPRika` is licensed under the BSD 3-Clause license.


===================


References
**********
1. Velez-Vega, C. & Gilson, M. K. Overcoming Dissipation in the Calculation of Standard Binding Free Energies by Ligand Extraction. J. Comput. Chem. 34, 2360–2371 (2013).
2. Fenley, A. T., Henriksen, N. M., Muddana, H. S. & Gilson, M. K. Bridging calorimetry and simulation through precise calculations of cucurbituril-guest binding enthalpies. J. Chem. Theory Comput. 10, 4069–4078 (2014).
3. Henriksen, N. M., Fenley, A. T. & Gilson, M. K. Computational Calorimetry: High-Precision Calculation of Host−Guest Binding Thermodynamics. J. Chem. Theory Comput. 11, 4377–4394 (2015).



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   Overview <self>
   install
   workflow


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

   build_docs
   api
   releasehistory
