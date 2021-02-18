Release History
===============

v1.1.0
------
**(02/17/2021)**

This release provides a number of enhancements to the paprika code and changes to the API.

* Adds support for writing APR restraints to ``Plumed`` and ``Colvars`` format.
* Refactoring of the ``AMBER`` simulation wrapper code and expanded to include wrappers for ``GROMACS`` and ``NAMD``.
* Namespace refactoring for the ``align``, ``dummy``, and ``tleap`` modules (now under ``paprika.build``).
* Namespace refactoring for the OpenFF-Evaluator related modules (now under ``paprika.evaluator``).
* Enhancements to the ``tleap`` class, which includes a wrapper around ``InterMol`` for converting ``AMBER`` files to other MD formats and support for solvating with the ``Bind3P`` water model.
* Migrate continuous integration (CI) to Github Actions (GHA).
* Minor bug fixes.


v1.0.4
------
**(08/18/2020)**

* Refactor codes for ``setup.py`` and ``analysis.py`` modules to allow a smoother integration with OpenFF-Evaluator
* OpenMM restraints code is now decoupled from ``setup.py`` and placed in a separate file ``restraints/openmm.py``
* Updated install requirements to Ambertools v20
* Added hydrogen mass repartitioning option to ``TLeap`` API module
* Minor bug fixes


v1.0.3
------
**(12/04/2019)**

Fixed API compatibility with the new pymbar v3.0.5 and the I/O for ``ref_structure`` in the ``static_DAT_restraint`` function


v1.0.2
------
**(11/24/2019)**

Fix parsing of ``CUDA_VISIBLE_DEVICES`` class property. This bugfix was intended to go into the ``v1.0.1`` release but was accidentally excluded.


v1.0.1
------
**(11/24/2019)**

Provides a more robust method to supply a random number seed to the ``compute_free_energy`` method.


v1.0.0
------
**(11/04/2019)**

* The API for ``paprika.setup`` and ``paprika.analyze`` are now converged.
* OpenMM support is proved by PropertyEstimator.
* The documentation has (mostly) been updated and expanded.


v0.1.1
------
**(07/30/2019)**

Minor code cleanup and fix reference values for test cases that were failing in the v0.1.0 release.


v0.1.0
------
**(07/27/2019)**

This release enables ``paprika`` to work with ``propertyestimator`` of the Open Force Field
Toolkit, and ``taproom`` to setup host-guest calculations with either GAFF or SMIRNOFF-based force fields using AMBER or
OpenMM as simulation backends.


v0.0.4
------
**(02/19/2019)**

Updates after merging Cookiecutter for Computational Molecular Sciences (CMS) Python Packages and updating tests to pass.


v0.0.3
------
**(05/16/2018)**

Initial code and was used to run `SAMPLing <https://github.com/slochower/SAMPLing>`_.