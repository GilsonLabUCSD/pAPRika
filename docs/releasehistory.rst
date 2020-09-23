Release History
===============



v1.0.4
------
* Refactor codes for ``setup.py`` and ``analysis.py`` modules to allow a smoother integration with OpenFF-Evaluator
* OpenMM restraints code is now decoupled from ``setup.py`` and placed in a separate file ``restraints/openmm.py``
* Updated install requirements to Ambertools v20
* Added hydrogen mass repartitioning option to ``TLeap`` API module
* Minor bug fixes

v1.0.3
------
Fixed API compatibility with the new pymbar v3.0.5 and the I/O for ``ref_structure`` in the ``static_DAT_restraint`` function

v1.0.2
------
Fix parsing of ``CUDA_VISIBLE_DEVICES`` class property. This bugfix was intended to go into the ``v1.0.1`` release but was accidentally excluded.

v1.0.1
------
Provides a more robust method to supply a random number seed to the ``compute_free_energy`` method.

v1.0.0
------
* The API for ``paprika.setup`` and ``paprika.analyze`` are now converged.
* OpenMM support is proved by PropertyEstimator.
* The documentation has (mostly) been updated and expanded.


v0.1.1
------
Minor code cleanup and fix reference values for test cases that were failing in the v0.1.0 release.


v0.1.0
------
This release enables ``paprika`` to work with ``propertyestimator`` of the Open Force Field
Toolkit, and ``taproom`` to setup host-guest calculations with either GAFF or SMIRNOFF-based force fields using AMBER or
OpenMM as simulation backends.


v0.0.4
------
Updates after merging Cookiecutter for Computational Molecular Sciences (CMS) Python Packages and updating tests to pass.


v0.0.3
------
Initial code and was used to run `SAMPLing <https://github.com/slochower/SAMPLing>`_.