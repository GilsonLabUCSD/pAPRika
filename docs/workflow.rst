*****
Usage
*****

Workflow
--------
todo: workflow-diagram (explain the caveats when using different MD engines).

Structure Preparation
---------------------


Defining Restraints
-------------------

In APR calculations we apply restraints on the host (or protein) and the guest molecules. The restraints can be grouped
into three different categories: (1) *static restraints*, (2) *changing restraints* and (3) *wall restraints*.

**(1) Static Restraints**

Static restraints do not change during the whole APR calculations and do not affect the free-energy. We apply static
restraints on the host (or protein) molecule to define a path for the guest molecule. The static restraints are a made up
of a combination of distance, angle and torsional (DAT) restraints based on the choice of anchor atoms. For host-guest
systems we need to define three anchor atoms ``[H1, H2, H3]`` and combined with three dummy atoms ``[D1, D2, D3]`` we
apply a total of six static restraints on the host molecule (three for the translation and three for orientation).

To generate a static restraint we use the function `static_DAT_restraints`. As an example, to apply a distance restraints
on ``D1`` and ``H1`` with a force constant of 5 kcal/mol/:math:`Å^2` we call

.. code :: python

    dist_static = static_DAT_restraint(
        restraint_mask_list = [D1, H1],
        num_window_list = windows,
        ref_structure = structure,
        force_constant = 5.0,
    )

where ``windows`` is a list of the length of each phase and ``structure`` is the structure file (usually a PDB).

**(2) Changing restraints**


**(3) Wall restraints**

Running a Simulation
--------------------


Analysis
--------
