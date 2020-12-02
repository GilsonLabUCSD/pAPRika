********
Workflow
********

This page provides a brief explanation of the workflow to perform APR calculations with `pAPRika`. For more detail users
are recommended to go through the tutorials, which details further on how to setup and run APR simulations from start to
finish.

.. figure :: _static/images/flowchart.png
   :figwidth: 500px
   :width: 170px
   :align: center

   Flowchart of the `paprika` workflow for a typical APR simulation.

Structure Preparation
---------------------
The starting structure for the APR simulation can be configured with `paprika`. The APR calculation is most efficient in
a rectangular box with the long axis parallel to the pulling axis (reduces the number of water molecules in the system).
To make this easy to set up, the ``Align`` module provides functions to shift and orient a structure. For example, we can
translate a structure to the origin and then orient the system to the :math:`z`-axis by running

.. code ::

    from paprika.build.align import *

    translated_structure = translate_to_origin(structure)
    aligned_structure = zalign(translated_structure, ":GST@C1", ":GST@C2")

To provide something to pull "against," we add dummy atoms that are fixed in place with strong positional restraints.
These dummy atoms can be added to a the host-guest structure using the ``Dummy Atoms`` module in `paprika`.

.. code ::

    from paprika.build.dummy import add_dummy

    structure = dummy.add_dummy(structure, residue_name="DM1", z=-6.0)
    structure = dummy.add_dummy(structure, residue_name="DM2", z=-9.0)
    structure = dummy.add_dummy(structure, residue_name="DM3", z=-11.2, y=2.2)
    structure.save("aligned_with_dummy.pdb", overwrite=True)

We will need the ``mol2`` and ``frcmod`` files for the dummy atoms, which we will need to generate the `AMBER` topology

.. code ::

   dummy.write_dummy_frcmod(filepath="complex/dummy.frcmod")
   dummy.write_dummy_mol2(residue_name="DM1", filepath="complex/dm1.mol2")
   dummy.write_dummy_mol2(residue_name="DM2", filepath="complex/dm2.mol2")
   dummy.write_dummy_mol2(residue_name="DM3", filepath="complex/dm3.mol2")

Finally, we can use the ``tleap`` wrapper to combine all of these components to generate the topology and coordinate files.

.. code ::

    from paprika.build.system import TLeap

    system = TLeap()
    system.output_prefix = "host-guest-dum"
    system.pbc_type = None
    system.neutralize = False

    system.template_lines = [
        "source leaprc.gaff",
        "HST = loadmol2 host.mol2",
        "GST = loadmol2 guest.mol2",
        "DM1 = loadmol2 dm1.mol2",
        "DM2 = loadmol2 dm2.mol2",
        "DM3 = loadmol2 dm3.mol2",
        "model = loadpdb aligned_with_dummy.pdb",
    ]
    system.build()


Defining Restraints
-------------------

.. figure :: _static/images/restraints.png
   :figwidth: 550px
   :align: center

In APR calculations we apply restraints on the host (or protein) and the guest molecules. The restraints can be grouped
into four categories: (1) *static restraints*, (2) *varying restraints*, (3) *wall restraints* and (4) *positional
restraints*.

**(1) Static Restraints**

Static restraints do not change during the whole APR process and do not affect the free energy. We apply static restraints
on the host (or protein) molecule to orient the host/protein degrees of freedom. The static restraints are composed of
distance, angle, and torsional (DAT) restraints based on the choice of anchor atoms. For host-guest systems, we need to
define three anchor atoms ``[H1,H2,H3]`` and combined with three dummy atoms ``[D1,D2,D3]``, we apply a total of six
static restraints on the host molecule (three for the translation and three for orientation).

To generate static restraints we use the function ``static_DAT_restraints``. As an example, to apply a distance restraint
on ``D1`` and ``H1`` with a force constant of 5 kcal/mol/:math:`Å^2` we call

.. code :: python

    from paprika.restraints import static_DAT_restraint

    dist_static = static_DAT_restraint(
        restraint_mask_list = [D1, H1],
        num_window_list = windows,  # list: [len(attach_lambda), len(pull_windows), len(release_lambda)]
        ref_structure = structure,  # Structure file (PDB)
        force_constant = 5.0,
    )

**(2) Varying Restraints**

As the name suggests, these restraints change during the APR process. During the `attach` and `release` phases, the force
constants of these restraints changes. In the `pull` phase, `varying restraints` can have their equilibrium position
change, and this can be used as the restraint to pull the guest molecule out of the host molecule.

To generate `varying restraints`, we use the ``DAT_restraint`` class. The code below shows a restraints `r` that starts
from 6.0 Å to 24 Å in the `pull` phase and stays restrained at 24 Å during the *release* phase.

.. code :: python

    from paprika.restraints import DAT_restraint

    r = DAT_restraint()
    r.mask1 = D1
    r.mask2 = G1
    r.topology = structure
    r.auto_apr = True
    r.continuous_apr = True

    r.attach["target"] = 6.0
    r.attach["fraction_list"] = attach_lambda
    r.attach["fc_final"] = 5.0

    r.pull["target_final"] = 24.0
    r.pull["num_windows"] = len(pull_windows)

    r.release["target"] = 24.0
    r.release["fraction_list"] = [1.0] * len(release_lambda)
    r.release["fc_final"] = 5.0

    r.initialize()

.. note ::

   The ``DAT_restraint`` class can also be used to apply conformational restraints on the host and/or guest molecule.
   For example, distance "jack" and dihedral restraints can be applied to cucurbiturils and cyclodextrins host molecules,
   respectively, to make the binding site more accessible.

**(3) Wall Restraints (optional)**

Wall restraints are half-harmonic potentials that is useful for preventing guest molecules from leaving the binding
site (for weak binding) or preventing the guest molecule from flipping during the attach phase. We still use the
``DAT_restraint`` class to generate the restraints but will use the ``custom_restraint_values`` method to generate
the half-harmonic potential.

.. note ::

   ``custom_restraint_values`` follows the *AMBER* NMR-restraint format, see Chapter 27 in the AMBER20 manual
   for more details.

Below is an example for generating a `"lower wall"` restraint that prevents the angle of ``[D1,G1,G2]`` from
decreasing below 91 degrees.

.. code :: python

    wall_orient = DAT_restraint()
    wall_orient.mask1 = D1
    wall_orient.mask2 = G1
    wall_orient.mask3 = G2
    wall_orient.topology = structure
    wall_orient.auto_apr = True
    wall_orient.continuous_apr = True

    wall_orient.attach["num_windows"] = attach_fractions
    wall_orient.attach["fc_initial"] = 200.0
    wall_orient.attach["fc_final"] = 200.0

    wall_orient.custom_restraint_values["r1"] = 91.0
    wall_orient.custom_restraint_values["r2"] = 0.0
    wall_orient.custom_restraint_values["rk2"] = 200.0
    wall_orient.custom_restraint_values["rk3"] = 0.0

    wall_orient.initialize()


**(4) Positional Restraints**

*Positional restraints* in APR simulations are applied to the dummy atoms. Together with *static restraints*, this
provides a laboratory frame of reference for the host-guest complex. Different MD programs handles `positional restraints`
differently. For example, in ``AMBER`` you can define positional restraints in the input configuration file using the
``ntr`` keyword (Chapter 19 in the AMBER20 manual). For other programs like ``GROMACS`` and ``NAMD`` that uses ``Plumed``,
*positional restraints* can be applied using the method ``add_dummy_atom_restraints()``.

.. note ::

   ``tleap`` may shift the coordinates of the system when it solvates the structure. Applying the *positional restraints*
   before the solvating the structure may lead to undesired errors during simulations. Therefore, special care needs to
   be taken when applying *positional restraints*. Take a look at tutorials `5 <tutorials/05-tutorial-cb6-but-plumed.ipynb>`_
   and `6 <tutorials/06-tutorial-cb6-but-gromacs.ipynb>`_ to see this distinction.


Running a Simulation
--------------------

`paprika` provides wrappers with the ``Simulate`` module for a number of MD engines enabling us to run the simulations
in python.

.. code :: python

   from paprika.simulate import AMBER

   simulation = AMBER()
   simulation.executable = "pmemd.cuda"
   simulation.path = "simulation"
   simulation.prefix = "equilibration"
   simulation.inpcrd = "minimize.rst7"
   simulation.ref = "host-guest-dum.rst7"
   simulation.topology = "host-guest-dum.prmtop"
   simulation.restraint_file = "disang.rest"

   simulation.config_pbc_md()

   # Positional restraints on dummy atoms
   simulation.cntrl["ntr"] = 1
   simulation.cntrl["restraint_wt"] = 50.0
   simulation.cntrl["restraintmask"] = "'@DUM'"

   print(f"Running equilibration in window {window}...")
   simulation.run()

Analysis
--------

Once the simulation is complete, the free energy can be obtained using the ``Analysis`` module, which will also
estimate the uncertainties.

.. code :: python

    from paprika.analysis import fe_calc

    free_energy = fe_calc()
    free_energy.prmtop = "host-guest-dum.prmtop"
    free_energy.trajectory = 'production.nc'
    free_energy.path = "windows"
    free_energy.restraint_list = guest_restraints
    free_energy.collect_data()
    free_energy.methods = ['ti-block']
    free_energy.ti_matrix = "full"
    free_energy.bootcycles = 1000
    free_energy.compute_free_energy()