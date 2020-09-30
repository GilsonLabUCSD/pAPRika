*****
Usage
*****

Basic Workflow
--------------

todo: workflow-diagram (explain the caveats when using different MD engines).


Structure Preparation
---------------------
The starting structure for APR calculations can be configured with *pAPRika*. The *Align* module provides functions to
shift and orient a structure. For example, we can translate a structure to the origin and the system to the :math:`z`-axis
by running

.. code ::

    from paprika.align import *

    translated_structure = translate_to_origin(structure)
    aligned_structure = zalign(translated_structure, ":GST@C1", ":GST@C2")

Dummy atoms are needed in APR calculations as anchor atoms that defines restraints (see the section below). *pAPRika*
provides utility functions to add dummy atoms to a structure.

.. code ::

    from paprika.dummy import add_dummy

    structure = dummy.add_dummy(structure, residue_name="DM1", z=-6.0)
    structure = dummy.add_dummy(structure, residue_name="DM2", z=-9.0)
    structure = dummy.add_dummy(structure, residue_name="DM3", z=-11.2, y=2.2)
    structure.save("aligned_with_dummy.pdb", overwrite=True)

*pAPRika* provides a wrapper to the ``tleap`` program from AmberTools and we can use this to generate the topology files.

.. code ::

    import paprika.tleap as tleap

    system = tleap.System()
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

In APR calculations we apply restraints on the host (or protein) and the guest molecules. The restraints can be grouped
into three different categories: (1) *static restraints*, (2) *changing restraints* and (3) *wall restraints*.

**(1) Static Restraints**

Static restraints do not change during the whole APR process and do not affect the free-energy. We apply static
restraints on the host (or protein) molecule to define a path for the guest molecule. The static restraints are a made
up of a combination of *distance*, *angle* and *torsional* (DAT) restraints based on the choice of anchor atoms. For
host-guest systems we need to define three anchor atoms ``[H1, H2, H3]`` and combined with three dummy atoms
``[D1, D2, D3]`` we apply a total of six static restraints on the host molecule (three for the translation and three
for orientation).

To generate a static restraint we use the function ``static_DAT_restraints``. As an example, to apply a distance restraints
on ``D1`` and ``H1`` with a force constant of 5 kcal/mol/:math:`Å^2` we call

.. code :: python

    from paprika.restraints import static_DAT_restraint

    dist_static = static_DAT_restraint(
        restraint_mask_list = [D1, H1],
        num_window_list = windows,  # list: [len(attach_lambda), len(pull_windows), len(release_lambda)]
        ref_structure = structure,  # Structure file (PDB)
        force_constant = 5.0,
    )


**(2) Changing Restraints**

As the name suggests, these restraints change during the APR process and can pull the guest molecule out of the host.
However, you can also use these to restraint conformations on the host and/or guest molecule. Unlike static restraints,
these restraints do affect the free-energy.

To generate changing restraints, we use the ``DAT_restraint`` class. The code below shows a changing restraints ``r``
that starts from 6.0 Å to 24 Å in the *pull* phase and stays restrained at 24 Å during the *release* phase.

.. code :: python

    from paprika.restraints import DAT_restraint

    r = DAT_restraint()
    r.mask1 = D1
    r.mask2 = G1
    r.topology = structure
    r.auto_apr = True
    r.continuous_apr = True

    r.attach["target"] = 6.0                     # Angstroms
    r.attach["fraction_list"] = attach_lambda    # Lambda values
    r.attach["fc_final"] = 5.0                   # kcal/mol/Angstroms**2

    r.pull["target_final"] = 24.0                # Angstroms
    r.pull["num_windows"] = len(pull_windows)    # Number of pull windows

    r.release["target"] = 24.0                   # Angstroms
    r.release["fraction_list"] = [1.0] * len(release_lambda)
    r.release["fc_final"] = 5.0                  # kcal/mol/Angstroms**2

    r.initialize()



**(3) Wall Restraints**

Wall restraints are half-harmonic potentials that is useful for preventing guest molecules from leaving the binding site
(for weak binding) or preventing the guest molecule from flipping during the attach phase. We still use the ``DAT_restraint``
class to generate wall restraints but we need to use the method ``custom_restraint_values`` to override build the half-harmonic
potential. Note: ``custom_restraint_values`` follows the *AMBER* NMR-restraint format (see the *AMBER* manual for more details).

Below is an example of generating a `"lower wall"` restraint that prevents the angle of ``[D1, G1, G2]`` from decreasing
below 91 degrees.

.. code :: python

    wall_orient = DAT_restraint()
    wall_orient.mask1 = D1
    wall_orient.mask2 = G1
    wall_orient.mask3 = G2
    wall_orient.topology = structure
    wall_orient.auto_apr = True
    wall_orient.continuous_apr = True

    wall_orient.custom_restraint_values["r1"] = 91.0
    wall_orient.custom_restraint_values["r2"] = 0.0
    wall_orient.custom_restraint_values["rk2"] = kwall
    wall_orient.custom_restraint_values["rk3"] = 0.0


Running a Simulation
--------------------

*pAPRika* provides wrappers for a few Molecular Dynamics (MD) engines and we can run the APR calculations in python

.. code :: python

    from paprika.simulate import Amber

    for window in window_list:
        simulation = Amber()
        simulation.executable = "pmemd"

        simulation.path = f"windows/{window}/"
        simulation.prefix = "production"

        simulation.inpcrd = "minimize.rst7"
        simulation.ref = "host-guest-dum.rst7"
        simulation.topology = "host-guest-dum.prmtop"
        simulation.restraint_file = "disang.rest"

        simulation.config_pbc_md()
        simulation.cntrl["ntr"] = 1
        simulation.cntrl["restraint_wt"] = 50.0
        simulation.cntrl["restraintmask"] = "'@DUM'"

        logging.info(f"Running production in window {window}...")
        simulation.run(overwrite=True)

Analysis
--------


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