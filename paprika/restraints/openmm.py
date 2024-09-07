"""A module aimed at applying restraints directly to OpenMM systems."""

import logging

import numpy as np

try:
    import openmm
    import openmm.unit as openmm_unit
except ImportError:
    import simtk.openmm as openmm
    import simtk.unit as openmm_unit

from typing import Optional, Union

import parmed as pmd
from openff.units import unit as openff_unit
from openff.units.openmm import to_openmm

from paprika.restraints import DAT_restraint

logger = logging.getLogger(__name__)
_PI_ = np.pi


def apply_positional_restraints(
    coordinate_path: str,
    system: openmm.System,
    atom_name: Optional[str] = "DUM",
    force_group: Optional[int] = 15,
    k_pos: Optional[Union[float, openmm_unit.Quantity, openff_unit.Quantity]] = 50.0,
):
    """A utility function which will add OpenMM harmonic positional restraints to
    any dummy atoms found within a system to restrain them to their initial
    positions.

    Parameters
    ----------
    coordinate_path : str
        The path to the coordinate file which the restraints will be applied to.
        This should contain either the host or the complex, the dummy atoms and
        solvent.
    system : :class:`openmm.System`
        The system object to add the positional restraints to.
    atom_name : str
        The name of the atom to restrain.
    force_group : int, optional
        The force group to add the positional restraints to.
    k_pos : float, openmm.unit.Quantity or openff.unit.Quantity, optional
        The force constant for restraining the dummy atoms (kcal/mol/Å^2 if float).
    """

    # noinspection PyTypeChecker
    structure: pmd.Structure = pmd.load_file(coordinate_path, structure=True)

    for atom in structure.atoms:
        if atom.name == atom_name:
            positional_restraint = openmm.CustomExternalForce(
                "k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
            )
            positional_restraint.addPerParticleParameter("k")
            positional_restraint.addPerParticleParameter("x0")
            positional_restraint.addPerParticleParameter("y0")
            positional_restraint.addPerParticleParameter("z0")

            # I haven't found a way to get this to use ParmEd's unit library here.
            # ParmEd correctly reports `atom.positions` as units of Ångstroms.
            # But then we can't access atom indices. Using `atom.xx` works for
            # coordinates, but is unitless.
            if isinstance(k_pos, float):
                k = k_pos * openmm_unit.kilocalories_per_mole / openmm_unit.angstroms**2
            elif isinstance(k_pos, openmm_unit.Quantity):
                k = k_pos
            elif isinstance(k_pos, openff_unit.Quantity):
                k = to_openmm(k_pos)

            # Convert position to nanometers
            x0 = 0.1 * atom.xx * openmm_unit.nanometers
            y0 = 0.1 * atom.xy * openmm_unit.nanometers
            z0 = 0.1 * atom.xz * openmm_unit.nanometers

            positional_restraint.addParticle(atom.idx, [k, x0, y0, z0])
            system.addForce(positional_restraint)
            positional_restraint.setForceGroup(force_group)


def apply_dat_restraint(
    system: openmm.System,
    restraint: DAT_restraint,
    phase: str,
    window_number: int,
    force_group: Optional[int] = None,
):
    """A utility function which takes in pAPRika DAT restraints and applies
    the restraints to an OpenMM System object.

    Parameters
    ----------
    system : :class:`openmm.System`
        The system object to add the positional restraints to.
    restraint : :class:`DAT_restraint`
        DAT_restraint to apply to `openmm.System`.
    phase : str
        Phase of calculation ("attach", "pull" or "release")
    window_number : int
        The corresponding window number of the current phase
    force_group : int, optional
        The force group to add the positional restraints to.

    """
    from paprika.restraints import BiasPotentialType, RestraintType
    from paprika.restraints.utils import get_bias_potential_type

    assert phase in {"attach", "pull", "release"}

    # Get bias potential type - `Harmonic`, `UpperWall`, `LowerWall`
    bias_type, restraint_values = get_bias_potential_type(
        restraint, phase, window_number, return_values=True
    )
    if bias_type is None:
        pass

    # Distance restraints
    if restraint.restraint_type == RestraintType.Distance:
        # colvar values
        r_0 = to_openmm(restraint.phase[phase]["targets"][window_number])
        k = to_openmm(restraint.phase[phase]["force_constants"][window_number])

        # Energy expression
        bond_energy_expression = "k_bond * (r - r_0)^2;"
        if bias_type == BiasPotentialType.UpperWall:
            bond_energy_expression = "step(r - r_0) * " + bond_energy_expression
            r_0 = to_openmm(restraint_values["r3"])
            k = to_openmm(restraint_values["rk3"])
        elif bias_type == BiasPotentialType.LowerWall:
            bond_energy_expression = "step(r_0 - r) * " + bond_energy_expression
            r_0 = to_openmm(restraint_values["r2"])
            k = to_openmm(restraint_values["rk2"])

        # Single particles
        if not restraint.group1 and not restraint.group2:
            bond_restraint = openmm.CustomBondForce(bond_energy_expression)
            bond_restraint.addPerBondParameter("k_bond")
            bond_restraint.addPerBondParameter("r_0")
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])

        # Group of particles
        else:
            bond_restraint = openmm.CustomCentroidBondForce(
                2, bond_energy_expression + "r = distance(g1, g2);"
            )
            bond_restraint.addPerBondParameter("k_bond")
            bond_restraint.addPerBondParameter("r_0")

            # Need to use geometrical instead of mass centers (all weights=1.0)
            g1 = bond_restraint.addGroup(
                restraint.index1, [1.0 for i in range(len(restraint.index1))]
            )
            g2 = bond_restraint.addGroup(
                restraint.index2, [1.0 for i in range(len(restraint.index2))]
            )

            bond_restraint.addBond([g1, g2], [k, r_0])

        system.addForce(bond_restraint)

        if force_group:
            bond_restraint.setForceGroup(force_group)

    # Angle restraints
    elif restraint.restraint_type == RestraintType.Angle:
        # colvar values
        theta_0 = to_openmm(restraint.phase[phase]["targets"][window_number])
        k = to_openmm(restraint.phase[phase]["force_constants"][window_number])

        # energy expression
        angle_energy_expression = "k_angle * (theta - theta_0)^2;"
        if bias_type == BiasPotentialType.UpperWall:
            angle_energy_expression = (
                "step(theta - theta_0) * "
            ) + angle_energy_expression
            theta_0 = to_openmm(restraint_values["r3"])
            k = to_openmm(restraint_values["rk3"])
        elif bias_type == BiasPotentialType.LowerWall:
            angle_energy_expression = (
                "step(theta_0 - theta) * "
            ) + angle_energy_expression
            theta_0 = to_openmm(restraint_values["r2"])
            k = to_openmm(restraint_values["rk2"])

        # Single particles
        if not restraint.group1 and not restraint.group2 and not restraint.group3:
            angle_restraint = openmm.CustomAngleForce(angle_energy_expression)
            angle_restraint.addPerAngleParameter("k_angle")
            angle_restraint.addPerAngleParameter("theta_0")
            angle_restraint.addAngle(
                restraint.index1[0],
                restraint.index2[0],
                restraint.index3[0],
                [k, theta_0],
            )

        # Group of particles
        else:
            angle_restraint = openmm.CustomCentroidBondForce(
                3, angle_energy_expression + "theta = angle(g1, g2, g3);"
            )
            angle_restraint.addPerBondParameter("k_angle")
            angle_restraint.addPerBondParameter("theta_0")

            # Need to use geometrical instead of mass centers (all weights=1.0)
            g1 = angle_restraint.addGroup(
                restraint.index1, [1.0 for i in range(len(restraint.index1))]
            )
            g2 = angle_restraint.addGroup(
                restraint.index2, [1.0 for i in range(len(restraint.index2))]
            )
            g3 = angle_restraint.addGroup(
                restraint.index3, [1.0 for i in range(len(restraint.index3))]
            )

            angle_restraint.addBond([g1, g2, g3], [k, theta_0])

        system.addForce(angle_restraint)

        if force_group:
            angle_restraint.setForceGroup(force_group)

    # Torsion restraints
    elif restraint.restraint_type == RestraintType.Torsion:
        # colvar values
        theta_0 = to_openmm(restraint.phase[phase]["targets"][window_number])
        k = to_openmm(restraint.phase[phase]["force_constants"][window_number])

        # energy expression
        dihedral_energy_expression = (
            f"k_torsion * min(min(abs(theta - theta_0), abs(theta - theta_0 + "
            f"2 * {_PI_})), abs(theta - theta_0 - 2 * {_PI_}))^2;"
        )
        if bias_type == BiasPotentialType.UpperWall:
            dihedral_energy_expression = (
                "step(theta - theta_0) * "
            ) + dihedral_energy_expression
            theta_0 = to_openmm(restraint_values["r3"])
            k = to_openmm(restraint_values["rk3"])
        elif bias_type == BiasPotentialType.LowerWall:
            dihedral_energy_expression = (
                "step(theta_0 - theta) * "
            ) + dihedral_energy_expression
            theta_0 = to_openmm(restraint_values["r2"])
            k = to_openmm(restraint_values["rk2"])

        # Single particles
        if (
            not restraint.group1
            and not restraint.group2
            and not restraint.group3
            and not restraint.group4
        ):
            dihedral_restraint = openmm.CustomTorsionForce(dihedral_energy_expression)
            dihedral_restraint.addPerTorsionParameter("k_torsion")
            dihedral_restraint.addPerTorsionParameter("theta_0")
            dihedral_restraint.addTorsion(
                restraint.index1[0],
                restraint.index2[0],
                restraint.index3[0],
                restraint.index4[0],
                [k, theta_0],
            )

        # Group of particles
        else:
            dihedral_restraint = openmm.CustomCentroidBondForce(
                4, dihedral_energy_expression + "theta = dihedral(g1,g2,g3,g4);"
            )
            dihedral_restraint.addPerBondParameter("k_torsion")
            dihedral_restraint.addPerBondParameter("theta_0")

            # Need to use geometrical instead of mass centers (all weights=1.0)
            g1 = dihedral_restraint.addGroup(
                restraint.index1, [1.0 for i in range(len(restraint.index1))]
            )
            g2 = dihedral_restraint.addGroup(
                restraint.index2, [1.0 for i in range(len(restraint.index2))]
            )
            g3 = dihedral_restraint.addGroup(
                restraint.index3, [1.0 for i in range(len(restraint.index3))]
            )
            g4 = dihedral_restraint.addGroup(
                restraint.index4, [1.0 for i in range(len(restraint.index4))]
            )

            dihedral_restraint.addBond([g1, g2, g3, g4], [k, theta_0])

        system.addForce(dihedral_restraint)

        if force_group:
            dihedral_restraint.setForceGroup(force_group)
