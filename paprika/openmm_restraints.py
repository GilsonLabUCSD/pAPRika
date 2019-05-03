import logging

import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

from paprika.restraints_utilities import parse_window

logger = logging.getLogger(__name__)


def add_restraint(restraint, window, system):
    """
    Apply a :class:`paprika.restraints.DAT_restraint` to an OpenMM System.

    Parameters
    ----------
    restraint : :class:`paprika.restraints.DAT_restraint`
        Restraint to add to the System
    window : str
        Simulation window, which indexes the restraint
    system : :class:`simtk.openmm.System`
        OpenMM System object to be modified

    Returns
    -------
    system  : :class:`simtk.openmm.System`
        OpenMM System object with restraint added

    """

    window, phase = parse_window(window)

    if None in {restraint.index3, restraint.index4}:
        logger.debug("Parsing bond restraint.")
        r_0 = restraint.phase[phase]["targets"][window] * unit.angstroms
        k = (
                restraint.phase[phase]["force_constants"][window]
                * unit.kilocalories_per_mole
                / unit.angstroms ** 2
        )

        if (
            restraint.index1 is not None
            and restraint.index2 is not None
            and 0 in {restraint.group1, restraint.group2}
        ):
            logger.info("Adding CustomBondForce")
            bond_restraint = mm.CustomBondForce("k * (r - r_0)^2")
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])

        else:
            logger.debug("Adding CustomCentroidBondForce")
            bond_restraint = mm.CustomCentroidBondForce("k * (r - r_0)^2")
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")
            bond_restraint.addBond(restraint.index1, restraint.index2, [k, r_0])

    # Make sure these are not called with `amber_index=True`.
    system.addForce(bond_restraint)

    return system


def setup_openmm_restraints(system, restraint, phase, window):
    """
    Add particle restraints with OpenMM.
    """

    # http://docs.openmm.org/7.1.0/api-c++/generated/OpenMM.CustomExternalForce.html
    # It's possible we might need to use `periodicdistance`.

    if (
        restraint.mask1 is not None
        and restraint.mask2 is not None
        and restraint.mask3 is None
        and restraint.mask4 is None
    ):

        if restraint.group1 is False and restraint.group2 is False:
            bond_restraint = mm.CustomBondForce("k * (r - r_0)^2")
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")

            r_0 = restraint.phase[phase]["targets"][window] * unit.angstrom
            k = (
                restraint.phase[phase]["force_constants"][window]
                * unit.kilocalorie_per_mole
                / unit.angstrom ** 2
            )
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])
            bond_restraint.setForceGroup(1)
            system.addForce(bond_restraint)
            logger.debug(
                "Added bond restraint between {} and {} with target value = "
                "{} and force constant = {}".format(
                    restraint.mask1, restraint.mask2, r_0, k
                )
            )
        elif restraint.group1 is True or restraint.group2 is True:
            # http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.openmm.CustomManyParticleForce.html
            # http://getyank.org/development/_modules/yank/restraints.html
            bond_restraint = mm.CustomCentroidBondForce(
                2, "k * (distance(g1, g2) - r_0)^2"
            )
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")

            r_0 = restraint.phase[phase]["targets"][window] * unit.angstrom
            k = (
                restraint.phase[phase]["force_constants"][window]
                * unit.kilocalorie_per_mole
                / unit.angstrom ** 2
            )
            g1 = bond_restraint.addGroup(restraint.index1)
            g2 = bond_restraint.addGroup(restraint.index2)
            bond_restraint.addBond([g1, g2], [k, r_0])
            bond_restraint.setForceGroup(1)
            system.addForce(bond_restraint)
            logger.debug(
                "Added bond restraint between {} and {} with target value = "
                "{} and force constant = {}".format(
                    restraint.mask1, restraint.mask2, r_0, k
                )
            )
    else:
        logger.error("Unable to add bond restraint...")
        logger.debug("restraint.index1 = {}".format(restraint.index1))
        logger.debug("restraint.index2 = {}".format(restraint.index2))
        raise Exception("Unable to add bond restraint...")

    if (
        restraint.mask1 is not None
        and restraint.mask2 is not None
        and restraint.mask3 is not None
        and restraint.mask4 is None
    ):
        if (
            restraint.group1 is not False
            and restraint.group2 is not False
            and restraint.group3 is not False
        ):
            logger.error("Unable to add a group angle restraint...")
            logger.debug("restraint.index1 = {}".format(restraint.index1))
            logger.debug("restraint.index2 = {}".format(restraint.index2))
            logger.debug("restraint.index3 = {}".format(restraint.index3))
            raise Exception("Unable to add a group angle restraint...")

        angle_restraint = mm.CustomAngleForce("k * (theta - theta_0)^2")
        angle_restraint.addPerAngleParameter("k")
        angle_restraint.addPerAngleParameter("theta_0")

        logger.debug(
            "Setting an angle restraint in degrees using a "
            "force constant in kcal per mol rad**2..."
        )
        theta_0 = restraint.phase[phase]["targets"][window] * unit.degrees
        k = (
            restraint.phase[phase]["force_constants"][window]
            * unit.kilocalorie_per_mole
            / unit.radian ** 2
        )
        angle_restraint.addAngle(
            restraint.index1[0], restraint.index2[0], restraint.index3[0], [k, theta_0]
        )
        system.addForce(angle_restraint)

    if (
        restraint.mask1 is not None
        and restraint.mask2 is not None
        and restraint.mask3 is not None
        and restraint.mask4 is not None
    ):
        if (
            restraint.group1 is not False
            and restraint.group2 is not False
            and restraint.group3 is not False
            and restraint.group4 is not False
        ):
            logger.error("Unable to add a group dihedral restraint...")
            logger.debug("restraint.index1 = {}".format(restraint.index1))
            logger.debug("restraint.index2 = {}".format(restraint.index2))
            logger.debug("restraint.index3 = {}".format(restraint.index3))
            logger.debug("restraint.index4 = {}".format(restraint.index4))
            raise Exception("Unable to add a group dihedral restraint...")

        dihedral_restraint = mm.CustomTorsionForce("k * (theta - theta_0)^2")
        dihedral_restraint.addPerTorsionParameter("k")
        dihedral_restraint.addPerTorsionParameter("theta_0")

        logger.debug(
            "Setting a torsion restraint in degrees using a "
            "force constant in kcal per mol rad**2..."
        )
        theta_0 = restraint.phase[phase]["targets"][window] * unit.degrees
        k = (
            restraint.phase[phase]["force_constants"][window]
            * unit.kilocalorie_per_mole
            / unit.radian ** 2
        )
        dihedral_restraint.addTorsion(
            restraint.index1[0],
            restraint.index2[0],
            restraint.index3[0],
            restraint.index4[0],
            [k, theta_0],
        )
        system.addForce(dihedral_restraint)

    return system
