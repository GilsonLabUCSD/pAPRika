import logging as log
import numpy as np
import parmed as pmd

try:
    import simtk.openmm as mm
    import simtk.unit as unit
except ImportError:
    pass


def check_restraints(restraint_list, create_window_list=False):
    """
    Do basic tests to ensure a list of DAT_restraints are consistent.
    We're gonna create the window list here too, because it needs the same code.
    """

    if all(restraint.continuous_apr is True for restraint in restraint_list):
        log.debug('All restraints are "continuous_apr" style.')
        all_continuous_apr = True
    elif all(restraint.continuous_apr is False for restraint in restraint_list):
        log.debug('All restraints are not "continuous_apr" style.')
        all_continuous_apr = False
    else:
        log.error("All restraints must have the same setting for .continuous_apr")
        # Should we do the following?
        raise Exception(
            "All restraints must have the same setting for .continuous_apr")

    window_list = []
    phases = ["attach", "pull", "release"]
    for phase in phases:
        win_counts = []
        for restraint in restraint_list:
            if restraint.phase[phase]["targets"] is not None:
                win_counts.append(len(restraint.phase[phase]["targets"]))
            else:
                win_counts.append(0)
        max_count = np.max(win_counts)
        if all(count == 0 or count == max_count for count in win_counts):
            if max_count > 0:
                if all_continuous_apr and phase == "attach":
                    max_count -= 1
                if max_count > 999:
                    log.info(
                        "Window name zero padding only applied for windows 0 - 999"
                    )
                # `continuous_apr` during attach means that the final attach window should be skipped and replaced
                # with `p000`.
                window_list += [
                    phase[0] + str("{:03.0f}".format(val))
                    for val in np.arange(0, max_count, 1)
                ]
                if all_continuous_apr and phase == "release":
                    max_count -= 1
                # `continuous_apr` during release means that `r000` should be skipped and replaced with the final
                # pull window.
                window_list += [
                    phase[0] + str("{:03.0f}".format(val))
                    for val in np.arange(1, max_count, 1)
                ]
        else:
            log.error(
                "Restraints have unequal number of windows during the {} phase.".format(
                    phase
                )
            )
            log.debug("Window counts for each restraint are as follows:")
            log.debug(win_counts)
            raise Exception(
                "Restraints have unequal number of windows during the {} phase.".format(
                    phase
                )
            )

    log.info("Restraints appear to be consistent")

    if create_window_list:
        return window_list


def create_window_list(restraint_list):
    """
    Create list of APR windows. Runs everything through check_restraints because
    we need to do that.
    """

    return check_restraints(restraint_list, create_window_list=True)


def amber_restraint_line(restraint, phase, window):
    """
    Return an AMBER restraint line a specific phase/window combination.
    For example:
    &rst iat= 3,109, r1= 0.0000, r2= 6.9665, r3= 6.9665, r4= 999.0000,
         rk2= 5.0000000, rk3= 5.0000000, &end
    Or:
    &rst iat= -1,-1, igr1=3,4,7,8,21,22,25,26,39,40,43,44,57,58,61,62,75,76,79,
    80,93,94,97,98, igr2=109,113,115,119,
    r1=     0.0000, r2=     5.9665, r3=     5.9665, r4=   999.0000,
    rk2=   5.0000000, rk3=   5.0000000, &end

    """
    if not restraint.index1:
        iat1 = " "
        raise Exception("There must be at least two atoms in a restraint.")
    elif not restraint.group1:
        iat1 = "{},".format(restraint.index1[0])
    else:
        iat1 = "-1,"
        igr1 = ""
        for index in restraint.index1:
            igr1 += "{},".format(index)

    if not restraint.index2:
        iat2 = " "
        raise Exception("There must be at least two atoms in a restraint.")
    elif not restraint.group2:
        iat2 = "{},".format(restraint.index2[0])
    else:
        iat2 = "-1,"
        igr2 = ""
        for index in restraint.index2:
            igr2 += "{},".format(index)

    if not restraint.index3:
        iat3 = ""
    elif not restraint.group3:
        iat3 = "{},".format(restraint.index3[0])
    else:
        iat3 = "-1,"
        igr3 = ""
        for index in restraint.index3:
            igr3 += "{},".format(index)

    if not restraint.index4:
        iat4 = ""
    elif not restraint.group4:
        iat4 = "{},".format(restraint.index4[0])
    else:
        iat4 = "-1,"
        igr4 = ""
        for index in restraint.index4:
            igr4 += "{},".format(index)

    # Set upper/lower bounds depending on whether distance, angle, or torsion
    lower_bound = 0.0
    upper_bound = 999.0

    if restraint.mask3 and not restraint.mask4:
        upper_bound = 180.0

    if restraint.mask3 and restraint.mask4:
        lower_bound = restraint.phase[phase]["targets"][window] - 180.0
        upper_bound = restraint.phase[phase]["targets"][window] + 180.0

    amber_restraint_values = {
        "r1": lower_bound,
        "r2": restraint.phase[phase]["targets"][window],
        "r3": restraint.phase[phase]["targets"][window],
        "r4": upper_bound,
        "rk2": restraint.phase[phase]["force_constants"][window],
        "rk3": restraint.phase[phase]["force_constants"][window],
    }

    for key, value in restraint.custom_restraint_values.items():
        if value is not None:
            log.debug("Overriding {} = {}".format(key, value))
            amber_restraint_values[key] = value

    # Prepare AMBER NMR-style restraint
    atoms = "".join([iat1, iat2, iat3, iat4])
    string = "&rst iat= {:16s} ".format(atoms)
    string += (
        " r1= {0:10.5f},".format(amber_restraint_values["r1"])
        + " r2= {0:10.5f},".format(amber_restraint_values["r2"])
        + " r3= {0:10.5f},".format(amber_restraint_values["r3"])
        + " r4= {0:10.5f},".format(amber_restraint_values["r4"])
        + " rk2= {0:10.5f},".format(amber_restraint_values["rk2"])
        + " rk3= {0:10.5f},".format(amber_restraint_values["rk3"])
    )

    if any([restraint.group1, restraint.group2, restraint.group3, restraint.group4]):
        string += "\n    "
        if restraint.group1:
            string += " igr1= {}".format(igr1)
        if restraint.group2:
            string += " igr2= {}".format(igr2)
        if restraint.group3:
            string += " igr3= {}".format(igr3)
        if restraint.group4:
            string += " igr4= {}".format(igr4)

    string += "  &end\n"
    return string


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
            bond_restraint.addBond(
                restraint.index1[0], restraint.index2[0], [k, r_0])
            bond_restraint.setForceGroup(1)
            system.addForce(bond_restraint)
            log.debug(
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
            log.debug(
                "Added bond restraint between {} and {} with target value = "
                "{} and force constant = {}".format(
                    restraint.mask1, restraint.mask2, r_0, k
                )
            )
    else:
        log.error("Unable to add bond restraint...")
        log.debug("restraint.index1 = {}".format(restraint.index1))
        log.debug("restraint.index2 = {}".format(restraint.index2))
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
            log.error("Unable to add a group angle restraint...")
            log.debug("restraint.index1 = {}".format(restraint.index1))
            log.debug("restraint.index2 = {}".format(restraint.index2))
            log.debug("restraint.index3 = {}".format(restraint.index3))
            raise Exception("Unable to add a group angle restraint...")

        angle_restraint = mm.CustomAngleForce("k * (theta - theta_0)^2")
        angle_restraint.addPerAngleParameter("k")
        angle_restraint.addPerAngleParameter("theta_0")

        log.debug(
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
            restraint.index1[0], restraint.index2[0], restraint.index3[0], [
                k, theta_0]
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
            log.error("Unable to add a group dihedral restraint...")
            log.debug("restraint.index1 = {}".format(restraint.index1))
            log.debug("restraint.index2 = {}".format(restraint.index2))
            log.debug("restraint.index3 = {}".format(restraint.index3))
            log.debug("restraint.index4 = {}".format(restraint.index4))
            raise Exception("Unable to add a group dihedral restraint...")

        dihedral_restraint = mm.CustomTorsionForce("k * (theta - theta_0)^2")
        dihedral_restraint.addPerTorsionParameter("k")
        dihedral_restraint.addPerTorsionParameter("theta_0")

        log.debug(
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


def clean_restraints_file(restraints, filename="restraints.in"):
    """
    Delete the restraints files for repeated testing.

    Parameters
    ----------
    restraints : object
    """

    log.warning("`clean_restraints_file()` needs to be tested.")
    for restraint in restraints:
        for window, _ in enumerate(restraint.phase["attach"]["force_constants"]):
            directory = "./windows/a{0:03d}".format(window)
            os.remove(directory + "/" + filename)

    for restraint in restraints:
        for window, _ in enumerate(restraint.phase["attach"]["targets"]):
            directory = "./windows/p{0:03d}".format(window)
            os.remove(directory + "/" + filename)


def to_json(restraint_list):
    for restraint in restraint_list:
        dictionary = restraint.__dict__
        # Get rid of the `topology` key if it is a ParmEd AmberParm, which is not natively
        # serializable (although this could be unpacked itself into another dictionary).
        print(dictionary)
        if isinstance(dictionary["topology"], pmd.amber.AmberParm):
            print("Removing pmd.AmberParm from json")
            dictionary["topology"] = None
        for phase in ["attach", "pull", "release"]:
            for array in ["force_constants", "targets"]:
                try:
                    dictionary["phase"][phase][array] = self.phase[array].tolist()
                except:
                    print("Could not convert {} to list".format(
                        self.phase[array]))
        return json.dumps(dictionary)


def from_json(filename):
    pass
