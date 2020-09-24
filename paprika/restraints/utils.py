import logging

import numpy as np

from paprika.utils import override_dict

logger = logging.getLogger(__name__)

_PI_ = np.pi


def parse_window(window):
    """
    Utility function to use a path to index a :class:`paprika.restraints.DAT_restraint` instance.

    Parameters
    ----------
    window : str
        A string representation of a particular simulation window

    Returns
    -------
    window_number : int
        The window number.
    phase : str
        The calculation phase.

    """
    if window[0] == "a":
        phase = "attach"
    elif window[0] == "p":
        phase = "pull"
    elif window[0] == "r":
        phase = "release"
    else:
        raise Exception("Cannot determine the phase for this restraint.")
    window_number = int(window[1:])

    return window_number, phase


def get_restraint_values(restraint, phase, window_number):
    """
    Extract the values of the restraints (Amber NMR-style) including positions
    and force constants. See the Amber documentation for further explanation
    on the NMR-style values.

    Parameters
    ----------
    restraint: :class:`paprika.restraints.DAT_restraint`
        Restraint object to extract the information from
    phase: str
        Current phase of the window.
    window_number: int
        Current window number.

    Returns
    -------
    restraint_values: dict
        Dictionary containing the Amber NMR-style values

    """
    lower_bound = 0.0
    upper_bound = 999.0

    if restraint.mask3 and not restraint.mask4:
        upper_bound = 180.0

    if restraint.mask3 and restraint.mask4:
        lower_bound = restraint.phase[phase]["targets"][window_number] - 180.0
        upper_bound = restraint.phase[phase]["targets"][window_number] + 180.0

    restraint_values = {
        "r1": lower_bound,
        "r2": restraint.phase[phase]["targets"][window_number],
        "r3": restraint.phase[phase]["targets"][window_number],
        "r4": upper_bound,
        "rk2": restraint.phase[phase]["force_constants"][window_number],
        "rk3": restraint.phase[phase]["force_constants"][window_number],
    }

    override_dict(restraint_values, restraint.custom_restraint_values)

    return restraint_values


def get_bias_potential_type(restraint, phase, window_number):
    """
    Function to determine the bias potential type for a particular restraint.
    The possible types of biases are: ``restraint``, ``upper_walls`` and ``lower_walls``.

    Parameters
    ----------
    restraint: :class:`paprika.restraints.DAT_restraint`
        The restraint to extract information from
    phase: str
        Current phase of the window
    window_number: int
        Current window number.

    Returns
    -------
    bias_type: str
        type of bias potential

    """

    bias_type = None

    amber_restraint_values = get_restraint_values(restraint, phase, window_number)

    if (
        amber_restraint_values["r2"] == amber_restraint_values["r3"]
        and amber_restraint_values["rk2"] == amber_restraint_values["rk3"]
    ):
        bias_type = "restraint"

    elif (
        amber_restraint_values["r2"] < amber_restraint_values["r3"]
        or amber_restraint_values["r2"] == 0.0
    ) and amber_restraint_values["rk2"] == amber_restraint_values["rk3"]:
        bias_type = "upper_walls"

    elif amber_restraint_values["r2"] == amber_restraint_values["r3"] and (
        amber_restraint_values["rk2"] < amber_restraint_values["rk3"]
        or amber_restraint_values["rk2"] == 0.0
    ):
        bias_type = "upper_walls"

    elif (
        amber_restraint_values["r2"] < amber_restraint_values["r3"]
        or amber_restraint_values["r2"] == 0.0
    ) and (
        amber_restraint_values["rk2"] <= amber_restraint_values["rk3"]
        or amber_restraint_values["rk2"] == 0.0
    ):
        bias_type = "upper_walls"

    elif (
        amber_restraint_values["r2"] > amber_restraint_values["r3"]
        or amber_restraint_values["r3"] == 0.0
    ) and amber_restraint_values["rk2"] == amber_restraint_values["rk3"]:
        bias_type = "lower_walls"

    elif amber_restraint_values["r2"] == amber_restraint_values["r3"] and (
        amber_restraint_values["rk2"] > amber_restraint_values["rk3"]
        or amber_restraint_values["rk3"] == 0.0
    ):
        bias_type = "lower_walls"

    elif (amber_restraint_values["r2"] < amber_restraint_values["r3"]) and (
        amber_restraint_values["rk2"] > amber_restraint_values["rk3"]
        or amber_restraint_values["rk3"] == 0.0
    ):
        bias_type = "lower_walls"

    if bias_type is None:
        raise Exception("Could not determine bias potential type from restraint.")

    return bias_type


def extract_guest_restraints(structure, restraints, guest_resname, dummy_prefix="DM"):
    """
    Utility function to extract the guest restraints from a list of restraints
    and return individual restraints in the form ``[r, theta, phi, alpha, beta, gamma]``.

    If there is no restraint applied to a particular reaction coordinate a
    ``None`` will be inserted.

    This function is useful for parsing guest restraints during analysis when
    computing ``ref_state_work``.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        The molecular structure of the system.
    restraints : list
        List of :class:`paprika.restraints.DAT_restraint` restraints.
    guest_resname : str
        Residue name of the guest molecule.
    dummy_prefix : str, optional, default="DM"
        The prefix for the dummy atoms residue name.

    Returns
    -------
    guest_restraints : list
        list of guest-specific :class:`paprika.restraints.DAT_restraint` restraints.

    Examples
    --------
        Extract guest restraints from a list of restraints and compute the work for reference state

        >>> guest_restraints = extract_guest_restraints(structure, restraints, "BEN")
        >>>
        >>> free_energy = analysis.fe_calc()
        >>> free_energy.compute_ref_state_work(guest_restraints)
        >>> print(free_energy["ref_state_work"])

    """
    guest_resname = guest_resname.upper()

    DM1 = f"{dummy_prefix}1"
    DM2 = f"{dummy_prefix}2"
    DM3 = f"{dummy_prefix}3"

    r = None
    theta = None
    phi = None
    alpha = None
    beta = None
    gamma = None

    for restraint in restraints:

        mask2_residue_name = structure[restraint.mask2].residues[0].name

        # Distance
        if (
            DM1 in restraint.mask1
            and guest_resname in mask2_residue_name
            and not restraint.mask3
            and not restraint.mask4
        ):
            r = restraint

        # Angle
        if restraint.mask3 and not restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name

            if (
                DM2 in restraint.mask1
                and DM1 in restraint.mask2
                and guest_resname in mask3_residue_name
            ):
                theta = restraint

            if (
                DM1 in restraint.mask1
                and guest_resname in mask2_residue_name
                and guest_resname in mask3_residue_name
            ):
                beta = restraint

        # Dihedral
        if restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name
            mask4_residue_name = structure[restraint.mask4].residues[0].name

            if (
                DM3 in restraint.mask1
                and DM2 in restraint.mask2
                and DM1 in restraint.mask3
                and guest_resname in mask4_residue_name
            ):
                phi = restraint

            if (
                DM2 in restraint.mask1
                and DM1 in restraint.mask2
                and guest_resname in mask3_residue_name
                and guest_resname in mask4_residue_name
            ):
                alpha = restraint

            if (
                DM1 in restraint.mask1
                and guest_resname in mask2_residue_name
                and guest_resname in mask3_residue_name
                and guest_resname in mask4_residue_name
            ):
                gamma = restraint

    guest_restraints = [r, theta, phi, alpha, beta, gamma]

    return guest_restraints


def restraints_from_ascii(filename):
    """
    Utility function to read in restraints from a simple ASCII file. This is
    useful when parsing restraints definition from VMD. Since you can mouse-click
    bonds, angles and dihedrals in VMD it can be faster to define these restraints.
    One scenario for this is when you want to apply conformational restraints on a
    host molecule based on a number of distance, angle and torsional restraints.

    The TCL script below is an example of extracting bonds selected in VMD to a file
    with a comma-separated format.

    .. code-block:: tcl

        set kbond 10.0
        set f [open "bond.dat" w]

        foreach bond [label list Bonds] {
            # Bond atom indices
            set i1 [lindex [split [lindex $bond 0]] 1]
            set i2 [lindex [split [lindex $bond 1]] 1]

            # Bond atom names
            set a1 [[atomselect top "index $i1"] get name]
            set a2 [[atomselect top "index $i2"] get name]

            # Bond atom resid
            set r1 [[atomselect top "index $i1"] get resname]
            set r2 [[atomselect top "index $i2"] get resname]

            # Bond value
            set dr [lindex [split [lindex $bond 2]] 1]

            # Print atom names
            puts $f ":${r1}@${a1},:${r2}@${a2},$dr,$kbond"
        }
        close $f

    The output from the TCL script sourced in VMD may look like

    .. code-block::

        :1@C12,:2@C4,12.5,5.0
        :1@O11,:1@C2,:1@C3,90.0,50.0
        :1@C12,:1@O13,:1@C14,1@C15,-121.16,6.0

    Parameters
    ----------
    filename : os.PathLike
        File name of template file.

    Returns
    -------
    restraints : dict
        A dictionary of restraints containing information of the atoms, target and spring constant.

    """
    restraints = {"atoms": [], "target": [], "k": [], "type": []}

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split(",")

                # Distance
                if len(line) == 4:
                    restraints["atoms"].append([line[0], line[1]])
                    restraints["target"].append(float(line[2]))
                    restraints["k"].append(float(line[3]))
                    restraints["type"].append("bond")

                # Angle
                elif len(line) == 5:
                    restraints["atoms"].append([line[0], line[1], line[2]])
                    restraints["target"].append(float(line[3]))
                    restraints["k"].append(float(line[4]))
                    restraints["type"].append("angle")

                # Torsion
                elif len(line) == 6:
                    restraints["atoms"].append([line[0], line[1], line[2], line[3]])
                    restraints["target"].append(float(line[4]))
                    restraints["k"].append(float(line[5]))
                    restraints["type"].append("dihedral")

                else:
                    logger.info(
                        "Restraint given is not a bond, angle or dihedral... skipping line."
                    )

    return restraints
