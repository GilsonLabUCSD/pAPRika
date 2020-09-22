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

    Return
    ------
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
        restraints to extract information from
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


def parse_restraints(static=None, host=None, guest=None, wall=None, symmetry=None, list_type='tuple'):
    """
    Utility function to parse restraints that is used when writing the
    restraints to file. If list_type='tuple' (default) the function simply
    returns a list of all the restraints, which is required for writing
    amber NMR restraint file (paprika.restraints.amber.amber_restraint_line).
    The option list_type='dict' is useful for parsing pAPRika restraints
    into a PLUMED-based restraint file.

    Parameters
    ----------
    static : list
        List of host static DAT_restraint()
    guest : list
        List of DAT_restraint() for guest
    host : list
        List of DAT_restraint() for host-conformation
    wall : list
        List of DAT_restraint() for guest-wall
    symmetry : list
        List of DAT_restraint() for guest-symmetry
    list_type : str
        Type of list to return (tuple or dict)

    Returns
    -------
    restraints_list : tuple/dict
        The list of restraints returned as a tuple or dictionary

    """

    if list_type is 'tuple':
        restraints_list = []
    elif list_type is 'dict':
        restraints_list = {}

    for restraint in ["static", "host", "guest", "wall", "symmetry"]:
        if eval(restraint) is not None:
            if len(eval(restraint)) != 0:
                if list_type is 'tuple':
                    restraints_list += eval(restraint)
                elif list_type is 'dict':
                    restraints_list[restraint] = eval(restraint)

    return restraints_list


def restraints_from_ascii(filename):
    """
    Utility function to read in restraints from a simple ASCII file. This is
    useful when parsing restraints definition from VMD (you can mouse click
    bond, angle and dihedral restraints and write a TCL script to print these
    to a file).

    Parameters
    ----------
    filename : str
        file name of template file.

    Returns
    -------
    restraints : dict
        dictionary of restraints containing information of the atoms, target and spring constant.

    Examples
    --------
    ASCII file should contain the atoms (2, 3 or 4), the equilibrium target and spring constant.

        :1@C12 :2@C4  12.5 5.0
        :1@O11 :1@C2 :1@C3  90.0 50.0
        :1@C12 :1@O13 :1@C14 :1@C15 -121.16 6.0

    """
    restraints = {'atoms': [], 'target': [], 'k': [], 'type': []}

    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split()

                if len(line) == 4:
                    restraints['atoms'].append([line[0], line[1]])
                    restraints['target'].append(float(line[2]))
                    restraints['k'].append(float(line[3]))
                    restraints['type'].append('bond')

                elif len(line) == 5:
                    restraints['atoms'].append([line[0], line[1], line[2]])
                    restraints['target'].append(float(line[3]))
                    restraints['k'].append(float(line[4]))
                    restraints['type'].append('angle')

                elif len(line) == 6:
                    restraints['atoms'].append([line[0], line[1], line[2], line[3]])
                    restraints['target'].append(float(line[4]))
                    restraints['k'].append(float(line[5]))
                    restraints['type'].append('dihedral')

                else:
                    print("Restraint given is not a bond, angle or dihedral... skipping line.")

    return restraints
