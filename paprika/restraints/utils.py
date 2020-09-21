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
