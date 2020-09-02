import logging
import numpy as np

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
    window : int
        The window number
    phase : str
        The calculation phase

    """
    if window[0] == "a":
        phase = "attach"
    elif window[0] == "p":
        phase = "pull"
    elif window[0] == "r":
        phase = "release"
    else:
        raise Exception("Cannot determine the phase for this restraint.")
    window = int(window[1:])

    return window, phase


def get_bias_potential_type(restraint, phase, window):
    """
    Function to determine the bias potential type for a particular restraint.
    The possible types of bias are: "harmonic", "upper_walls" and "lower_walls"

    Parameters
    ----------
    restraint: DAT_restraints
        restraints to extract information from
    phase: str
        current phase of the window
    window: int
        window number

    Returns
    -------
    bias_type: str
        type of bias potential

    """
    bias_type = None

    # Method borrowed from amber_restraint_line module
    amber_restraint_values = {
        "r2": restraint.phase[phase]["targets"][window],
        "r3": restraint.phase[phase]["targets"][window],
        "rk2": restraint.phase[phase]["force_constants"][window],
        "rk3": restraint.phase[phase]["force_constants"][window],
    }
    for key, value in restraint.custom_restraint_values.items():
        if value is not None:
            logger.debug("Overriding {} = {}".format(key, value))
            amber_restraint_values[key] = value

    if amber_restraint_values["r2"] == amber_restraint_values["r3"] and \
            amber_restraint_values["rk2"] == amber_restraint_values["rk3"]:
        bias_type = "harmonic"

    if amber_restraint_values["r2"] < amber_restraint_values["r3"] and \
            amber_restraint_values["rk2"] == amber_restraint_values["rk3"]:
        bias_type = "upper_walls"

    if amber_restraint_values["r2"] == amber_restraint_values["r3"] and \
            amber_restraint_values["rk2"] < amber_restraint_values["rk3"]:
        bias_type = "upper_walls"

    elif amber_restraint_values["r2"] == amber_restraint_values["r3"] and \
            amber_restraint_values["rk2"] > amber_restraint_values["rk3"]:
        bias_type = "lower_walls"

    if bias_type is None:
        raise Exception("Could not determine bias potential type from restraint.")

    return bias_type
