import logging

from paprika.restraints.utils import parse_window
from paprika.utils import override_dict

logger = logging.getLogger(__name__)


def amber_restraint_line(restraint, window):
    """Writes an AMBER NMR-style restraint line for a specific window.

    Parameters
    ----------
    restraint: :class:`paprika.restraints.DAT_restraint`
        The pAPRika restraint to be used.
    window: str
        The calculation window that will be used to index the restraint values.

    Returns
    -------
    string: str
        A string that can be written to a file that can be read by AMBER.

    Examples
    --------

    .. code-block::

        &rst iat= 3,109, r1= 0.0000, r2= 6.9665, r3= 6.9665, r4= 999.0000,
             rk2= 5.0000000, rk3= 5.0000000, &end

    Or:

    .. code-block::

        &rst iat= -1,-1, igr1=3,4,7,8,21,22,25,26,39,40,43,44,57,58,61,62,75,76,79,
        80,93,94,97,98, igr2=109,113,115,119,
        r1=     0.0000, r2=     5.9665, r3=     5.9665, r4=   999.0000,
        rk2=   5.0000000, rk3=   5.0000000, &end

    """

    window, phase = parse_window(window)
    if (
        restraint.phase[phase]["force_constants"] is None
        and restraint.phase[phase]["targets"] is None
    ):
        return ""

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
    override_dict(amber_restraint_values, restraint.custom_restraint_values)

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
