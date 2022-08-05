import logging

from openff.units import unit as openff_unit

from paprika.restraints.utils import get_restraint_values, parse_window

logger = logging.getLogger(__name__)


def amber_restraint_line(restraint, window):
    """Writes an `AMBER` NMR-style restraint line for a specific window.

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

    if not restraint.group1:
        iat1 = "{},".format(restraint.index1[0])
    else:
        iat1 = "-1,"
        igr1 = ""
        for index in restraint.index1:
            igr1 += "{},".format(index)

    if not restraint.group2:
        iat2 = "{},".format(restraint.index2[0])
    else:
        iat2 = "-1,"
        igr2 = ""
        for index in restraint.index2:
            igr2 += "{},".format(index)

    iat3 = ""
    if restraint.index3 and not restraint.group3:
        iat3 = "{},".format(restraint.index3[0])
    elif restraint.group3:
        iat3 = "-1,"
        igr3 = ""
        for index in restraint.index3:
            igr3 += "{},".format(index)

    iat4 = ""
    if restraint.index4 and not restraint.group4:
        iat4 = "{},".format(restraint.index4[0])
    elif restraint.group4:
        iat4 = "-1,"
        igr4 = ""
        for index in restraint.index4:
            igr4 += "{},".format(index)

    # Restraint values - Amber NMR-style
    amber_restraint_values = get_restraint_values(restraint, phase, window)

    # Amber units
    energy_unit = openff_unit.kcal / openff_unit.mole
    target_unit = (
        openff_unit.angstrom
        if restraint.restraint_type == "distance"
        else openff_unit.degrees
    )
    force_constant_unit = energy_unit / target_unit**2
    if not restraint.restraint_type == "distance":
        force_constant_unit = energy_unit / openff_unit.radians**2

    # Prepare AMBER NMR-style restraint
    atoms = "".join([iat1, iat2, iat3, iat4])
    string = "&rst iat= {:16s} ".format(atoms)
    string += (
        " r1= {0:10.5f},".format(amber_restraint_values["r1"].to(target_unit).magnitude)
        + " r2= {0:10.5f},".format(
            amber_restraint_values["r2"].to(target_unit).magnitude
        )
        + " r3= {0:10.5f},".format(
            amber_restraint_values["r3"].to(target_unit).magnitude
        )
        + " r4= {0:10.5f},".format(
            amber_restraint_values["r4"].to(target_unit).magnitude
        )
        + " rk2= {0:10.5f},".format(
            amber_restraint_values["rk2"].to(force_constant_unit).magnitude
        )
        + " rk3= {0:10.5f},".format(
            amber_restraint_values["rk3"].to(force_constant_unit).magnitude
        )
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
