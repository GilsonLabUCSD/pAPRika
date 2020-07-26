import logging

from paprika.restraints.utils import parse_window, restraint_to_colvar

logger = logging.getLogger(__name__)

PI = 3.1415926535


def colvar_module_file(file, restraints, window, legacy_k=True):
    """
    Writes a COLVAR Module file for a specific APR window. The COLVAR
    Module is supported in the MD codes NAMD and LAMMPS.

    Giacomo Fiorin, Michael L. Klein and Jérôme Hénin. Using collectie
    variables to drive molecular dynamics simulations. Mol. Phys. 111,
    3345 (2013)

    Parameters
    ----------
    file: class '_io.TextIOWrapper'
        The file object handle to save the plumed file.
    restraints: dict
        The pAPRika restraints in dictionary form.
    window: str
        The APR window that will be used to index the restraint values.
    legacy_k : bool
        Option to specify whether the restraints defined in static_DAT_restraint()
        and DAT_restraint() uses legacy k. Old MD codes like AMBER and CHARMM
        requires the user to multiply the force constant by 1/2 beforehand. New MD
        codes like GROMACS and NAMD requires the user to set the force constant
        without the 1/2 factor. NOTE: The COLVAR Module uses the latter for the
        spring constant.

    Examples
    --------
    Below is an example of converting restraints defined by static_DAT_restraint()
    and DAT_restraint() to a COLVAR Module format. We can use the utility function
    paprika.restraints.utils.parse_restraints to convert the list of restraints to
    a dictionary.

        >>> for window in window_list:
        >>>     with open(f"windows/{window}/plumed.dat", "w") as file:
        >>>      restraints =  parse_restraints(
        >>>         static=static_restraints,
        >>>         host=conformational_restraints,
        >>>         guest=guest_restraints,
        >>>         list_type='dict',
        >>>     )
        >>>     colvar_module_file(file, restraints, window)

    """

    window, phase = parse_window(window)
    file.write("ColvarsTrajFrequency    5000\n")
    file.write("ColvarsRestartFrequency 50000\n")

    if "static" in restraints.keys():
        colvar = restraint_to_colvar(restraints["static"], phase, window, legacy_k=legacy_k, radians=False,)
        write_colvar_module(file, colvar, "static")

    if "host" in restraints.keys():
        colvar = restraint_to_colvar(restraints["host"], phase, window, legacy_k=legacy_k, radians=False)
        write_colvar_module(file, colvar, "host")

    if "guest" in restraints.keys():
        colvar = restraint_to_colvar(restraints["guest"], phase, window, legacy_k=legacy_k, radians=False)
        write_colvar_module(file, colvar, "guest")

    if "wall" in restraints.keys():
        colvar = restraint_to_colvar(restraints["wall"], phase, window, legacy_k=legacy_k, radians=False)
        write_colvar_module(file, colvar, "wall")

    if "symmetry" in restraints.keys():
        colvar = restraint_to_colvar(restraints["symmetry"], phase, window, legacy_k=legacy_k, radians=False)
        write_colvar_module(file, colvar, "symmetry")


def write_colvar_module(file, colvar, label, convert_kangle=True):
    """
    Write collective variable and restraints to file.

    Parameters
    ----------
    file : class '_io.TextIOWrapper'
        The file object handle to save the plumed file.
    colvar : dict
        Dictionary containing information about the collective variable.
    label : str
        Restraint type for naming purposes.
    convert_kangle : bool
        Convert angle force constant from kcal/mol/rad^2 to kcal/mol/deg^2.

    Examples
    --------
    # static restraints

    """
    conversion = 1.0
    if convert_kangle:
        conversion = (PI / 180) ** 2

    file.write(f"# {label} restraints\n")
    arg = ""
    at = ""
    kappa = ""

    # Write COLVAR definition
    for ndx in range(colvar["ncolvar"]):
        colvar_name = f"{label[0]}_{ndx + 1}"
        arg += f"{colvar_name} "
        at += f"{colvar['AT'][ndx]:0.4f} "
        colvar_def = ""

        if colvar["type"][ndx] == "BOND":
            colvar_def = f"colvar {{\n" \
                         f"\tname {colvar_name}\n" \
                         f"\tdistance {{\n" \
                         f"\t\tforceNoPBC yes\n" \
                         f"\t\tgroup1 {{ atomNumbers {{ {colvar['atoms'][ndx][0]} }} }}\n" \
                         f"\t\tgroup2 {{ atomNumbers {{ {colvar['atoms'][ndx][1]} }} }}\n" \
                         f"\t}}\n}}\n"
            kappa += f"{colvar['KAPPA'][ndx] * colvar['factor']:0.5f} "

        elif colvar["type"][ndx] == "ANGLE":
            colvar_def = f"colvar {{\n" \
                         f"\tname {colvar_name}\n" \
                         f"\tangle {{\n" \
                         f"\t\tforceNoPBC yes\n" \
                         f"\t\tgroup1 {{ atomNumbers {{ {colvar['atoms'][ndx][0]} }} }}\n" \
                         f"\t\tgroup2 {{ atomNumbers {{ {colvar['atoms'][ndx][1]} }} }}\n" \
                         f"\t\tgroup3 {{ atomNumbers {{ {colvar['atoms'][ndx][2]} }} }}\n" \
                         f"\t}}\n}}\n"
            kappa += f"{colvar['KAPPA'][ndx] * colvar['factor'] * conversion:0.5f} "

        elif colvar["type"][ndx] == "TORSION":
            colvar_def = f"colvar {{\n" \
                         f"\tname {colvar_name}\n" \
                         f"\tdihedral {{\n" \
                         f"\t\tforceNOPBC yes\n" \
                         f"\t\tgroup1 {{ atomNumbers {{ {colvar['atoms'][ndx][0]} }} }}\n" \
                         f"\t\tgroup2 {{ atomNumbers {{ {colvar['atoms'][ndx][1]} }} }}\n" \
                         f"\t\tgroup3 {{ atomNumbers {{ {colvar['atoms'][ndx][2]} }} }}\n" \
                         f"\t\tgroup4 {{ atomNumbers {{ {colvar['atoms'][ndx][3]} }} }}\n" \
                         f"\t}}\n}}\n"
            kappa += f"{colvar['KAPPA'][ndx] * colvar['factor'] * conversion:0.5f} "

        file.write(colvar_def)

    # Write Harmonic restraints
    if label == "wall" or label == "symmetry":
        bias = "harmonicWalls"
        target = "upperWalls"
        force = "upperWallConstant"
    else:
        bias = "harmonic"
        target = "centers"
        force = "forceConstant"

    harmonic = f"{bias} {{\n" \
               f"\tcolvars {arg}\n" \
               f"\t{target} {at}\n" \
               f"\t{force} {kappa}\n" \
               f"}}\n"

    file.write(harmonic)
