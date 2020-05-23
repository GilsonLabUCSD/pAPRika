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


def restraints_from_ascii(filename):
    """
    Utility function to read in restraints from a simple ASCII file.

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


def extract_guest_restraints(structure, guest_resname, restraints):
    """
    Utility function to extract the guest restraints from a list of restraints
    and return individual restraints in the form:
        [r, theta, phi, alpha, beta, gamma]

    If there is no restraint applied to a particular reaction coordinate
    a `None` will be inserted.

    This function is useful for parsing guest restraints in analysis when
    computing `ref_state_work`.

    Parameters
    ----------
    structure : parmed.Structure
        parmed structure of the system.
    guest_resname : str
        Residue name of the guest molecule.
    restraints : list
        list of restraints.

    Returns
    -------
    list
        list of guest-specific DAT_restraint().

    """
    guest_resname = guest_resname.upper()

    r = None
    theta = None
    phi = None
    alpha = None
    beta = None
    gamma = None

    for restraint in restraints:

        mask2_residue_name = structure[restraint.mask2].residues[0].name

        # Distance
        if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and not restraint.mask3 and not \
                restraint.mask4:
            r = restraint

        # Angle
        if restraint.mask3 and not restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name

            if "DM2" in restraint.mask1 and "DM1" in restraint.mask2 and guest_resname in mask3_residue_name:
                theta = restraint

            if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and guest_resname in \
                    mask3_residue_name:
                beta = restraint

        # Dihedral
        if restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name
            mask4_residue_name = structure[restraint.mask4].residues[0].name

            if "DM3" in restraint.mask1 and "DM2" in restraint.mask2 and "DM1" in restraint.mask3 and guest_resname \
                    in mask4_residue_name:
                phi = restraint

            if "DM2" in restraint.mask1 and "DM1" in restraint.mask2 and guest_resname in mask3_residue_name \
                    and guest_resname in mask4_residue_name:
                alpha = restraint

            if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and guest_resname in \
                    mask3_residue_name and guest_resname in mask4_residue_name:
                gamma = restraint

    return [r, theta, phi, alpha, beta, gamma]
