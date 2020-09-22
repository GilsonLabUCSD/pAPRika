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
