import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

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

def add_restraint(restraint, window, system):
    """
    Apply a :class:`paprika.restraints.DAT_restraint` to an OpenMM System.

    Parameters
    ----------
    restraint : :class:`paprika.restraints.DAT_restraint`
        Restraint to add to the System
    window : str
        Simulation window, which indexes the restraint
    system : :class:`simtk.openmm.System`
        OpenMM System object to be modified

    Returns
    -------
    system  : :class:`simtk.openmm.System`
        OpenMM System object with restraint added

    """

    window, phase = parse_window(window)

    if restraint.index1 an

    # Do we need special group handling code here?

    bond_restraint = mm.CustomBondForce('k * (r - r_0)^2')
    bond_restraint.addPerBondParameter('k')
    bond_restraint.addPerBondParameter('r_0')

    r_0 = restraint.phase[phase]["targets"][window] * unit.angstroms
    k = restraint.phase[phase]["force_constants"][window] * unit.kilocalories_per_mole / unit.angstroms ** 2
    # Make sure these are not called with `amber_index=True`.
    bond_restraint.addBond(restraint.index1, restraint.index2, [k, r_0])
    system.addForce(bond_restraint)

    return system