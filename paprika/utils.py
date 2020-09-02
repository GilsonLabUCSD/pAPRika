import logging
import os as os
import shutil
from datetime import datetime
from functools import lru_cache

import parmed as pmd
import pytraj as pt
from parmed.structure import Structure as ParmedStructureClass

logger = logging.getLogger(__name__)


def return_parmed_structure(filename):
    """
    Return a structure object from a filename.

    Parameters
    ----------
    filename : str
        The name of the file to load

    Returns
    -------
    structure : :class:`parmed.structure.Structure`

    """
    # `parmed` can read both PDBs and
    # .inpcrd/.prmtop files with the same function call.
    try:
        structure = pmd.load_file(filename)
        logger.info("Loaded {}...".format(filename))
    except IOError:
        logger.error("Unable to load file: {}".format(filename))
    return structure


@lru_cache(maxsize=32)
def index_from_mask(structure, mask, amber_index=False):
    """
    Return the atom indicies for a given selection mask.

    Parameters
    ----------
    structure : `class`:`parmed.structure.Structure`
        The structure that contains the atoms
    mask : str
        The atom mask
    amber_index : bool
        If true, 1 will be added to the returned indices

    Returns
    -------
    indices : int
        Atom index or indices corresponding to the mask

    """
    if amber_index:
        index_offset = 1
    else:
        index_offset = 0
    if isinstance(structure, str):
        structure = return_parmed_structure(structure)
    elif isinstance(structure, ParmedStructureClass):
        pass
    else:
        raise Exception(
            "index_from_mask does not support the type associated with structure:"
            + type(structure)
        )
    # http://parmed.github.io/ParmEd/html/api/parmed/parmed.amber.mask.html?highlight=mask#module-parmed.amber.mask
    indices = [
        i + index_offset for i in pmd.amber.mask.AmberMask(structure, mask).Selected()
    ]
    logger.debug("There are {} atoms in the mask {}  ...".format(len(indices), mask))
    return indices


def make_window_dirs(
        window_list, stash_existing=False, path="./", window_dir_name="windows"
):
    """
    Make a series of windows to hold simulation data.

    Parameters
    ----------
    window_list : list
        List of simulation windows. The names in this list will be used for the folders.
    stash_existing : bool
        Whether to move an existing windows directory to a backup with the current time
    path :
        Root path for the directories
    window_dir_name :
        Name for the top level directory

    Returns
    -------

    """

    win_dir = os.path.join(path, window_dir_name)

    if stash_existing and os.path.isdir(win_dir):
        stash_dir = os.path.join(
            path, window_dir_name + "_{:%Y.%m.%d_%H.%M.%S}".format(datetime.now())
        )
        shutil.move(win_dir, stash_dir)

    for window in window_list:
        window_path = os.path.join(win_dir, window)
        if not os.path.exists(window_path):
            os.makedirs(window_path)


def strip_prmtop(prmtop, mask=":WAT,:Na+,:Cl-"):
    """Strip residues from a structure and write a new parameter file. This could probably also be done with ParmEd.

    Parameters:
    ----------
    prmtop : {str}
        Existing parameter file
    mask : {str}, optional
        The list of atom masks to strip (the default is [':WAT,:Na+,:Cl-'])

    Returns:
    -------
    stripped.topology
        The stripped topology that can be used to read stripped trajectories with `pytraj`

    """

    structure = pt.load_topology(os.path.normpath(prmtop))
    stripped = pt.strip(mask, structure)
    # stripped_name = os.path.join(os.path.splitext(prmtop)[0], '-stripped', os.path.splitext(prmtop)[1])
    # stripped.save(filename=stripped_name)
    # logger.debug('Stripping {} from parameter file and writing {}...'.format(mask, stripped_name))
    return stripped


def parse_mden(file):
    """
    Return energies from an AMBER `mden` file.

    Parameters
    ----------
    file : str
        Name of output file

    Returns
    -------
    energies : dict
        A dictionary containing VDW, electrostatic, bond, angle, dihedral, V14, E14, and total energy.

    """

    vdw, ele, bnd, ang, dih, v14, e14 = [], [], [], [], [], [], []

    with open(file, "r") as f:
        for line in f.readlines()[10:]:
            words = line.rstrip().split()
            if words[0] == "L6":
                vdw.append(float(words[3]))
                ele.append(float(words[4]))
            elif words[0] == "L7":
                bnd.append(float(words[2]))
                ang.append(float(words[3]))
                dih.append(float(words[4]))
            elif words[0] == "L8":
                v14.append(float(words[1]))
                e14.append(float(words[2]))

    energies = {
        "Bond": bnd,
        "Angle": ang,
        "Dihedral": dih,
        "V14": v14,
        "E14": e14,
        "VDW": vdw,
        "Ele": ele,
        "Total": [sum(x) for x in zip(bnd, ang, dih, v14, e14, vdw, ele)],
    }

    return energies


def parse_mdout(file):
    """
    Return energies from an AMBER `mdout` file.

    Parameters
    ----------
    file : str
        Name of output file

    Returns
    -------
    energies : dict
        A dictionary containing VDW, electrostatic, bond, angle, dihedral, V14, E14, and total energy.

    """

    vdw, ele, bnd, ang, dih, v14, e14 = [], [], [], [], [], [], []
    restraint = []

    with open(file, "r") as f:
        for line in f.readlines():
            words = line.rstrip().split()
            if len(words) > 1:
                if "BOND" in words[0]:
                    bnd.append(float(words[2]))
                    ang.append(float(words[5]))
                    dih.append(float(words[8]))
                if "VDWAALS" in words[0]:
                    vdw.append(float(words[2]))
                    ele.append(float(words[5]))
                if "1-4" in words[0]:
                    v14.append(float(words[3]))
                    e14.append(float(words[7]))
                    restraint.append(float(words[10]))

    energies = {
        "Bond": bnd,
        "Angle": ang,
        "Dihedral": dih,
        "V14": v14,
        "E14": e14,
        "VDW": vdw,
        "Ele": ele,
        "Restraint": restraint,
        "Total": [sum(x) for x in zip(bnd, ang, dih, v14, e14, vdw, ele)],
    }

    return energies


def extract_dummy_atoms(structure, resname=None, serial=True):
    """
    Extract information about dummy atoms from a parmed structure and
    returns the information as a dictionary.

    Parameters
    ----------
    structure : :class:`parmed.structure.Structure`
        The parmed structure object we want to extract from
    resname : list
        List of residue name for the dummy atoms (default: ["DM1", "DM2", "DM3"])
    serial : bool
        get indices in serial (starts from 1) or index (starts from 0)

    Returns
    -------
    dummy_atoms : dict

    Output example
    --------------

        main keys: {'DM1', 'DM2', 'DM3'}
        sub keys: 'pos'      - cartesian coordinates (x,y,z)
                  'idx'      - atom indices
                  'idx_type' - type of atom index (serial or index)
                  'mass'     - mass of dummy atom
    """

    if resname is None:
        resname = ["DM1", "DM2", "DM3"]

    dummy_atoms = {name: {} for name in resname}

    for dummy_atom in resname:
        residue = f":{dummy_atom}"
        dummy_atoms[dummy_atom]["pos"] = structure[residue].coordinates[0]
        dummy_atoms[dummy_atom]["mass"] = [atom.mass for atom in structure[residue].atoms][0]
        dummy_atoms[dummy_atom]["idx"] = index_from_mask(structure, residue, amber_index=serial)[0]
        dummy_atoms[dummy_atom]["idx_type"] = "serial" if serial else "index"

    return dummy_atoms
