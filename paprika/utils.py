import logging
import os as os
import shutil
from datetime import datetime
from functools import lru_cache

import numpy as np
import parmed as pmd
import pytraj as pt
from openff.units import unit as openff_unit
from parmed.structure import Structure as ParmedStructureClass

logger = logging.getLogger(__name__)


def get_key(dct, value):
    """
    Get dictionary key given the value.

    .. note ::
        This function will return a list of keys if there are more than one key with the same value in the order they
        are found.

    Parameters
    ----------
    dct: dict
        Python dictionary to get the key from.
    value: str
        The value to match to the dictionary.

    Returns
    -------
    key: str
        The dictionary key that matched the value.


    """
    key = [key for key in dct if (dct[key] == value)]

    if len(key) > 1:
        logger.warning(
            "There more than one key with the same value. Please check if this is the desired output."
        )

    return key


def get_dict_without_keys(dct, *keys):
    """
    Returns a copy of a dictionary without specific keys.

    Parameters
    ----------
    dct: dict
        A python dictionary.
    keys: int or str
        The keys of the dictionary to negate.

    Returns
    -------
    new_dict: dict
        A new dictionary without *keys.

    """
    # https://stackoverflow.com/a/41422467
    new_dict = dict(filter(lambda key_value: key_value[0] not in keys, dct.items()))

    return new_dict


def override_dict(dct, custom):
    """Overrides dictionary values from that of a custom dictionary.

    Parameters
    ----------
    dct: dict
        Python dictionary to override.
    custom: dict
        A custom dictionary which will overwrite the original dictionary.

    """
    for key, value in custom.items():
        if value is not None:
            logger.debug("Overriding {} = {}".format(key, value))
            dct[key] = value


def return_parmed_structure(filename):
    """
    Return a structure object from a filename.

    Parameters
    ----------
    filename : os.PathLike
        The name of the file to load.

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
    structure : :class:`parmed.structure.Structure`
        The structure that contains the atoms
    mask : str
        The atom mask.
    amber_index : bool, optional, default=False
        If true, 1 will be added to the returned indices.

    Returns
    -------
    indices : int
        Atom index or indices corresponding to the mask.

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
    stash_existing : bool, optional, default=False
        Whether to move an existing windows directory to a backup with the current time.
    path : os.PathLike, optional, default='./'
        Root path for the directories.
    window_dir_name : os.PathLike, optional, default='windows'
        Name for the top level directory.
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

    Parameters
    ----------
    prmtop : os.PathLike
        Existing Amber parameter file.
    mask : str, optional, default=':WAT,:Na+,:Cl-'
        The list of atom masks to strip.

    Returns
    -------
    stripped: :class:`pytraj.topology`
        The stripped topology that can be used to read stripped trajectories with :class:`pytraj`.

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
    file : os.PathLike
        Name of the Amber output file.

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
    Return energies from an AMBER ``mdout` file.

    Parameters
    ----------
    file : os.PathLike
        Name of Amber output file

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


def is_file_and_not_empty(file_path):
    """Util function to check if a file both exists at the specified ``path`` and is not empty.

    Parameters
    ----------
    file_path: os.PathLike
        The file path to check.

    Returns
    -------
    bool
        That a file both exists at the specified ``path`` and is not empty.
    """
    # This function is copied from OpenFF-Evaluator
    return os.path.isfile(file_path) and (os.path.getsize(file_path) != 0)


def check_unit(variable, base_unit):
    """
    Run a check if variable is a float or openff_unit.Quantity. If float or int,
    assign the unit of ``base_unit``, if openff_unit.Quantity, check the dimensionality.

    Parameters
    ----------
    variable: int or float or openff_unit.Quantity
    base_unit: openff_unit.Quantity
    """
    quantity = variable
    if isinstance(variable, float) or isinstance(variable, int):
        quantity = openff_unit.Quantity(variable, units=base_unit)
    elif isinstance(variable, openff_unit.Quantity):
        assert variable.dimensionality == base_unit.dimensionality
    elif isinstance(variable, list):
        if all(isinstance(x, float) for x in variable):
            quantity = openff_unit.Quantity(variable, units=base_unit)
        elif all(isinstance(x, openff_unit.Quantity) for x in variable):
            assert all(
                variable[x].dimensionality == base_unit.dimensionality for x in variable
            )
        else:
            raise KeyError(
                "Please make my life easier by either specifying a list of all float or all openff.unit.Quantity."
            )
    elif isinstance(variable, np.ndarray):
        quantity = openff_unit.Quantity(variable, units=base_unit)
    else:
        raise KeyError("``variable`` should be a float or openff.unit.Quantity.")

    return quantity
