import logging
import os as os
import shutil
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
import pytraj as pt
from datetime import datetime

logger = logging.getLogger(__name__)

def return_parmed_structure(filename):
    """
    Return structure object from file name.
    """
    # `parmed` can read both PDBs and
    # .inpcrd/.prmtop files with the same function call.
    try:
        structure = pmd.load_file(filename)
        logger.info("Loaded {}...".format(filename))
    except BaseException:
        logger.error("Unable to load file: {}".format(filename))
    return structure


def index_from_mask(structure, mask, amber_index=False):
    """
    Return the atom indicies for a given mask.
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
    Make a series of directories to hold the simulation setup files
    and the data. Here we could check if the directories already exist and prompt
    the user or quit or do something else.
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


def decompose_openmm_energy(
    simulation, groups=[0, 1], names=["non-restraint", "restraint"]
):
    """Return individual energy components.
    """

    energies = dict()
    # Get the total potential energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
    energies.update({"total": energy})

    for index, group in enumerate(groups):
        state = simulation.context.getState(getEnergy=True, groups={group})
        energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
        energies.update({names[index]: energy})

    return energies


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
