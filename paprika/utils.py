import logging as log
import os as os
import subprocess as sp
import shutil
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from paprika import align
from datetime import datetime

global HAS_OPENMM
try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit
    log.debug('OpenMM support: Yes')
    HAS_OPENMM = True
except ImportError:
    log.debug('OpenMM support: No')
    HAS_OPENMM = False


def check_for_leap_log(path='./'):
    """Check if `leap.log` exists, and if so, delete so the current run doesn't append."""
    filename = 'leap.log'
    try:
        os.remove(path + filename)
        log.debug('Deleted existing leap.log file...')
    except OSError:
        pass

def return_parmed_structure(filename):
    """
    Return structure object from file name.
    """
    # `parmed` can read both PDBs and
    # .inpcrd/.prmtop files with the same function call.
    try:
        structure = pmd.load_file(filename)
        log.info('Loaded {}...'.format(filename))
    except:
        log.error('Unable to load file: {}'.format(filename))
    return structure


def index_from_mask(structure, mask, amber_index=False):
    """
    Return the atom indicies for a given mask.
    """
    if amber_index:
        index_offset = 1
    else:
        index_offset = 0
    if type(structure) is str:
        structure = return_parmed_structure(structure)
    elif type(structure) is ParmedStructureClass:
        pass
    else:
        raise Exception('index_from_mask does not support the type associated with structure:'+type(structure))
    # http://parmed.github.io/ParmEd/html/api/parmed/parmed.amber.mask.html?highlight=mask#module-parmed.amber.mask
    indices = [i + index_offset for i in pmd.amber.mask.AmberMask(structure, mask).Selected()]
    log.debug('There are {} atoms in the mask {}  ...'.format(len(indices), mask))
    return indices


def make_window_dirs(window_list, stash_existing=False):
    """
    Make a series of directories to hold the simulation setup files
    and the data. Here we could check if the directories already exist and prompt
    the user or quit or do something else.
    """

    cwd = os.getcwd()
    win_dir = cwd+'/windows'

    if stash_existing and os.path.isdir(win_dir):
        stash_dir = cwd+"/windows_{:%Y.%m.%d_%H.%M.%S}".format(datetime.now())
        shutil.move(win_dir, stash_dir)

    for window in window_list:
        if not os.path.exists(win_dir+'/'+window):
            os.makedirs(win_dir+'/'+window)


def amber_to_pdb(topology, coordinates):
    """
    In certain cases, we may need to convert AMBER-formatted files to
    PDB format, with proper CONECT records. `cpptraj` seems to be the write
    tool for the job at the moment.
    """

    pdb_output = '.'.join(toplogy.split('.')[0:-1]) + '.pdb'
    pdb_input = '.'.join(toplogy.split('.')[0:-1]) + '.in'
    log.info('Converting AMBER coordinates and topology to PDB format.')
    log.debug('Calling `cpptraj`...')
    with open(pdb_input, 'w') as file:
        file.write('parm {}\n'.format(toplogy))
        file.write('trajin {}\n'.format(coordinates))
        file.write('trajout {} conect\n'.format(pdb))
    sp.check_call(['cpptraj', '-i', pdb_input])


def decompose_openmm_energy(simulation, groups=[0, 1], names=['non-restraint', 'restraint']):
    """Return individual energy components.
    """

    energies = dict()
    # Get the total potential energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
    energies.update({'total': energy})

    for index, group in enumerate(groups):
        state = simulation.context.getState(getEnergy=True, groups={group})
        energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
        energies.update({names[index]: energy})

    return energies
