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
        raise Exception('index_from_mask does not support the type associated with structure:' + type(structure))
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
    win_dir = cwd + '/windows'

    if stash_existing and os.path.isdir(win_dir):
        stash_dir = cwd + "/windows_{:%Y.%m.%d_%H.%M.%S}".format(datetime.now())
        shutil.move(win_dir, stash_dir)

    for window in window_list:
        if not os.path.exists(win_dir + '/' + window):
            os.makedirs(win_dir + '/' + window)


def create_pdb_with_conect(amber_inpcrd, amber_prmtop, output_pdb, path='./'):
    """
    Create a PDB file containing CONECT records.
    `cpptraj` must be in your PATH.
    Parameters
    ----------
    amber_inpcrd : str
        AMBER (or other) coordinates
    amber_prmtop : str
        AMBER (or other) parameters
    output_pdb : str
        Output PDB file name
    path : str
        Directory for input and output files
    """
    cpptraj = \
        '''
    parm {}
    trajin {}
    trajout {} conect
    '''.format(amber_prmtop, amber_inpcrd, output_pdb)

    cpptraj_input = output_pdb + '.in'
    cpptraj_output = output_pdb + '.out'

    with open(path + cpptraj_input, 'w') as file:
        file.write(cpptraj)
    with open(path + cpptraj_output, 'w') as file:
        p = sp.Popen(['cpptraj', '-i', cpptraj_input], cwd=path, stdout=file, stderr=file)
        output, error = p.communicate()
    if p.returncode == 0:
        log.debug('PDB file written by cpptraj.')
    else:
        log.error('Error returned by cpptraj.')
        log.error(f'Output: {output}')
        log.error(f'Error: {error}')
        p = sp.Popen(['cat', cpptraj_output], cwd=path, stdout=sp.PIPE)
        for line in p.stdout:
            log.error(line.decode("utf-8").strip(), )


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
