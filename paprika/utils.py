import logging as log
import os as os
import subprocess as sp

import parmed as pmd
from parmed.structure import Structure as RefStructureForCmp
from paprika import align

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


def index_from_mask(input_structure, mask, amber):
    """
    Return the atom indicies for a given mask.
    """
    if amber:
        index_offset = 1
    else:
        index_offset = 0
    if type(input_structure) is str:
        structure = align.return_structure(input_structure)
    elif type(input_structure) is RefStructureForCmp:
        structure = input_structure
    else:
        raise Exception('index_from_mask does not support the type associated with input_structure')
    # http://parmed.github.io/ParmEd/html/api/parmed/parmed.amber.mask.html?highlight=mask#module-parmed.amber.mask
    indices = [i + index_offset for i in pmd.amber.mask.AmberMask(structure, mask).Selected()]
    log.debug('There are {} atoms in the mask {}  ...'.format(len(indices), mask))
    return indices


def make_window_dirs(window_list):
    """
    Make a series of directories to hold the simulation setup files
    and the data. Here we could check if the directories already exist and prompt
    the user or quit or do something else.
    """

    for window in window_list:
        # It seems unprudent to use '.' for CWD here.
        directory = os.getcwd()
        if not os.path.exists(directory + '/windows/' + window):
            os.makedirs(directory + '/windows/' + window)


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
