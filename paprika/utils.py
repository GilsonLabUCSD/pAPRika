import logging as log
import os as os
import subprocess as sp

import parmed as pmd
from paprika import align


def index_from_mask(structure_file, mask, index_offset=1):
    """
    Return the atom indicies for a given mask.
    The index_offset keyword sets the index offset, commonly 0 or 1.
    """

    structure = align.return_structure(structure_file)
    # http://parmed.github.io/ParmEd/html/api/parmed/parmed.amber.mask.html?highlight=mask#module-parmed.amber.mask
    indices = [
        i + index_offset for i in pmd.amber.mask.AmberMask(structure, mask).Selected()]
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
