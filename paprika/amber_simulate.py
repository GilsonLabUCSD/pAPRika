import logging as log
import subprocess as sp
from collections import OrderedDict

import numpy as np
import parmed as pmd
from paprika.restraints import *


def _amber_write_input_file(filename, dictionary, title='Input.'):
    """
    Process a dictionary to generate an AMBER input file.

    Parameters
    ---------
    filename : string
        Name of the AMBER input file (to be written)
    dictionary : dict
        Python dictionary containing `&cntrl` keywords
    title : str
        Title of the AMBER input file
    """

    log.debug('Writing {}'.format(filename))
    with open(filename, 'w') as f:
        f.write("{}\n".format(title))
        for namelist in ['cntrl', 'ewald']:
            if dictionary[namelist] is not None:
                f.write(" &{}\n".format(namelist))
                for key, val in dictionary[namelist].items():
                    if val is not None:
                        f.write("  {:15s} {:s},\n".format(
                            key + ' =', str(val)))
                f.write(" /\n")
        if dictionary['cntrl']['nmropt'] == 1:
            if dictionary['wt'] is not None:
                for line in dictionary['wt']:
                    f.write(" " + line + "\n")
            f.write(" &wt type = 'END', /\n")
            if dictionary['DISANG'] is not None:
                f.write("DISANG = {}\n".format(dictionary['DISANG']))
                f.write("LISTOUT = POUT\n\n")
        if dictionary['GROUP'] is not None:
            f.write("{:s}".format(dictionary['GROUP']))



class AMBER_GB_simulation(phase=None, window=None):
    """
    """

    def __init__(self):

        # Setup simulation directory and files
        self.path = './'
        self.amber_executable = 'pmemd'
        self.topology = self.path + 'prmtop'
        self.restraints = self.path + 'restraints.in'

        self.min = {}
        self.min['prefix'] = 'minimize'
        self.min['input'] = self.path + self.min['prefix'] + '.in'
        self.min['output'] = self.path + self.min['prefix'] + '.out'
        self.min['restart'] = self.path + self.min['prefix'] + '.rst'
        self.min['info'] = self.path + self.min['prefix'] + '.inf'
        self.min['coords'] = self.path + self.min['prefix'] + '.inpcrd'

        self.min['cntrl'] = OrderedDict()
        # Minimize
        self.min['cntrl']['imin'] = 1
        # Coordinates will be read
        self.min['cntrl']['ntx'] = 1
        # Not a restart
        self.min['cntrl']['irest'] = 0
        # Maximum number of cycles
        self.min['cntrl']['maxcyc'] = 5000
        # Number of steps until switching from SD to CG
        self.min['cntrl']['ncyc'] = 1000
        # Number of steps to output human-readable energy
        self.min['cntrl']['ntpr'] = 100
        # Output format
        self.min['cntrl']['ntxo'] = 1
        # Whether all interactions are calculated
        self.min['cntrl']['ntf'] = 1
        # Whether SHAKE is performed
        self.min['cntrl']['ntc'] = 1
        # Nonbonded cutoff
        self.min['cntrl']['cut'] = 999.0
        # Flag for generalized Born
        self.min['cntrl']['igb'] = 1
        self.min['cntrl']['ntr'] = None
        self.min['cntrl']['restraint_wt'] = None
        self.min['cntrl']['restraintmask'] = None
        self.min['cntrl']['nmropt'] = 1
        self.min['cntrl']['pencut'] = -1
        self.min['ewald'] = None
        self.min['wt'] = None    # or []
        self.min['DISANG'] = self.restraints
        self.min['GROUP'] = None    # or []

        self.md = {}
        self.md['prefix'] = 'md'
        self.md['input'] = self.path + self.md['prefix'] + '.in'
        self.md['output'] = self.path + self.md['prefix'] + '.out'
        self.md['restart'] = self.path + self.md['prefix'] + '.rst'
        self.md['info'] = self.path + self.md['prefix'] + '.inf'
        self.md['coords'] = self.path + self.md['prefix'] + '.inpcrd'
        self.md['ref_coords'] = None

        self.md['cntrl'] = OrderedDict()
        self.md['cntrl']['imin'] = 0
        self.md['cntrl']['ntx'] = 1
        self.md['cntrl']['irest'] = 0
        self.md['cntrl']['dt'] = 0.002
        self.md['cntrl']['nstlim'] = 5000
        self.md['cntrl']['ntpr'] = 500
        self.md['cntrl']['ntwe'] = 0
        self.md['cntrl']['ntwr'] = 5000
        self.md['cntrl']['ntxo'] = 1
        self.md['cntrl']['ntwx'] = 500
        self.md['cntrl']['ioutfm'] = 1
        self.md['cntrl']['ntf'] = 2
        self.md['cntrl']['ntc'] = 2
        self.md['cntrl']['igb'] = 1
        self.md['cntrl']['cut'] = 999.0
        self.md['cntrl']['ntt'] = 3
        self.md['cntrl']['gamma_ln'] = 1.0
        self.md['cntrl']['ig'] = -1
        self.md['cntrl']['temp0'] = 298.15
        self.md['cntrl']['ntr'] = None
        self.md['cntrl']['restraint_wt'] = None
        self.md['cntrl']['restraintmask'] = None
        self.md['cntrl']['nmropt'] = 1
        self.md['cntrl']['pencut'] = -1

        self.md['ewald'] = None
        self.md['wt'] = None
        self.md['DISANG'] = self.restraints
        self.md['GROUP'] = None


    def minimize_with_amber(self, soft=False):
        """
        Minimize the system.

        If soft=True, slowly turn on non-bonded interactions during minimization
        so that the restraints get enforced first.

        """

        if soft:
            self.min['wt'] = [
                "&wt type = 'NB', istep1=0, istep2=2000, value1 = 0.0, value2=0.0, IINC=50, /",
                "&wt type = 'NB', istep1=2000, istep2=4000, value1 = 0.0, value2=1.0, IINC=50, /"]
        _amber_write_input_file(self.min['input'], self.min, title='GB Minimization.')

        log.info('Running GB Minimization at {}'.format(self.path))
        sp.call("{6} -O -p {0} -c {1} -ref {1} -i {2}\
                 -o {3} -r {4} -inf {5}"
                .format(self.toplogy, self.min['coords'], self.min['input'],
                        self.min['output'], self.min['restart'], self.min['info'],
                        self.amber_executable), cwd=self.path, shell=True)
        log.debug('TODO: Does the above command need `shell=True` ...?')
        log.debug('TODO: Catch errors here...')

    def run_md_with_amber(self):
        """
        Run MD with AMBER.
        """

        _amber_write_input_file(self.md['input'], self.md, title='MD.')
        log.info('Running AMBER MD at {}'.format(self.path))

        if self.md['ref_coords']:
            log.warning('Check the specificiation of reference coordinates is correct.')
            execution_string = "{6} -O -p {0} -c {1} -ref {7} -i {2}\
                     -o {3} -r {4} -inf {5}".format(self.toplogy, self.md['coords'], self.md['input'],
                                                    self.md['output'], self.md['restart'],
                                                    self.md['restart'], self.amber_executable,
                                                    self.md['ref_coords'])

            log.debug(execution_string)
            sp.call(execution_string, cwd=self.path)

        else:
            execution_string = "{6} -O -p {0} -c {1} -i {2}\
                     -o {3} -r {4} -inf {5}".format(self.toplogy, self.md['coords'], self.md['input'],
                                                    self.md['output'], self.md['restart'],
                                                    self.md['restart'], self.amber_executable)

            log.debug(execution_string)
            sp.call(execution_string, cwd=self.path)
        log.info('Minimization completed...')


