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

class MutableDict(dict):
    """ If prefix gets updated, updated all associated files """

    suffixes = {
        'input':        '.in',
        'inpcrd':       '.inpcrd',
        'ref':          '.inpcrd',
        'output':       '.out',
        'restart':      '.rst',
        'mdinfo':       '.mdinfo',
        'mdcrd':        '.nc',
    }

    def __setitem__(self,key,value):
        if key == 'prefix':
            dict.__setitem__(self, 'prefix', value)
            for filetype, suff in self.suffixes.items():
                if dict.__contains__(self, filetype):
                    dict.__setitem__(self, filetype, value+suff)
        else:
            dict.__setitem__(self,key,value)


class Simulation(object):
    """
    AMBER simulation class.
    """

    def __init__(self):

        ### Setup simulation directory and files
        self.path = '.' # Assume everything will be created/executed in this path
        self.amber_executable = 'pmemd'
        self.phase = None   
        self.window = None
        self.topology = 'prmtop'
        self.restraint_file = 'restraints.in'


        ### Minimization Settings
        self._min = MutableDict({})
        self._min['prefix'] = 'minimize'
        self._min['input'] = self._min['prefix'] + '.in'
        self._min['inpcrd'] = self._min['prefix'] + '.inpcrd'
        #self._min['ref'] = None  # We're gonna leave this unset
        self._min['output'] = self._min['prefix'] + '.out'
        self._min['restart'] = self._min['prefix'] + '.rst'
        self._min['mdinfo'] = self._min['prefix'] + '.mdinfo'

        self.min['cntrl'] = OrderedDict()
        self.min['cntrl']['imin'] = 1
        self.min['cntrl']['ntx'] = 1
        self.min['cntrl']['irest'] = 0
        self.min['cntrl']['maxcyc'] = 5000
        self.min['cntrl']['ncyc'] = 1000
        self.min['cntrl']['ntpr'] = 100
        self.min['cntrl']['ntxo'] = 1
        self.min['cntrl']['ntf'] = 1
        self.min['cntrl']['ntc'] = 1
        self.min['cntrl']['cut'] = 999.0
        self.min['cntrl']['ntb'] = 0
        self.min['cntrl']['igb'] = 1
        self.min['cntrl']['ntr'] = None
        self.min['cntrl']['restraint_wt'] = None
        self.min['cntrl']['restraintmask'] = None
        self.min['cntrl']['nmropt'] = 1
        self.min['cntrl']['pencut'] = -1

        self.min['ewald'] = None
        self.min['wt'] = None    # or []
        self.min['DISANG'] = self.restraint_file
        self.min['GROUP'] = None    # or []


        ### MD settings
        self._md = MutableDict({})
        self._md['prefix'] = 'md'
        self._md['input'] = self._md['prefix'] + '.in'
        self._md['inpcrd'] = self._md['prefix'] + '.inpcrd'
        #self._md['ref'] = None  # Leave unset
        self._md['output'] = self._md['prefix'] + '.out'
        self._md['restart'] = self._md['prefix'] + '.rst'
        self._md['mdcrd'] = self._md['prefix'] + '.nc'
        self._md['mdinfo'] = self._md['prefix'] + '.mdinfo'
        self._md['coords'] = self._md['prefix'] + '.inpcrd'

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
        self.md['DISANG'] = self.restraint_file
        self.md['GROUP'] = None

    ### Refresh the min/md dicts if prefix changes
    @property
    def min(self):
        return self._min
    @property
    def md(self):
        return self._md

    def config_ntp(self):
        """
        Configure input settings to default NTP.
        """

        self.min['cntrl']['ntb'] = 1
        self.min['cntrl']['cut'] = 8.0
    
    def minimize(self, soft=False):
        """
        Minimize the system.

        If soft=True, slowly turn on non-bonded interactions during minimization
        so that the restraints get enforced first.

        """

        # These settings hardcoded at the moment ... possibly expose for editing in the future
        if soft:
            # Set a burn in value that is 25% of the way between ncyc and maxcyc
            burn_in = int(float(self.min['cntrl']['ncyc']) + 0.20*(float(self.min['cntrl']['maxcyc']) - float(self.min['cntrl']['ncyc'])))
            # If the burn_in value is nuts, then just put it to zero
            if burn_in < 0 or burn_in >= self.min['cntrl']['maxcyc']:
                burn_in = 0
            # Set an end_soft value that is 75% of way between ncyc and maxcyc
            end_soft = int(float(self.min['cntrl']['ncyc']) + 0.60*(float(self.min['cntrl']['maxcyc']) - float(self.min['cntrl']['ncyc'])))
            self.min['wt'] = [
                "&wt type = 'NB', istep1=0, istep2={:.0f}, value1 = 0.0, value2=0.0, IINC=50, /".format(burn_in),
                "&wt type = 'NB', istep1={:.0f}, istep2={:.0f}, value1 = 0.0, value2=1.0, IINC=50, /".format(burn_in,end_soft)]

        _amber_write_input_file(self.path+'/'+self.min['input'], self.min, title='GB Minimization.')

        log.info('Running GB Minimization at {}'.format(self.path))

        if 'ref' in self.min:
            exec_list = [
                self.amber_executable, '-O', '-p', self.topology, '-ref',  self.min['ref'],
                '-c', self.min['inpcrd'], '-i', self.min['input'], '-o', self.min['output'],
                '-r', self.min['restart'], '-inf', self.min['mdinfo']
            ]
        else:
            exec_list = [
                self.amber_executable, '-O', '-p', self.topology,
                '-c', self.min['inpcrd'], '-i', self.min['input'], '-o', self.min['output'],
                '-r', self.min['restart'], '-inf', self.min['mdinfo']
            ]


        log.debug('Exec line: '+' '.join(exec_list))
        sp.call(exec_list, cwd=self.path)
        log.debug('TODO: Catch errors here...')
        log.info('MD completed...')


    def run_md(self):
        """
        Run MD with AMBER.
        """

        _amber_write_input_file(self.path+'/'+self.md['input'], self.md, title='MD.')

        log.info('Running AMBER MD at {}'.format(self.path))

        if 'ref' in self.md:
            exec_list = [
                self.amber_executable, '-O', '-p', self.topology, '-ref',  self.md['ref'],
                '-c', self.md['inpcrd'], '-i', self.md['input'], '-o', self.md['output'],
                '-r', self.md['restart'], '-x', self.md['mdcrd'], '-inf', self.md['mdinfo']
            ]
        else:
            exec_list = [
                self.amber_executable, '-O', '-p', self.topology,
                '-c', self.md['inpcrd'], '-i', self.md['input'], '-o', self.md['output'],
                '-r', self.md['restart'], '-x', self.md['mdcrd'], '-inf', self.md['mdinfo']
            ]

        log.debug('Exec Line: '+' '.join(exec_list))
        sp.call(exec_list, cwd=self.path)
        log.debug('TODO: Catch errors here...')
        log.info('MD completed...')


