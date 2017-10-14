import logging as log
import numpy as np
import os as os
import parmed as pmd
import subprocess as sp

from collections import OrderedDict


def _write_input_file(filename, input_dict, title='Input.'):
    """
    Process a dictionary to generate a pmemd input file.
    """

    log.debug('Writing {}'.format(filename))
    with open(filename, 'w') as f:
        f.write("{}\n".format(title))
        for namelist in ['cntrl','ewald']:
            if input_dict[namelist] is not None:
                f.write(" &{}\n".format(namelist))
                for key,val in input_dict[namelist].items():
                    if val is not None:
                        f.write("  {:15s} {:s},\n".format(key+' =',str(val)))
                f.write(" /\n")
        if input_dict['cntrl']['nmropt'] == 1:          ### Need to check if this exists first!!
            if input_dict['wt'] is not None:
                for line in  input_dict['wt']:
                    f.write(" "+line+"\n")
            f.write(" &wt type = 'END', /\n")
            if input_dict['DISANG'] is not None:
                f.write("DISANG = {}\n".format(input_dict['DISANG']))
                f.write("LISTOUT = POUT\n\n")
        if input_dict['GROUP'] is not None:
            f.write("{:s}".format(input_dict['GROUP']))



class pme_simulation():
    """
    All PME/PBC simulation methods.
    """

    def __init__(self):
        self.todo = 'todo'


class gb_simulation():
    """
    All GB simulation methods.
    """

    def __init__(self):

        
        self.path = '.'
        self.prmtop = 'prmtop'
        self.inpcrd = 'inpcrd'
        self.job = {}

        # See AMBER manual for definitions
        # I'm not really pleased with the format here, but I want to keep
        # the order of the keys, so I have to go line by line for the
        # OrderedDict.  If I make the dict all at once, it won't be in order.
        # This format should allow flexible editing of the input files

        # Minimization input file settings
        self.job['minimization'] = {}
        self.job['minimization']['cntrl'] =                     OrderedDict()
        self.job['minimization']['cntrl']['imin'] =             1
        self.job['minimization']['cntrl']['ntx'] =              1
        self.job['minimization']['cntrl']['irest'] =            0
        self.job['minimization']['cntrl']['maxcyc'] =           5000
        self.job['minimization']['cntrl']['ncyc'] =             1000
        self.job['minimization']['cntrl']['ntpr'] =             100
        self.job['minimization']['cntrl']['ntxo'] =             1
        self.job['minimization']['cntrl']['ntf'] =              1
        self.job['minimization']['cntrl']['ntc'] =              1
        self.job['minimization']['cntrl']['cut'] =              999.0
        self.job['minimization']['cntrl']['igb'] =              1
        self.job['minimization']['cntrl']['ntr'] =              None
        self.job['minimization']['cntrl']['restraint_wt'] =     None
        self.job['minimization']['cntrl']['restraintmask'] =    None
        self.job['minimization']['cntrl']['nmropt'] =           1
        self.job['minimization']['cntrl']['pencut'] =           -1
        self.job['minimization']['ewald'] =                     None    # or OrderedDict() (Add the other namelists?)
        self.job['minimization']['wt'] =                        None    # or []
        self.job['minimization']['DISANG'] =                    'restraints.in'
        self.job['minimization']['GROUP'] =                     None    # or []

        # Production input file settings
        self.job['md'] = {}
        self.job['md']['cntrl'] =                       OrderedDict()
        self.job['md']['cntrl']['imin'] =               0
        self.job['md']['cntrl']['ntx'] =                1
        self.job['md']['cntrl']['irest'] =              0
        self.job['md']['cntrl']['dt'] =                 0.002
        self.job['md']['cntrl']['nstlim'] =             5000
        self.job['md']['cntrl']['ntpr'] =               500
        self.job['md']['cntrl']['ntwe'] =               0
        self.job['md']['cntrl']['ntwr'] =               5000
        self.job['md']['cntrl']['ntxo'] =               1
        self.job['md']['cntrl']['ntwx'] =               500
        self.job['md']['cntrl']['ioutfm'] =             1
        self.job['md']['cntrl']['ntf'] =                2
        self.job['md']['cntrl']['ntc'] =                2
        self.job['md']['cntrl']['igb'] =                1
        self.job['md']['cntrl']['cut'] =                999.0
        self.job['md']['cntrl']['ntt'] =                3
        self.job['md']['cntrl']['gamma_ln'] =           1.0
        self.job['md']['cntrl']['ig'] =                 -1
        self.job['md']['cntrl']['temp0'] =              298.15
        self.job['md']['cntrl']['ntr'] =                None
        self.job['md']['cntrl']['restraint_wt'] =       None
        self.job['md']['cntrl']['restraintmask'] =      None
        self.job['md']['cntrl']['nmropt'] =             1
        self.job['md']['cntrl']['pencut'] =             -1
        self.job['md']['ewald'] =                       None    # or OrderedDict() (Add the other namelists?)
        self.job['md']['wt'] =                          None
        self.job['md']['DISANG'] =                      'restraints.in'
        self.job['md']['GROUP'] =                       None

        #self.print_time = 1.0 #ps.  Figure out a way to update dicts if this is reset by user. 
        # offer way for user to just provide the input file?

    def minimize(self, Soft=False):
        """
        Minimize the system.

            If Soft=True, slowly turn on non-bonded interactions during minimization
            so that the restraints get enforced first.

        """

        if Soft:
            self.job['minimization']['wt'] = [ "&wt type = 'NB', istep1=0, istep2=2000, value1 = 0.0, value2=0.0, IINC=50, /",
                                                "&wt type = 'NB', istep1=2000, istep2=4000, value1 = 0.0, value2=1.0, IINC=50, /" ]

        _write_input_file(self.path+'/gb_minimize.in', self.job['minimization'], title='GB Minimization.')

        ### ToDo: Use nice error catching method dave figured out.
        log.debug('Running GB Minimization at {}'.format(self.path))
        sp.call("pmemd -O -p {0} -c {1} -ref {1} -i gb_minimize.in\
                 -o mdout.gb_minimize -r rst.gb_minimize -inf mdinfo.gb_minimize".format(self.prmtop,self.inpcrd), cwd=self.path, shell=True)


    def runmd(self, input_file, input_crds, output_suffix, ref_crds=None):
        """
        Run MD.
        """

        _write_input_file(self.path+'/'+input_file, self.job['md'], title='MD.')

        if ref_crds is None:
            ref_crds = input_crds

        ### ToDo: Use nice error catching method dave figured out.
        ### ToDo: Is this the best way to do file naming?
        log.debug('Running MD at {}'.format(self.path))
        sp.call("pmemd -O -p {0} -c {1} -ref {2} -i {3} -o mdout.{4} -r rst.{4} -x traj.{4} -e mden.{4} -inf mdinfo.{4}\
                ".format(self.prmtop,self.inpcrd,ref_crds,input_file,output_suffix), cwd=self.path, shell=True)

