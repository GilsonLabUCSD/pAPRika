import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd

class pme_simulation():
    """
    All PME/PBC simulation methods.  Minimization, equilibration, production.
    """

    def __init__(self):

        self.minimization_settings = {
            'imin':             1,
            'ntx':              1,
            'irest':            0,
            'maxcyc':           20000,
            'ncyc':             5000,
            'ntpr':             50,
            'ntxo':             1,
            'ntf':              1,
            'ntc':              1,
            'ntb':              1,
            'cut':              9.0,
            'ntr':              0,
            'restraint_wt':     None,
            'restraintmask':    None,
            'nmropt':           1,
            'pencut':           -1,
            '&namelist_string': None,  # For example: " &ewald\n  eedmeth=5, /\n"
            'wt_string':        None,  # For example: '&wt type='TEMP0', istep1=0, istep2=500, value1=10.0, value2=300.0, /'
            'DISANG':           'restraints.in',
            'LISTOUT':          'POUT'
            'GROUP_string':     None,  # "GROUP\nTITLE\nRES 1 15\nEND"
        }


class gb_simulation():
    """
    All GB simulation methods.
    """

    def __init__(self):

        self.minimization_settings = {
            ### See AMBER manual for definitions
            'imin':             1,
            'ntx':              1,
            'irest':            0,
            'maxcyc':           20000,
            'ncyc':             5000,
            'ntpr':             50,
            'ntxo':             1,
            'ntf':              1,
            'ntc':              1,
            'ntb':              0,
            'cut':              999.0,
            'igb':              1,
            'ntr':              0,
            'restraint_wt':     None,
            'restraintmask':    None,
            'nmropt':           1,
            'pencut':           -1,
            '&namelist_string': None,  # For example: " &ewald\n  eedmeth=5, /\n"
            'wt_string':        None,  # For example: '&wt type='TEMP0', istep1=0, istep2=500, value1=10.0, value2=300.0, /'
            'DISANG':           'restraints.in',
            'LISTOUT':          'POUT'
            'GROUP_string':     None,  # "Group Title\5.0\nRES 1 15\nEND\nEND\n"
        }

        self.equilibration = {
            ### See AMBER manual for definitions
            'imin':             0,
            'ntx':              1,
            'irest':            0,
            'nstlim':           500000,
            '
            'ntpr':             500,
            'ntwr':             500000,
            'ntxo':             1,
            'ntwx':             500,
            'ioutfm':           1,
            'ntf':              2,
            'ntc':              2,
            'ntb':              0,
            'cut':              999.0,
            'igb':              1,
            'ntt':              3,
            'gamma_ln':         1.0,
            'temp0':            298.15,
            'ntr':              0,
            'restraint_wt':     None,
            'restraintmask':    None,
            'nmropt':           1,
            'pencut':           -1,
            '&namelist_string': None,  # For example: " &ewald\n  eedmeth=5, /\n"
            'wt_string':        None,  # For example: '&wt type='TEMP0', istep1=0, istep2=500, value1=10.0, value2=300.0, /'
            'DISANG':           'restraints.in',
            'LISTOUT':          'POUT'
            'GROUP_string':     None,  # "Group Title\5.0\nRES 1 15\nEND\nEND\n"
        }
             

        self.print_time = 1.0 #ps.  Controls 

