import numpy as np
import subprocess as sp
import logging as log
import parmed as pmd
from paprika import align

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')


# Page 379 of Amber 7 manual has mask specifications
class DAT_restraint(object):
    """
    Distance or angle or torsion restraints.
    """

    def __init__(self):
        self.structure_file = None
        self.mask1 = None
        self.mask2 = None
        self.mask3 = None
        self.mask4 = None
        self.index1 = None
        self.index2 = None
        self.index3 = None
        self.index4 = None

        self.attach =  {'target':            None,
                        'force_initial':     None,
                        'force_final':       None,
                        'force_increment':   None
                       }
        self.pull =    {'target_initial':    None,
                        'target_final':      None,
                        'target_increment':  None,
                        'force_initial':     None,
                        'force_final':       None,
                        'force_increment':   None
                       }
        self.release = {'target':            None,
                        'force_initial':     None,
                        'force_final':       None,
                        'force_increment':   None
                       }

        self.attach_windows  = []
        self.pull_windows    = []
        self.release_windows = []

    def index_from_mask(self, mask):
        """
        Return the atom indicies, given a mask.
        """
        structure = align.return_structure(self.structure_file)
        return [i.idx for i in structure[mask]]

    def initialize(self):
        """
        If the user hasn't already declared the windows list, we will
        construct one from the initial, final, and increment values.
        """
        if not self.attach_windows and self.attach['force_initial']:
            self.attach_windows = np.arange(self.attach['force_initial'],
                                            self.attach['force_final'],
                                            self.attach['force_increment'])
        if not self.pull_windows and self.pull['force_final']:
            self.pull_windows = np.arange(self.pull['force_initial'],
                                          self.pull['force_final'],
                                          self.pull['force_increment'])
        if not self.release_windows and self.release['force_final']:
            self.release_windows = np.arange(self.release['force_initial'],
                                            self.release['force_final'],
                                            self.release['force_increment'])
        self.phase = {'attach':  self.attach_windows,
                      'pull':    self.pull_windows,
                      'release': self.release_windows
                     }

        self.index1 = self.index_from_mask(self.mask1)
        self.index2 = self.index_from_mask(self.mask2)
        if self.mask3:
            self.index3 = self.index_from_mask(self.mask3)
        else:
            self.index3 = None
        if self.mask4:
            self.index4 = self.index_from_mask(self.mask4)
        else:
            self.index4 = None


def return_restraint_line(restraint, phase, window, group=False):
    """
    Write the restraints to a file for each window.
    For example:
    &rst iat= 3,109, r1= 0.0000, r2= 6.9665, r3= 6.9665, r4= 999.0000,
         rk2= 5.0000000, rk3= 5.0000000, &end
    """

    string = '&rst ' + \
             '\tiat = {}, {}, {}, {},'.format(restraint.index1, restraint.index2, \
                                              restraint.index3, restraint.index4) + \
             '\tr1  = {0:4.4f}'.format(0) + \
             '\tr2  = {0:4.4f}'.format(restraint.phase[phase][window]) + \
             '\tr3  = {0:4.4f}'.format(restraint.phase[phase][window]) + \
             '\tr4  = {0:4.4f}'.format(999) + \
             '\trk2 = {0:4.4f}'.format(restraint.phase[phase] + \
             '\trk3 = {0:4.4f}'.format(restraint.phase[phase] + \
             '\t&end'

    if group:
        '''
        # Group Distance Restraint
        &rst iat= -1,-1, igr1=3,4,7,8,21,22,25,26,39,40,43,44,57,58,61,62,75,76,79,80,93,94,97,98, igr2=109,113,115,119, r1=     0.0000, r2=     5.9665, r3=     5.9665, r4=   999.0000, rk2=   5.0000000, rk3=   5.0000000, &end
        '''
        raise Exception('Not yet!')
    print(string)


# Initialize a distance restraint that acts on :BUT and :CB6...
this = DAT_restraint()
this.structure_file = '../test/cb6-but/cb6-but.pdb'
this.mask1 = ':BUT'
this.mask2 = ':CB6'

# User specifies the windows directly...
# this.pull_windows = [0, 1.1, 2.4, 3.6]
# User specifies initial, final, and increment, and we build the list...
this.pull = {'target_initial'   : 0,  # angstroms (by definition, I guess)
             'target_final'     : 10, # angstroms
             'target_increment' : 1,  # angstroms
             'force_initial'    : 5,  # kcal per mol
             'force_final'      : 5,  # kcal per mol
             'force_increment'  : 1   # kcal per mol (but shouldn't matter)
            }
this.initialize()
return_restraint_line(this, phase='pull', window=0)

# Generate windows for each
# for number, window in enumerate(this.pull_windows):
    # sp.call(['mkdir -p windows/p{}'.format(number)], shell=True)

# In each window, write all restraints
