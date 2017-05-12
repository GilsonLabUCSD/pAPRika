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

        self.attach_forces, self.attach_targets  = [], []
        self.pull_forces, self.pull_targets    = [], []
        self.release_forces, self.release_targets = [], []

    def index_from_mask(self, mask):
        """
        Return the atom indicies, given a mask.
        """
        structure = align.return_structure(self.structure_file)
        # Use a `view` of the whole object to avoid reindexing each selection.
        masked = structure.view[mask]
        log.debug('There are {} atoms in the mask...'.format(len(masked.atoms)))
        indices = [i.idx for i in masked]
        # log.debug('Found indices {} for mask {}...'.format(indices, mask))
        return indices

    def initialize(self):
        """
        If the user hasn't already declared the windows list, we will
        construct one from the initial, final, and increment values.
        """
        if not self.attach_forces and self.attach['force_final']:
            self.attach_forces = np.arange(self.attach['force_initial'],
                                           self.attach['force_final'],
                                           self.attach['force_increment'])
            # Make a list of targets as long as the attachment windows.
            self.attach_targets = [self.attach['target']] * len(self.attach_forces)
        if not self.pull_forces and self.pull['force_final']:
            # In the pulling phase, the target distance also changes.
            self.pull_targets = np.arange(self.pull['target_initial'],
                                          self.pull['target_final'],
                                          self.pull['target_increment'])

            if self.pull['force_initial'] != self.pull['force_final']:
                self.pull_forces = np.arange(self.pull['force_initial'],
                                             self.pull['force_final'],
                                             self.pull['force_increment'])
            else:
                self.pull_forces = [self.pull['force_initial']] * len(self.pull_targets)
        if not self.release_forces and self.release['force_final']:
            self.release_forces = np.arange(self.release['force_initial'],
                                            self.release['force_final'],
                                            self.release['force_increment'])
            # Make a list of targets as long as the release windows.
            self.release_targets = [self.release['target']] * len(self.release_forces)

        self.phase = {'attach':  {'forces' : self.attach_forces,
                                  'targets' : self.attach_targets
                                 },
                      'pull':    {'forces' : self.pull_forces,
                                  'targets' : self.pull_targets
                                 },
                      'release': {'forces' : self.release_forces,
                                  'targets' : self.release_targets
                                 }
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
    Or:
    &rst iat= -1,-1, igr1=3,4,7,8,21,22,25,26,39,40,43,44,57,58,61,62,75,76,79,
    80,93,94,97,98, igr2=109,113,115,119,
    r1=     0.0000, r2=     5.9665, r3=     5.9665, r4=   999.0000,
    rk2=   5.0000000, rk3=   5.0000000, &end

    """

    # This is not very elegant, but it seems to do the trick!
    if not restraint.index1:
        iat1 = ' '
        raise Exception('There must be at least two atoms in a restraint.')
    elif len(restraint.index1) == 1:
        group1 = False
        iat1 = '{} '.format(restraint.index1[0])
    else:
        group1 = True
        iat1 = '-1 '
        igr1 = ''
        for index in restraint.index1:
            igr1 += '{}, '.format(index)

    if not restraint.index2:
        iat2 = ' '
        raise Exception('There must be at least two atoms in a restraint.')
    elif len(restraint.index2) == 1:
        group2 = False
        iat2 = '{} '.format(restraint.index2[0])
    else:
        group2 = True
        iat2 = '-1 '
        igr2 = ''
        for index in restraint.index2:
            igr2 += '{}, '.format(index)

    if not restraint.index3:
        iat3 = ' '
        group3 = False
    elif len(restraint.index3) == 1:
        group3 = False
        iat3 = '{} '.format(restraint.index3[0])
    else:
        group3 = True
        iat3 = '-1 '
        igr3 = ''
        for index in restraint.index3:
            igr3 += '{}, '.format(index)

    if not restraint.index4:
        iat4 = ' '
        group4 = False
    elif len(restraint.index4) == 1:
        group4 = False
        iat4 = '{} '.format(restraint.index4[0])
    else:
        group4 = True
        iat4 = '-1 '
        igr4 = ''
        for index in restraint.index4:
            igr4 += '{}, '.format(index)

    string = '&rst ' + \
            '\tiat = {}, {}, {}, {},'.format(iat1, \
                                             iat2, \
                                             iat3, \
                                             iat4)
    if group1:
        string += 'igr1 = {}'.format(igr1)
    if group2:
        string += 'igr2 = {}'.format(igr2)
    if group3:
        string += 'igr3 = {}'.format(igr3)
    if group4:
        string += 'igr4 = {}'.format(igr4)
    string += \
            '\tr1  = {0:4.4f}'.format(0) + \
            '\tr2  = {0:4.4f}'.format(restraint.phase[phase]['targets'][window]) + \
            '\tr3  = {0:4.4f}'.format(restraint.phase[phase]['targets'][window]) + \
            '\tr4  = {0:4.4f}'.format(999) + \
            '\trk2 = {0:4.4f}'.format(restraint.phase[phase]['forces'][window]) + \
            '\trk3 = {0:4.4f}'.format(restraint.phase[phase]['forces'][window]) + \
            '\t&end'
    return string


# Initialize a distance restraint that acts on :BUT and :CB6...
this = DAT_restraint()
this.structure_file = '../test/cb6-but/cb6-but.pdb'
this.mask1 = ':BUT@C3'
this.mask2 = ':CB6@C10'

this.pull = {'target_initial'   : 0,  # angstroms (by definition, I guess)
             'target_final'     : 10, # angstroms
             'target_increment' : 1,  # angstroms
             'force_initial'    : 5,  # kcal per mol
             'force_final'      : 5,  # kcal per mol
             'force_increment'  : 1   # kcal per mol (but shouldn't matter)
            }
this.initialize()
line = return_restraint_line(this, phase='pull', window=5)
print(line)

# Initialize a distance restraint that acts on :BUT and :CB6...
this = DAT_restraint()
this.structure_file = '../test/cb6-but/cb6-but.pdb'
this.mask1 = ':BUT'
this.mask2 = ':CB6'

this.pull = {'target_initial'   : 4,  # angstroms (by definition, I guess)
             'target_final'     : 12, # angstroms
             'target_increment' : 1,  # angstroms
             'force_initial'    : 7,  # kcal per mol
             'force_final'      : 23,  # kcal per mol
             'force_increment'  : 1   # kcal per mol (but shouldn't matter)
            }
this.initialize()
line = return_restraint_line(this, phase='pull', window=5)
print(line)