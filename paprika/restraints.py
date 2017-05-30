import numpy as np
import subprocess as sp
import logging as log
import os as os
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

        self.attach =  {'target_initial':    None, # Percent of force constant
                        'target_final':      None,
                        'target_increment':  None,
                        'force_constant':    None,
                        'targets':           None, # List to hold actual target values
                        'forces':            None  # List to hold actual force values
                       }
        self.pull =    {'target_initial':    None, # Distance for restraint
                        'target_final':      None,
                        'target_increment':  None,
                        'force_constant' :   None,
                        'targets':           None,  # List to hold actual target values
                        'forces':            None   # List to hold actual force values
                       }
        self.release = {'target':            None,
                        'force_initial':     None,
                        'force_final':       None,
                        'force_increment':   None
                       }

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
        if not self.attach['targets']:
            log.debug('Building attach targets from fractions...')
            self.attach['targets'] = np.arange(self.attach['target_initial'],
                                               self.attach['target_final'] +
                                               self.attach['target_increment'],
                                               self.attach['target_increment'])
        if not self.attach['forces']:
            log.debug('Building attach force targets from fractions...')
            self.attach['forces'] = [self.attach['force_constant'] * i
                                     for i in self.attach['targets']]

        if not self.pull['targets']:
            log.debug('Building pull targets from fractions...')
            self.pull['targets'] = np.arange(self.pull['target_initial'],
                                             self.pull['target_final'] +
                                             self.pull['target_increment'],
                                             self.pull['target_increment'])
        if not self.pull['forces']:
            log.debug('Building pull force targets from fractions...')
            self.pull['forces'] = [self.pull['force_constant'] * i
                                     for i in self.pull['targets']]



        if not self.attach_targets and self.attach['force']:
            log.debug('Building attach force targets from a range...')
            # Make a list of targets as long as the attachment windows.
            # I'm not sure the following line is working...
            log.warn('Investigate!')
            # self.attach_targets = [self.attach['target']] * len(self.attach_forces)

        if not self.pull_forces and self.pull['force_final']:
            # In the pulling phase, the target distance also changes.
            self.pull_targets = np.arange(self.pull['target_initial'],
                                          self.pull['target_final'] +
                                          self.pull['target_increment'],
                                          self.pull['target_increment'])

            if self.pull['force_initial'] != self.pull['force_final']:
                self.pull_forces = np.arange(self.pull['force_initial'],
                                             self.pull['force_final'] +
                                             self.pull['force_increment'],
                                             self.pull['force_increment'])
            else:
                self.pull_forces = [self.pull['force_initial']] * len(self.pull_targets)
        if not self.release_forces and self.release['force_final']:
            self.release_forces = np.arange(self.release['force_initial'],
                                            self.release['force_final'] +
                                            self.release['force_increment'],
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
        # Right?
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
            '\tr1  = {0:4.4f},'.format(0) + \
            '\tr2  = {0:4.4f},'.format(restraint.phase[phase]['targets'][window]) + \
            '\tr3  = {0:4.4f},'.format(restraint.phase[phase]['targets'][window]) + \
            '\tr4  = {0:4.4f},'.format(999) + \
            '\trk2 = {0:4.4f},'.format(restraint.phase[phase]['forces'][window]) + \
            '\trk3 = {0:4.4f},'.format(restraint.phase[phase]['forces'][window]) + \
            '\t&end'
    return string

def make_directories(restraint):
    """
    Make a series of directories to hold the simulation setup files
    and the data. This function should probably end up in a separate
    file eventually.
    """
    # Here we could check if the directories already exist and prompt
    # the user or quit or do something else.
    # If exist_ok is False (the default), an OSError is raised if the target directory already
    # exists.

    log.debug('We ought to make sure somewhere that all restraints have'
              ' the same number of windows. Here I am creating directories based '
              ' on the number of windows in a single restraint.')
    for window, _ in enumerate(restraint.attach_forces):
        os.makedirs('./windows/a{0:03d}'.format(window), exist_ok=True)
    for window, _ in enumerate(restraint.pull_forces):
        os.makedirs('./windows/p{0:03d}'.format(window), exist_ok=True)

def write_restraints_file(restraints, filename='restraints.in'):
    """
    Take all the restraints and write them to a file in each window.
    """
    for restraint in restraints:
        for window, _ in enumerate(restraint.attach_forces):
            line = return_restraint_line(restraint, phase='attach', window=window)
            # Using append mode is crucial for multiple restraints.
            directory = './windows/a{0:03d}'.format(window)
            f = open(directory + '/' + filename, 'a')
            f.write(line + '\n')
            f.close()

    for restraint in restraints:
        for window, _ in enumerate(restraint.pull_forces):
            line = return_restraint_line(restraint, phase='pull', window=window)
            # Using append mode is crucial for multiple restraints.
            directory = './windows/p{0:03d}'.format(window)
            f = open(directory + '/' + filename, 'a')
            f.write(line + '\n')
            f.close()

def clean_restraints_file(restraints, filename='restraints.in'):
    """
    Delete the restraints files for repeated testing.
    """
    for restraint in restraints:
        for window, _ in enumerate(restraint.attach_forces):

            directory = './windows/a{0:03d}'.format(window)
            os.remove(directory + '/' + filename)

    for restraint in restraints:
        for window, _ in enumerate(restraint.pull_forces):

            directory = './windows/p{0:03d}'.format(window)
            os.remove(directory + '/' + filename)
