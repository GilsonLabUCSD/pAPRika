import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd
from paprika import align

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(format='%(asctime)s %(message)s',
                datefmt='%Y-%m-%d %I:%M:%S %p')


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

        self.attach = {
            'target_values':    None,  # If the user wants to specify a list...
            'target_initial':   None,  # Percent of force constant
            'target_final':     None,
            'target_increment': None,
            'force_values':     None,  # If the user wants to specify a list...
            'force_initial':    None,
            'force_final':      None,
            'force_increment':  None
        }

        self.pull = {
            'target_values':    None,
            'target_initial':   None,
            'target_final':     None,
            'target_increment': None,
            'force_values':     None,
            'force_initial':    None,
            'force_final':      None,
            'force_increment':  None
        }

        self.release = {
            'target_values':    None,
            'target_initial':   None,
            'target_final':     None,
            'target_increment': None,
            'force_values':     None,
            'force_initial':    None,
            'force_final':      None,
            'force_increment':  None
        }
        log.info("1. Specify a list of values like `restraint.attach['target_values']`")
        log.info("2. Specify a range of values using `restraint.attach['target_initial]`, "
                 "`restraint.attach['target_final']` and "
                 "`restraint.attach['target_increment']` ")
        log.info("3. Specify a single value like `restraint.attach['force_final']`")

    def initialize(self):
        """
        If the user hasn't already declared the windows list, we will
        construct one from the initial, final, and increment values.
        """

        # Create the container to hold the actual target and force values that
        # are computed from the class lists. This is a roundabout way of keeping track
        # of the data we need to calculate the work done by the restraints while allowing
        # flexibility in how the user specifies the targets and force constants using the
        # __init__ method. There is probably a more elegant way of keeping track of this, but
        # this should work for now.

        self.phase = {
            'attach':  {
                'force_constants': None,
                'targets':         None
            },
            'pull':    {
                'force_constants': None,
                'targets':         None
            },
            'release': {
                'force_constants': None,
                'targets':         None
            }
        }
        # ------------------------------------ ATTACH ------------------------------------ #
        # Either the user specifies a direct list of target values during attachment...
        if self.attach['target_values']:
            log.debug('User directly specified list of target values for attachment...')
            self.phase['attach']['targets'] = self.attach['target_values']
        # ... or we build a list of targets using a range supplied.
        elif (self.attach['target_initial'] is not None) and \
             (self.attach['target_increment'] is not None) and \
             (self.attach['target_final'] is not None):
            log.debug('Building attachment targets from fractions...')
            self.phase['attach']['targets'] = np.arange(self.attach['target_initial'],
                                                        self.attach['target_final'] +
                                                        self.attach['target_increment'],
                                                        self.attach['target_increment'])
        # After we have the list of targets, we calculate the force constant in each window during
        # attachment. Either the user specifies the list of force constants
        # directly...
        if self.attach['force_values']:
            log.debug('User directly specified list of force constants for attachment...')
            self.phase['attach']['force_constants'] = self.attach['force_values']
        # ...or we find the force constant in each window by multiplying the target fraction
        # by the final force constant.
        elif self.attach['force_final']:
            log.debug('Building attachment force constants from fractions...')
            self.phase['attach']['force_constants'] = [self.attach['force_final'] * percent
                                                       for percent in self.phase['attach']['targets']]

        # ------------------------------------ RELEASE ------------------------------------ #
        # Since release is the inverse of attach, the process is identical.
        # Either the user specifies a direct list of target values during release...
        if self.release['target_values']:
            log.debug('User directly specified list of target values for release...')
            self.phase['release']['targets'] = self.release['target_values']
        # ... or we build a list of targets using a range supplied.
        elif (self.release['target_initial'] is not None) and \
             (self.release['target_increment'] is not None) and \
             (self.release['target_final'] is not None):
            log.debug('Building release targets from fractions...')
            self.phase['release']['targets'] = np.arange(self.release['target_initial'],
                                                         self.release['target_final'] +
                                                         self.release['target_increment'],
                                                         self.release['target_increment'])
        # After we have the list of targets, we calculate the force constant in each window during
        # release. Either the user specifies the list of force constants
        # directly...
        if self.release['force_values']:
            log.debug('User directly specified list of force constants for release...')
            self.phase['release']['force_constants'] = self.release['force_values']
        # ...or we find the force constant in each window by multiplying the target fraction
        # by the final force constant.
        elif self.release['force_final']:
            log.debug('Building release force constants from fractions...')
            self.phase['release']['force_constants'] = [self.release['force_final'] * percent
                                                        for percent in self.phase['release']['targets']]

        # ------------------------------------ PULL ------------------------------------ #
        # Either the user specifies a direct list of target values for puling...
        if self.pull['target_values']:
            log.debug('User directly specified a list of target values for pulling...')
            self.phase['pull']['targets'] = self.pull['target_values']
        # ...or we build a list of targets using a range supplied.
        elif (self.pull['target_initial'] is not None) and \
             (self.pull['target_increment'] is not None) and \
             (self.pull['target_final'] is not None):
            log.debug('Building pull targets from fractions...')
            self.phase['pull']['targets'] = np.arange(self.pull['target_initial'],
                                                      self.pull['target_final'] +
                                                      self.pull['target_increment'],
                                                      self.pull['target_increment'])

        # After we have the list of targets, we calculate the force constant in each window during
        # release. Either the user specifies the list of force constants
        # directly...
        if self.pull['force_values']:
            log.debug('User directly specified list of force constants for pull...')
            self.phase['pull']['force_constants'] = self.pull['force_values']
        # ...or we find the force constant in each window by multiplying the target fraction
        # by the final force constant.
        elif self.pull['force_final']:
            log.debug('Building pull force constants from fractions...')
            self.phase['pull']['force_constants'] = [self.pull['force_final'] * percent
                                                     for percent in self.phase['pull']['targets']]

        # ---------------------------------- ATOM MASKS ---------------------------------- #
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
             '\tiat = {}, {}, {}, {},'.format(iat1,
                                              iat2,
                                              iat3,
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
        '\trk2 = {0:4.4f},'.format(restraint.phase[phase]['force_constants'][window]) + \
        '\trk3 = {0:4.4f},'.format(restraint.phase[phase]['force_constants'][window]) + \
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
            line = return_restraint_line(
                    restraint, phase='attach', window=window)
            # Using append mode is crucial for multiple restraints.
            directory = './windows/a{0:03d}'.format(window)
            f = open(directory + '/' + filename, 'a')
            f.write(line + '\n')
            f.close()

    for restraint in restraints:
        for window, _ in enumerate(restraint.pull_forces):
            line = return_restraint_line(
                    restraint, phase='pull', window=window)
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
