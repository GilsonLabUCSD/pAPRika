import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd
import weakref as weakref
from collections import defaultdict
from paprika import align

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(format='%(asctime)s %(message)s',
                datefmt='%Y-%m-%d %I:%M:%S %p')

# https://stackoverflow.com/questions/328851/printing-all-instances-of-a-class
class KeepRefs(object):
    __refs__ = defaultdict(list)

    def __init__(self):
        self.__refs__[self.__class__].append(weakref.ref(self))

    @classmethod
    def get_instances(cls):
        for inst_ref in cls.__refs__[cls]:
            inst = inst_ref()
            if inst is not None:
                yield inst


class DAT_restraint(KeepRefs):
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

        # For attach and release, in most cases the target distance, angle, or torsion value does not change
        # throughout the windows, but we setup the restraints with flexibility.
        self.attach = {
            'target':             None, # The target value for the restraint (mandatory)
            'fc_initial':         None, # The initial force constant (optional)
            'fc_final':           None, # The final force constant (optional)
            'num_windows':        None, # The number of windows (optional)
            'fc_increment':       None, # The force constant increment (optional)
            'fraction_increment': None, # The percentage of the force constant increment (optional)
            'fraction_list':      None, # The list of force constant percentages (optional)
            'fc_list':            None  # The list of force constants (will be created if not given)
        }
        # For the pull phase, the target distance changes and usually the force constant does not. But we will be
        # general and allow the force constant to change as well.
        self.pull = {
            'fc':                 None, # The force constant for the restraint (mandatory)
            'target_initial':     None, # The initial target value (optional)
            'target_final':       None, # The final target value (optional)
            'num_windows':        None, # The number of windows (optional)
            'target_increment':   None, # The target value increment (optional)
            'fraction_increment': None, # The percentage of the target value increment (optional)
            'fraction_list':      None, # The list of target value percentages (optional)
            'target_list':        None  # The list of target values (will be created if not given)
        }

        self.release = {
            'target':             None, # The target value for the restraint (mandatory)
            'fc_initial':         None, # The initial force constant (optional)
            'fc_final':           None, # The final force constant (optional)
            'num_windows':        None, # The number of windows (optional)
            'fc_increment':       None, # The force constant increment (optional)
            'fraction_increment': None, # The percentage of the force constant increment (optional)
            'fraction_list':      None, # The list of force constant percentages (optional)
            'fc_list':            None  # The list of force constants (will be created if not given)
        }
        super(DAT_restraint, self).__init__()

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
        log.debug('Calculating attach targets and force constants...')

        if self.attach['num_windows'] is not None  and  self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['attach']['force_constants'] = np.linspace(self.attach['fc_initial'],
                                                                      self.attach['fc_final'],
                                                                      self.attach['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['attach']['force_constants'] = np.linspace(0.0, self.attach['fc_final'],
                                                                      self.attach['num_windows'])
            self.phase['attach']['targets'] = [self.attach['target']] * self.attach['num_windows']

        elif self.attach['fc_increment'] is not None  and  self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['attach']['force_constants'] = np.arange(self.attach['fc_initial'],
                                                                    self.attach['fc_final'] +
                                                                    self.attach['fc_increment'],
                                                                    self.attach['fc_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['attach']['force_constants'] = np.arange(0.0,
                                                                    self.attach['fc_final'] +
                                                                    self.attach['fc_increment'],
                                                                    self.attach['fc_increment'])
            self.phase['attach']['targets'] = [self.attach['target']] * len(self.phase['attach']['force_constants'])

        elif self.attach['fraction_list'] is not None  and  self.attach['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['attach']['force_constants'] = [fraction * self.attach['fc_final'] for
                                                       fraction in self.attach['fraction_list']]
            self.phase['attach']['targets'] = [self.attach['target']] * len(self.phase['attach']['force_constants'])
            
        elif self.attach['fraction_increment'] is not None  and  self.attach['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(0, 1.0 + self.attach['fraction_increment'], self.attach['fraction_increment'])
            self.phase['attach']['force_constants'] = [fraction * self.attach['fc_final'] for fraction in fractions]
            self.phase['attach']['targets'] = [self.attach['target']] * len(self.phase['attach']['force_constants'])

        elif self.attach['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['attach']['force_constants'] = self.attach['fc_list']
            self.phase['attach']['targets'] = [self.attach['target']] * len(self.phase['attach']['force_constants'])

        elif all(v is None for k, v in self.attach.items()):
            log.debug('No restraint info set for this phase! Skipping ...')

        else: 
            log.error('ERROR: Restraint input did not match one of the supported methods ...')
            log.error('Input:')
            for k, v in self.attach.items():
                #if v is not None: <-- NMH: Just print 'em all right?
                log.error('{} = {}'.format(k, v))
            ### PROBABLY SHOULD RAISE AN EXCEPTION HERE

        # ------------------------------------ PULL ------------------------------------ #
        log.debug('Calculating pull targets and force constants...')

        if self.pull['num_windows'] is not None  and  self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['pull']['targets'] = np.linspace(self.pull['target_initial'],
                                                            self.pull['target_final'],
                                                            self.pull['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['pull']['targets'] = np.linspace(0.0, self.pull['target_final'],
                                                            self.pull['num_windows'])
            self.phase['pull']['force_constants'] = [self.pull['fc']] * self.pull['num_windows']

        elif self.pull['target_increment'] is not None  and  self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['pull']['targets'] = np.arange(self.pull['target_initial'],
                                                          self.pull['target_final'] +
                                                          self.pull['target_increment'],
                                                          self.pull['target_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['pull']['targets'] = np.arange(0.0, self.pull['target_final'] +
                                                          self.pull['target_increment'],
                                                          self.pull['target_increment'])
            self.phase['pull']['force_constants'] = [self.pull['fc']] * len(self.phase['pull']['targets'])

        elif self.pull['fraction_list'] is not None  and  self.pull['target_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['pull']['targets'] = [fraction * self.pull['target_final'] for
                                             fraction in self.pull['fraction_list']]
            self.phase['pull']['force_constants'] = [self.pull['fc']] * len(self.phase['pull']['targets'])

        elif self.pull['fraction_increment'] is not None  and  self.pull['target_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(0, 1.0 + self.pull['fraction_increment'], self.pull['fraction_increment'])
            self.phase['pull']['targets'] = [fraction * self.pull['target_final'] for fraction in fractions]
            self.phase['pull']['force_constants'] = [self.pull['fc']] * len(self.phase['pull']['targets'])

        elif self.pull['target_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['pull']['targets'] = self.pull['fc_list']
            self.phase['pull']['force_constants'] = [self.pull['fc']] * len(self.phase['pull']['targets'])

        elif all(v is None for k, v in self.pull.items()):
            log.debug('No restraint info set for this phase! Skipping ...')

        else:
            log.error('ERROR: Restraint input did not match one of the supported methods ...')
            log.error('Input:')
            for k, v in self.pull.items():
                #if v is not None: <-- NMH: Just print 'em all right?
                log.error('{} = {}'.format(k, v))
            ### PROBABLY SHOULD RAISE AN EXCEPTION HERE

        # ------------------------------------ RELEASE ------------------------------------ #
        log.debug('Calculating release targets and force constants...')

        if self.release['num_windows'] is not None  and  self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['release']['force_constants'] = np.linspace(self.release['fc_initial'],
                                                                       self.release['fc_final'],
                                                                       self.release['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['release']['force_constants'] = np.linspace(0.0, self.release['fc_final'],
                                                                       self.release['num_windows'])
            self.phase['release']['targets'] = [self.release['target']] * self.release['num_windows']

        elif self.release['fc_increment'] is not None  and  self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['release']['force_constants'] = np.arange(self.release['fc_initial'],
                                                                     self.release['fc_final'] +
                                                                     self.release['fc_increment'],
                                                                     self.release['fc_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['release']['force_constants'] = np.arange(0.0,
                                                                     self.release['fc_final'] +
                                                                     self.release['fc_increment'],
                                                                     self.release['fc_increment'])
            self.phase['release']['targets'] = [self.release['target']] * len(self.phase['release']['force_constants'])

        elif self.release['fraction_list'] is not None  and  self.release['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['release']['force_constants'] = [fraction * self.release['fc_final'] for
                                                        fraction in self.release['fraction_list']]
            self.phase['release']['targets'] = [self.release['target']] * len(self.phase['release']['force_constants'])

        elif self.release['fraction_increment'] is not None  and  self.release['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(0, 1.0 + self.release['fraction_increment'], self.release['fraction_increment'])
            self.phase['release']['force_constants'] = [fraction * self.release['fc_final'] for fraction in fractions]
            self.phase['release']['targets'] = [self.release['target']] * len(self.phase['release']['force_constants'])

        elif self.release['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['release']['force_constants'] = self.release['fc_list']
            self.phase['release']['targets'] = [self.release['target']] * len(self.phase['release']['force_constants'])

        elif all(v is None for k, v in self.release.items()):
            log.debug('No restraint info set for this phase! Skipping ...')

        else:
            log.error('ERROR: Restraint input did not match one of the supported methods ...')
            log.error('Input:')
            for k, v in self.release.items():
                #if v is not None: <-- NMH: Just print 'em all right?
                log.error('{} = {}'.format(k, v))
            ### PROBABLY SHOULD RAISE AN EXCEPTION HERE

        # ----------------------------------- DEBUG -------------------------------------- #

        for phase in 'attach pull release'.split():
            if self.phase[phase]['targets'] is not None:
                log.debug('Number of {} windows = {}'.format(phase,len(self.phase[phase]['targets'])))
            else:
                log.debug('This restraint will be skipped in the {} phase'.format(phase))

        # ---------------------------------- ATOM MASKS ---------------------------------- #
        log.debug('Assigning atom indices ...')
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

### NMMH

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
        iat1 = '{},'.format(restraint.index1[0])
    else:
        group1 = True
        iat1 = '-1,'
        igr1 = ''
        for index in restraint.index1:
            igr1 += '{},'.format(index)

    if not restraint.index2:
        iat2 = ' '
        raise Exception('There must be at least two atoms in a restraint.')
    elif len(restraint.index2) == 1:
        group2 = False
        iat2 = '{},'.format(restraint.index2[0])
    else:
        group2 = True
        iat2 = '-1,'
        igr2 = ''
        for index in restraint.index2:
            igr2 += '{},'.format(index)

    if not restraint.index3:
        iat3 = ''
        group3 = False
    elif len(restraint.index3) == 1:
        group3 = False
        iat3 = '{},'.format(restraint.index3[0])
    else:
        group3 = True
        iat3 = '-1,'
        igr3 = ''
        for index in restraint.index3:
            igr3 += '{},'.format(index)

    if not restraint.index4:
        iat4 = ''
        group4 = False
    elif len(restraint.index4) == 1:
        group4 = False
        iat4 = '{},'.format(restraint.index4[0])
    else:
        group4 = True
        iat4 = '-1,'
        igr4 = ''
        for index in restraint.index4:
            igr4 += '{},'.format(index)

    string = '&rst ' + \
             '\tiat = {}{}{}{}'.format(iat1,iat2,iat3,iat4)
    if group1:
        string += ' igr1 = {}'.format(igr1)
    if group2:
        string += ' igr2 = {}'.format(igr2)
    if group3:
        string += ' igr3 = {}'.format(igr3)
    if group4:
        string += ' igr4 = {}'.format(igr4)
    ### If angle or torsion, need to get correct endpoints, not 0 - 999
    string += \
        '\tr1 = {0:4.4f},'.format(0) + \
        '\tr2 = {0:4.4f},'.format(restraint.phase[phase]['targets'][window]) + \
        '\tr3 = {0:4.4f},'.format(restraint.phase[phase]['targets'][window]) + \
        '\tr4 = {0:4.4f},'.format(999) + \
        '\trk2 = {0:4.4f},'.format(restraint.phase[phase]['force_constants'][window]) + \
        '\trk3 = {0:4.4f},'.format(restraint.phase[phase]['force_constants'][window]) + \
        '\t&end'
    return string




def write_restraints_file(restraints, filename='restraints.in'):
    """
    Take all the restraints and write them to a file in each window.
    """
    ### NMH: Should we just do `for restraint in DAT_restraint.get_instances():`?
    for restraint in restraints:
        for window, _ in enumerate(restraint.phase['attach']['force_constants']):
            line = return_restraint_line(
                    restraint, phase='attach', window=window)
            # Using append mode is crucial for multiple restraints.
            directory = './windows/a{0:03d}'.format(window)
            f = open(directory + '/' + filename, 'a')
            f.write(line + '\n')
            f.close()

    for restraint in restraints:
        for window, _ in enumerate(restraint.phase['pull']['targets']):
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
        for window, _ in enumerate(restraint.phase['attach']['force_constants']):
            directory = './windows/a{0:03d}'.format(window)
            os.remove(directory + '/' + filename)

    for restraint in restraints:
        for window, _ in enumerate(restraint.phase['attach']['targets']):
            directory = './windows/p{0:03d}'.format(window)
            os.remove(directory + '/' + filename)

def error_checking():
    pass
