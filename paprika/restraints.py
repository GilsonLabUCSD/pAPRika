import logging as log
import numpy as np
import os as os
import sys as sys
import parmed as pmd
import subprocess as sp
import weakref as weakref
import pytraj as pt
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from collections import defaultdict
from paprika import utils

try:
    import simtk.openmm as mm
    import simtk.unit as unit
except:
    pass


class DAT_restraint(object):
    """
    Distance or angle or torsion restraints on atoms in the simulation.
    """

    def __init__(self):

        self.topology = None
        self.mask1 = None
        self.mask2 = None
        self.mask3 = None
        self.mask4 = None
        self.index1 = None
        self.index2 = None
        self.index3 = None
        self.index4 = None

        # In the case of a non-harmonic restraint, the pre-calculated values can be overridden with ones from this dictionary.
        self.custom_restraint_values = {}

        self.auto_apr = False  # If True, sets some pull and release values automatically.
        # If True, the first window of pull is re-used as last window of attach and the last window of pull is re-used as first window of release.
        self.continuous_apr = True
        self.amber_index = False

        self.attach = {
            'target': None,  # The target value for the restraint (mandatory)
            'fc_initial': None,  # The initial force constant (optional)
            'fc_final': None,  # The final force constant (optional)
            'num_windows': None,  # The number of windows (optional)
            'fc_increment': None,  # The force constant increment (optional)
            'fraction_increment': None,  # The percentage of the force constant increment (optional)
            'fraction_list': None,  # The list of force constant percentages (optional)
            'fc_list': None  # The list of force constants (will be created if not given)
        }

        self.pull = {
            'fc': None,  # The force constant for the restraint (mandatory)
            'target_initial': None,  # The initial target value (optional)
            'target_final': None,  # The final target value (optional)
            'num_windows': None,  # The number of windows (optional)
            'target_increment': None,  # The target value increment (optional)
            'fraction_increment': None,  # The percentage of the target value increment (optional)
            'fraction_list': None,  # The list of target value percentages (optional)
            'target_list': None  # The list of target values (will be created if not given)
        }

        self.release = {
            'target': None,  # The target value for the restraint (mandatory)
            'fc_initial': None,  # The initial force constant (optional)
            'fc_final': None,  # The final force constant (optional)
            'num_windows': None,  # The number of windows (optional)
            'fc_increment': None,  # The force constant increment (optional)
            'fraction_increment': None,  # The percentage of the force constant increment (optional)
            'fraction_list': None,  # The list of force constant percentages (optional)
            'fc_list': None  # The list of force constants (will be created if not
            # given)
        }

    def _calc_meth(self, phase, rdict, meth):
        """ Return the appropriate list of force_constants and targets depending on the method """

        force_constants = None
        targets = None

        # Attach/Release, Force Constant Method 1
        if phase in ('a', 'r') and meth == '1':
            force_constants = np.linspace(rdict['fc_initial'], rdict['fc_final'], rdict['num_windows'])

        # Attach/Release, Force Constant Method 1a
        elif phase in ('a', 'r') and meth == '1a':
            force_constants = np.linspace(0.0, rdict['fc_final'], rdict['num_windows'])

        # Attach/Release, Force Constant Method 2
        elif phase in ('a', 'r') and meth == '2':
            force_constants = np.arange(rdict['fc_initial'], rdict['fc_final'] + rdict['fc_increment'],
                                        rdict['fc_increment'])

        # Attach/Release, Force Constant Method 2a
        elif phase in ('a', 'r') and meth == '2a':
            force_constants = np.arange(0.0, rdict['fc_final'] + rdict['fc_increment'], rdict['fc_increment'])

        # Attach/Release, Force Constant Method 3
        elif phase in ('a', 'r') and meth == '3':
            force_constants = np.asarray([fraction * rdict['fc_final'] for fraction in rdict['fraction_list']])

        # Attach/Release, Force Constant Method 4
        elif phase in ('a', 'r') and meth == '4':
            fractions = np.arange(0, 1.0 + rdict['fraction_increment'], rdict['fraction_increment'])
            force_constants = np.asarray([fraction * rdict['fc_final'] for fraction in fractions])

        # Attach/Release, Force Constant Method 5
        elif phase in ('a', 'r') and meth == '5':
            force_constants = np.asarray(rdict['fc_list'])

        # Attach/Release, Target Method
        if phase in ('a', 'r'):
            targets = np.asarray([rdict['target']] * len(force_constants))

        # Pull, Target Method 1
        if phase == 'p' and meth == '1':
            targets = np.linspace(rdict['target_initial'], rdict['target_final'], rdict['num_windows'])

        # Pull, Target Method 1a
        elif phase == 'p' and meth == '1a':
            targets = np.linspace(0.0, rdict['target_final'], rdict['num_windows'])

        # Pull, Target Method 2
        elif phase == 'p' and meth == '2':
            targets = np.arange(rdict['target_initial'], rdict['target_final'] + rdict['target_increment'],
                                rdict['target_increment'])

        # Pull, Target Method 2a
        elif phase == 'p' and meth == '2a':
            targets = np.arange(0.0, rdict['target_final'] + rdict['target_increment'], rdict['target_increment'])

        # Pull, Target Method 3
        elif phase == 'p' and meth == '3':
            targets = np.asarray([fraction * rdict['target_final'] for fraction in rdict['fraction_list']])

        # Pull, Target Method 4
        elif phase == 'p' and meth == '4':
            fractions = np.arange(0, 1.0 + rdict['fraction_increment'], rdict['fraction_increment'])
            targets = np.asarray([fraction * rdict['target_final'] for fraction in fractions])

        # Pull, Target Method 5
        elif phase == 'p' and meth == '5':
            targets = np.asarray(rdict['target_list'])

        # Pull, Force Constant Method
        if phase == 'p':
            force_constants = np.asarray([rdict['fc']] * len(targets))

        if force_constants is None and targets is None:
            log.error('Unsupported Phase/Method: {} / {}'.format(phase, meth))
            raise Exception('Unexpected phase/method combination passed to _calc_meth')

        return force_constants, targets

    def initialize(self):
        """
        Depending on which dict values are provided for each phase, a different method will
        be used to determine the list of force_constants and targets (below).

        For Attach/Release, a `target` value is required and the method is determined if the
        following dict values are not `None`:
            Method 1:   num_windows, fc_initial, fc_final
            Method 1a:  num_windows, fc_final
            Method 2:   fc_increment, fc_initial, fc_final
            Method 2a:  fc_increment, fc_final
            Method 3:   fraction_list, fc_final
            Method 4:   fraction_increment, fc_final
            Method 5:   fc_list

        For Pull, a `fc` value is required and the method is determined if the
        following dict values are not `None`:
            Method 1:   num_windows, target_initial, target_final
            Method 1a:  num_windows, target_final
            Method 2:   target_increment, target_initial, target_final
            Method 2a:  target_increment, target_final
            Method 3:   fraction_list, target_final
            Method 4:   fraction_increment, target_final
            Method 5:   target_list
        """

        ### Setup. These are the lists that will be most used by other modules
        self.phase = {
            'attach': {
                'force_constants': None,
                'targets': None
            },
            'pull': {
                'force_constants': None,
                'targets': None
            },
            'release': {
                'force_constants': None,
                'targets': None
            }
        }
        # ------------------------------------ ATTACH ------------------------------------ #
        log.debug('Calculating attach targets and force constants...')

        ### Temporary variables to improve readability
        force_constants = None
        targets = None

        if self.attach['num_windows'] is not None and self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Attach, Method #1')
                force_constants, targets = self._calc_meth('a', self.attach, '1')
            else:
                ### METHOD 1a ###
                log.debug('Attach, Method #1a')
                force_constants, targets = self._calc_meth('a', self.attach, '1a')

        elif self.attach['fc_increment'] is not None and self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Attach, Method #2')
                force_constants, targets = self._calc_meth('a', self.attach, '2')
            else:
                ### METHOD 2a ###
                log.debug('Attach, Method #2a')
                force_constants, targets = self._calc_meth('a', self.attach, '2a')

        elif self.attach['fraction_list'] is not None and self.attach['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Attach, Method #3')
            force_constants, targets = self._calc_meth('a', self.attach, '3')

        elif self.attach['fraction_increment'] is not None and self.attach['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Attach, Method #4')
            force_constants, targets = self._calc_meth('a', self.attach, '4')

        elif self.attach['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Attach, Method #5')
            force_constants, targets = self._calc_meth('a', self.attach, '5')

        elif all(v is None for k, v in self.attach.items()):
            log.debug('No restraint info set for the attach phase! Skipping...')

        else:
            log.error('Attach restraint input did not match one of the supported methods...')
            for k, v in self.attach.items():
                log.debug('{} = {}'.format(k, v))
            raise Exception('Attach restraint input did not match one of the supported methods...')

        if force_constants is not None and targets is not None:
            self.phase['attach']['force_constants'] = force_constants
            self.phase['attach']['targets'] = targets

        # ------------------------------------ PULL ------------------------------------ #
        log.debug('Calculating pull targets and force constants...')

        force_constants = None
        targets = None

        if self.auto_apr and self.pull['target_final'] is not None:
            self.pull['fc'] = self.phase['attach']['force_constants'][-1]
            self.pull['target_initial'] = self.phase['attach']['targets'][-1]

        if self.pull['num_windows'] is not None and self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Pull, Method #1')
                force_constants, targets = self._calc_meth('p', self.pull, '1')
            else:
                ### METHOD 1a ###
                log.debug('Pull, Method #1a')
                force_constants, targets = self._calc_meth('p', self.pull, '1a')

        elif self.pull['target_increment'] is not None and self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Pull, Method #2')
                force_constants, targets = self._calc_meth('p', self.pull, '2')
            else:
                ### METHOD 2a ###
                log.debug('Pull, Method #2a')
                force_constants, targets = self._calc_meth('p', self.pull, '2a')

        elif self.pull['fraction_list'] is not None and self.pull['target_final'] is not None:
            ### METHOD 3 ###
            log.debug('Pull, Method #3')
            force_constants, targets = self._calc_meth('p', self.pull, '3')

        elif self.pull['fraction_increment'] is not None and self.pull['target_final'] is not None:
            ### METHOD 4 ###
            log.debug('Pull, Method #4')
            force_constants, targets = self._calc_meth('p', self.pull, '4')

        elif self.pull['target_list'] is not None:
            ### METHOD 5 ###
            log.debug('Pull, Method #5')
            force_constants, targets = self._calc_meth('p', self.pull, '5')

        elif all(v is None for k, v in self.pull.items()):
            log.debug('No restraint info set for the pull phase! Skipping...')

        else:
            log.error('Pull restraint input did not match one of the supported methods...')
            for k, v in self.pull.items():
                log.debug('{} = {}'.format(k, v))
            raise Exception('Pull restraint input did not match one of the supported methods...')

        if force_constants is not None and targets is not None:
            self.phase['pull']['force_constants'] = force_constants
            self.phase['pull']['targets'] = targets

        # ------------------------------------ RELEASE ------------------------------------ #
        log.debug('Calculating release targets and force constants...')

        force_constants = None
        targets = None

        # I don't want auto_apr to make release restraints, unless I'm sure the user wants them.
        # I'm gonna assume that specifying self.attach['fc_final'] indicates you want it,
        # although this weakens the whole purpose of auto_apr
        if self.auto_apr and self.release['fc_final'] is not None:
            self.release['target'] = self.phase['pull']['targets'][-1]
            for key in [
                    'fc_final', 'fc_initial', 'num_windows', 'fc_increment', 'fraction_increment', 'fraction_list',
                    'fc_list'
            ]:
                if self.attach[key] is not None and self.release[key] is None:
                    self.release[key] = self.attach[key]

        if self.release['num_windows'] is not None and self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Release, Method #1')
                force_constants, targets = self._calc_meth('r', self.release, '1')
            else:
                ### METHOD 1a ###
                log.debug('Release, Method #1a')
                force_constants, targets = self._calc_meth('r', self.release, '1a')

        elif self.release['fc_increment'] is not None and self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Release, Method #2')
                force_constants, targets = self._calc_meth('r', self.release, '2')
            else:
                ### METHOD 2a ###
                log.debug('Release, Method #2a')
                force_constants, targets = self._calc_meth('r', self.release, '2a')

        elif self.release['fraction_list'] is not None and self.release['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Release, Method #3')
            force_constants, targets = self._calc_meth('r', self.release, '3')

        elif self.release['fraction_increment'] is not None and self.release['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Release, Method #4')
            force_constants, targets = self._calc_meth('r', self.release, '4')

        elif self.release['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Release, Method #5')
            force_constants, targets = self._calc_meth('r', self.release, '5')

        elif all(v is None for k, v in self.release.items()):
            log.debug('No restraint info set for the release phase! Skipping...')

        else:
            log.error('Release restraint input did not match one of the supported methods...')
            for k, v in self.release.items():
                log.debug('{} = {}'.format(k, v))
            raise Exception('Release restraint input did not match one of the supported methods...')

        if force_constants is not None and targets is not None:
            self.phase['release']['force_constants'] = force_constants
            self.phase['release']['targets'] = targets

        # ----------------------------------- WINDOWS ------------------------------------ #

        for phase in ['attach', 'pull', 'release']:
            if self.phase[phase]['targets'] is not None:
                window_count = len(self.phase[phase]['targets'])
                #DAT_restraint.window_counts[phase].append(window_count)
                log.debug('Number of {} windows = {}'.format(phase, window_count))
            else:
                #DAT_restraint.window_counts[phase].append(None)
                log.debug('This restraint will be skipped in the {} phase'.format(phase))

        # ---------------------------------- ATOM MASKS ---------------------------------- #
        log.debug('Assigning atom indices...')
        self.index1 = utils.index_from_mask(self.topology, self.mask1, self.amber_index)
        self.index2 = utils.index_from_mask(self.topology, self.mask2, self.amber_index)
        if self.mask3:
            self.index3 = utils.index_from_mask(self.topology, self.mask3, self.amber_index)
        else:
            self.index3 = None
        if self.mask4:
            self.index4 = utils.index_from_mask(self.topology, self.mask4, self.amber_index)
        else:
            self.index4 = None
        # If any `index` has more than one atom, mark it as a group restraint.
        #print('Masks:',self.mask1, self.mask2, self.mask3, self.mask4)
        #print('index:',self.index1, self.index2, self.index3, self.index4)
        if self.mask1 and len(self.index1) > 1:
            self.group1 = True
        else:
            self.group1 = False
        if self.mask2 and len(self.index2) > 1:
            self.group2 = True
        else:
            self.group2 = False
        if self.mask3 and len(self.index3) > 1:
            self.group3 = True
        else:
            self.group3 = False
        if self.mask4 and len(self.index4) > 1:
            self.group4 = True
        else:
            self.group4 = False
        #print('index:',self.group1, self.group2, self.group3, self.group4)


def static_DAT_restraint(restraint_mask_list,
                         num_window_list,
                         ref_structure,
                         force_constant,
                         continuous_apr=True,
                         amber_index=False):
    """ Create a static restraint """

    # Setup reference structure
    if isinstance(ref_structure, str):
        ref_structure = utils.return_parmed_structure(ref_structure)
    elif isinstance(ref_structure, pmd.structure.Structure):
        pass
    else:
        raise Exception(
            'static_DAT_restraint does not support the type associated with ref_structure:' + type(ref_structure))
    ref_traj = pt.load_parmed(ref_structure, traj=True)

    # Check num_window_list
    if len(num_window_list) != 3:
        raise Exception(
            'The num_window_list needs to contain three integers corresponding to the number of windows in the attach, pull, and release phase, respectively'
        )

    # Setup restraint
    rest = DAT_restraint()
    rest.continuous_apr = continuous_apr
    rest.amber_index = amber_index
    rest.topology = ref_structure
    rest.mask1 = restraint_mask_list[0]
    rest.mask2 = restraint_mask_list[1]
    if len(restraint_mask_list) >= 3:
        rest.mask3 = restraint_mask_list[2]
    if len(restraint_mask_list) == 4:
        rest.mask4 = restraint_mask_list[3]

    # Target value
    mask_string = ' '.join(restraint_mask_list)
    if len(restraint_mask_list) == 2:
        # Distance restraint
        target = pt.distance(ref_traj, mask_string)[0]
    elif len(restraint_mask_list) == 3:
        # Angle restraint
        target = pt.angle(ref_traj, mask_string)[0]
    elif len(restraint_mask_list) == 4:
        # Dihedral restraint
        target = pt.dihedral(ref_traj, mask_string)[0]
    else:
        raise Exception(
            'The number of masks (' + str(len(restraint_mask_list)) +
            ') in restraint_mask_list is not 2, 3, or 4 and thus is not one of the supported types: distance, angle, dihedral'
        )

    # Attach phase
    if num_window_list[0] is not None and num_window_list[0] != 0:
        rest.attach['target'] = target
        rest.attach['fc_initial'] = force_constant
        rest.attach['fc_final'] = force_constant
        rest.attach['num_windows'] = num_window_list[0]

    # Pull phase
    if num_window_list[1] is not None and num_window_list[1] != 0:
        rest.pull['fc'] = force_constant
        rest.pull['target_initial'] = target
        rest.pull['target_final'] = target
        rest.pull['num_windows'] = num_window_list[1]

    # Release phase
    if num_window_list[2] is not None and num_window_list[2] != 0:
        rest.attach['target'] = target
        rest.attach['fc_initial'] = force_constant
        rest.attach['fc_final'] = force_constant
        rest.attach['num_windows'] = num_window_list[2]

    rest.initialize()

    return rest


def check_restraints(restraint_list, create_window_list=False):
    """
    Do basic tests to ensure a list of DAT_restraints are consistent.
    We're gonna create the window list here too, because it needs the same code.
    """

    if all(restraint.continuous_apr is True for restraint in restraint_list):
        log.debug('All restraints are "continuous_apr" style.')
        all_continuous_apr = True
    elif all(restraint.continuous_apr is False for restraint in restraint_list):
        log.debug('All restraints are not "continuous_apr" style.')
        all_continuous_apr = False
    else:
        log.error('All restraints must have the same setting for .continuous_apr')
        ### Should we do the following?
        raise Exception('All restraints must have the same setting for .continuous_apr')

    window_list = []
    phases = ['attach', 'pull', 'release']
    for phase in phases:
        win_counts = []
        for restraint in restraint_list:
            if restraint.phase[phase]['targets'] is not None:
                win_counts.append(len(restraint.phase[phase]['targets']))
            else:
                win_counts.append(0)
        max_count = np.max(win_counts)
        if all(count == 0 or count == max_count for count in win_counts):
            if max_count > 0:
                if all_continuous_apr and phase in ('attach', 'release'):
                    max_count -= 1
                if max_count > 999:
                    log.info('Window name zero padding only applied for windows 0 - 999')
                window_list += [phase[0] + str('{:03.0f}'.format(val)) for val in np.arange(0, max_count, 1)]
        else:
            log.error('Restraints have unequal number of windows during the {} phase.'.format(phase))
            log.debug('Window counts for each restraint are as follows:')
            log.debug(win_counts)
            raise Exception('Restraints have unequal number of windows during the {} phase.'.format(phase))

    log.info('Restraints appear to be consistent')

    if create_window_list:
        return window_list


def create_window_list(restraint_list):
    """
    Create list of APR windows. Runs everything through check_restraints because
    we need to do that.
    """

    return check_restraints(restraint_list, create_window_list=True)


def amber_restraint_line(restraint, phase, window):
    """
    Return an AMBER restraint line a specific phase/window combination.
    For example:
    &rst iat= 3,109, r1= 0.0000, r2= 6.9665, r3= 6.9665, r4= 999.0000,
         rk2= 5.0000000, rk3= 5.0000000, &end
    Or:
    &rst iat= -1,-1, igr1=3,4,7,8,21,22,25,26,39,40,43,44,57,58,61,62,75,76,79,
    80,93,94,97,98, igr2=109,113,115,119,
    r1=     0.0000, r2=     5.9665, r3=     5.9665, r4=   999.0000,
    rk2=   5.0000000, rk3=   5.0000000, &end

    """
    if not restraint.index1:
        iat1 = ' '
        raise Exception('There must be at least two atoms in a restraint.')
    elif not restraint.group1:
        iat1 = '{},'.format(restraint.index1[0])
    else:
        iat1 = '-1,'
        igr1 = ''
        for index in restraint.index1:
            igr1 += '{},'.format(index)

    if not restraint.index2:
        iat2 = ' '
        raise Exception('There must be at least two atoms in a restraint.')
    elif not restraint.group2:
        iat2 = '{},'.format(restraint.index2[0])
    else:
        iat2 = '-1,'
        igr2 = ''
        for index in restraint.index2:
            igr2 += '{},'.format(index)

    if not restraint.index3:
        iat3 = ''
    elif not restraint.group3:
        iat3 = '{},'.format(restraint.index3[0])
    else:
        iat3 = '-1,'
        igr3 = ''
        for index in restraint.index3:
            igr3 += '{},'.format(index)

    if not restraint.index4:
        iat4 = ''
    elif not restraint.group4:
        iat4 = '{},'.format(restraint.index4[0])
    else:
        iat4 = '-1,'
        igr4 = ''
        for index in restraint.index4:
            igr4 += '{},'.format(index)

    # Set upper/lower bounds depending on whether distance, angle, or torsion
    lower_bound = 0.0
    upper_bound = 999.0

    if restraint.mask3 and not restraint.mask4:
        upper_bound = 180.0

    if restraint.mask3 and restraint.mask4:
        lower_bound = restraint.phase[phase]['targets'][window] - 180.0
        upper_bound = restraint.phase[phase]['targets'][window] + 180.0

    amber_restraint_values = {
        'r1': lower_bound,
        'r2': restraint.phase[phase]['targets'][window],
        'r3': restraint.phase[phase]['targets'][window],
        'r4': upper_bound,
        'rk2': restraint.phase[phase]['force_constants'][window],
        'rk3': restraint.phase[phase]['force_constants'][window]
    }

    for key, value in restraint.custom_restraint_values.items():
        if value is not None:
            log.debug('Overriding {} = {}'.format(key, value))
            amber_restraint_values[key] = value

    # Prepare AMBER NMR-style restraint
    atoms = ''.join([iat1, iat2, iat3, iat4])
    string = '&rst iat = {:16s} '.format(atoms)
    string += \
         ' r1 = {0:10.5f},'.format(amber_restraint_values['r1']) + \
         ' r2 = {0:10.5f},'.format(amber_restraint_values['r2']) + \
         ' r3 = {0:10.5f},'.format(amber_restraint_values['r3']) + \
         ' r4 = {0:10.5f},'.format(amber_restraint_values['r4']) + \
        ' rk2 = {0:10.5f},'.format(amber_restraint_values['rk2']) + \
        ' rk3 = {0:10.5f},'.format(amber_restraint_values['rk3'])

    if any([restraint.group1, restraint.group2, restraint.group3, restraint.group4]):
        string += '\n    '
        if restraint.group1:
            string += ' igr1 = {}'.format(igr1)
        if restraint.group2:
            string += ' igr2 = {}'.format(igr2)
        if restraint.group3:
            string += ' igr3 = {}'.format(igr3)
        if restraint.group4:
            string += ' igr4 = {}'.format(igr4)

    string += "  &end\n"
    return string


def setup_openmm_restraints(system, restraint, phase, window):
    """
    Add particle restraints with OpenMM.
    """

    # http://docs.openmm.org/7.1.0/api-c++/generated/OpenMM.CustomExternalForce.html
    # It's possible we might need to use `periodicdistance`.

    if restraint.mask1 is not None and \
                    restraint.mask2 is not None and \
                    restraint.mask3 is None and \
                    restraint.mask4 is None:

        if restraint.group1 is False and restraint.group2 is False:
            bond_restraint = mm.CustomBondForce('k * (r - r_0)^2')
            bond_restraint.addPerBondParameter('k')
            bond_restraint.addPerBondParameter('r_0')

            r_0 = restraint.phase[phase]['targets'][window] * unit.angstrom
            k = restraint.phase[phase]['force_constants'][window] * unit.kilocalorie_per_mole / unit.angstrom**2
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])
            bond_restraint.setForceGroup(1)
            system.addForce(bond_restraint)
            log.debug('Added bond restraint between {} and {} with target value = '
                      '{} and force constant = {}'.format(restraint.mask1, restraint.mask2, r_0, k))
        elif restraint.group1 is True or restraint.group2 is True:
            # http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.openmm.CustomManyParticleForce.html
            # http://getyank.org/development/_modules/yank/restraints.html
            bond_restraint = mm.CustomCentroidBondForce(2, 'k * (distance(g1, g2) - r_0)^2')
            bond_restraint.addPerBondParameter('k')
            bond_restraint.addPerBondParameter('r_0')

            r_0 = restraint.phase[phase]['targets'][window] * unit.angstrom
            k = restraint.phase[phase]['force_constants'][window] * unit.kilocalorie_per_mole / unit.angstrom**2
            g1 = bond_restraint.addGroup(restraint.index1)
            g2 = bond_restraint.addGroup(restraint.index2)
            bond_restraint.addBond([g1, g2], [k, r_0])
            bond_restraint.setForceGroup(1)
            system.addForce(bond_restraint)
            log.debug('Added bond restraint between {} and {} with target value = '
                      '{} and force constant = {}'.format(restraint.mask1, restraint.mask2, r_0, k))
    else:
        log.error('Unable to add bond restraint...')
        log.debug('restraint.index1 = {}'.format(restraint.index1))
        log.debug('restraint.index2 = {}'.format(restraint.index2))
        raise Exception('Unable to add bond restraint...')

    if restraint.mask1 is not None and \
       restraint.mask2 is not None and \
       restraint.mask3 is not None and \
       restraint.mask4 is None:
        if restraint.group1 is not False and \
           restraint.group2 is not False and \
           restraint.group3 is not False:
            log.error('Unable to add a group angle restraint...')
            log.debug('restraint.index1 = {}'.format(restraint.index1))
            log.debug('restraint.index2 = {}'.format(restraint.index2))
            log.debug('restraint.index3 = {}'.format(restraint.index3))
            raise Exception('Unable to add a group angle restraint...')

        angle_restraint = mm.CustomAngleForce('k * (theta - theta_0)^2')
        angle_restraint.addPerAngleParameter('k')
        angle_restraint.addPerAngleParameter('theta_0')

        log.debug('Setting an angle restraint in degrees using a ' 'force constant in kcal per mol rad**2...')
        theta_0 = restraint.phase[phase]['targets'][window] * unit.degrees
        k = restraint.phase[phase]['force_constants'][window] * \
            unit.kilocalorie_per_mole / unit.radian ** 2
        angle_restraint.addAngle(restraint.index1[0], restraint.index2[0], restraint.index3[0], [k, theta_0])
        system.addForce(angle_restraint)

    if restraint.mask1 is not None and \
       restraint.mask2 is not None and \
       restraint.mask3 is not None and \
       restraint.mask4 is not None:
        if restraint.group1 is not False and \
           restraint.group2 is not False and \
           restraint.group3 is not False and \
           restraint.group4 is not False:
            log.error('Unable to add a group dihedral restraint...')
            log.debug('restraint.index1 = {}'.format(restraint.index1))
            log.debug('restraint.index2 = {}'.format(restraint.index2))
            log.debug('restraint.index3 = {}'.format(restraint.index3))
            log.debug('restraint.index4 = {}'.format(restraint.index4))
            raise Exception('Unable to add a group dihedral restraint...')

        dihedral_restraint = mm.CustomTorsionForce('k * (theta - theta_0)^2')
        dihedral_restraint.addPerTorsionParameter('k')
        dihedral_restraint.addPerTorsionParameter('theta_0')

        log.debug('Setting a torsion restraint in degrees using a ' 'force constant in kcal per mol rad**2...')
        theta_0 = restraint.phase[phase]['targets'][window] * unit.degrees
        k = restraint.phase[phase]['force_constants'][window] * \
            unit.kilocalorie_per_mole / unit.radian ** 2
        dihedral_restraint.addTorsion(restraint.index1[0], restraint.index2[0], restraint.index3[0],
                                      restraint.index4[0], [k, theta_0])
        system.addForce(dihedral_restraint)

    return system


def clean_restraints_file(restraints, filename='restraints.in'):
    """
    Delete the restraints files for repeated testing.

    Parameters
    ----------
    restraints : object
    """

    log.warning('`clean_restraints_file()` needs to be tested.')
    for restraint in restraints:
        for window, _ in enumerate(restraint.phase['attach']['force_constants']):
            directory = './windows/a{0:03d}'.format(window)
            os.remove(directory + '/' + filename)

    for restraint in restraints:
        for window, _ in enumerate(restraint.phase['attach']['targets']):
            directory = './windows/p{0:03d}'.format(window)
            os.remove(directory + '/' + filename)


def error_checking():
    print('Error checking needs to be implemented.')
    pass
