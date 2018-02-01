import logging as log
import numpy as np
import os as os
import sys as sys
import parmed as pmd
import subprocess as sp
import weakref as weakref

from collections import defaultdict
from paprika import utils

try:
    import simtk.openmm as mm
    import simtk.unit as unit
except:
    pass


class KeepRefs(object):
    # https://stackoverflow.com/questions/328851/printing-all-instances-of-a-class
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
    Distance or angle or torsion restraints on atoms in the simulation.
    """
    # Global lists to keep track of restraints and window counts.
    # This is helpful because I can do `max(DAT_restraint.window_counts['attach'])`
    # to get the expected number of attach windows. Some individual restraints may
    # not have attach windows, but this ensures I can check all of them.
    restraint_list = []
    window_counts = {'attach': [], 'pull': [], 'release': []}

    def __init__(self):
        # https://stackoverflow.com/questions/3484019/python-list-to-store-class-instance
        DAT_restraint.restraint_list.append(self)

        self.structure_file = None
        self.mask1 = None
        self.mask2 = None
        self.mask3 = None
        self.mask4 = None
        self.index1 = None
        self.index2 = None
        self.index3 = None
        self.index4 = None
        self.auto_apr = False  # If True, sets some pull and release values automatically.
        # If True, the first window of pull is re-used as last window of attach and the last window of pull is re-used as first window of release.
        self.continuous_apr = True

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
        super(DAT_restraint, self).__init__()

    def initialize(self):
        """
        If the user hasn't already declared the windows list, we will
        construct one from the initial, final, and increment values.
        """
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

        if self.attach['num_windows'] is not None and self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['attach']['force_constants'] = np.linspace(
                    self.attach['fc_initial'], self.attach['fc_final'], self.attach['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['attach']['force_constants'] = np.linspace(0.0, self.attach['fc_final'],
                                                                      self.attach['num_windows'])
            self.phase['attach']['targets'] = np.asarray([self.attach['target']] * self.attach['num_windows'])

        elif self.attach['fc_increment'] is not None and self.attach['fc_final'] is not None:
            if self.attach['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['attach']['force_constants'] = np.arange(
                    self.attach['fc_initial'], self.attach['fc_final'] + self.attach['fc_increment'],
                    self.attach['fc_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['attach']['force_constants'] = np.arange(0.0,
                                                                    self.attach['fc_final'] +
                                                                    self.attach['fc_increment'],
                                                                    self.attach['fc_increment'])
            self.phase['attach']['targets'] = np.asarray([self.attach['target']
                ] * len(self.phase['attach']['force_constants']))

        elif self.attach['fraction_list'] is not None and self.attach['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['attach']['force_constants'] = np.asarray([fraction * self.attach['fc_final'] for
                                                       fraction in self.attach['fraction_list']])
            self.phase['attach']['targets'] = np.asarray([self.attach['target']
                ] * len(self.phase['attach']['force_constants']))

        elif self.attach['fraction_increment'] is not None and self.attach['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(
                0, 1.0 + self.attach['fraction_increment'],
                self.attach['fraction_increment'])
            self.phase['attach']['force_constants'] = np.asarray([fraction * self.attach['fc_final']
                for fraction in fractions])
            self.phase['attach']['targets'] = np.asarray([self.attach['target']
                ] * len(self.phase['attach']['force_constants']))

        elif self.attach['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['attach']['force_constants'] = np.asarray(self.attach['fc_list'])
            self.phase['attach']['targets'] = np.asarray([self.attach['target']
                ] * len(self.phase['attach']['force_constants']))

        elif all(v is None for k, v in self.attach.items()):
            log.debug('No restraint info set for this phase! Skipping...')

        else:
            log.error('Restraint input did not match one of the supported methods...')
            for k, v in self.attach.items():
                log.debug('{} = {}'.format(k, v))
            sys.exit(1)
        # ------------------------------------ PULL ------------------------------------ #
        log.debug('Calculating pull targets and force constants...')

        if self.auto_apr:
            self.pull['fc'] = self.phase['attach']['force_constants'][-1]
            self.pull['target_initial'] = self.phase['attach']['targets'][-1]

        if self.pull['num_windows'] is not None and self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['pull']['targets'] = np.linspace(self.pull['target_initial'], self.pull['target_final'],
                                                            self.pull['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['pull']['targets'] = np.linspace(0.0, self.pull['target_final'],
                                                            self.pull['num_windows'])
            self.phase['pull']['force_constants'] = np.asarray([self.pull['fc']] * self.pull['num_windows'])

        elif self.pull['target_increment'] is not None and self.pull['target_final'] is not None:
            if self.pull['target_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['pull']['targets'] = np.arange(self.pull['target_initial'],
                                                          self.pull['target_final'] + self.pull['target_increment'],
                                                          self.pull['target_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['pull']['targets'] = np.arange(0.0, self.pull['target_final'] +
                                                          self.pull['target_increment'],
                                                          self.pull['target_increment'])
            self.phase['pull']['force_constants'] = np.asarray([self.pull['fc']] * len(self.phase['pull']['targets']))

        elif self.pull['fraction_list'] is not None and self.pull['target_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['pull']['targets'] = np.asarray([fraction * self.pull['target_final'] for
                                             fraction in self.pull['fraction_list']])
            self.phase['pull']['force_constants'] = np.asarray([
                self.pull['fc']] * len(self.phase['pull']['targets']))

        elif self.pull['fraction_increment'] is not None and self.pull['target_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(
                0, 1.0 + self.pull['fraction_increment'],
                self.pull['fraction_increment'])
            self.phase['pull']['targets'] = np.asarray([fraction * self.pull['target_final']
                for fraction in fractions])
            self.phase['pull']['force_constants'] = np.asarray([
                self.pull['fc']] * len(self.phase['pull']['targets']))

        elif self.pull['target_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['pull']['targets'] = np.asarray(self.pull['target_list'])
            self.phase['pull']['force_constants'] = np.asarray([
                self.pull['fc']] * len(self.phase['pull']['targets']))

        elif all(v is None for k, v in self.pull.items()):
            log.debug('No restraint info set for this phase! Skipping...')

        else:
            log.error('Restraint input did not match one of the supported methods...')
            for k, v in self.pull.items():
                log.debug('{} = {}'.format(k, v))
            sys.exit(1)

        # ------------------------------------ RELEASE ------------------------------------ #
        log.debug('Calculating release targets and force constants...')

        if self.auto_apr:
            self.release['target'] = self.phase['pull']['targets'][-1]
            for key in ['fc_final', 'fc_initial', 'num_windows', 'fc_increment', 'fraction_increment',
                        'fraction_list', 'fc_list']:
                if self.attach[key] is not None and self.release[key] is None:
                    self.release[key] = self.attach[key]

        if self.release['num_windows'] is not None and self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 1 ###
                log.debug('Method #1')
                self.phase['release']['force_constants'] = np.linspace(
                    self.release['fc_initial'], self.release['fc_final'], self.release['num_windows'])
            else:
                ### METHOD 1a ###
                log.debug('Method #1a')
                self.phase['release']['force_constants'] = np.linspace(0.0, self.release['fc_final'],
                                                                       self.release['num_windows'])
            self.phase['release']['targets'] = np.asarray([
                self.release['target']] * self.release['num_windows'])

        elif self.release['fc_increment'] is not None and self.release['fc_final'] is not None:
            if self.release['fc_initial'] is not None:
                ### METHOD 2 ###
                log.debug('Method #2')
                self.phase['release']['force_constants'] = np.arange(
                    self.release['fc_initial'], self.release['fc_final'] + self.release['fc_increment'],
                    self.release['fc_increment'])
            else:
                ### METHOD 2a ###
                log.debug('Method #2a')
                self.phase['release']['force_constants'] = np.arange(0.0,
                                                                     self.release['fc_final'] +
                                                                     self.release['fc_increment'],
                                                                     self.release['fc_increment'])
            self.phase['release']['targets'] = np.asarray([self.release['target']
                ] * len(self.phase['release']['force_constants']))

        elif self.release['fraction_list'] is not None and self.release['fc_final'] is not None:
            ### METHOD 3 ###
            log.debug('Method #3')
            self.phase['release']['force_constants'] = np.asarray([fraction * self.release['fc_final'] for
                                                        fraction in self.release['fraction_list']])
            self.phase['release']['targets'] = np.asarray([self.release['target']
                ] * len(self.phase['release']['force_constants']))

        elif self.release['fraction_increment'] is not None and self.release['fc_final'] is not None:
            ### METHOD 4 ###
            log.debug('Method #4')
            fractions = np.arange(
                0, 1.0 + self.release['fraction_increment'],
                self.release['fraction_increment'])
            self.phase['release']['force_constants'] = np.asarray([
                fraction * self.release['fc_final'] for fraction in fractions])
            self.phase['release']['targets'] = np.asarray([self.release['target']] * len(self.phase['release']['force_constants']))

        elif self.release['fc_list'] is not None:
            ### METHOD 5 ###
            log.debug('Method #5')
            self.phase['release']['force_constants'] = np.asarray(self.release['fc_list'])
            self.phase['release']['targets'] = np.asarray([self.release['target']
                ] * len(self.phase['release']['force_constants']))

        elif all(v is None for k, v in self.release.items()):
            log.debug('No restraint info set for this phase! Skipping...')

        else:
            log.error('Restraint input did not match one of the supported methods...')
            for k, v in self.release.items():
                log.debug('{} = {}'.format(k, v))
            sys.exit(1)

        # ----------------------------------- WINDOWS ------------------------------------ #

        for phase in ['attach', 'pull', 'release']:
            if self.phase[phase]['targets'] is not None:
                window_count = len(self.phase[phase]['targets'])
                DAT_restraint.window_counts[phase].append(window_count)
                log.debug('Number of {} windows = {}'.format(phase, window_count))
            else:
                DAT_restraint.window_counts[phase].append(None)
                log.debug('This restraint will be skipped in the {} phase'.format(phase))

        # ---------------------------------- ATOM MASKS ---------------------------------- #
        log.debug('Assigning atom indices...')
        self.index1 = utils.index_from_mask(self.structure_file, self.mask1)
        self.index2 = utils.index_from_mask(self.structure_file, self.mask2)
        if self.mask3:
            self.index3 = utils.index_from_mask(self.structure_file, self.mask3)
        else:
            self.index3 = None
        if self.mask4:
            self.index4 = utils.index_from_mask(self.structure_file, self.mask4)
        else:
            self.index4 = None
        # If any `index` has more than one atom, mark it as a group restraint.
        if self.mask1 and len(self.index1) == 1:
            self.group1 = False
        else:
            self.group1 = True
        if self.mask2 and len(self.index2) == 1:
            self.group2 = False
        else:
            self.group2 = True
        if self.mask3 and len(self.index3) == 1:
            self.group3 = False
        else:
            self.group3 = True
        if self.mask4 and len(self.index4) == 1:
            self.group4 = False
        else:
            self.group4 = True


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
                win_counts.append(None)
        max_count = np.max(win_counts)
        if all(count is None or count == max_count for count in win_counts):
            if all_continuous_apr and phase in ('attach', 'release'):
                max_count -= 1
            if max_count > 999:
                log.info('Window name zero padding only applied for windows 0 - 999')
            window_list += [phase[0] + str('{:03.0f}'.format(val)) for val in np.arange(0,max_count,1)]
        else:
            log.error('Restraints have unequal number of windows during the {} phase.'.format(phase))
            log.debug('Window counts for each restraint are as follows:')
            log.debug(DAT_restraint.window_counts[phase])
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


def amber_restraint_line(restraint, phase, window, group=False):
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

    if restraint.group3 and not restraint.group4:
        upper_bound = 180.0

    if restraint.group3 and restraint.group4:
        lower_bound = restraint.phase[phase]['targets'][window] - 180.0
        upper_bound = restraint.phase[phase]['targets'][window] + 180.0

    # Prepare AMBER NMR-style restraint
    string = '&rst iat = {:6s}{:6s}{:6s}{:6s} '.format(iat1, iat2, iat3, iat4)
    string += \
         ' r1 = {0:10.5f},'.format(lower_bound) + \
         ' r2 = {0:10.5f},'.format(restraint.phase[phase]['targets'][window]) + \
         ' r3 = {0:10.5f},'.format(restraint.phase[phase]['targets'][window]) + \
         ' r4 = {0:10.5f},'.format(upper_bound) + \
        ' rk2 = {0:10.5f},'.format(restraint.phase[phase]['force_constants'][window]) + \
        ' rk3 = {0:10.5f},'.format(restraint.phase[phase]['force_constants'][window])

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

    string += '  &end'
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
            bond_restraint = mm.CustomBondForce('k * (r - r_0)**2')
            bond_restraint.addPerBondParameter('k')
            bond_restraint.addPerBondParameter('r_0')

            r_0 = restraint.phase[phase]['targets'][window] * \
                  0.1 * unit.nanometers
            k = restraint.phase[phase]['force_constants'][window] / \
                0.239 / 0.01 * unit.kilojoules_per_mole / unit.nanometers ** 2
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])
            system.addForce(bond_restraint)
            log.debug('Added bond restraint between {} and {} with target value = '
                      '{} and force constant = {}'.format(restraint.mask1, restraint.mask2, r_0, k))
        elif restraint.group1 is True or restraint.group2 is True:
            # http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.openmm.CustomManyParticleForce.html
            # http://getyank.org/development/_modules/yank/restraints.html
            bond_restraint = mm.CustomCentroidBondForce(2, 'k * (distance(g1, g2) - r_0)^2')
            bond_restraint.addPerBondParameter('k')
            bond_restraint.addPerBondParameter('r_0')

            r_0 = restraint.phase[phase]['targets'][window] * \
                  0.1 * unit.nanometers
            k = restraint.phase[phase]['force_constants'][window] / \
                0.239 / 0.01 * unit.kilojoules_per_mole / unit.nanometers ** 2
            g1 = bond_restraint.addGroup(restraint.index1)
            g2 = bond_restraint.addGroup(restraint.index2)
            bond_restraint.addBond([g1, g2], [k, r_0])
            system.addForce(bond_restraint)
            log.debug('Added bond restraint between {} and {} with target value = '
                      '{} and force constant = {}'.format(restraint.mask1, restraint.mask2, r_0, k))
    else:
        log.error('Unable to add bond restraint...')
        log.debug('restraint.index1 = {}'.format(restraint.index1))
        log.debug('restraint.index2 = {}'.format(restraint.index2))
        sys.exit(1)

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
            sys.exit(1)

        angle_restraint = mm.CustomAngleForce('k * (theta - theta_0)**2')
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
            sys.exit(1)

        dihedral_restraint = mm.CustomTorsionForce('k * (theta - theta_0)**2')
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
