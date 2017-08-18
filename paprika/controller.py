import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd
from paprika import align
from paprika.restraints import DAT_restraint


def check_restraint_windows():
    """
    Check that all restraints have the same number of attach, pull, and
    release windows.
    """

    phases = ['attach', 'pull', 'release']
    windows = {}
    for phase in phases:
        windows[phase] = []
    # For each restraint, for each phase, if it exists, count the windows and add to a list
    for restraint in DAT_restraint.get_instances():
        for phase in phases:
            if restraint.phase[phase]['targets'] is not None:
                windows[phase].append(len(restraint.phase[phase]['targets']))

    # Take the window count for the first restraint, add up how many total restraints 
    # have the same number, does that total equal the restraint total?
    for phase in phases:
        if windows[phase]:
            if not windows[phase].count(windows[phase][0]) == len(windows[phase]):
                log.error('Restraints have unequal number of windows during the {} phase.'.format(phase))
                log.info(windows[phase])
            log.debug('Passed {} check. There are {} windows.'.format(phase,windows[phase][0]))
        else:
            log.debug('There are no {} windows.'.format(phase))


def make_directories(restraint):
    """
    Make a series of directories to hold the simulation setup files
    and the data. Here we could check if the directories already exist and prompt
    the user or quit or do something else.
    """

    for phase in ['attach', 'pull', 'release']:
        if restraint.phase[phase]['targets'] is not None:
            for window in range(len(restraint.phase[phase]['targets'])):
                if not os.path.exists('./windows/{0:1s}{1:03d}'.format(phase[0:1],window)):
                    os.makedirs('./windows/{0:1s}{1:03d}'.format(phase[0:1],window))


