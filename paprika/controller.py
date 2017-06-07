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
    attach_windows = []
    pull_windows = []
    release_windows = []
    for restraint in DAT_restraint.get_instances():
        attach_windows.append(len(restraint.phase['attach']['targets']))
        pull_windows.append(len(restraint.phase['pull']['targets']))
        release_windows.append(len(restraint.phase['release']['targets']))

    if not attach_windows.count(attach_windows[0]) == len(attach_windows):
        log.error('Restraints have unequal number of windows during the attachment phase.')
        log.info(attach_windows)
    if not pull_windows.count(pull_windows[0]) == len(pull_windows):
        log.error('Restraints have unequal number of windows during the pulling phase.')
        log.info(pull_windows)
    if not release_windows.count(release_windows[0]) == len(release_windows):
        log.error('Restraints have unequal number of windows during the release phase.')
        log.info(release_windows)

def make_directories(num_attach_windows, num_pull_windows, num_release_windows):
    """
    Make a series of directories to hold the simulation setup files
    and the data. Here we could check if the directories already exist and prompt
    the user or quit or do something else.
    """

    for window in range(num_attach_windows):
        if not os.path.exists('./windows/a{0:03d}'.format(window)):
            os.makedirs('./windows/a{0:03d}'.format(window))
    for window in range(num_pull_windows):
        if not os.path.exists('./windows/p{0:03d}'.format(window)):
            os.makedirs('./windows/p{0:03d}'.format(window))
    for window in range(num_release_windows):
        if not os.path.exists('./windows/r{0:03d}'.format(window)):
            os.makedirs('./windows/r{0:03d}'.format(window))