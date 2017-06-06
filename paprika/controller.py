import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd
from paprika import align


def make_directories(restraint):
    """
    Make a series of directories to hold the simulation setup files
    and the data. This function should probably end up in a separate
    file eventually.
    """
    # Here we could check if the directories already exist and prompt
    # the user or quit or do something else.

    log.debug('We ought to make sure somewhere that all restraints have'
              ' the same number of windows. Here I am creating directories based '
              ' on the number of windows in a single restraint.')
    for window, _ in enumerate(restraint.phase['attach']['force_constants']):
        if not os.path.exists('./windows/a{0:03d}'.format(window)):
            os.makedirs('./windows/a{0:03d}'.format(window))
    for window, _ in enumerate(restraint.phase['pull']['targets']):
        if not os.path.exists('./windows/p{0:03d}'.format(window)):
            os.makedirs('./windows/p{0:03d}'.format(window))
