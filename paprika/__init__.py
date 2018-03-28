"""
The pAPRika package sets up and performs attach-pull-release calculations.
"""

import logging as log
import subprocess as sp
import re as re

from paprika.version import find_version

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

__version__ = find_version()
