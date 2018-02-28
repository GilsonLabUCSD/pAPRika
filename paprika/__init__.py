"""
The pAPRika package sets up and performs attach-pull-release calculations.
"""

import logging as log
import subprocess as sp
import re as re

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(
    format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

__version__ = '0.0.3'
try:
    # Try to use git to find current commit.
    p = sp.Popen(
        ["git", "describe", "--always"], stdout=sp.PIPE, stderr=sp.PIPE)
    output, error = p.communicate()
    git_describe = output.decode("utf-8").strip()
    __version__ = re.sub('-g[0-9a-f]*$', '', git_describe)
except sp.CalledProcessError as e:
    print('Git error...')
    pass
