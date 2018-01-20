"""
The pAPRika package sets up and performs attach-pull-release calculations.
"""

from ._version import get_versions
import logging as log
import parmed as pmd

__version__ = get_versions()['version']
del get_versions

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(
    format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

__all__ = ["align", "build", "restraints", "simulate"]
