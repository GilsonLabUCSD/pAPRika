"""
The pAPRika package sets up and performs attach-pull-release calculations.
"""

import logging as log
import parmed as pmd

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(
    format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

__all__ = ["align", "build", "restraints", "simulate"]
