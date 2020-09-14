"""
pAPRika
Advanced toolkit for binding free energy calculations
"""
# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Handle versioneer
from ._version import get_versions

# Add imports here
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

import logging

logger = logging.getLogger(__name__)

try:
    from simtk import openmm

    from paprika.setup import Setup

    setup = Setup
except ImportError as e:
    logging.info("OpenMM not found.")
    logging.info("`paprika.setup()` requires OpenMM.")
    setup = None

from paprika.analyze import Analyze

analyze = Analyze

if setup is None:
    __all__ = ["setup", "analyze"]
else:
    __all__ = ["analyze"]
