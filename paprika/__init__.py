"""
pAPRika
Advanced toolkit for binding free energy calculations
"""
# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
import logging

from paprika import log

# Handle versioneer
from ._version import get_versions

logger = logging.getLogger(__name__)
log.config_root_logger(verbose=False)

versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
