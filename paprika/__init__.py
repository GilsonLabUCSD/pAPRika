"""
pAPRika
Advanced toolkit for binding free energy calculations
"""
# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
import logging
logger = logging.getLogger(__name__)
from paprika import log
log.config_root_logger(verbose=False)

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
