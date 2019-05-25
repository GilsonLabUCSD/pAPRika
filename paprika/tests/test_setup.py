"""
Tests the setup wrapper of paprika.
"""

import logging

import pytest

import paprika
from paprika import log

log.config_root_logger(verbose=True)
logger = logging.getLogger(__name__)

def test_setup():
    """ Test that we can load setup YAML files. """
    paprika.setup(host="cb6", guest="but")