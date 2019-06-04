"""
Tests the setup wrapper of paprika.
"""

import logging

import paprika
from paprika import log

log.config_root_logger(verbose=True)
logger = logging.getLogger(__name__)


def test_setup():
    """ Test that we can load setup YAML files. """
    setup_object = paprika.setup(host="cb6", guest="but")
    print(setup_object.desolvated_window_paths)


if __name__ == "__main__":
    test_setup()
