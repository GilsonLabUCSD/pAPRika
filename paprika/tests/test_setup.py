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

def test_generate_gaff():
    """ Test that we can load setup YAML files. """
    setup_object = paprika.setup(host="cb6", guest="but", generate_gaff_files=True)
    print(setup_object.desolvated_window_paths)


if __name__ == "__main__":
    test_generate_gaff()
