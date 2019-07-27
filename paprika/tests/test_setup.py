"""
Tests the setup wrapper of paprika.
"""

import logging
import pytest

import paprika
from paprika import log

log.config_root_logger(verbose=True)
logger = logging.getLogger(__name__)

try:
    import taproom
    taproom_not_found = False
except ImportError as e:
    taproom_not_found = True

taproom_status = pytest.mark.skipif(
    taproom_not_found, reason="Benchmarks not installed"
)


@taproom_status
def test_single_orientation():
    """ Test that we can load setup YAML files. """
    setup_object = paprika.setup(host="cb6", guest="but")
    print(setup_object.desolvated_window_paths)

@taproom_status
def test_double_orientation():
    """ Test that we can load setup YAML files. """
    setup_object = paprika.setup(host="bcd", guest="hex", guest_orientation="p")
    print(setup_object.desolvated_window_paths)
    setup_object = paprika.setup(host="bcd", guest="hex", guest_orientation="s")
    print(setup_object.desolvated_window_paths)

@taproom_status
def test_generate_gaff():
    """ Test that we can load setup YAML files. """
    setup_object = paprika.setup(host="cb6", guest="but", generate_gaff_files=True)
    print(setup_object.desolvated_window_paths)

@taproom_status
def test_release():
    setup_object = paprika.setup(host="bcd", guest=None)
    print(setup_object.desolvated_window_paths)