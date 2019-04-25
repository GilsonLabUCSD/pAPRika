"""
Test utility functions
"""

import os

import parmed as pmd

from paprika.utils import *


def test_mkdirs():
    """ Test that we can make directories for the windows. """
    window_list = ["a000", "p014", "r089"]
    make_window_dirs(window_list, path="tmp")
    for window in window_list:
        assert os.path.exists(os.path.join("tmp", "windows", window))


def test_strip_prmtop():
    """ Test that we can remove items from structures. """
    cb6_only = strip_prmtop(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.prmtop"),
        mask=":BUT",
    )
    residues = [r.name for r in cb6_only.residues]
    assert "CB6" in residues
    assert "BUT" not in residues
