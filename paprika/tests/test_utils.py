"""
Test utility functions
"""

import os
import shutil

import pytest

from paprika.utils import is_file_and_not_empty, make_window_dirs, strip_prmtop


@pytest.fixture
def clean_files(directory=os.path.join(os.path.dirname(__file__), "tmp")):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


def test_mkdirs():
    """Test that we can make directories for the windows."""
    window_list = ["a000", "p014", "r089"]
    make_window_dirs(window_list, path="tmp")
    for window in window_list:
        assert os.path.exists(os.path.join("tmp", "windows", window))


def test_strip_prmtop():
    """Test that we can remove items from structures."""
    cb6_only = strip_prmtop(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.prmtop"),
        mask=":BUT",
    )
    residues = [r.name for r in cb6_only.residues]
    assert "CB6" in residues
    assert "BUT" not in residues


def test_is_file_and_not_empty(clean_files):
    # Check for existing file
    but_frcmod = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.frcmod")
    )

    assert os.path.isfile(but_frcmod)
    assert os.path.getsize(but_frcmod) != 0
    assert is_file_and_not_empty(but_frcmod)

    # Check for emtpy
    path = os.path.join(os.path.dirname(__file__), "tmp", "testfile1")
    with open(path, "w") as f:
        pass

    assert os.path.isfile(path)
    assert os.path.getsize(path) == 0
    assert not is_file_and_not_empty(path)

    # Check file after writing
    path = os.path.join(os.path.dirname(__file__), "tmp", "testfile2")
    with open(path, "w") as f:
        f.write("Hello World\n")

    assert os.path.isfile(path)
    assert os.path.getsize(path) != 0
    assert is_file_and_not_empty(path)
