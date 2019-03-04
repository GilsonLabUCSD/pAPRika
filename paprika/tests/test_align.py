"""
Tests the alignment of residues to the z axis.
"""

import os

from paprika.align import *


def test_center_mask():
    """ Test that the first mask is centered """
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    test_coordinates = check_coordinates(aligned_cb6, ":CB6")
    assert np.allclose(test_coordinates, np.zeros(3))


def test_alignment_after_offset():
    """ Test that molecule is properly aligned after random offset. """
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    random_coordinates = np.random.randint(10) * np.random.rand(1, 3)
    cb6_offset = offset_structure(cb6, random_coordinates)
    aligned_cb6 = zalign(cb6_offset, ":CB6", ":BUT")
    test_coordinates = check_coordinates(aligned_cb6, ":CB6")
    assert np.allclose(test_coordinates, np.zeros(3))


def test_theta_after_alignment():
    """ Test that molecule is properly aligned after random offset. """
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    assert get_theta(aligned_cb6, ":CB6", ":BUT", axis="z") == 0
