"""
Tests the alignment of residues to the z axis.
"""

import os

import numpy as np
import parmed as pmd
import pytest

from paprika.align import (
    check_coordinates,
    get_theta,
    offset_structure,
    translate_to_origin,
    zalign,
)


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
    assert (
        pytest.approx(get_theta(aligned_cb6, ":CB6", ":BUT", axis="x"), 0.001) == 1.5708
    )
    assert (
        pytest.approx(get_theta(aligned_cb6, ":CB6", ":BUT", axis="y"), 0.001) == 1.5708
    )


def test_translate_to_origin():
    """ Test that molecule is properly aligned after translated to the origin."""
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb"),
        structure=True,
    )

    # Translate molecule to origin
    translated_cb6 = translate_to_origin(cb6)
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], 0.001) == 0.0
    assert pytest.approx(centroid[1], 0.001) == 0.0
    assert pytest.approx(centroid[2], 0.001) == 0.0

    # Shift then translate only in the z-axis
    cb6_offset = offset_structure(cb6, np.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension="z")
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], 0.001) != 0.0
    assert pytest.approx(centroid[1], 0.001) != 0.0
    assert pytest.approx(centroid[2], 0.001) == 0.0

    # Randomly shift then translate only in the x- and z-axis
    cb6_offset = offset_structure(cb6, np.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension=[1, 0, 1])
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], 0.001) == 0.0
    assert pytest.approx(centroid[1], 0.001) != 0.0
    assert pytest.approx(centroid[2], 0.001) == 0.0
