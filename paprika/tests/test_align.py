"""
Tests the alignment of residues to the z axis.
"""

import os

import numpy as np
import parmed as pmd
import pytest

from paprika.build.align import (
    align_principal_axes,
    check_coordinates,
    get_principal_axis_vector,
    get_theta,
    offset_structure,
    rotate_around_axis,
    translate_to_origin,
    zalign,
)


def test_center_mask():
    """Test that the first mask is centered"""
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    test_coordinates = check_coordinates(aligned_cb6, ":CB6")
    assert np.allclose(test_coordinates, np.zeros(3))


def test_alignment_after_offset():
    """Test that molecule is properly aligned after random offset."""
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    random_coordinates = np.random.randint(10) * np.random.rand(1, 3)
    cb6_offset = offset_structure(cb6, random_coordinates)
    aligned_cb6 = zalign(cb6_offset, ":CB6", ":BUT")
    test_coordinates = check_coordinates(aligned_cb6, ":CB6")
    assert np.allclose(test_coordinates, np.zeros(3))


def test_theta_after_alignment():
    """Test that molecule is properly aligned after random offset."""
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    assert get_theta(aligned_cb6, ":CB6", ":BUT", axis="z") == 0
    assert (
        pytest.approx(get_theta(aligned_cb6, ":CB6", ":BUT", axis="x"), abs=1e-3)
        == 1.5708
    )
    assert (
        pytest.approx(get_theta(aligned_cb6, ":CB6", ":BUT", axis="y"), abs=1e-3)
        == 1.5708
    )


def test_translate_to_origin():
    """Test that molecule is properly aligned after translated to the origin."""
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb"),
        structure=True,
    )

    # Translate molecule to origin
    translated_cb6 = translate_to_origin(cb6)
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) == 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0

    # Shift then translate only in the z-axis
    cb6_offset = offset_structure(cb6, np.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension="z")
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) != 0.0
    assert pytest.approx(centroid[1], abs=1e-3) != 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0

    # Randomly shift then translate only in the x- and z-axis
    cb6_offset = offset_structure(cb6, np.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension=[1, 0, 1])
    coordinates = translated_cb6.coordinates
    masses = np.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) != 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0


def test_get_principal_axis():
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    principal_axis = get_principal_axis_vector(cb6, principal_axis=1)

    assert pytest.approx(principal_axis[0], abs=1e-1) == 0.0
    assert pytest.approx(principal_axis[1], abs=1e-1) == 0.0
    assert pytest.approx(principal_axis[2], abs=1e-1) == 1.0


def test_align_principal_axes():
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_aligned = align_principal_axes(
        cb6, atom_mask=":CB6", principal_axis=1, axis="y"
    )

    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="z") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 90.0

    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="y") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 0.0

    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="x") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 90.0


def test_rotate_around_axis():
    # Cartesian axes
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_aligned = rotate_around_axis(cb6, axis="z", angle=90.0)
    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="z") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 0.0

    cb6_aligned = rotate_around_axis(cb6_aligned, axis="x", angle=90.0)
    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="x") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 90.0

    cb6_aligned = rotate_around_axis(cb6_aligned, axis="y", angle=90.0)
    angle = get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="y") * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 180.0

    # Arbitrary axes
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_axis = rotate_around_axis(cb6, axis=[0, 0, 1], angle=90.0)
    angle = get_theta(cb6_axis, ":BUT@C", ":BUT@C3", axis=[0, 0, 1]) * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 0.0

    cb6_axis = rotate_around_axis(cb6, axis=[1, 1, 1], angle=90.0)
    angle = get_theta(cb6_axis, ":BUT@C", ":BUT@C3", axis=[1, 1, 1]) * 180 / np.pi
    assert pytest.approx(angle, abs=1e-1) == 54.7


def test_check_coordinates():
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    com = check_coordinates(cb6, mask=":BUT")
    assert pytest.approx(com[0], abs=1e-3) == 0.0
    assert pytest.approx(com[1], abs=1e-3) == 0.0
    assert pytest.approx(com[2], abs=1e-1) == 1.9
