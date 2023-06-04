"""
Tests the alignment of residues to the z axis.
"""

import os

import numpy
import parmed
import pytest
from openff.units import unit as openff_unit

from paprika.build.align import (
    align_principal_axes,
    get_centroid,
    get_principal_axis_vector,
    get_theta,
    rotate_around_axis,
    shift_structure,
    translate_to_origin,
    zalign,
)


def test_center_mask():
    """Test that the first mask is centered."""
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT", weight="mass")
    test_coordinates = get_centroid(aligned_cb6, ":CB6", weight="mass")
    assert numpy.allclose(test_coordinates, numpy.zeros(3))
    test_coordinates = get_centroid(aligned_cb6, ":CB6", weight="geo")
    assert numpy.allclose(
        test_coordinates, numpy.array([-0.0002163, 0.00113288, -0.00072443])
    )


def test_alignment_after_offset():
    """Test that molecule is properly aligned after random offset."""
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    random_coordinates = numpy.random.randint(10) * numpy.random.rand(1, 3)
    cb6_offset = shift_structure(cb6, random_coordinates)
    aligned_cb6 = zalign(cb6_offset, ":CB6", ":BUT")
    test_coordinates = get_centroid(aligned_cb6, ":CB6", weight="mass")
    assert numpy.allclose(test_coordinates, numpy.zeros(3))
    test_coordinates = get_centroid(aligned_cb6, ":CB6", weight="geo")
    assert numpy.allclose(
        test_coordinates, numpy.array([-0.0002163, 0.00113288, -0.00072443])
    )


def test_theta_after_alignment():
    """Test that molecule is properly aligned after random offset."""
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    assert (
        get_theta(aligned_cb6, ":CB6", ":BUT", axis="z")
        .to(openff_unit.radians)
        .magnitude
        == 0
    )
    assert (
        pytest.approx(
            get_theta(aligned_cb6, ":CB6", ":BUT", axis="x")
            .to(openff_unit.radians)
            .magnitude,
            abs=1e-3,
        )
        == 1.5708
    )
    assert (
        pytest.approx(
            get_theta(aligned_cb6, ":CB6", ":BUT", axis="y")
            .to(openff_unit.radians)
            .magnitude,
            abs=1e-3,
        )
        == 1.5708
    )


def test_translate_to_origin():
    """Test that molecule is properly aligned after translated to the origin."""
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb"),
        structure=True,
    )

    # Translate molecule to origin
    translated_cb6 = translate_to_origin(cb6)
    coordinates = translated_cb6.coordinates
    masses = numpy.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = parmed.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) == 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0

    # Shift then translate only in the z-axis
    cb6_offset = shift_structure(cb6, numpy.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension="z")
    coordinates = translated_cb6.coordinates
    masses = numpy.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = parmed.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) != 0.0
    assert pytest.approx(centroid[1], abs=1e-3) != 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0

    # Randomly shift then translate only in the x- and z-axis
    cb6_offset = shift_structure(cb6, numpy.array([3, 5, 10]))
    translated_cb6 = translate_to_origin(cb6_offset, dimension=[1, 0, 1])
    coordinates = translated_cb6.coordinates
    masses = numpy.asarray([atom.mass for atom in translated_cb6.atoms])
    centroid = parmed.geometry.center_of_mass(coordinates, masses)

    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) != 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 0.0


def test_get_principal_axis():
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    principal_axis = get_principal_axis_vector(cb6, principal_axis=1)

    assert pytest.approx(principal_axis[0], abs=1e-1) == 0.0
    assert pytest.approx(principal_axis[1], abs=1e-1) == 0.0
    assert pytest.approx(principal_axis[2], abs=1e-1) == 1.0


def test_align_principal_axes():
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_aligned = align_principal_axes(
        cb6, atom_mask=":CB6", principal_axis=1, axis="y"
    )

    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="z")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 90.0

    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="y")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 0.0

    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="x")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 90.0


def test_rotate_around_axis():
    # Cartesian axes
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_aligned = rotate_around_axis(cb6, axis="z", angle=90.0 * openff_unit.degrees)
    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="z")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 0.0

    cb6_aligned = rotate_around_axis(
        cb6_aligned, axis="x", angle=90.0 * openff_unit.degrees
    )
    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="x")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 90.0

    cb6_aligned = (
        rotate_around_axis(cb6_aligned, axis="y", angle=90.0) * openff_unit.degrees
    )
    angle = (
        get_theta(cb6_aligned, ":BUT@C", ":BUT@C3", axis="y")
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 180.0

    # Arbitrary axes
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    cb6_axis = rotate_around_axis(cb6, axis=[0, 0, 1], angle=90.0 * openff_unit.degrees)
    angle = (
        get_theta(cb6_axis, ":BUT@C", ":BUT@C3", axis=[0, 0, 1])
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 0.0

    cb6_axis = rotate_around_axis(cb6, axis=[1, 1, 1], angle=90.0 * openff_unit.degrees)
    angle = (
        get_theta(cb6_axis, ":BUT@C", ":BUT@C3", axis=[1, 1, 1])
        .to(openff_unit.degrees)
        .magnitude
    )
    assert pytest.approx(angle, abs=1e-1) == 54.7


def test_get_centroid():
    cb6 = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    centroid = get_centroid(cb6, atom_mask=":BUT", weight="mass")
    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) == 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 1.918

    centroid = get_centroid(cb6, atom_mask=":BUT", weight="geo")
    assert pytest.approx(centroid[0], abs=1e-3) == 0.0
    assert pytest.approx(centroid[1], abs=1e-3) == 0.0
    assert pytest.approx(centroid[2], abs=1e-3) == 1.918
