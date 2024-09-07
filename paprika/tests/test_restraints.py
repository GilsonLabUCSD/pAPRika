"""
Tests the restraints utilities.
"""

import logging
import os

import numpy

try:
    import openmm
    import openmm.app as app
    import openmm.unit as openmm_unit
except ImportError:
    import simtk.openmm as openmm
    import simtk.openmm.app as app
    import simtk.unit as openmm_unit

import parmed
import pytest
from openff.units import unit as openff_unit

from paprika.restraints import create_window_list
from paprika.restraints.openmm import apply_dat_restraint, apply_positional_restraints
from paprika.restraints.restraints import (
    BiasPotentialType,
    DAT_restraint,
    RestraintType,
)
from paprika.restraints.utils import (
    extract_guest_restraints,
    get_bias_potential_type,
    get_restraint_values,
    restraints_from_ascii,
)
from paprika.tests.test_tleap import clean_files

logger = logging.getLogger(__name__)


def test_DAT_restraint():
    # Method 1
    logger.info("### Testing restraint 1, Method 1")
    rest1 = DAT_restraint()
    rest1.amber_index = True
    rest1.continuous_apr = False
    rest1.auto_apr = False
    rest1.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest1.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest1.mask2 = ":BUT@C3"
    rest1.attach["target"] = 3.0
    rest1.attach["num_windows"] = 4
    rest1.attach["fc_initial"] = 0.0
    rest1.attach["fc_final"] = 3.0
    rest1.pull["fc"] = rest1.attach["fc_final"]
    rest1.pull["num_windows"] = 4
    rest1.pull["target_initial"] = rest1.attach["target"]
    rest1.pull["target_final"] = 6.0
    rest1.release["target"] = rest1.pull["target_final"]
    rest1.release["num_windows"] = rest1.attach["num_windows"]
    rest1.release["fc_initial"] = rest1.attach["fc_initial"]
    rest1.release["fc_final"] = rest1.attach["fc_final"]
    rest1.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / target_units**2
    assert rest1.restraint_type == RestraintType.Distance
    assert rest1.index1 == [13, 31, 49, 67, 85, 103]
    assert rest1.index2 == [119]
    assert rest1.index3 is None
    assert rest1.index4 is None
    assert numpy.allclose(
        rest1.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert numpy.allclose(
        rest1.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([3.0, 3.0, 3.0, 3.0]),
    )
    assert numpy.allclose(
        rest1.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([3.0, 3.0, 3.0, 3.0]),
    )
    assert numpy.allclose(
        rest1.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([3.0, 4.0, 5.0, 6.0]),
    )
    assert numpy.allclose(
        rest1.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert numpy.allclose(
        rest1.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([6.0, 6.0, 6.0, 6.0]),
    )
    window_list = create_window_list([rest1])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 1a
    logger.info("### Testing restraint 2, Method 1a")
    rest2 = DAT_restraint()
    rest2.amber_index = True
    rest2.continuous_apr = False
    rest2.auto_apr = False
    rest2.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest2.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest2.mask2 = ":BUT@C3"
    rest2.mask3 = ":BUT@C"
    rest2.attach["target"] = 180.0
    rest2.attach["num_windows"] = 4
    rest2.attach["fc_final"] = 75.0
    rest2.pull["fc"] = rest2.attach["fc_final"]
    rest2.pull["num_windows"] = 4
    rest2.pull["target_final"] = 180.0
    rest2.release["target"] = rest2.pull["target_final"]
    rest2.release["num_windows"] = rest2.attach["num_windows"]
    rest2.release["fc_final"] = rest2.attach["fc_final"]
    rest2.initialize()

    target_units = openff_unit.degrees
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.radians**2
    assert rest2.restraint_type == RestraintType.Angle
    assert rest2.index1 == [13, 31, 49, 67, 85, 103]
    assert rest2.index2 == [119]
    assert rest2.index3 == [109]
    assert rest2.index4 is None
    assert numpy.allclose(
        rest2.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest2.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([180.0, 180.0, 180.0, 180.0]),
    )
    assert numpy.allclose(
        rest2.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([75.0, 75.0, 75.0, 75.0]),
    )
    assert numpy.allclose(
        rest2.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 60.0, 120.0, 180.0]),
    )
    assert numpy.allclose(
        rest2.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest2.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([180.0, 180.0, 180.0, 180.0]),
    )
    window_list = create_window_list([rest2])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 2 (Note auto_apr = True)
    logger.info("### Testing restraint 3, Method 2")
    rest3 = DAT_restraint()
    rest3.amber_index = True
    rest3.continuous_apr = False
    rest3.auto_apr = True
    rest3.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest3.mask1 = ":CB6@O2"
    rest3.mask2 = ":CB6@O"
    rest3.mask3 = ":BUT@C3"
    rest3.mask4 = ":BUT@C"
    rest3.attach["target"] = 90.0
    rest3.attach["fc_increment"] = 25.0
    rest3.attach["fc_initial"] = 0.0
    rest3.attach["fc_final"] = 75.0
    rest3.pull["target_increment"] = 1.0
    rest3.pull["target_final"] = 93.0
    rest3.release["fc_final"] = 75.0
    rest3.initialize()

    target_units = openff_unit.degrees
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.radians**2
    assert rest3.restraint_type == RestraintType.Torsion
    assert rest3.index1 == [31]
    assert rest3.index2 == [13]
    assert rest3.index3 == [119]
    assert rest3.index4 == [109]
    assert numpy.allclose(
        rest3.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest3.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([90.0, 90.0, 90.0, 90.0]),
    )
    assert numpy.allclose(
        rest3.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([75.0, 75.0, 75.0, 75.0]),
    )
    assert numpy.allclose(
        rest3.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([90.0, 91.0, 92.0, 93.0]),
    )
    assert numpy.allclose(
        rest3.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest3.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([93.0, 93.0, 93.0, 93.0]),
    )
    window_list = create_window_list([rest3])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 2a
    logger.info("### Testing restraint 4, Method 2a")
    rest4 = DAT_restraint()
    rest4.amber_index = True
    rest4.continuous_apr = False
    rest4.auto_apr = False
    rest4.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest4.mask1 = ":CB6@O2"
    rest4.mask2 = ":CB6@O"
    rest4.mask3 = ":BUT@C3"
    rest4.mask4 = ":BUT@C"
    rest4.attach["target"] = 0.0
    rest4.attach["fc_increment"] = 25.0
    rest4.attach["fc_final"] = 75.0
    rest4.pull["fc"] = 75.0
    rest4.pull["target_increment"] = 1.0
    rest4.pull["target_final"] = 3.0
    rest4.release["target"] = 3.0
    rest4.release["fc_increment"] = 25.0
    rest4.release["fc_final"] = 75.0
    rest4.initialize()

    target_units = openff_unit.degrees
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.radians**2
    assert rest4.index1 == [31]
    assert rest4.index2 == [13]
    assert rest4.index3 == [119]
    assert rest4.index4 == [109]
    assert numpy.allclose(
        rest4.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest4.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0, 0.0]),
    )
    assert numpy.allclose(
        rest4.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([75.0, 75.0, 75.0, 75.0]),
    )
    assert numpy.allclose(
        rest4.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert numpy.allclose(
        rest4.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 25.0, 50.0, 75.0]),
    )
    assert numpy.allclose(
        rest4.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([3.0, 3.0, 3.0, 3.0]),
    )
    window_list = create_window_list([rest4])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 3
    logger.info("### Testing restraint 5, Method 3")
    rest5 = DAT_restraint()
    rest5.amber_index = True
    rest5.continuous_apr = False
    rest5.auto_apr = False
    rest5.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest5.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest5.mask2 = ":BUT@C*"
    rest5.attach["target"] = 0.0
    rest5.attach["fraction_list"] = [0.0, 0.2, 0.5, 1.0]
    rest5.attach["fc_final"] = 5.0
    rest5.pull["fc"] = rest5.attach["fc_final"]
    rest5.pull["fraction_list"] = [0.0, 0.5, 1.0]
    rest5.pull["target_final"] = 1.0
    rest5.release["target"] = rest5.pull["target_final"]
    rest5.release["fraction_list"] = [0.0, 0.3, 0.6, 1.0]
    rest5.release["fc_final"] = rest5.attach["fc_final"]
    rest5.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest5.index1 == [13, 31, 49, 67, 85, 103]
    assert rest5.index2 == [109, 113, 115, 119]
    assert rest5.index3 is None
    assert rest5.index4 is None
    assert numpy.allclose(
        rest5.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.0, 2.5, 5.0]),
    )
    assert numpy.allclose(
        rest5.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0, 0.0]),
    )
    assert numpy.allclose(
        rest5.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([5.0, 5.0, 5.0]),
    )
    assert numpy.allclose(
        rest5.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.5, 1.0]),
    )
    assert numpy.allclose(
        rest5.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.5, 3.0, 5.0]),
    )
    assert numpy.allclose(
        rest5.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([1.0, 1.0, 1.0, 1.0]),
    )
    window_list = create_window_list([rest5])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 4
    logger.info("### Testing restraint 6, Method 4")
    rest6 = DAT_restraint()
    rest6.amber_index = True
    rest6.continuous_apr = False
    rest6.auto_apr = False
    rest6.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest6.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest6.mask2 = ":BUT@C*"
    rest6.attach["target"] = 0.0
    rest6.attach["fraction_increment"] = 0.25
    rest6.attach["fc_final"] = 5.0
    rest6.pull["fc"] = rest6.attach["fc_final"]
    rest6.pull["fraction_increment"] = 0.5
    rest6.pull["target_final"] = 1.0
    rest6.release["target"] = rest6.pull["target_final"]
    rest6.release["fraction_increment"] = 0.33
    rest6.release["fc_final"] = rest6.attach["fc_final"]
    rest6.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest6.index1 == [13, 31, 49, 67, 85, 103]
    assert rest6.index2 == [109, 113, 115, 119]
    assert rest6.index3 is None
    assert rest6.index4 is None
    assert numpy.allclose(
        rest6.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.25, 2.5, 3.75, 5.0]),
    )
    assert numpy.allclose(
        rest6.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0, 0.0, 0.0]),
    )
    assert numpy.allclose(
        rest6.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([5.0, 5.0, 5.0]),
    )
    assert numpy.allclose(
        rest6.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.5, 1.0]),
    )
    assert numpy.allclose(
        rest6.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.65, 3.3, 5.0]),
    )
    assert numpy.allclose(
        rest6.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([1.0, 1.0, 1.0, 1.0]),
    )
    window_list = create_window_list([rest6])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "a004",
        "p000",
        "p001",
        "p002",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 5 (Note continuous_apr = True)
    logger.info("### Testing restraint 7, Method 5")
    rest7 = DAT_restraint()
    rest7.amber_index = True
    rest7.continuous_apr = True
    rest7.auto_apr = False
    rest7.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest7.mask1 = ":1@O,O1,:BUT@H1"
    rest7.mask2 = ":CB6@N"
    rest7.attach["target"] = 0.0
    rest7.attach["fc_list"] = [0.0, 0.5, 1.0, 2.0]
    rest7.pull["fc"] = 2.0
    rest7.pull["target_list"] = [0.0, 0.5, 1.0, 1.5]
    rest7.release["target"] = 1.5
    rest7.release["fc_list"] = [0.0, 0.66, 1.2, 2.0]
    rest7.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest7.index1 == [13, 14, 111]
    assert rest7.index2 == [3]
    assert rest7.index3 is None
    assert rest7.index4 is None
    assert numpy.allclose(
        rest7.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 0.5, 1.0, 2.0]),
    )
    assert numpy.allclose(
        rest7.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0, 0.0]),
    )
    assert numpy.allclose(
        rest7.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([2.0, 2.0, 2.0, 2.0]),
    )
    assert numpy.allclose(
        rest7.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.5, 1.0, 1.5]),
    )
    assert numpy.allclose(
        rest7.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 0.66, 1.2, 2.0]),
    )
    assert numpy.allclose(
        rest7.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([1.5, 1.5, 1.5, 1.5]),
    )
    window_list = create_window_list([rest7])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "p000",
        "p001",
        "p002",
        "p003",
        "r001",
        "r002",
        "r003",
    ]

    # Just Attach
    logger.info("### Testing restraint 8, just attach")
    rest8 = DAT_restraint()
    rest8.amber_index = True
    rest8.continuous_apr = False
    rest8.auto_apr = False
    rest8.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest8.mask1 = ":CB6@O"
    rest8.mask2 = ":BUT@C3"
    rest8.attach["target"] = 0.0
    rest8.attach["num_windows"] = 4
    rest8.attach["fc_initial"] = 0.0
    rest8.attach["fc_final"] = 3.0
    rest8.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest8.index1 == [13]
    assert rest8.index2 == [119]
    assert rest8.index3 is None
    assert rest8.index4 is None
    assert numpy.allclose(
        rest8.phase["attach"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert numpy.allclose(
        rest8.phase["attach"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0, 0.0]),
    )
    assert rest8.phase["pull"]["force_constants"] is None
    assert rest8.phase["pull"]["targets"] is None
    assert rest8.phase["release"]["force_constants"] is None
    assert rest8.phase["release"]["targets"] is None
    window_list = create_window_list([rest8])
    assert window_list == ["a000", "a001", "a002", "a003"]

    # Just Pull
    logger.info("### Testing restraint 9, just pull")
    rest9 = DAT_restraint()
    rest9.amber_index = True
    rest9.continuous_apr = False
    rest9.auto_apr = False
    rest9.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest9.mask1 = ":CB6@O"
    rest9.mask2 = ":BUT@C3"
    rest9.pull["fc"] = 3.0
    rest9.pull["num_windows"] = 4
    rest9.pull["target_initial"] = 0.0
    rest9.pull["target_final"] = 3.0
    rest9.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest9.index1 == [13]
    assert rest9.index2 == [119]
    assert rest9.index3 is None
    assert rest9.index4 is None
    assert rest9.phase["attach"]["force_constants"] is None
    assert rest9.phase["attach"]["targets"] is None
    assert numpy.allclose(
        rest9.phase["pull"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([3.0, 3.0, 3.0, 3.0]),
    )
    assert numpy.allclose(
        rest9.phase["pull"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert rest9.phase["release"]["force_constants"] is None
    assert rest9.phase["release"]["targets"] is None
    window_list = create_window_list([rest9])
    assert window_list == ["p000", "p001", "p002", "p003"]

    # Just Release
    logger.info("### Testing restraint 10, just release")
    rest10 = DAT_restraint()
    rest10.amber_index = True
    rest10.continuous_apr = False
    rest10.auto_apr = False
    rest10.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest10.mask1 = ":CB6@O"
    rest10.mask2 = ":BUT@C3"
    rest10.release["target"] = 0.0
    rest10.release["num_windows"] = 3
    rest10.release["fc_initial"] = 0.0
    rest10.release["fc_final"] = 2.0
    rest10.initialize()

    target_units = openff_unit.angstrom
    force_constant_units = openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    assert rest10.index1 == [13]
    assert rest10.index2 == [119]
    assert rest10.index3 is None
    assert rest10.index4 is None
    assert rest10.phase["attach"]["force_constants"] is None
    assert rest10.phase["attach"]["targets"] is None
    assert rest10.phase["pull"]["force_constants"] is None
    assert rest10.phase["pull"]["targets"] is None
    assert numpy.allclose(
        rest10.phase["release"]["force_constants"].to(force_constant_units).magnitude,
        numpy.array([0.0, 1.0, 2.0]),
    )
    assert numpy.allclose(
        rest10.phase["release"]["targets"].to(target_units).magnitude,
        numpy.array([0.0, 0.0, 0.0]),
    )
    window_list = create_window_list([rest10])
    assert window_list == ["r000", "r001", "r002"]

    # Test inconsistent continuous_apr:
    with pytest.raises(Exception):
        window_list = create_window_list([rest7, rest8])

    # Test inconsistent windows:
    with pytest.raises(Exception):
        window_list = create_window_list([rest1, rest10])


def test_get_restraint_values():
    # Test Harmonic restraint
    attach_fractions = numpy.linspace(0, 1.0, 25)
    initial_distance = 2.65
    pull_distances = numpy.linspace(0 + initial_distance, 16.0 + initial_distance, 40)

    restraint = DAT_restraint()
    restraint.continuous_apr = True
    restraint.amber_index = True
    restraint.topology = os.path.join(
        os.path.dirname(__file__), "../data/k-cl/k-cl.pdb"
    )
    restraint.mask1 = "@K+"
    restraint.mask2 = "@Cl-"

    restraint.attach["target"] = initial_distance
    restraint.attach["fraction_list"] = attach_fractions
    restraint.attach["fc_final"] = 10.0

    restraint.pull["fc"] = restraint.attach["fc_final"]
    restraint.pull["target_list"] = pull_distances

    restraint.initialize()

    restraint_values = get_restraint_values(restraint, "attach", 0)
    assert restraint_values["r1"].magnitude == 0.0
    assert restraint_values["r2"].magnitude == 2.65
    assert restraint_values["r3"].magnitude == 2.65
    assert restraint_values["r4"].magnitude == 999.0
    assert restraint_values["rk2"].magnitude == 0.0
    assert restraint_values["rk3"].magnitude == 0.0

    restraint_values = get_restraint_values(restraint, "pull", 0)
    assert restraint_values["r1"].magnitude == 0.0
    assert restraint_values["r2"].magnitude == 2.65
    assert restraint_values["r3"].magnitude == 2.65
    assert restraint_values["r4"].magnitude == 999.0
    assert restraint_values["rk2"].magnitude == 10.0
    assert restraint_values["rk3"].magnitude == 10.0

    # Test custom values
    wall = DAT_restraint()
    wall.auto_apr = False
    wall.amber_index = True
    wall.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    wall.mask1 = "@K+"
    wall.mask2 = "@Cl-"

    wall.attach["fc_initial"] = 1.0
    wall.attach["fc_final"] = 1.0

    wall.custom_restraint_values["rk2"] = 1.0
    wall.custom_restraint_values["rk3"] = 1.0
    wall.custom_restraint_values["r2"] = 0.0
    wall.custom_restraint_values["r3"] = 3.5

    wall.attach["target"] = 3.5
    wall.attach["num_windows"] = len(attach_fractions)

    wall.initialize()

    restraint_values = get_restraint_values(wall, "attach", 0)
    assert restraint_values["r1"].magnitude == 0.0
    assert restraint_values["r2"].magnitude == 0.0
    assert restraint_values["r3"].magnitude == 3.5
    assert restraint_values["r4"].magnitude == 999.0
    assert restraint_values["rk2"].magnitude == 1.0
    assert restraint_values["rk3"].magnitude == 1.0

    wall = DAT_restraint()
    wall.auto_apr = False
    wall.amber_index = True
    wall.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    wall.mask1 = "@K+"
    wall.mask2 = "@Cl-"

    wall.attach["fc_initial"] = 1.0
    wall.attach["fc_final"] = 1.0

    wall.custom_restraint_values["rk2"] = 1.0
    wall.custom_restraint_values["rk3"] = 1.0
    wall.custom_restraint_values["r2"] = 3.5
    wall.custom_restraint_values["r3"] = 0.0

    wall.attach["target"] = 3.5
    wall.attach["num_windows"] = len(attach_fractions)

    wall.initialize()

    restraint_values = get_restraint_values(wall, "attach", 0)
    assert restraint_values["r1"].magnitude == 0.0
    assert restraint_values["r2"].magnitude == 3.5
    assert restraint_values["r3"].magnitude == 0.0
    assert restraint_values["r4"].magnitude == 999.0
    assert restraint_values["rk2"].magnitude == 1.0
    assert restraint_values["rk3"].magnitude == 1.0

    restraint_values = get_restraint_values(wall, "pull", 0)
    assert restraint_values is None


def test_get_bias_potential_type():
    # Test Harmonic restraint
    attach_fractions = numpy.linspace(0, 1.0, 25)
    initial_distance = 2.65
    pull_distances = numpy.linspace(0 + initial_distance, 16.0 + initial_distance, 40)

    restraint = DAT_restraint()
    restraint.continuous_apr = True
    restraint.amber_index = True
    restraint.topology = os.path.join(
        os.path.dirname(__file__), "../data/k-cl/k-cl.pdb"
    )
    restraint.mask1 = "@K+"
    restraint.mask2 = "@Cl-"

    restraint.attach["target"] = initial_distance
    restraint.attach["fraction_list"] = attach_fractions
    restraint.attach["fc_final"] = 10.0

    restraint.pull["fc"] = restraint.attach["fc_final"]
    restraint.pull["target_list"] = pull_distances

    restraint.initialize()

    assert get_bias_potential_type(restraint, "attach", 0) == BiasPotentialType.Harmonic
    assert get_bias_potential_type(restraint, "pull", 0) == BiasPotentialType.Harmonic
    assert get_bias_potential_type(restraint, "release", 0) is None

    # Test upper wall restraint (1)
    upper = DAT_restraint()
    upper.auto_apr = False
    upper.amber_index = True
    upper.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    upper.mask1 = "@K+"
    upper.mask2 = "@Cl-"

    upper.attach["fc_initial"] = 1.0
    upper.attach["fc_final"] = 1.0

    upper.custom_restraint_values["rk2"] = 1.0
    upper.custom_restraint_values["rk3"] = 1.0
    upper.custom_restraint_values["r2"] = 0.0
    upper.custom_restraint_values["r3"] = 3.5

    upper.attach["target"] = 3.5
    upper.attach["num_windows"] = len(attach_fractions)

    upper.initialize()

    assert get_bias_potential_type(upper, "attach", 0) == BiasPotentialType.UpperWall
    assert get_bias_potential_type(upper, "attach", 1) == BiasPotentialType.UpperWall

    # Test upper wall restraint (2)
    upper = DAT_restraint()
    upper.auto_apr = False
    upper.amber_index = True
    upper.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    upper.mask1 = "@K+"
    upper.mask2 = "@Cl-"

    upper.attach["fc_initial"] = 1.0
    upper.attach["fc_final"] = 1.0

    upper.custom_restraint_values["rk2"] = 0.0
    upper.custom_restraint_values["rk3"] = 1.0
    upper.custom_restraint_values["r2"] = 3.5
    upper.custom_restraint_values["r3"] = 3.5

    upper.attach["target"] = 3.5
    upper.attach["num_windows"] = len(attach_fractions)

    upper.initialize()

    assert get_bias_potential_type(upper, "attach", 0) == BiasPotentialType.UpperWall
    assert get_bias_potential_type(upper, "attach", 1) == BiasPotentialType.UpperWall

    # Test upper wall restraint (3)
    upper = DAT_restraint()
    upper.auto_apr = False
    upper.amber_index = True
    upper.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    upper.mask1 = "@K+"
    upper.mask2 = "@Cl-"

    upper.attach["fc_initial"] = 1.0
    upper.attach["fc_final"] = 1.0

    upper.custom_restraint_values["rk2"] = 0.0
    upper.custom_restraint_values["rk3"] = 1.0
    upper.custom_restraint_values["r2"] = 0.0
    upper.custom_restraint_values["r3"] = 3.5

    upper.attach["target"] = 3.5
    upper.attach["num_windows"] = len(attach_fractions)

    upper.initialize()

    assert get_bias_potential_type(upper, "attach", 0) == BiasPotentialType.UpperWall
    assert get_bias_potential_type(upper, "attach", 1) == BiasPotentialType.UpperWall

    # Test lower wall restraint (1)
    lower = DAT_restraint()
    lower.auto_apr = False
    lower.amber_index = True
    lower.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    lower.mask1 = "@K+"
    lower.mask2 = "@Cl-"

    lower.attach["fc_initial"] = 1.0
    lower.attach["fc_final"] = 1.0

    lower.custom_restraint_values["rk2"] = 1.0
    lower.custom_restraint_values["rk3"] = 1.0
    lower.custom_restraint_values["r2"] = 3.5
    lower.custom_restraint_values["r3"] = 0.0

    lower.attach["target"] = 3.5
    lower.attach["num_windows"] = len(attach_fractions)

    lower.initialize()

    assert get_bias_potential_type(lower, "attach", 0) == BiasPotentialType.LowerWall
    assert get_bias_potential_type(lower, "attach", 1) == BiasPotentialType.LowerWall

    # Test lower wall restraint (2)
    lower = DAT_restraint()
    lower.auto_apr = False
    lower.amber_index = True
    lower.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    lower.mask1 = "@K+"
    lower.mask2 = "@Cl-"

    lower.attach["fc_initial"] = 1.0
    lower.attach["fc_final"] = 1.0

    lower.custom_restraint_values["rk2"] = 1.0
    lower.custom_restraint_values["rk3"] = 0.0
    lower.custom_restraint_values["r2"] = 3.5
    lower.custom_restraint_values["r3"] = 3.5

    lower.attach["target"] = 3.5
    lower.attach["num_windows"] = len(attach_fractions)

    lower.initialize()

    assert get_bias_potential_type(lower, "attach", 0) == BiasPotentialType.LowerWall
    assert get_bias_potential_type(lower, "attach", 1) == BiasPotentialType.LowerWall

    # Test lower wall restraint (3)
    lower = DAT_restraint()
    lower.auto_apr = False
    lower.amber_index = True
    lower.topology = os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb")
    lower.mask1 = "@K+"
    lower.mask2 = "@Cl-"

    lower.attach["fc_initial"] = 1.0
    lower.attach["fc_final"] = 1.0

    lower.custom_restraint_values["rk2"] = 1.0
    lower.custom_restraint_values["rk3"] = 0.0
    lower.custom_restraint_values["r2"] = 3.5
    lower.custom_restraint_values["r3"] = 6.5

    lower.attach["target"] = 3.5
    lower.attach["num_windows"] = len(attach_fractions)

    lower.initialize()

    assert get_bias_potential_type(lower, "attach", 0) == BiasPotentialType.LowerWall
    assert get_bias_potential_type(lower, "attach", 1) == BiasPotentialType.LowerWall


def test_extract_guest_restraints():
    """Test extract_guest_restraints to make sure it is extracting the correct restraints."""
    restraints = []

    # Guest - r
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM1"
    r.mask2 = ":BUT@C3"
    r.attach["target"] = 6.0
    r.attach["num_windows"] = 4
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 5.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 24.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    # Guest - theta
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM1"
    r.mask2 = ":BUT@C3"
    r.mask3 = ":BUT@C"
    r.attach["target"] = 180.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 180.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    # Guest - beta
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM2"
    r.mask2 = ":DM1"
    r.mask3 = ":BUT@C3"
    r.attach["target"] = 180.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 180.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    structure = parmed.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        structure=True,
    )
    guest_restraints = extract_guest_restraints(
        structure, restraints, "BUT", return_type="list"
    )

    assert guest_restraints[0] is not None
    assert guest_restraints[1] is not None
    assert guest_restraints[2] is None
    assert guest_restraints[3] is None
    assert guest_restraints[4] is not None
    assert guest_restraints[5] is None

    # Guest - phi
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM3"
    r.mask2 = ":DM2"
    r.mask3 = ":DM1"
    r.mask4 = ":BUT@C3"
    r.attach["target"] = 120.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 120.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    guest_restraints = extract_guest_restraints(
        structure, restraints, "BUT", return_type="list"
    )

    assert guest_restraints[0] is not None
    assert guest_restraints[1] is not None
    assert guest_restraints[2] is not None
    assert guest_restraints[3] is None
    assert guest_restraints[4] is not None
    assert guest_restraints[5] is None

    # Guest - alpha
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM2"
    r.mask2 = ":DM1"
    r.mask3 = ":BUT@C3"
    r.mask4 = ":BUT@C"
    r.attach["target"] = 50.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 50.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    guest_restraints = extract_guest_restraints(
        structure, restraints, "BUT", return_type="list"
    )

    assert guest_restraints[0] is not None
    assert guest_restraints[1] is not None
    assert guest_restraints[2] is not None
    assert guest_restraints[3] is not None
    assert guest_restraints[4] is not None
    assert guest_restraints[5] is None

    # Guest - gamma
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM1"
    r.mask2 = ":BUT@C3"
    r.mask3 = ":BUT@C"
    r.mask4 = ":BUT@C2"
    r.attach["target"] = 50.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 50.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    restraints.append(r)

    # Test list
    guest_restraints = extract_guest_restraints(
        structure, restraints, "BUT", return_type="list"
    )
    assert isinstance(guest_restraints, list)

    assert guest_restraints[0] is not None
    assert guest_restraints[1] is not None
    assert guest_restraints[2] is not None
    assert guest_restraints[3] is not None
    assert guest_restraints[4] is not None
    assert guest_restraints[5] is not None

    # Test dict
    guest_restraints = extract_guest_restraints(
        structure, restraints, "BUT", return_type="dict"
    )
    assert isinstance(guest_restraints, dict)

    assert guest_restraints["r"] is not None
    assert guest_restraints["theta"] is not None
    assert guest_restraints["phi"] is not None
    assert guest_restraints["alpha"] is not None
    assert guest_restraints["beta"] is not None
    assert guest_restraints["gamma"] is not None


def test_restraints_output_modules(clean_files):
    """Test restraints output modules (colvars, plumed, openmm, amber)."""
    import paprika.build.dummy as dummy
    from paprika.restraints.amber import amber_restraint_line
    from paprika.restraints.colvars import Colvars
    from paprika.restraints.plumed import Plumed
    from paprika.restraints.utils import parse_window

    window = "p000"
    window_number, phase = parse_window(window)
    os.makedirs(os.path.join("tmp", window))

    # Test OpenMM restraint modules
    prmtop_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.prmtop")
    )
    inpcrd_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.rst7")
    )

    # Add Dummy atoms
    structure = parmed.load_file(prmtop_path, inpcrd_path, structure=True)

    dummy_atoms = {
        0: {"resname": "DM1", "x": 0.0, "y": 0.0, "z": -6.0},
        1: {"resname": "DM2", "x": 0.0, "y": 0.0, "z": -9.0},
        2: {"resname": "DM3", "x": 0.0, "y": 2.2, "z": -11.2},
    }

    # Create DAT restraints
    guest_restraints = []

    # Guest - r
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM1"
    r.mask2 = ":BUT@C3"
    r.attach["target"] = 6.0
    r.attach["num_windows"] = 4
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 5.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 24.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    guest_restraints.append(r)

    # Guest - theta
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM1"
    r.mask2 = ":BUT@C3"
    r.mask3 = ":BUT@C"
    r.attach["target"] = 180.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 180.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    guest_restraints.append(r)

    # Guest - beta
    r = DAT_restraint()
    r.amber_index = True
    r.continuous_apr = True
    r.auto_apr = True
    r.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"
    )
    r.mask1 = ":DM2"
    r.mask2 = ":DM1"
    r.mask3 = ":BUT@C3"
    r.attach["target"] = 180.0
    r.attach["num_windows"] = 15
    r.attach["fc_initial"] = 0.0
    r.attach["fc_final"] = 100.0
    r.pull["fc"] = r.attach["fc_final"]
    r.pull["num_windows"] = 46
    r.pull["target_initial"] = r.attach["target"]
    r.pull["target_final"] = 180.0
    r.release["target"] = r.pull["target_final"]
    r.release["num_windows"] = r.attach["num_windows"]
    r.release["fc_initial"] = r.attach["fc_initial"]
    r.release["fc_final"] = r.attach["fc_final"]
    r.initialize()
    guest_restraints.append(r)

    # Create OpenMM System
    system = structure.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )

    # Create restraints for OpenMM system
    k_pos = 50.0 * openmm_unit.kilocalories_per_mole / openmm_unit.angstrom**2
    apply_positional_restraints(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb"),
        system,
        k_pos=k_pos,
    )
    for restraint in guest_restraints:
        apply_dat_restraint(system, restraint, phase=phase, window_number=window_number)

    positional_restraints = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.CustomExternalForce)
    ]
    DAT_restraint_list = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.CustomBondForce)
        or isinstance(force, openmm.CustomAngleForce)
        or isinstance(force, openmm.CustomTorsionForce)
    ]

    # Test dummy atom positional restraint
    for i, force in enumerate(positional_restraints):
        particle, parameters = force.getParticleParameters(0)
        assert pytest.approx(parameters[0], abs=1e-3) == k_pos.value_in_unit(
            openmm_unit.kilojoule_per_mole / openmm_unit.nanometer**2
        )
        assert pytest.approx(parameters[1], abs=1e-3) == dummy_atoms[i]["x"] / 10
        assert pytest.approx(parameters[2], abs=1e-3) == dummy_atoms[i]["y"] / 10
        assert pytest.approx(parameters[3], abs=1e-3) == dummy_atoms[i]["z"] / 10

    # Test Amber NMR-style restraints
    atom1, atom2, parameters = DAT_restraint_list[0].getBondParameters(0)
    assert pytest.approx(parameters[0]) == 2092.0
    assert pytest.approx(parameters[1]) == 0.6

    atom1, atom2, atom3, parameters = DAT_restraint_list[1].getAngleParameters(0)
    assert pytest.approx(parameters[0]) == 418.4
    assert pytest.approx(parameters[1]) == numpy.pi

    atom1, atom2, atom3, parameters = DAT_restraint_list[2].getAngleParameters(0)
    assert pytest.approx(parameters[0]) == 418.4
    assert pytest.approx(parameters[1]) == numpy.pi

    # Test Amber restraints
    r_string = amber_restraint_line(guest_restraints[0], window)
    assert r_string.split()[2].split(",")[0] == "123"
    assert r_string.split()[2].split(",")[1] == "119"
    assert float(r_string.split()[4].split(",")[0]) == 0.0
    assert float(r_string.split()[6].split(",")[0]) == 6.0
    assert float(r_string.split()[8].split(",")[0]) == 6.0
    assert float(r_string.split()[10].split(",")[0]) == 999.0
    assert float(r_string.split()[12].split(",")[0]) == 5.0
    assert float(r_string.split()[14].split(",")[0]) == 5.0

    theta_string = amber_restraint_line(guest_restraints[1], window)
    assert theta_string.split()[2].split(",")[0] == "123"
    assert theta_string.split()[2].split(",")[1] == "119"
    assert theta_string.split()[2].split(",")[2] == "109"
    assert float(theta_string.split()[4].split(",")[0]) == 0.0
    assert float(theta_string.split()[6].split(",")[0]) == 180.0
    assert float(theta_string.split()[8].split(",")[0]) == 180.0
    assert float(theta_string.split()[10].split(",")[0]) == 180.0
    assert float(theta_string.split()[12].split(",")[0]) == 100.0
    assert float(theta_string.split()[14].split(",")[0]) == 100.0

    # Test Plumed output
    plumed = Plumed()
    plumed.path = "tmp"
    plumed.file_name = "plumed.dat"
    plumed.window_list = [window]
    plumed.restraint_list = guest_restraints
    plumed.dump_to_file()
    plumed.add_dummy_atom_restraints(structure, window)

    with open(os.path.join(plumed.path, window, "plumed.dat"), "r") as f:
        plumed_string = f.readlines()

    for iline, line in enumerate(plumed_string):
        if "DISTANCE" in line:
            restraint_line = line.split()
            if restraint_line[0] == ":c1":
                assert restraint_line[2] == "ATOMS=123,119"
            elif restraint_line[1] == ":c2":
                assert restraint_line[2] == "ATOMS=123,119,109"
            elif restraint_line[1] == ":c3":
                assert restraint_line[2] == "ATOMS=124,123,119"

        if line.startswith("RESTRAINT ..."):
            line = plumed_string[iline + 2]
            restraint_line = line.split(",")
            assert float(restraint_line[2]) == dummy_atoms[0]["z"]
            assert float(restraint_line[5]) == dummy_atoms[1]["z"]
            assert float(restraint_line[7]) == dummy_atoms[2]["y"]
            assert float(restraint_line[8]) == dummy_atoms[2]["z"]

        elif line.startswith("RESTRAINT"):
            restraint_line = line.split()
            if "c1" in restraint_line[1]:
                assert float(restraint_line[2].split("=")[1]) == 6.0
                assert float(restraint_line[3].split("=")[1]) == 10.0
            elif "c2" in restraint_line[1]:
                assert float(restraint_line[2].split("=")[1]) == 3.1416
                assert float(restraint_line[3].split("=")[1]) == 200.0
            elif "c3" in restraint_line[1]:
                assert float(restraint_line[2].split("=")[1]) == 3.1416
                assert float(restraint_line[3].split("=")[1]) == 200.0

    # Test Colvar output
    colvar = Colvars()
    colvar.path = "tmp"
    colvar.file_name = "colvars.dat"
    colvar.window_list = [window]
    colvar.restraint_list = guest_restraints
    colvar.dump_to_file()

    f = open(os.path.join(colvar.path, window, "colvars.dat"), "r+")
    for line in f:
        if "name c1" in line:
            line = f.readline()
            line = f.readline()
            line = f.readline()
            assert line.strip() == "group1 { atomNumbers 123 }"
            line = f.readline()
            assert line.strip() == "group2 { atomNumbers 119 }"

        elif "name c2" in line:
            line = f.readline()
            line = f.readline()
            line = f.readline()
            assert line.strip() == "group1 { atomNumbers 123 }"
            line = f.readline()
            assert line.strip() == "group2 { atomNumbers 119 }"
            line = f.readline()
            assert line.strip() == "group3 { atomNumbers 109 }"

        elif "name c3" in line:
            line = f.readline()
            line = f.readline()
            line = f.readline()
            assert line.strip() == "group1 { atomNumbers 124 }"
            line = f.readline()
            assert line.strip() == "group2 { atomNumbers 123 }"
            line = f.readline()
            assert line.strip() == "group3 { atomNumbers 119 }"

        elif "colvars c1" in line:
            line = f.readline()
            assert line.strip() == "centers 6.0000"
            line = f.readline()
            assert line.strip() == "forceConstant 10.0000"

        elif "colvars c2" in line:
            line = f.readline()
            assert line.strip() == "centers 180.0000"
            line = f.readline()
            assert line.strip() == "forceConstant 0.0609"

        elif "colvars c3" in line:
            line = f.readline()
            assert line.strip() == "centers 180.0000"
            line = f.readline()
            assert line.strip() == "forceConstant 0.0609"


def test_openmm_centroid_and_wall(clean_files):
    # Create OpenMM System
    prmtop_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.prmtop")
    )
    prmtop = app.AmberPrmtopFile(prmtop_path)
    system = prmtop.createSystem(
        nonbondedCutoff=app.NoCutoff,
        constraints=app.HBonds,
    )

    logger.info("### Testing OpenMM distance centroid selection.")

    distance_target = 3.0 * openff_unit.angstrom
    k_initial = 0.0 * openff_unit.kcal / openff_unit.mol / openff_unit.angstrom**2
    k_final = 3.0 * openff_unit.kcal / openff_unit.mol / openff_unit.angstrom**2

    rest1 = DAT_restraint()
    rest1.amber_index = False
    rest1.continuous_apr = False
    rest1.auto_apr = False
    rest1.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest1.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest1.mask2 = ":BUT@C"

    rest1.attach["target"] = distance_target
    rest1.attach["num_windows"] = 4
    rest1.attach["fc_initial"] = k_initial
    rest1.attach["fc_final"] = k_final
    rest1.initialize()

    assert rest1.index1 == [12, 30, 48, 66, 84, 102]
    assert rest1.index2 == [108]

    apply_dat_restraint(system, rest1, "attach", 3, force_group=1)

    bond_force = [force for force in system.getForces() if force.getForceGroup() == 1][
        0
    ]
    assert bond_force.getNumBonds() == 1
    assert bond_force.getNumGroups() == 2
    assert bond_force.getNumPerBondParameters() == 2

    bond_energy_expression = "k_bond * (r - r_0)^2;r = distance(g1, g2);"
    assert bond_force.getEnergyFunction() == bond_energy_expression

    groups, params = bond_force.getBondParameters(0)
    assert groups == (0, 1)
    assert (
        pytest.approx(params[0], abs=1e-3)
        == k_final.to(
            openff_unit.kJ / openff_unit.mol / openff_unit.nanometer**2
        ).magnitude
    )
    assert (
        pytest.approx(params[1], abs=1e-3)
        == distance_target.to(openff_unit.nanometer).magnitude
    )

    particles, weights = bond_force.getGroupParameters(0)
    assert particles == (12, 30, 48, 66, 84, 102)
    assert weights == (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    particles, weights = bond_force.getGroupParameters(1)
    assert particles == (108,)
    assert weights == (1.0,)

    logger.info("### Testing OpenMM distance upper wall selection.")

    rest1.custom_restraint_values["r2"] = 0.0
    rest1.custom_restraint_values["rk2"] = 0.0
    rest1.initialize()

    apply_dat_restraint(system, rest1, "attach", 3, force_group=2)
    bond_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 2
    ][0]
    assert "step(r - r_0)" in bond_wall_force.getEnergyFunction()

    logger.info("### Testing OpenMM distance lower wall selection.")

    rest1.custom_restraint_values = {}
    rest1.custom_restraint_values["r3"] = 0.0
    rest1.custom_restraint_values["r4"] = 0.0
    rest1.initialize()

    apply_dat_restraint(system, rest1, "attach", 3, force_group=3)
    bond_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 3
    ][0]
    assert "step(r_0 - r)" in bond_wall_force.getEnergyFunction()

    logger.info("### Testing OpenMM angle centroid selection.")

    angle_target = 180.0 * openff_unit.degrees
    k_initial = 0.0 * openff_unit.kcal / openff_unit.mol / openff_unit.radians**2
    k_final = 100.0 * openff_unit.kcal / openff_unit.mol / openff_unit.radians**2

    rest2 = DAT_restraint()
    rest2.amber_index = False
    rest2.continuous_apr = False
    rest2.auto_apr = False
    rest2.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest2.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest2.mask2 = ":BUT@C"
    rest2.mask3 = ":BUT@C3"

    rest2.attach["target"] = angle_target
    rest2.attach["num_windows"] = 4
    rest2.attach["fc_initial"] = k_initial
    rest2.attach["fc_final"] = k_final
    rest2.initialize()

    assert rest2.index1 == [12, 30, 48, 66, 84, 102]
    assert rest2.index2 == [108]
    assert rest2.index3 == [118]

    apply_dat_restraint(system, rest2, "attach", 3, force_group=4)

    angle_force = [force for force in system.getForces() if force.getForceGroup() == 4][
        0
    ]
    assert angle_force.getNumBonds() == 1
    assert angle_force.getNumGroups() == 3
    assert angle_force.getNumPerBondParameters() == 2

    angle_energy_expression = "k_angle * (theta - theta_0)^2;theta = angle(g1, g2, g3);"
    assert angle_force.getEnergyFunction() == angle_energy_expression

    groups, params = angle_force.getBondParameters(0)
    assert groups == (0, 1, 2)
    assert (
        pytest.approx(params[0], abs=1e-3)
        == k_final.to(
            openff_unit.kJ / openff_unit.mol / openff_unit.radians**2
        ).magnitude
    )
    assert (
        pytest.approx(params[1], abs=1e-3)
        == angle_target.to(openff_unit.radians).magnitude
    )

    particles, weights = angle_force.getGroupParameters(0)
    assert particles == (12, 30, 48, 66, 84, 102)
    assert weights == (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    particles, weights = angle_force.getGroupParameters(1)
    assert particles == (108,)
    assert weights == (1.0,)

    particles, weights = angle_force.getGroupParameters(2)
    assert particles == (118,)
    assert weights == (1.0,)

    logger.info("### Testing OpenMM angle upper wall selection.")

    rest2.custom_restraint_values["r2"] = 0.0
    rest2.custom_restraint_values["rk2"] = 0.0
    rest2.initialize()

    apply_dat_restraint(system, rest2, "attach", 3, force_group=5)
    angle_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 5
    ][0]
    assert "step(theta - theta_0)" in angle_wall_force.getEnergyFunction()

    logger.info("### Testing OpenMM angle lower wall selection.")

    rest2.custom_restraint_values = {}
    rest2.custom_restraint_values["r3"] = 0.0
    rest2.custom_restraint_values["r4"] = 0.0
    rest2.initialize()

    apply_dat_restraint(system, rest2, "attach", 3, force_group=6)
    angle_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 6
    ][0]
    assert "step(theta_0 - theta)" in angle_wall_force.getEnergyFunction()

    logger.info("### Testing OpenMM dihedral centroid selection.")

    rest3 = DAT_restraint()
    rest3.amber_index = False
    rest3.continuous_apr = False
    rest3.auto_apr = False
    rest3.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest3.mask1 = ":CB6@C"
    rest3.mask2 = ":CB6@O,O2,O4,O6,O8,O10"
    rest3.mask3 = ":BUT@C"
    rest3.mask4 = ":BUT@C3"

    rest3.attach["target"] = angle_target
    rest3.attach["num_windows"] = 4
    rest3.attach["fc_initial"] = k_initial
    rest3.attach["fc_final"] = k_final
    rest3.initialize()

    assert rest3.index1 == [0]
    assert rest3.index2 == [12, 30, 48, 66, 84, 102]
    assert rest3.index3 == [108]
    assert rest3.index4 == [118]

    apply_dat_restraint(system, rest3, "attach", 3, force_group=7)

    dihedral_force = [
        force for force in system.getForces() if force.getForceGroup() == 7
    ][0]
    assert dihedral_force.getNumBonds() == 1
    assert dihedral_force.getNumGroups() == 4
    assert dihedral_force.getNumPerBondParameters() == 2

    _PI_ = numpy.pi
    dihedral_energy_expression = (
        f"k_torsion * min(min(abs(theta - theta_0), abs(theta - theta_0 + "
        f"2 * {_PI_})), abs(theta - theta_0 - 2 * {_PI_}))^2;"
        "theta = dihedral(g1,g2,g3,g4);"
    )
    assert dihedral_force.getEnergyFunction() == dihedral_energy_expression

    groups, params = dihedral_force.getBondParameters(0)
    assert groups == (0, 1, 2, 3)
    assert (
        pytest.approx(params[0], abs=1e-3)
        == k_final.to(
            openff_unit.kJ / openff_unit.mol / openff_unit.radians**2
        ).magnitude
    )
    assert (
        pytest.approx(params[1], abs=1e-3)
        == angle_target.to(openff_unit.radians).magnitude
    )

    particles, weights = dihedral_force.getGroupParameters(0)
    assert particles == (0,)
    assert weights == (1.0,)

    particles, weights = dihedral_force.getGroupParameters(1)
    assert particles == (12, 30, 48, 66, 84, 102)
    assert weights == (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    particles, weights = dihedral_force.getGroupParameters(2)
    assert particles == (108,)
    assert weights == (1.0,)

    particles, weights = dihedral_force.getGroupParameters(3)
    assert particles == (118,)
    assert weights == (1.0,)

    logger.info("### Testing OpenMM dihedral upper wall selection.")

    rest3.custom_restraint_values["r2"] = 0.0
    rest3.custom_restraint_values["rk2"] = 0.0
    rest3.initialize()

    apply_dat_restraint(system, rest3, "attach", 3, force_group=8)
    dihedral_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 8
    ][0]
    assert "step(theta - theta_0)" in dihedral_wall_force.getEnergyFunction()

    logger.info("### Testing OpenMM dihedral lower wall selection.")

    rest3.custom_restraint_values = {}
    rest3.custom_restraint_values["r3"] = 0.0
    rest3.custom_restraint_values["r4"] = 0.0
    rest3.initialize()

    apply_dat_restraint(system, rest3, "attach", 3, force_group=9)
    dihedral_wall_force = [
        force for force in system.getForces() if force.getForceGroup() == 9
    ][0]
    assert "step(theta_0 - theta)" in dihedral_wall_force.getEnergyFunction()


def test_restraints_from_ascii(clean_files):
    code_block = [
        ":1@C12,:2@C4,12.5,5.0",
        ":1@O11,:1@C2,:1@C3,90.0,50.0",
        ":1@C12,:1@O13,:1@C14,1@C15,-121.16,6.0",
    ]
    tmp_file = "tmp/restraints.txt"
    with open(tmp_file, "w") as f:
        for line in code_block:
            f.writelines(f"{line}\n")

    restraints_dict = restraints_from_ascii(tmp_file)

    energy_unit = openff_unit.kcal / openff_unit.mol
    distance_unit = openff_unit.angstrom
    angle_unit = openff_unit.degree
    kdist_unit = energy_unit / distance_unit**2
    kangle_unit = energy_unit / openff_unit.radian**2

    # Distance
    assert restraints_dict["atoms"][0][0] == ":1@C12"
    assert restraints_dict["atoms"][0][1] == ":2@C4"
    assert restraints_dict["target"][0] == 12.5 * distance_unit
    assert restraints_dict["k"][0] == 5.0 * kdist_unit
    assert restraints_dict["type"][0] == RestraintType.Distance

    # Angle
    assert restraints_dict["atoms"][1][0] == ":1@O11"
    assert restraints_dict["atoms"][1][1] == ":1@C2"
    assert restraints_dict["atoms"][1][2] == ":1@C3"
    assert restraints_dict["target"][1] == 90.0 * angle_unit
    assert restraints_dict["k"][1] == 50.0 * kangle_unit
    assert restraints_dict["type"][1] == RestraintType.Angle

    # Torsion
    assert restraints_dict["atoms"][2][0] == ":1@C12"
    assert restraints_dict["atoms"][2][1] == ":1@O13"
    assert restraints_dict["atoms"][2][2] == ":1@C14"
    assert restraints_dict["atoms"][2][3] == "1@C15"
    assert restraints_dict["target"][2] == -121.16 * angle_unit
    assert restraints_dict["k"][2] == 6.0 * kangle_unit
    assert restraints_dict["type"][2] == RestraintType.Torsion
