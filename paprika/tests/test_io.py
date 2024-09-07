"""
Test that we can save and load restraints as JSON.
"""

import os
import shutil

import pytest

from paprika.io import (
    load_restraints,
    load_trajectory,
    read_restraint_data,
    save_restraints,
)
from paprika.restraints import DAT_restraint


@pytest.fixture(scope="function", autouse=True)
def clean_files(directory="tmp"):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


def test_save_and_load_single_restraint(clean_files):
    """Test we can save a simple restraint"""
    rest = DAT_restraint()
    rest.amber_index = True
    rest.continuous_apr = False
    rest.auto_apr = False
    rest.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    rest.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest.mask2 = ":BUT@C3"
    rest.attach["target"] = 3.0
    rest.attach["num_windows"] = 4
    rest.attach["fc_initial"] = 0.0
    rest.attach["fc_final"] = 3.0
    rest.pull["fc"] = rest.attach["fc_final"]
    rest.pull["num_windows"] = 4
    rest.pull["target_initial"] = rest.attach["target"]
    rest.pull["target_final"] = 6.0
    rest.release["target"] = rest.pull["target_final"]
    rest.release["num_windows"] = rest.attach["num_windows"]
    rest.release["fc_initial"] = rest.attach["fc_initial"]
    rest.release["fc_final"] = rest.attach["fc_final"]
    rest.initialize()
    save_restraints([rest], os.path.join("tmp", "rest.json"))
    restraints = load_restraints(os.path.join("tmp", "rest.json"))
    assert rest == restraints[0]


def test_save_and_load_list_restraints(clean_files):
    """Test we can save and load a list of restraints"""
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

    save_restraints([rest1, rest2], os.path.join("tmp", "restraints.json"))
    restraints = load_restraints(os.path.join("tmp", "restraints.json"))
    assert rest1 == restraints[0]
    assert rest2 == restraints[1]


def test_load_trajectory_and_restraints():
    window = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/a000")
    trajectory_name = "md.nc"
    topology = "../vac.prmtop"
    trajectory = load_trajectory(
        window, trajectory_name, topology, single_topology=False
    )
    assert trajectory is not None

    restraint = DAT_restraint()
    restraint.amber_index = True
    restraint.continuous_apr = False
    restraint.auto_apr = False
    restraint.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
    )
    restraint.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    restraint.mask2 = ":BUT@C3"
    restraint.attach["target"] = 3.0
    restraint.attach["num_windows"] = 4
    restraint.attach["fc_initial"] = 0.0
    restraint.attach["fc_final"] = 3.0
    restraint.initialize()

    value = read_restraint_data(trajectory, restraint)

    assert value is not None
