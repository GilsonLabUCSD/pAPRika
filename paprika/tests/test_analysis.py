import logging
import os
import shutil

import numpy as np
import parmed as pmd
import pytest
from pytest import approx

from paprika import analysis, log, restraints

log.config_root_logger(verbose=True)
logger = logging.getLogger(__name__)


@pytest.fixture(scope="module", autouse=True)
def clean_files(directory="tmp"):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


@pytest.fixture(scope="module", autouse=True)
def setup_free_energy_calculation():
    input_pdb = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )

    # Distance restraint
    rest1 = restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = input_pdb
    rest1.mask1 = ":CB6@O"
    rest1.mask2 = ":BUT@C1"
    rest1.attach["target"] = 4.5
    rest1.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest1.attach["fc_final"] = 5.0
    rest1.pull["fc"] = rest1.attach["fc_final"]
    rest1.pull["target_initial"] = rest1.attach["target"]
    rest1.pull["target_final"] = 18.5
    rest1.pull["num_windows"] = 19
    rest1.initialize()

    # Angle restraint
    rest2 = restraints.DAT_restraint()
    rest2.continuous_apr = True
    rest2.amber_index = True
    rest2.topology = input_pdb
    rest2.mask1 = ":CB6@O1"
    rest2.mask2 = ":CB6@O"
    rest2.mask3 = ":BUT@C1"
    rest2.attach["target"] = 8.0
    rest2.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest2.attach["fc_final"] = 50.0
    rest2.pull["fc"] = rest2.attach["fc_final"]
    rest2.pull["target_initial"] = rest2.attach["target"]
    rest2.pull["target_final"] = rest2.attach["target"]
    rest2.pull["num_windows"] = 19
    rest2.initialize()

    # Dihedral restraint
    rest3 = restraints.DAT_restraint()
    rest3.continuous_apr = True
    rest3.amber_index = True
    rest3.topology = input_pdb
    rest3.mask1 = ":CB6@O11"
    rest3.mask2 = ":CB6@O1"
    rest3.mask3 = ":CB6@O"
    rest3.mask4 = ":BUT@C1"
    rest3.attach["target"] = -60.0
    rest3.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest3.attach["fc_final"] = 50.0
    rest3.pull["fc"] = rest3.attach["fc_final"]
    rest3.pull["target_initial"] = rest3.attach["target"]
    rest3.pull["target_final"] = rest3.attach["target"]
    rest3.pull["num_windows"] = 19
    rest3.initialize()

    # Create window directories
    restraints.restraints.create_window_list([rest1, rest2, rest3])

    seed = 12345

    fecalc = analysis.fe_calc()
    fecalc.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but-apr/vac.prmtop"
    )
    fecalc.trajectory = "*.nc"
    fecalc.path = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/")
    fecalc.restraint_list = [rest1, rest2, rest3]
    fecalc.methods = ["mbar-block", "ti-block", "mbar-autoc"]
    fecalc.bootcycles = 100
    fecalc.ti_matrix = "diagonal"
    fecalc.compute_largest_neighbor = True
    fecalc.compute_roi = True
    fecalc.collect_data(single_topology=True)
    fecalc.compute_free_energy(seed=seed)
    fecalc.compute_ref_state_work([rest1, rest2, rest3, None, None, None])

    return fecalc


def test_setup(clean_files, setup_free_energy_calculation):
    pass


def test_mbar_block(clean_files, setup_free_energy_calculation):
    results = setup_free_energy_calculation.results
    method = "mbar-block"

    # Test mbar-block free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = [13.267731176, 0.16892084090, -2.1791430735, 0.93638948302]
    assert reference_values == approx(test_vals, abs=0.01)

    # Test attach mbar-block largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"]
    reference_values = np.array([0.0198918, 0.0451676, 0.0564517, 0.1079282, 0.1079282])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull mbar-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"]
    reference_values = np.array(
        [
            0.2053769,
            0.2053769,
            0.1617423,
            0.1747668,
            0.5255023,
            0.5255023,
            0.1149945,
            0.1707901,
            0.2129136,
            0.2129136,
            0.1942189,
            0.1768906,
            0.1997338,
            0.1997338,
            0.2014766,
            0.2014766,
            0.1470727,
            0.1442517,
            0.1434395,
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_ti_block(clean_files, setup_free_energy_calculation):
    results = setup_free_energy_calculation.results

    method = "ti-block"

    # Test ti-block free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = np.array([13.35, 0.26, -1.85, 0.78])
    assert reference_values == approx(test_vals, abs=0.01)

    # ROI only runs during TI.

    # Test attach ti-block largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"]
    reference_values = np.array([0.03, 0.07, 0.10, 0.18, 0.18])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull ti-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"]

    reference_values = np.array(
        [
            0.33156402,
            0.33156402,
            0.22515133,
            0.2219127,
            0.2219127,
            0.1311959,
            0.13514015,
            0.15078472,
            0.15078472,
            0.12448228,
            0.10678047,
            0.10678047,
            0.10157904,
            0.14122943,
            0.16608568,
            0.16608568,
            0.14718857,
            0.14090383,
            0.11005729,
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_reference_state_work(clean_files, setup_free_energy_calculation):
    results = setup_free_energy_calculation.results
    assert np.isclose(-4.34372240, results["ref_state_work"].magnitude)


def test_temperature(clean_files):
    input_pdb = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-dum.pdb")
    )

    # Distance restraint
    rest1 = restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = input_pdb
    rest1.mask1 = ":DM1"
    rest1.mask2 = ":BUT@C"
    rest1.attach["target"] = 6.0
    rest1.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest1.attach["fc_final"] = 5.0
    rest1.pull["fc"] = rest1.attach["fc_final"]
    rest1.pull["target_initial"] = rest1.attach["target"]
    rest1.pull["target_final"] = 24.0
    rest1.pull["num_windows"] = 19
    rest1.initialize()

    # theta restraint
    rest2 = restraints.DAT_restraint()
    rest2.continuous_apr = True
    rest2.amber_index = True
    rest2.topology = input_pdb
    rest2.mask1 = ":DM2"
    rest2.mask2 = ":DM1"
    rest2.mask3 = ":BUT@C"
    rest2.attach["target"] = 180.0
    rest2.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest2.attach["fc_final"] = 100.0
    rest2.pull["fc"] = rest2.attach["fc_final"]
    rest2.pull["target_initial"] = rest2.attach["target"]
    rest2.pull["target_final"] = rest2.attach["target"]
    rest2.pull["num_windows"] = 19
    rest2.initialize()

    # beta restraint
    rest3 = restraints.DAT_restraint()
    rest3.continuous_apr = True
    rest3.amber_index = True
    rest3.topology = input_pdb
    rest3.mask1 = ":DM1"
    rest3.mask2 = ":BUT@C"
    rest3.mask3 = ":BUT@C3"
    rest3.attach["target"] = 180.0
    rest3.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest3.attach["fc_final"] = 100.0
    rest3.pull["fc"] = rest3.attach["fc_final"]
    rest3.pull["target_initial"] = rest3.attach["target"]
    rest3.pull["target_final"] = rest3.attach["target"]
    rest3.pull["num_windows"] = 19
    rest3.initialize()

    fecalc = analysis.fe_calc()
    fecalc.temperature = 298.15
    assert pytest.approx(fecalc.beta.magnitude, abs=1e-3) == 1.68780
    fecalc.compute_ref_state_work([rest1, rest2, None, None, rest3, None])
    assert (
        pytest.approx(fecalc.results["ref_state_work"].magnitude, abs=1e-3) == -7.14151
    )

    fecalc = analysis.fe_calc()
    fecalc.temperature = 328.15
    assert pytest.approx(fecalc.beta.magnitude, abs=1e-3) == 1.53285
    fecalc.compute_ref_state_work([rest1, rest2, None, None, rest3, None])
    assert (
        pytest.approx(fecalc.results["ref_state_work"].magnitude, abs=1e-3) == -7.70392
    )

    fecalc = analysis.fe_calc()
    fecalc.temperature = 278.15
    assert pytest.approx(fecalc.beta.magnitude, abs=1e-3) == 1.80839
    fecalc.compute_ref_state_work([rest1, rest2, None, None, rest3, None])
    assert (
        pytest.approx(fecalc.results["ref_state_work"].magnitude, abs=1e-3) == -6.75834
    )
