import logging
import os
import shutil

import numpy as np
import parmed as pmd
import pytest
from pytest import approx

from paprika import analysis, log, restraints
from paprika.analysis import utils
from paprika.utils import is_file_and_not_empty

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

    seed = 12345

    fecalc = analysis.fe_calc()
    fecalc.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but-apr/vac.prmtop"
    )
    fecalc.trajectory = "*.nc"
    fecalc.path = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/")
    fecalc.restraint_list = [rest1, rest2, rest3]
    fecalc.methods = ["ti-block", "mbar-block", "mbar-autoc", "mbar-boot"]
    fecalc.boot_cycles = 100
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
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = np.array([0.0198918, 0.0451676, 0.0564517, 0.1079282, 0.1079282])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull mbar-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude
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
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = np.array([0.03, 0.07, 0.10, 0.18, 0.18])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull ti-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude

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


def test_save_and_loading(clean_files, setup_free_energy_calculation):
    # Save FE results to file
    setup_free_energy_calculation.save_results("tmp/results.json")
    assert is_file_and_not_empty("tmp/results.json") is True
    setup_free_energy_calculation.save_results("tmp/results.json", overwrite=True)
    assert is_file_and_not_empty("tmp/results.json") is True

    # Load results
    fe_calc = analysis.fe_calc()
    fe_calc.load_results("tmp/results.json")
    assert np.isclose(-4.34372240, fe_calc.results["ref_state_work"].magnitude)

    # Save Simulation data to file
    setup_free_energy_calculation.save_data("tmp/simulation_data.json")
    assert is_file_and_not_empty("tmp/simulation_data.json") is True
    setup_free_energy_calculation.save_data("tmp/simulation_data.json", overwrite=True)
    assert is_file_and_not_empty("tmp/simulation_data.json") is True

    # Load simulation data
    fe_calc.load_data("tmp/simulation_data.json")
    assert fe_calc.simulation_data is not None
    assert fe_calc.changing_restraints is not None
    assert fe_calc.orders is not None


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


def test_bootstrap():
    """Test the utility modules in `analysis`"""

    # Test regression statistics
    x = np.linspace(0, 10, 11)
    y = np.linspace(0, 10, 11)
    stats = analysis.summarize_statistics(x, y)

    assert all(stats == np.array([1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0]))

    # Test Regression bootstrap
    x_sem = np.ones_like(x)
    y_sem = np.ones_like(y)

    np.random.seed(0)
    results = analysis.regression_bootstrap(x, x_sem, y, y_sem, cycles=1)
    compare = {
        "slope": 0.7703030693742714,
        "intercept": 0.7777706430286471,
        "R": 0.9184956779750857,
        "R**2": 0.8436343104589124,
        "RMSE": 1.637675285465762,
        "MSE": -0.40918918963339473,
        "MUE": 1.2814477376354574,
        "Tau": 0.8545454545454545,
    }

    for stat in results["mean"]:
        assert pytest.approx(results["mean"][stat], abs=1e-3) == compare[stat]

    # Test dG Bootstrap
    np.random.seed(0)
    results = analysis.dG_bootstrap(-12, 2, -12, 2, cycles=1, with_uncertainty=True)
    assert pytest.approx(results["mean"], abs=1e-3) == -11.205587976469898
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -11.20558798
    assert pytest.approx(results["ci"][1], abs=1e-3) == -11.20558798

    results = analysis.dG_bootstrap(-12, 2, -12, 2, cycles=1, with_uncertainty=False)
    assert pytest.approx(results["mean"], abs=1e-3) == -0.4106792724182964
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -0.41067927
    assert pytest.approx(results["ci"][1], abs=1e-3) == -0.41067927

    # Test dH Bootstrap
    np.random.seed(0)
    results = analysis.dH_bootstrap(
        -15,
        2,
        -15,
        2,
        -10,
        2,
        -10,
        2,
        cycles=1,
        with_uncertainty=True,
    )
    assert pytest.approx(results["mean"], abs=1e-3) == -11.50986102184404
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -11.5098610
    assert pytest.approx(results["ci"][1], abs=1e-3) == -11.5098610


def test_utils():
    """Test the utility modules in `analysis`"""
    assert utils.get_factors(10) == [1, 2, 5, 10]
    assert utils.get_nearest_max(100) == 90
    np.random.seed(0)
    results = utils.get_block_sem(np.random.normal(10.0, 2.0, 100))
    assert pytest.approx(results, abs=1e-3) == 0.3916368835724714
    assert utils.get_subsampled_indices(10, 2.0) == [0, 2, 4, 6, 8]
