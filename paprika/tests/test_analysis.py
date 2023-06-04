import logging
import os
import shutil
from copy import deepcopy

import numpy
import parmed
import pytest
from openff.units import unit as openff_unit
from pytest import approx

from paprika import analysis, log, restraints
from paprika.analysis import utils
from paprika.utils import is_file_and_not_empty

log.config_root_logger(verbose=True)
logger = logging.getLogger(__name__)


random_seed = 12345


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
    input_pdb = parmed.load_file(
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

    # Create Analysis Instance
    fecalc = analysis.fe_calc()
    fecalc.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but-apr/vac.prmtop"
    )
    fecalc.trajectory = "*.nc"
    fecalc.path = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/")
    fecalc.restraint_list = [rest1, rest2, rest3]
    fecalc.boot_cycles = 100
    fecalc.ti_matrix = "diagonal"
    fecalc.compute_largest_neighbor = True
    fecalc.compute_roi = True
    fecalc.conservative_subsample = False
    fecalc.exact_sem_each_ti_fraction = False
    fecalc.fractions = [0.2, 0.4, 0.6, 0.8, 1.0]
    fecalc.energy_unit = openff_unit.kcal / openff_unit.mole
    fecalc.distance_unit = openff_unit.angstrom
    fecalc.angle_unit = openff_unit.degrees
    fecalc.temperature_unit = openff_unit.kelvin
    fecalc.collect_data(single_topology=True)
    fecalc.compute_ref_state_work([rest1, rest2, rest3, None, None, None])

    return fecalc


def test_setup(clean_files, setup_free_energy_calculation):
    assert setup_free_energy_calculation.temperature_unit == openff_unit.kelvin
    assert setup_free_energy_calculation.distance_unit == openff_unit.angstrom
    assert setup_free_energy_calculation.angle_unit == openff_unit.degrees
    assert (
        setup_free_energy_calculation.energy_unit == openff_unit.kcal / openff_unit.mole
    )
    assert setup_free_energy_calculation.fractions == [0.2, 0.4, 0.6, 0.8, 1.0]
    results = deepcopy(setup_free_energy_calculation.results)
    setup_free_energy_calculation.results = results
    assert setup_free_energy_calculation.results == results
    assert setup_free_energy_calculation.exact_sem_each_ti_fraction is False
    assert setup_free_energy_calculation.conservative_subsample is False

    # Test save and load results -- JSON
    results = deepcopy(setup_free_energy_calculation.results)
    setup_free_energy_calculation.save_results("tmp/results.json", overwrite=True)
    assert is_file_and_not_empty("tmp/results.json")
    setup_free_energy_calculation.results = {}
    setup_free_energy_calculation.load_results("tmp/results.json")
    assert len(setup_free_energy_calculation.results) == len(results)

    # Test save and load simulation data -- JSON
    setup_free_energy_calculation.save_simulation_data_to_json(
        "tmp/simulation.json", overwrite=True
    )
    assert is_file_and_not_empty("tmp/simulation.json")
    setup_free_energy_calculation.changing_restraints = None
    setup_free_energy_calculation.orders = None
    setup_free_energy_calculation.simulation_data = None
    setup_free_energy_calculation.load_simulation_data_from_json("tmp/simulation.json")
    assert setup_free_energy_calculation.changing_restraints is not None
    assert setup_free_energy_calculation.orders is not None
    assert setup_free_energy_calculation.simulation_data is not None


def test_mbar_block(clean_files, setup_free_energy_calculation):
    method = "mbar-block"

    # Estimate FE with `mbar-block`
    setup_free_energy_calculation.methods = [method]
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test `mbar-block` free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = [13.2677, 0.1689, -2.1791, 0.9364]
    assert reference_values == approx(test_vals, abs=0.01)

    # Test attach `mbar-block` largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.0199, 0.0452, 0.0565, 0.1079, 0.1079])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull `mbar-block` largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array(
        [
            # fmt: off
            0.20538, 0.20538, 0.16174, 0.17477, 0.52550,
            0.52550, 0.11499, 0.17079, 0.21291, 0.21291,
            0.19422, 0.17689, 0.19973, 0.19973, 0.20148,
            0.20148, 0.14707, 0.14425, 0.14344
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_mbar_autoc(clean_files, setup_free_energy_calculation):
    method = "mbar-autoc"

    setup_free_energy_calculation.methods = [method]
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test `mbar-autoc` free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = [13.26773, 0.0808, -2.1791, 0.6969]
    assert reference_values == approx(test_vals, abs=0.01)

    # Test attach `mbar-autoc` largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.0199, 0.0259, 0.0336, 0.0383, 0.0383])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull `mbar-autoc` largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude
    numpy.savetxt("tmp/test.txt", test_vals, fmt="%.5f")
    reference_values = numpy.array(
        [
            # fmt: off
            0.10274, 0.11361, 0.13074, 0.14136, 0.38928,
            0.38928, 0.11215, 0.12951, 0.13145, 0.13145,
            0.13113, 0.13113, 0.13751, 0.14186, 0.14186,
            0.12847, 0.13550, 0.13550, 0.13404
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_ti_block(clean_files, setup_free_energy_calculation):
    method = "ti-block"

    # Estimate FE with exact sem
    setup_free_energy_calculation.methods = [method]
    setup_free_energy_calculation.exact_sem_each_ti_fraction = True
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test ti-block free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = numpy.array([13.31, 0.25, -1.62, 0.87])
    assert reference_values == approx(test_vals, abs=0.01)

    # ROI only runs during TI.

    # Test attach ti-block largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.03, 0.07, 0.10, 0.18, 0.18])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull ti-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude

    reference_values = numpy.array(
        [
            # fmt: off
            0.33156, 0.33156, 0.21509, 0.22191, 0.22191,
            0.10746, 0.13514, 0.15078, 0.15078, 0.15518,
            0.10678, 0.10678, 0.10158, 0.14123, 0.16609,
            0.16609, 0.14719, 0.14090, 0.11006
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)

    # Estimate FE without exact sem
    setup_free_energy_calculation.exact_sem_each_ti_fraction = False
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test ti-block free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = numpy.array([13.31, 0.25, -1.62, 0.87])
    assert reference_values == approx(test_vals, abs=0.01)

    # ROI only runs during TI.

    # Test attach ti-block largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.03, 0.07, 0.10, 0.18, 0.18])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull ti-block largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude

    reference_values = numpy.array(
        [
            # fmt: off
            0.33156, 0.33156, 0.21509, 0.22191, 0.22191,
            0.10746, 0.13514, 0.15078, 0.15078, 0.15518,
            0.10678, 0.10678, 0.10157, 0.14122, 0.16608,
            0.16608, 0.14718, 0.14090, 0.11005
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_ti_nocor(clean_files, setup_free_energy_calculation):
    method = "ti-nocor"

    # Estimate FE with exact sem
    setup_free_energy_calculation.methods = [method]
    setup_free_energy_calculation.exact_sem_each_ti_fraction = True
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test `ti-nocor` free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = numpy.array([13.34, 0.09, -1.71, 0.56])
    assert reference_values == approx(test_vals, abs=0.01)

    # ROI only runs during TI.

    # Test attach `ti-nocor` largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.014, 0.036, 0.055, 0.066, 0.066])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull `ti-nocor` largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude

    reference_values = numpy.array(
        [
            # fmt: off
            0.09992, 0.09992, 0.08598, 0.10866, 0.10866,
            0.08704, 0.10761, 0.10761, 0.11872, 0.11872,
            0.10678, 0.10678, 0.10158, 0.09112, 0.10342,
            0.10342, 0.10319, 0.10608, 0.11005
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)

    # No-exact sem
    setup_free_energy_calculation.exact_sem_each_ti_fraction = False
    setup_free_energy_calculation.compute_free_energy(seed=random_seed)
    results = setup_free_energy_calculation.results

    # Test `ti-nocor` free energies and uncertainties
    test_vals = [
        results["attach"][method]["fe"].magnitude,
        results["attach"][method]["sem"].magnitude,
        results["pull"][method]["fe"].magnitude,
        results["pull"][method]["sem"].magnitude,
    ]
    reference_values = numpy.array([13.34, 0.098, -1.71, 0.56])
    assert reference_values == approx(test_vals, abs=0.01)

    # ROI only runs during TI.

    # Test attach `ti-nocor` largest_neighbor values
    test_vals = results["attach"][method]["largest_neighbor"].magnitude
    reference_values = numpy.array([0.0138, 0.0362, 0.055, 0.0657, 0.0657])
    assert reference_values == approx(test_vals, abs=0.01)

    # Test pull `ti-nocor` largest_neighbor values
    test_vals = results["pull"][method]["largest_neighbor"].magnitude

    reference_values = numpy.array(
        [
            # fmt: off
            0.09992, 0.09992, 0.08598, 0.10867, 0.10866,
            0.08704, 0.10761, 0.10761, 0.11872, 0.11872,
            0.10678, 0.10678, 0.10158, 0.09112, 0.10342,
            0.10342, 0.10319, 0.10608, 0.11006
            # fmt: on
        ]
    )
    assert reference_values == approx(test_vals, abs=0.01)


def test_reference_state_work(clean_files, setup_free_energy_calculation):
    results = setup_free_energy_calculation.results
    assert numpy.isclose(-4.34372, results["ref_state_work"].magnitude)


def test_save_and_loading(clean_files, setup_free_energy_calculation):
    # Save FE results to file
    setup_free_energy_calculation.save_results("tmp/results.json")
    assert is_file_and_not_empty("tmp/results.json") is True
    setup_free_energy_calculation.save_results("tmp/results.json", overwrite=True)
    assert is_file_and_not_empty("tmp/results.json") is True

    # Load results
    fe_calc = analysis.fe_calc()
    fe_calc.load_results("tmp/results.json")
    assert numpy.isclose(-4.34372240, fe_calc.results["ref_state_work"].magnitude)

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
    input_pdb = parmed.load_file(
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
    x = numpy.linspace(0, 10, 11)
    y = numpy.linspace(0, 10, 11)
    stats = analysis.summarize_statistics(x, y)

    assert all(stats == numpy.array([1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0]))

    # Test Regression bootstrap
    x_sem = numpy.ones_like(x)
    y_sem = numpy.ones_like(y)

    numpy.random.seed(0)
    results = analysis.regression_bootstrap(x, x_sem, y, y_sem, cycles=1)
    # fmt: off
    compare = {
        "slope": 0.77030,
        "intercept": 0.77778,
        "R": 0.91850,
        "R**2": 0.84363,
        "RMSE": 1.63768,
        "MSE": -0.40919,
        "MUE": 1.28145,
        "Tau": 0.85454,
    }
    # fmt: on

    for stat in results["mean"]:
        assert pytest.approx(results["mean"][stat], abs=1e-3) == compare[stat]

    # Test dG Bootstrap
    numpy.random.seed(0)
    results = analysis.dG_bootstrap(-12, 2, -12, 2, cycles=1, with_uncertainty=True)
    assert pytest.approx(results["mean"], abs=1e-3) == -11.20559
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -11.20559
    assert pytest.approx(results["ci"][1], abs=1e-3) == -11.20559

    results = analysis.dG_bootstrap(-12, 2, -12, 2, cycles=1, with_uncertainty=False)
    assert pytest.approx(results["mean"], abs=1e-3) == -0.41068
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -0.41068
    assert pytest.approx(results["ci"][1], abs=1e-3) == -0.41068

    # Test dH Bootstrap
    numpy.random.seed(0)
    # fmt: off
    results = analysis.dH_bootstrap(
        -15, 2, -15, 2, -10, 2, -10, 2,
        cycles=1,
        with_uncertainty=True,
    )
    # fmt: on
    assert pytest.approx(results["mean"], abs=1e-3) == -11.50986
    assert pytest.approx(results["sem"], abs=1e-3) == 0.0
    assert pytest.approx(results["ci"][0], abs=1e-3) == -11.50986
    assert pytest.approx(results["ci"][1], abs=1e-3) == -11.50986


def test_utils():
    """Test the utility modules in `analysis`"""
    assert utils.get_factors(10) == [1, 2, 5, 10]
    assert utils.get_nearest_max(100) == 90

    numpy.random.seed(0)
    results = utils.get_block_sem(numpy.random.normal(10.0, 2.0, 100))
    assert pytest.approx(results, abs=1e-3) == 0.39164
    assert utils.get_subsampled_indices(10, 2.0) == [0, 2, 4, 6, 8]
