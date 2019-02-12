import parmed as pmd
import numpy as np
import os
import shutil

from paprika import restraints
from paprika import analysis

from paprika.tests import addons

import pytest


@pytest.fixture(scope="function", autouse=True)
def clean_files(directory="tmp"):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


def test_fe_calc(clean_files):

    inputpdb = pmd.load_file(os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb"))

    # Distance restraint
    rest1 = restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = inputpdb
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
    rest2.topology = inputpdb
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
    rest3.topology = inputpdb
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
    restraints.create_window_list([rest1, rest2, rest3])

    # Phase abbreviations
    phase_dict = {"a": "attach", "p": "pull", "r": "release"}

    np.random.seed(12345)

    fecalc = analysis.fe_calc()
    fecalc.prmtop = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/vac.prmtop")
    fecalc.trajectory = "*.nc"
    fecalc.path = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/")
    fecalc.restraint_list = [rest1, rest2, rest3]
    fecalc.methods = ["mbar-block", "ti-block"]
    fecalc.bootcycles = 100
    fecalc.ti_matrix = "diagonal"
    fecalc.compute_largest_neighbor = True
    fecalc.compute_roi = True
    fecalc.collect_data(single_prmtop=True)
    fecalc.compute_free_energy()

    ##################
    # Test mbar-block
    ##################

    method = "mbar-block"
    # Test mbar-block free energies and uncertainties
    test_vals = [
        fecalc.results["attach"][method]["fe"],
        fecalc.results["attach"][method]["sem"],
        fecalc.results["pull"][method]["fe"],
        fecalc.results["pull"][method]["sem"],
    ]
    ref_vals = [13.267731176, 0.16892084090, -2.1791430735, 0.93638948302]
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i], rtol=0.0, atol=0.00001)

    # Test attach mbar-block largest_neighbor values
    test_vals = fecalc.results["attach"][method]["largest_neighbor"]
    ref_vals = np.array([0.0198918, 0.0451676, 0.0564517, 0.1079282, 0.1079282])
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i], rtol=0.0, atol=0.00001)

    # Test pull mbar-block largest_neighbor values
    test_vals = fecalc.results["pull"][method]["largest_neighbor"]
    ref_vals = np.array(
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
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i], rtol=0.0, atol=0.00001)

    ##################
    # Test ti-block
    ##################

    method = "ti-block"

    # Test ti-block free energies and uncertainties
    test_vals = [
        fecalc.results["attach"][method]["fe"],
        fecalc.results["attach"][method]["sem"],
        fecalc.results["pull"][method]["fe"],
        fecalc.results["pull"][method]["sem"],
    ]
    ref_vals = np.array([13.35823327, 0.25563407, -1.77873259, 0.95162741])
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i])

    # Test attach roi values
    test_vals = fecalc.results["attach"][method]["roi"]
    ref_vals = np.array(
        [-0.00027503, -0.00020536, 0.00000691, -0.00062548, -0.00029785]
    )
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i])

    # Test pull roi values
    test_vals = fecalc.results["pull"][method]["roi"]
    ref_vals = np.array(
        [
            -0.00191275,
            -0.00047460,
            -0.00284036,
            -0.00144093,
            -0.00199245,
            -0.00160141,
            -0.00071836,
            -0.00188979,
            -0.00162371,
            -0.00108091,
            -0.00197727,
            -0.00223339,
            -0.00144556,
            -0.00272537,
            -0.00098154,
            -0.00071218,
            -0.00073691,
            -0.00195884,
            -0.00155617,
        ]
    )
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i])

    # Test attach ti-block largest_neighbor values
    test_vals = fecalc.results["attach"][method]["largest_neighbor"]
    ref_vals = np.array([0.02809982, 0.07439878, 0.09719278, 0.16782417, 0.16782417])
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i])

    # Test pull ti-block largest_neighbor values
    test_vals = fecalc.results["pull"][method]["largest_neighbor"]
    ref_vals = np.array(
        [
            0.31296449,
            0.31296449,
            0.23297017,
            0.24211187,
            0.24211187,
            0.13217678,
            0.11879266,
            0.15423428,
            0.15423428,
            0.12861334,
            0.11276903,
            0.11276903,
            0.11150722,
            0.13483801,
            0.17510575,
            0.17510575,
            0.14725404,
            0.12972769,
            0.11803089,
        ]
    )
    for i in range(len(test_vals)):
        assert np.isclose(ref_vals[i], test_vals[i])

    #######################
    # Test ref_state_work
    #######################

    # Test reference state calculation
    fecalc.compute_ref_state_work([rest1, rest2, rest3, None, None, None])
    assert np.isclose(-4.34372240, fecalc.results["ref_state_work"])


def reprint_values(results, method):
    """
    Hack method to reprint values when I rearrange things and mess up the order of random
    number generation.
    """
    reprint = np.concatenate(
        (
            np.array(
                [
                    results["attach"][method]["fe"],
                    results["attach"][method]["sem"],
                    results["pull"][method]["fe"],
                    results["pull"][method]["sem"],
                ]
            ),
            np.array([9999]),
            results["attach"][method]["roi"],
            np.array([9999]),
            results["pull"][method]["roi"],
            np.array([9999]),
            results["attach"][method]["largest_neighbor"],
            np.array([9999]),
            results["pull"][method]["largest_neighbor"],
            np.array([9999]),
        )
    )
    for val in reprint:
        if val == 9999:
            print("")
        else:
            # Comment this to make python 2.7 pass.  This is just for redoing the
            # the ref_vals anyway
            #            print("{:.8f}, ".format(val), end='')
            pass
