import os
import shutil

import parmed as pmd
import pytest

from paprika import restraints
from paprika.amber import Simulation
from paprika.restraints import amber
from paprika.restraints.restraints import create_window_list
from paprika.tests import addons
from paprika.utils import parse_mden, parse_mdout


@pytest.fixture(scope="function", autouse=True)
def clean_files(directory="tmp"):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


@addons.using_sander
# @addons.using_pmemd_cuda
def test_amber_single_window_gbmin(clean_files):
    # Distance restraint
    restraint = restraints.DAT_restraint()
    restraint.continuous_apr = True
    restraint.amber_index = True
    restraint.topology = os.path.join(
        os.path.dirname(__file__), "../data/k-cl/k-cl.pdb"
    )
    restraint.mask1 = ":K+"
    restraint.mask2 = ":Cl-"
    restraint.attach["target"] = 4.5
    restraint.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    restraint.attach["fc_final"] = 5.0
    restraint.pull["fc"] = restraint.attach["fc_final"]
    restraint.pull["target_initial"] = restraint.attach["target"]
    restraint.pull["target_final"] = 18.5
    restraint.pull["num_windows"] = 5
    restraint.initialize()

    windows_directory = os.path.join("tmp", "k-cl", "windows")
    window_list = create_window_list([restraint])

    for window in window_list:
        os.makedirs(os.path.join(windows_directory, window))
        with open(
            os.path.join(windows_directory, window, "restraints.in"), "a"
        ) as file:
            string = amber.amber_restraint_line(restraint, window)
            file.write(string)

    for window in window_list:
        if window[0] == "a":
            structure = pmd.load_file(
                os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl-sol.prmtop"),
                os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl-sol.rst7"),
                structure=True,
            )
            for atom in structure.atoms:
                if atom.name == "Cl-":
                    atom.xz = 2.65
            structure.save(
                os.path.join(windows_directory, window, "k-cl.prmtop"), overwrite=True
            )
            structure.save(
                os.path.join(windows_directory, window, "k-cl.rst7"), overwrite=True
            )
            structure.save(
                os.path.join(windows_directory, window, "k-cl.pdb"), overwrite=True
            )

        elif window[0] == "p":
            structure = pmd.load_file(
                os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl-sol.prmtop"),
                os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl-sol.rst7"),
                structure=True,
            )
            target_difference = (
                restraint.phase["pull"]["targets"][int(window[1:])]
                - restraint.phase["pull"]["targets"][0]
            )

            for atom in structure.atoms:
                if atom.name == "Cl-":
                    atom.xz += target_difference
            structure.save(
                os.path.join(windows_directory, window, "k-cl.prmtop"), overwrite=True
            )
            structure.save(
                os.path.join(windows_directory, window, "k-cl.rst7"), overwrite=True
            )
            structure.save(
                os.path.join(windows_directory, window, "k-cl.pdb"), overwrite=True
            )

    gbsim = Simulation()
    gbsim.path = os.path.join("tmp", "k-cl", "windows", "a003")
    gbsim.executable = "sander"
    gbsim.topology = "k-cl.prmtop"
    gbsim.prefix = "minimize"
    gbsim.inpcrd = "k-cl.rst7"
    gbsim.config_gb_min()
    gbsim.cntrl["maxcyc"] = 1
    gbsim.cntrl["ncyc"] = 1
    gbsim.run()

    gbsim.config_gb_md()
    gbsim.prefix = "md"
    gbsim.inpcrd = "minimize.rst7"
    gbsim.cntrl["nstlim"] = 1
    gbsim.cntrl["ntwe"] = 1
    gbsim.cntrl["ntpr"] = 1
    gbsim.cntrl["ig"] = 777
    gbsim.run()

    mden = parse_mden(os.path.join("tmp", "k-cl", "windows", "a003", "md.mden"))

    assert pytest.approx(mden["Bond"][0]) == 0
    assert pytest.approx(mden["Angle"][0]) == 0
    assert pytest.approx(mden["Dihedral"][0]) == 0
    assert pytest.approx(mden["V14"][0]) == 0
    assert pytest.approx(mden["E14"][0]) == 0

    assert pytest.approx(mden["VDW"][0], 0.1) == 25956.13225
    assert pytest.approx(mden["Ele"][0], 0.1) == -18828.99631
    assert pytest.approx(mden["Total"][0], 0.1) == 7127.13594


def test_amber_minimization(clean_files):
    simulation = Simulation()
    simulation.path = os.path.join("tmp")

    shutil.copy(
        os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.prmtop"), "tmp"
    )
    shutil.copy(
        os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.rst7"), "tmp"
    )

    simulation.executable = "sander"
    simulation.restraint_file = None

    simulation.prefix = "minimize"
    simulation.topology = "k-cl.prmtop"
    simulation.inpcrd = "k-cl.rst7"

    simulation.config_gb_min()
    # Turn off GB for now.
    simulation.cntrl["igb"] = 0
    simulation.cntrl["ntb"] = 0

    simulation.run()

    mdout = parse_mdout(os.path.join("tmp", "minimize.out"))

    assert pytest.approx(mdout["Bond"][-1]) == 0
    assert pytest.approx(mdout["Angle"][-1]) == 0
    assert pytest.approx(mdout["Dihedral"][-1]) == 0
    assert pytest.approx(mdout["V14"][-1]) == 0
    assert pytest.approx(mdout["E14"][-1]) == 0

    assert pytest.approx(mdout["VDW"][0], 0.1) == 6.5734
    assert pytest.approx(mdout["Ele"][0], 0.1) == -211.7616
