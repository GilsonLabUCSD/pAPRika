import os
import shutil

import parmed as pmd
import pytest

from paprika import restraints
from paprika.restraints import amber
from paprika.restraints.restraints import create_window_list
from paprika.simulate import AMBER, GROMACS, NAMD
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

    gbsim = AMBER()
    gbsim.path = os.path.join("tmp", "k-cl", "windows", "a003")
    gbsim.executable = "sander"
    gbsim.topology = "k-cl.prmtop"
    gbsim.prefix = "minimize"
    gbsim.coordinates = "k-cl.rst7"
    gbsim.config_gb_min()
    gbsim.cntrl["maxcyc"] = 1
    gbsim.cntrl["ncyc"] = 1
    gbsim.run()

    gbsim.config_gb_md()
    gbsim.prefix = "md"
    gbsim.coordinates = "minimize.rst7"
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
    simulation = AMBER()
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
    simulation.coordinates = "k-cl.rst7"

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


def test_gromacs_config(clean_files):
    # Test PBC
    simulation = GROMACS()
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_pbc_md(
        ensemble=GROMACS.Ensemble.NPT,
        integrator=GROMACS.Integrator.VelocityVerlet,
        thermostat=GROMACS.Thermostat.NoseHoover,
        barostat=GROMACS.Barostat.ParrinelloRahman,
    )
    simulation.control["nsteps"] = 1250000
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".mdp"))
    lines = f.readlines()

    # Check if critical keys are written
    assert "nsteps" in "".join(lines)
    assert "tcoupl" in "".join(lines)
    assert "pcoupl" in "".join(lines)

    # Check if the values are correct
    for line in lines:
        if line.startswith("nsteps"):
            assert int(line.split()[-1]) == 1250000

        elif line.startswith("tcoupl"):
            assert line.split()[-1] == "nose-hoover"

        elif line.startswith("pcoupl") and not line.startswith("pcoupltype"):
            assert line.split()[-1] == "Parrinello-Rahman"

        elif line.startswith("pcoupltype"):
            assert line.split()[-1] == "isotropic"

        elif line.startswith("integrator"):
            assert line.split()[-1] == "md-vv"

    # Test for tc-groups
    simulation.tc_groups = ["CB8", "AMT", "HOH"]
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".mdp"))
    lines = f.readlines()

    for line in lines:
        if line.startswith("tc-grps"):
            assert line.split()[2] == "CB8"
            assert line.split()[3] == "AMT"
            assert line.split()[4] == "HOH"

        elif line.startswith("tau-t"):
            assert len(line.split()) == 5

        elif line.startswith("ref-t"):
            assert len(line.split()) == 5

    # Test vacuum
    simulation = GROMACS()
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_vac_md(
        integrator=GROMACS.Integrator.LangevinDynamics,
        thermostat=GROMACS.Thermostat.VelocityRescaling,
    )
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".mdp"))
    lines = f.readlines()

    assert "tcoupl" not in "".join(lines)

    for line in lines:
        if line.startswith("integrator"):
            assert line.split()[-1] == "sd"

        elif line.startswith("tau-t"):
            assert float(line.split()[-1]) == 0.1

        elif line.startswith("rcoulomb"):
            assert float(line.split()[-1]) == 333.3

        elif line.startswith("rvdw"):
            assert float(line.split()[-1]) == 333.3

        elif line.startswith("DispCorr"):
            assert line.split()[-1] == "no"

    # Test min
    simulation = GROMACS()
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_pbc_min(
        optimizer=GROMACS.Optimizer.ConjugateGradient,
    )
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".mdp"))
    lines = f.readlines()

    for line in lines:
        if line.startswith("integrator"):
            assert line.split()[-1] == "cg"


def test_namd_config(clean_files):
    # Test PBC
    simulation = NAMD()
    simulation.topology = "k-cl.prmtop"
    simulation.coordinates = "k-cl.rst7"
    simulation.checkpoint = "equilibration"
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_pbc_md(
        ensemble=NAMD.Ensemble.NPT,
        thermostat=NAMD.Thermostat.LoweAnderson,
        barostat=NAMD.Barostat.Berendsen,
    )
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".conf"))
    lines = f.readlines()

    assert "ambercoor" in "".join(lines)
    assert "parmfile" in "".join(lines)
    assert "readexclusions" in "".join(lines)
    assert "scnb" in "".join(lines)
    assert "bincoordinates" in "".join(lines)
    assert "binvelocities" in "".join(lines)
    assert "extendedSystem" in "".join(lines)
    assert "loweAndersen" in "".join(lines)

    for line in lines:
        if line.startswith("amber") and not line.startswith("ambercoor"):
            assert line.split()[-1] == "yes"

        elif line.startswith("readexclusions"):
            assert line.split()[-1] == "yes"

        elif line.startswith("scnb"):
            assert float(line.split()[-1]) == 2.0

        elif line.startswith("loweAndersenCutoff"):
            assert float(line.split()[-1]) == 2.7

        elif line.startswith("BerendsenPressureCompressibility"):
            assert float(line.split()[-1]) == 4.57e-5

        elif line.startswith("bincoordinates"):
            assert line.split()[-1] == simulation.checkpoint+".coor"

        elif line.startswith("binvelocities"):
            assert line.split()[-1] == simulation.checkpoint+".vel"

        elif line.startswith("extendedSystem"):
            assert line.split()[-1] == simulation.checkpoint+".xsc"

    # Test vac
    simulation = NAMD()
    simulation.topology = "k-cl.prmtop"
    simulation.coordinates = "k-cl.rst7"
    simulation.checkpoint = "equilibration"
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_vac_md(
        thermostat=NAMD.Thermostat.Langevin,
    )
    simulation.control["run"] = 500000
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".conf"))
    lines = f.readlines()

    assert "GBIS" in "".join(lines)

    for line in lines:
        if line.startswith("PME"):
            assert line.split()[-1] == "off"

        elif line.startswith("GBIS"):
            assert line.split()[-1] == "off"

        elif line.startswith("cutoff"):
            assert float(line.split()[-1]) == 999.0

        elif line.startswith("run"):
            assert int(line.split()[-1]) == 500000

    # Test GB
    simulation = NAMD()
    simulation.topology = "k-cl.prmtop"
    simulation.coordinates = "k-cl.rst7"
    simulation.checkpoint = "equilibration"
    simulation.path = os.path.join("tmp")
    simulation.prefix = "test"
    simulation.config_gb_md(
        thermostat=NAMD.Thermostat.Langevin,
    )
    simulation.implicit["ionConcentration"] = 0.15
    simulation._write_input_file()

    f = open(os.path.join(simulation.path, simulation.prefix + ".conf"))
    lines = f.readlines()

    assert "GBIS" in "".join(lines)
    assert "solventDielectric" in "".join(lines)
    assert "ionConcentration" in "".join(lines)

    for line in lines:
        if line.startswith("PME"):
            assert line.split()[-1] == "off"

        elif line.startswith("GBIS"):
            assert line.split()[-1] == "on"

        elif line.startswith("ionConcentration"):
            assert float(line.split()[-1]) == 0.15

        elif line.startswith("cutoff"):
            assert float(line.split()[-1]) == 999.0

        elif line.startswith("run"):
            assert int(line.split()[-1]) == 5000
