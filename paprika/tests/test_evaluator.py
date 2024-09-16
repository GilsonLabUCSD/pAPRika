"""
Tests evaluator modules.
"""

import logging
import os
import shutil

import numpy as np
import parmed as pmd
import pytest
import pytraj as pt
from openff.units import unit as openff_unit

from paprika.evaluator import Analyze, Setup
from paprika.evaluator.amber import generate_gaff
from paprika.restraints import DAT_restraint

logger = logging.getLogger(__name__)


@pytest.fixture
def clean_files(directory=os.path.join(os.path.dirname(__file__), "tmp")):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


def test_evaluator_setup_structure(clean_files):
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    G1 = ":BUT@C"
    G2 = ":BUT@C3"

    # Test prepare_complex_structure
    host_guest_pdb = os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    guest_atom_indices = []
    structure = pmd.load_file(host_guest_pdb, structure=True)
    for atom in structure.topology.atoms():
        if atom.residue.name == "BUT":
            guest_atom_indices.append(atom.index)

    host_guest_structure_initial = Setup.prepare_complex_structure(
        host_guest_pdb,
        guest_atom_indices,
        f"{G1} {G2}",
        24.0,
        0,
        46,
    )

    assert all(host_guest_structure_initial[G1].coordinates[0] == [0.0, 0.0, 0.0])

    host_guest_structure_final = Setup.prepare_complex_structure(
        host_guest_pdb,
        guest_atom_indices,
        f"{G1} {G2}",
        24.0,
        45,
        46,
    )
    assert all(host_guest_structure_final[G1].coordinates[0] == [0.0, 0.0, 24.0])

    cG1 = host_guest_structure_final[G1].coordinates[0]
    cG2 = host_guest_structure_final[G2].coordinates[0]
    vec = cG2 - cG1
    axis = np.array([0, 0, 1])
    theta = np.arccos(np.dot(vec, axis) / (np.linalg.norm(vec) * np.linalg.norm(axis)))
    assert theta == 0.0

    # Test prepare_host_structure
    structure = pmd.load_file(host_guest_pdb, structure=True)
    structure[":CB6"].save(os.path.join(temporary_directory, "cb6.pdb"))
    host_pdb = os.path.join(temporary_directory, "cb6.pdb")
    host_atom_indices = []
    for atom in structure.topology.atoms():
        if atom.residue.name == "CB6":
            host_atom_indices.append(atom.index)

    host_structure = Setup.prepare_host_structure(
        host_pdb,
        host_atom_indices,
    )
    center_of_mass = pmd.geometry.center_of_mass(
        host_structure.coordinates, masses=np.ones(len(host_structure.coordinates))
    )
    assert pytest.approx(center_of_mass[0], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[1], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[2], abs=1e-3) == 0.0

    inertia_tensor = np.dot(
        host_structure.coordinates.transpose(), host_structure.coordinates
    )
    eig_val, eig_vec = np.linalg.eig(inertia_tensor)
    assert pytest.approx(eig_vec[0, -1], abs=1e-3) == 0.0
    assert pytest.approx(eig_vec[1, -1], abs=1e-3) == 0.0
    assert pytest.approx(eig_vec[2, -1], abs=1e-3) == 1.0

    # Test add dummy
    Setup.add_dummy_atoms_to_structure(
        host_structure,
        [
            np.array([0, 0, 0]),
            np.array([0, 0, -3.0]),
            np.array([0, 2.2, -5.2]),
        ],
        np.zeros(3),
    )
    dummy_atoms = []
    for atom in host_structure.topology.atoms():
        if atom.name == "DUM":
            dummy_atoms.append(atom)

    assert len(dummy_atoms) == 3
    assert pytest.approx(host_structure[":DM1"].coordinates[0][2], abs=1e-3) == 0.0
    assert pytest.approx(host_structure[":DM2"].coordinates[0][2], abs=1e-3) == -3.0
    assert pytest.approx(host_structure[":DM3"].coordinates[0][1], abs=1e-3) == 2.2
    assert pytest.approx(host_structure[":DM3"].coordinates[0][2], abs=1e-3) == -5.2


def test_evaluator_analyze(clean_files):
    input_pdb = os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    structure = pmd.load_file(input_pdb, structure=True)

    guest_atom_indices = []
    for atom in structure.topology.atoms():
        if atom.residue.name == "BUT":
            guest_atom_indices.append(atom.index)

    host_guest_structure = Setup.prepare_complex_structure(
        input_pdb,
        guest_atom_indices,
        ":BUT@C :BUT@C3",
        24.0,
        0,
        46,
    )
    Setup.add_dummy_atoms_to_structure(
        host_guest_structure,
        [
            np.array([0, 0, 0]),
            np.array([0, 0, -3.0]),
            np.array([0, 2.2, -5.2]),
        ],
        np.zeros(3),
    )

    # Distance restraint
    rest1 = DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = host_guest_structure
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

    # Angle 1 restraint
    rest2 = DAT_restraint()
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

    # Angle 2
    rest3 = DAT_restraint()
    rest3.continuous_apr = True
    rest3.amber_index = True
    rest3.topology = input_pdb
    rest3.mask1 = ":DM1"
    rest3.mask2 = ":BUT@C"
    rest3.mask3 = ":BUT@C3"
    rest3.attach["target"] = 180.0
    rest3.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest3.attach["fc_final"] = 100.0
    rest3.pull["fc"] = rest2.attach["fc_final"]
    rest3.pull["target_initial"] = rest2.attach["target"]
    rest3.pull["target_final"] = rest2.attach["target"]
    rest3.pull["num_windows"] = 19
    rest3.initialize()

    temperature = 298.15
    guest_restraints = [rest1, rest2, rest3]
    ref_state_work = Analyze.compute_ref_state_work(temperature, guest_restraints)
    assert (
        pytest.approx(
            ref_state_work.to(openff_unit.kcal / openff_unit.mole).magnitude, abs=1e-3
        )
        == -7.14151
    )

    fe_sym = Analyze.symmetry_correction(n_microstates=1, temperature=298.15)
    assert fe_sym == 0.0
    fe_sym = Analyze.symmetry_correction(n_microstates=2, temperature=298.15)
    assert pytest.approx(fe_sym, abs=1e-3) == -0.410679


def test_evaluator_gaff(clean_files):
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    butane = os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.mol2")
    resname = "BUT"

    gaff_version = "gaff"
    generate_gaff(
        butane,
        resname,
        output_name="but",
        gaff_version=gaff_version,
        directory_path=temporary_directory,
    )

    structure = pt.iterload(
        os.path.join(temporary_directory, f"but.{gaff_version}.mol2")
    )

    # fmt: off
    butane_atom_type = ["c3", "hc", "hc", "hc", "c3", "hc", "c3", "hc", "hc", "hc", "c3", "hc", "hc", "hc"]
    # fmt: on

    residue_names = []
    for i, atom in enumerate(structure.topology.atoms):
        if atom.resname not in residue_names:
            residue_names.append(atom.resname)

        assert atom.type == butane_atom_type[i]

    assert len(residue_names) == 1
    assert residue_names[0] == resname
