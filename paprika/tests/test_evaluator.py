"""
Tests evaluator modules.
"""
import logging
import os
import shutil

import numpy
import parmed
import pytest
import pytraj
import yaml
from openff.units import unit as openff_unit

from paprika.evaluator import Analyze, Setup
from paprika.evaluator.amber import generate_gaff
from paprika.restraints import DAT_restraint
from paprika.taproom.taproom import read_yaml_schema
from paprika.taproom.utils import convert_string_to_quantity, de_alias

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


@pytest.fixture()
def complex_file():
    complex_pdb = os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")

    butane_molecule = []
    structure = parmed.load_file(complex_pdb, structure=True)
    for atom in structure.topology.atoms():
        if atom.residue.name == "BUT":
            butane_molecule.append(atom.index)

    G1 = ":BUT@C"
    G2 = ":BUT@C3"

    host_guest = Setup.prepare_complex_structure(
        complex_pdb,
        butane_molecule,
        f"{G1} {G2}",
        24.0,
        0,
        46,
    )

    Setup.add_dummy_atoms_to_structure(
        host_guest,
        [
            numpy.array([0, 0, 0]),
            numpy.array([0, 0, -3.0]),
            numpy.array([0, 2.2, -5.2]),
        ],
        numpy.zeros(3),
    )

    return host_guest


@pytest.fixture(scope="module")
def yaml_restraint_schema():
    yaml_file = """name: bam
structure: bam.mol2
complex: a-bam.pdb
net_charge: +1e
aliases:
    - D1: :DM1
    - D2: :DM2
    - D3: :DM3
    - G1: :BAM@C4
    - G2: :BAM@N1
restraints:
  guest:
    - restraint:
        atoms: D1 G1
        attach:
          # During the 'attach' phase, the `force_constant` argument is the
          # final force constant.
          force_constant: 5.0 * kilocalorie / mole / angstrom**2
          target: 6.0 * angstrom
        pull:
          # During the 'pull' phase, the `target` argument is the final value of
          # the restraint.
          force_constant: 5.0 * kilocalorie / mole / angstrom**2
          target: 24.0 * angstrom
    - restraint:
        atoms: D2 D1 G1
        attach:
          force_constant: 100.0 * kilocalorie / mole / radians**2
          target: 180.0 * degrees
        pull:
          force_constant: 100.0 * kilocalorie / mole / radians**2
          target: 180.0 * degrees
    - restraint:
        atoms: D1 G1 G2
        attach:
          force_constant: 100.0 * kilocalorie / mole / radians**2
          target: 180.0 * degrees
        pull:
          force_constant: 100.0 * kilocalorie / mole / radians**2
          target: 180.0 * degrees

  wall_restraints:
    - restraint:
        atoms: ":1@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":2@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":3@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":4@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":5@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":6@O2 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 9.3 * angstrom
    - restraint:
        atoms: ":1@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom
    - restraint:
        atoms: ":2@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom
    - restraint:
        atoms: ":3@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom
    - restraint:
        atoms: ":4@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom
    - restraint:
        atoms: ":5@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom
    - restraint:
        atoms: ":6@O6 G1"
        force_constant: 50.0 * kilocalorie / mole / angstrom**2
        target: 11.3 * angstrom

symmetry_correction:
  restraints:
    - restraint:
        atoms: D2 G1 G2
        force_constant: 200.0 * kilocalorie / mole / radian**2
        target: 91 * degrees
  # Do not attempt to automatically correct for the symmetry restraint by adding -RT \ln (microstates).
  # Instead, we will apply the symmetry restraint, which locks in a particular binding orientation, and then
  # perform separate calculations.
  microstates: 1"""
    return yaml_file


@pytest.fixture(scope="module")
def restraints_schema():
    schema = {
        "static": [
            {
                "atoms": ":DM1 :CB6@O",
                "force_constant": 5.0
                * openff_unit.kcal
                / openff_unit.mol
                / openff_unit.angstrom**2,
            }
        ],
        "conformational": [
            {
                "atoms": ":CB6@O :CB6@O2 :CB6@O4 :CB6@O6",
                "force_constant": 6.0
                * openff_unit.kcal
                / openff_unit.mol
                / openff_unit.radians**2,
                "target": 104.3 * openff_unit.degrees,
            }
        ],
        "symmetry": [
            {
                "atoms": ":DM2 :BUT@C :BUT@C3",
                "force_constant": 50.0
                * openff_unit.kcal
                / openff_unit.mol
                / openff_unit.radians**2,
                "target": 11.0 * openff_unit.degrees,
            }
        ],
        "wall": [
            {
                "atoms": ":CB6@O :BUT@C",
                "force_constant": 50.0
                * openff_unit.kcal
                / openff_unit.mol
                / openff_unit.angstrom**2,
                "target": 11.0 * openff_unit.angstrom,
            }
        ],
        "guest": [
            {
                "atoms": ":DM1 :BUT@C",
                "attach": {
                    "force_constant": 5.0
                    * openff_unit.kcal
                    / openff_unit.mol
                    / openff_unit.angstrom**2,
                    "target": 6.0 * openff_unit.angstrom,
                },
                "pull": {
                    "force_constant": 5.0
                    * openff_unit.kcal
                    / openff_unit.mol
                    / openff_unit.angstrom**2,
                    "target": 24.0 * openff_unit.angstrom,
                },
            }
        ],
    }
    return schema


def test_taproom_yaml(clean_files, yaml_restraint_schema):
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    k_dist = 5.0 * openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    k_wall = 50.0 * openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    k_angle = 100.0 * openff_unit.kcal / openff_unit.mole / openff_unit.radian**2
    r_initial = 6.0 * openff_unit.angstrom
    r_final = 24.0 * openff_unit.angstrom
    angle = 180.0 * openff_unit.degrees

    # Write yaml string to file
    with open(f"{temporary_directory}/guest.yaml", "w") as f:
        f.write(yaml_restraint_schema)

    # Check de_alias
    with open(f"{temporary_directory}/guest.yaml", "r") as f:
        yaml_data = yaml.safe_load(f)

    if "aliases" in yaml_data.keys():
        yaml_data = de_alias(yaml_data)

    distance_restraint = yaml_data["restraints"]["guest"][0]["restraint"]
    theta_restraint = yaml_data["restraints"]["guest"][1]["restraint"]
    beta_restraint = yaml_data["restraints"]["guest"][2]["restraint"]
    symmetry_restraint = yaml_data["symmetry_correction"]["restraints"][0]["restraint"]

    assert distance_restraint["atoms"] == ":DM1 :BAM@C4"
    assert theta_restraint["atoms"] == ":DM2 :DM1 :BAM@C4"
    assert beta_restraint["atoms"] == ":DM1 :BAM@C4 :BAM@N1"
    assert symmetry_restraint["atoms"] == ":DM2 :BAM@C4 :BAM@N1"

    # Check convert string to quantity
    convert_string_to_quantity(yaml_data)

    assert distance_restraint["attach"]["force_constant"] == k_dist
    assert theta_restraint["attach"]["force_constant"] == k_angle
    assert beta_restraint["attach"]["force_constant"] == k_angle

    assert distance_restraint["pull"]["force_constant"] == k_dist
    assert theta_restraint["pull"]["force_constant"] == k_angle
    assert beta_restraint["pull"]["force_constant"] == k_angle

    assert distance_restraint["attach"]["target"] == r_initial
    assert distance_restraint["pull"]["target"] == r_final

    assert theta_restraint["attach"]["target"] == angle
    assert theta_restraint["pull"]["target"] == angle
    assert beta_restraint["attach"]["target"] == angle
    assert beta_restraint["pull"]["target"] == angle

    for i, restraint in enumerate(yaml_data["restraints"]["wall_restraints"]):
        assert restraint["restraint"]["force_constant"] == k_wall
        if i < 6:
            assert restraint["restraint"]["target"] == 9.3 * openff_unit.angstrom
        else:
            assert restraint["restraint"]["target"] == 11.3 * openff_unit.angstrom

    # Check full conversion
    guest_spec = read_yaml_schema(f"{temporary_directory}/guest.yaml")

    distance_restraint = guest_spec["restraints"]["guest"][0]["restraint"]
    theta_restraint = guest_spec["restraints"]["guest"][1]["restraint"]
    beta_restraint = guest_spec["restraints"]["guest"][2]["restraint"]

    assert distance_restraint["attach"]["force_constant"] == k_dist
    assert theta_restraint["attach"]["force_constant"] == k_angle
    assert beta_restraint["attach"]["force_constant"] == k_angle

    assert distance_restraint["pull"]["force_constant"] == k_dist
    assert theta_restraint["pull"]["force_constant"] == k_angle
    assert beta_restraint["pull"]["force_constant"] == k_angle

    assert distance_restraint["attach"]["target"] == r_initial
    assert distance_restraint["pull"]["target"] == r_final

    assert theta_restraint["attach"]["target"] == angle
    assert theta_restraint["pull"]["target"] == angle
    assert beta_restraint["attach"]["target"] == angle
    assert beta_restraint["pull"]["target"] == angle

    for i, restraint in enumerate(guest_spec["restraints"]["wall_restraints"]):
        assert restraint["restraint"]["force_constant"] == k_wall
        if i < 6:
            assert restraint["restraint"]["target"] == 9.3 * openff_unit.angstrom
        else:
            assert restraint["restraint"]["target"] == 11.3 * openff_unit.angstrom


def test_evaluator_setup_structure(clean_files):
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    G1 = ":BUT@C"
    G2 = ":BUT@C3"

    # Test prepare_complex_structure
    host_guest_pdb = os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    guest_atom_indices = []
    structure = parmed.load_file(host_guest_pdb, structure=True)
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
    axis = numpy.array([0, 0, 1])
    theta = numpy.arccos(
        numpy.dot(vec, axis) / (numpy.linalg.norm(vec) * numpy.linalg.norm(axis))
    )
    assert theta == 0.0

    # Test prepare_host_structure
    structure = parmed.load_file(host_guest_pdb, structure=True)
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
    center_of_mass = parmed.geometry.center_of_mass(
        host_structure.coordinates, masses=numpy.ones(len(host_structure.coordinates))
    )
    assert pytest.approx(center_of_mass[0], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[1], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[2], abs=1e-3) == 0.0

    host_structure = Setup.prepare_host_structure(
        host_pdb,
        host_atom_indices=None,
    )
    center_of_mass = parmed.geometry.center_of_mass(
        host_structure.coordinates, masses=numpy.ones(len(host_structure.coordinates))
    )
    assert pytest.approx(center_of_mass[0], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[1], abs=1e-3) == 0.0
    assert pytest.approx(center_of_mass[2], abs=1e-3) == 0.0

    inertia_tensor = numpy.dot(
        host_structure.coordinates.transpose(), host_structure.coordinates
    )
    eig_val, eig_vec = numpy.linalg.eig(inertia_tensor)
    assert pytest.approx(eig_vec[0, -1], abs=1e-3) == 0.0
    assert pytest.approx(eig_vec[1, -1], abs=1e-3) == 0.0
    assert pytest.approx(eig_vec[2, -1], abs=1e-3) == 1.0

    # Test add dummy
    Setup.add_dummy_atoms_to_structure(
        host_structure,
        [
            numpy.array([0, 0, 0]),
            numpy.array([0, 0, -3.0]),
            numpy.array([0, 2.2, -5.2]),
        ],
        numpy.zeros(3),
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


def test_evaluator_setup_restraints(clean_files, complex_file, restraints_schema):
    complex_file_path = os.path.join(os.path.dirname(__file__), "complex.pdb")
    complex_file.save(complex_file_path, overwrite=True)

    attach_lambdas = list(numpy.linspace(0, 1, 15))
    n_attach = 15
    n_pull = 46
    n_release = 15

    static_restraints = Setup.build_static_restraints(
        complex_file_path,
        n_attach,
        n_pull,
        n_release,
        restraints_schema["static"],
        use_amber_indices=False,
    )
    assert static_restraints[0].mask1 == ":DM1"
    assert static_restraints[0].mask2 == ":CB6@O"
    for k in static_restraints[0].phase["attach"]["force_constants"]:
        assert k == restraints_schema["static"][0]["force_constant"]

    conformational_restraint = Setup.build_conformational_restraints(
        complex_file_path,
        attach_lambdas,
        None,
        None,
        restraints_schema["conformational"],
        use_amber_indices=False,
    )
    assert conformational_restraint[0].mask1 == ":CB6@O"
    assert conformational_restraint[0].mask2 == ":CB6@O2"
    assert conformational_restraint[0].mask3 == ":CB6@O4"
    assert conformational_restraint[0].mask4 == ":CB6@O6"
    for i, (k, target) in enumerate(
        zip(
            conformational_restraint[0].phase["attach"]["force_constants"],
            conformational_restraint[0].phase["attach"]["targets"],
        )
    ):
        assert (
            k
            == attach_lambdas[i]
            * restraints_schema["conformational"][0]["force_constant"]
        )
        assert target == restraints_schema["conformational"][0]["target"]

    symmetry_restraints = Setup.build_symmetry_restraints(
        complex_file_path,
        n_attach,
        restraints_schema["symmetry"],
        use_amber_indices=False,
    )
    assert symmetry_restraints[0].mask1 == ":DM2"
    assert symmetry_restraints[0].mask2 == ":BUT@C"
    assert symmetry_restraints[0].mask3 == ":BUT@C3"
    assert symmetry_restraints[0].custom_restraint_values["r1"] is None
    assert symmetry_restraints[0].custom_restraint_values["r2"] is None
    assert (
        symmetry_restraints[0].custom_restraint_values["r3"]
        == 0.0 * openff_unit.degrees
    )
    assert (
        symmetry_restraints[0].custom_restraint_values["r4"]
        == 0.0 * openff_unit.degrees
    )
    assert (
        symmetry_restraints[0].custom_restraint_values["rk2"]
        == restraints_schema["symmetry"][0]["force_constant"]
    )
    assert (
        symmetry_restraints[0].custom_restraint_values["rk3"]
        == restraints_schema["symmetry"][0]["force_constant"]
    )

    wall_restraints = Setup.build_wall_restraints(
        complex_file_path,
        n_attach,
        restraints_schema["wall"],
        use_amber_indices=False,
    )
    assert wall_restraints[0].mask1 == ":CB6@O"
    assert wall_restraints[0].mask2 == ":BUT@C"
    assert (
        wall_restraints[0].custom_restraint_values["r1"] == 0.0 * openff_unit.angstrom
    )
    assert (
        wall_restraints[0].custom_restraint_values["r2"] == 0.0 * openff_unit.angstrom
    )
    assert wall_restraints[0].custom_restraint_values["r3"] is None
    assert wall_restraints[0].custom_restraint_values["r4"] is None
    assert (
        wall_restraints[0].custom_restraint_values["rk2"]
        == restraints_schema["wall"][0]["force_constant"]
    )
    assert (
        wall_restraints[0].custom_restraint_values["rk3"]
        == restraints_schema["wall"][0]["force_constant"]
    )

    guest_restraints = Setup.build_guest_restraints(
        complex_file_path,
        attach_lambdas,
        None,
        restraints_schema["guest"],
        use_amber_indices=False,
    )
    assert guest_restraints[0].mask1 == ":DM1"
    assert guest_restraints[0].mask2 == ":BUT@C"
    for i, (k, target) in enumerate(
        zip(
            guest_restraints[0].phase["attach"]["force_constants"],
            guest_restraints[0].phase["attach"]["targets"],
        )
    ):
        assert (
            k
            == attach_lambdas[i]
            * restraints_schema["guest"][0]["attach"]["force_constant"]
        )
        assert target == restraints_schema["guest"][0]["attach"]["target"]


def test_evaluator_analyze(clean_files):
    # ---------------------------------------------------------------- #
    # Configure system and restraints
    # ---------------------------------------------------------------- #
    input_pdb = os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    structure = parmed.load_file(input_pdb, structure=True)

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
            numpy.array([0, 0, 0]),
            numpy.array([0, 0, -3.0]),
            numpy.array([0, 2.2, -5.2]),
        ],
        numpy.zeros(3),
    )

    angle = 180 * openff_unit.degrees
    k_angle = 100 * openff_unit.kcal / openff_unit.mole / openff_unit.radians**2
    r_initial = 6.0 * openff_unit.angstrom
    r_final = 24.0 * openff_unit.angstrom
    k_r = 5.0 * openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
    attach_fractions = [0.00, 0.04, 0.181, 0.496, 1.000]

    # Distance restraint
    rest1 = DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = host_guest_structure
    rest1.mask1 = ":DM1"
    rest1.mask2 = ":BUT@C"
    rest1.attach["target"] = r_initial
    rest1.attach["fraction_list"] = attach_fractions
    rest1.attach["fc_final"] = k_r
    rest1.pull["fc"] = rest1.attach["fc_final"]
    rest1.pull["target_initial"] = rest1.attach["target"]
    rest1.pull["target_final"] = r_final
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
    rest2.attach["target"] = angle
    rest2.attach["fraction_list"] = attach_fractions
    rest2.attach["fc_final"] = k_angle
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
    rest3.attach["target"] = angle
    rest3.attach["fraction_list"] = attach_fractions
    rest3.attach["fc_final"] = k_angle
    rest3.pull["fc"] = rest2.attach["fc_final"]
    rest3.pull["target_initial"] = rest2.attach["target"]
    rest3.pull["target_final"] = rest2.attach["target"]
    rest3.pull["num_windows"] = 19
    rest3.initialize()

    # Angle 2
    rest4 = DAT_restraint()
    rest4.continuous_apr = True
    rest4.amber_index = True
    rest4.topology = input_pdb
    rest4.mask1 = ":DM1"
    rest4.mask2 = ":BUT@C"
    rest4.mask3 = ":BUT@C2"
    rest4.mask4 = ":BUT@C3"
    rest4.attach["target"] = angle
    rest4.attach["fraction_list"] = attach_fractions
    rest4.attach["fc_final"] = k_angle
    rest4.pull["fc"] = rest2.attach["fc_final"]
    rest4.pull["target_initial"] = rest2.attach["target"]
    rest4.pull["target_final"] = rest2.attach["target"]
    rest4.pull["num_windows"] = 19
    rest4.initialize()

    temperature = 298.15 * openff_unit.kelvin
    guest_restraints = [rest1, rest2, rest3]

    # ---------------------------------------------------------------- #
    # Test reference state work
    # ---------------------------------------------------------------- #
    ref_state_work = Analyze.compute_ref_state_work(temperature, [rest1, rest2, rest3])
    assert (
        pytest.approx(
            ref_state_work.to(openff_unit.kcal / openff_unit.mole).magnitude, abs=1e-3
        )
        == -7.14151
    )
    ref_state_work = Analyze.compute_ref_state_work(
        temperature, [rest1, rest2, rest3, rest4]
    )
    assert (
        pytest.approx(
            ref_state_work.to(openff_unit.kcal / openff_unit.mole).magnitude, abs=1e-3
        )
        == -9.41062
    )

    # ---------------------------------------------------------------- #
    # Test symmetry correction
    # ---------------------------------------------------------------- #
    fe_sym = Analyze.symmetry_correction(n_microstates=1, temperature=298.15)
    assert fe_sym == 0.0
    fe_sym = Analyze.symmetry_correction(n_microstates=2, temperature=298.15)
    assert pytest.approx(fe_sym, abs=1e-3) == -0.410679

    # ---------------------------------------------------------------- #
    # Test APR phase
    # ---------------------------------------------------------------- #
    guest_restraints[0].mask1 = ":CB6@O"
    guest_restraints[1].mask1 = ":CB6@O2"
    guest_restraints[1].mask2 = ":CB6@O"
    guest_restraints[2].mask1 = ":CB6@O"
    for restraint in guest_restraints:
        restraint.initialize()

    apr_data_path = os.path.join(os.path.dirname(__file__), "../data/cb6-but-apr/")
    tmp_path = os.path.join(os.path.dirname(__file__), "tmp")
    window_list = ["a000", "a001", "a002", "a003"] + [f"p{i:03}" for i in range(19)]
    for window in window_list:
        shutil.copytree(f"{apr_data_path}/{window}", f"{tmp_path}/{window}")
        shutil.copy(f"{apr_data_path}/vac.pdb", f"{tmp_path}/{window}/vac.pdb")

    results = Analyze.compute_phase_free_energy(
        phase="attach",
        restraints=guest_restraints,
        windows_directory=tmp_path,
        topology_name="vac.pdb",
        trajectory_mask="*.nc",
        analysis_method="ti-block",
    )

    # loose comparison due to short trajectory
    assert pytest.approx(results["attach"]["ti-block"]["fe"].magnitude, abs=2) == 946
    assert pytest.approx(results["attach"]["ti-block"]["sem"].magnitude, abs=2) == 10


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

    structure = pytraj.iterload(
        os.path.join(temporary_directory, f"but.{gaff_version}.mol2")
    )

    # fmt: off
    butane_atom_type = [
        "c3", "hc", "hc", "hc", "c3", "hc", "c3",
        "hc", "hc", "hc", "c3", "hc", "hc", "hc",
    ]
    # fmt: on

    residue_names = []
    for i, atom in enumerate(structure.topology.atoms):
        if atom.resname not in residue_names:
            residue_names.append(atom.resname)

        assert atom.type == butane_atom_type[i]

    assert len(residue_names) == 1
    assert residue_names[0] == resname
