"""
Tests tleap tools.
"""

import logging
import os
import random as random
import shutil
import subprocess as sp

import numpy as np

try:
    import openmm
    import openmm.unit as openmm_unit
except ImportError:
    import simtk.openmm as openmm
    import simtk.unit as openmm_unit

import parmed as pmd
import pytest

from paprika.build.align import zalign
from paprika.build.system import TLeap
from paprika.build.system.utils import ANGSTROM_CUBED_TO_LITERS, PBCBox
from paprika.utils import is_file_and_not_empty

logger = logging.getLogger(__name__)


@pytest.fixture
def clean_files(directory=os.path.join(os.path.dirname(__file__), "tmp")):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


@pytest.mark.slow
def test_solvation_simple(clean_files):
    """Test that we can solvate CB6-BUT using default settings."""
    waters = np.random.randint(100, 10000)
    logger.debug("Trying {} waters with default settings...".format(waters))
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.build()
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == waters


@pytest.mark.parametrize("shape", [PBCBox.octahedral, PBCBox.cubic])
def test_solvation_shapes(shape, clean_files):
    """Test that we can solvate CB6-BUT with a truncated octahedron."""
    waters = np.random.randint(1000, 10000)
    logger.debug("Trying {} waters in a truncated octahedron...".format(waters))
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.pbc_type = shape
    sys.build()
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == waters


@pytest.mark.parametrize("water_model", ["tip4p", "opc"])
def test_solvation_water_model(water_model, clean_files):
    """Test that we can solvate CB6-BUT with a truncated octahedron."""
    waters = np.random.randint(1000, 10000)
    logger.debug("Trying {} waters in a truncated octahedron...".format(waters))
    sys = TLeap()
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.pbc_type = PBCBox.cubic
    sys.set_water_model(water_model)
    sys.template_lines = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff",
        "loadamberparams ../paprika/data/cb6-but/cb6.frcmod",
        "CB6 = loadmol2 ../paprika/data/cb6-but/cb6.mol2",
        "loadamberparams ../paprika/data/cb6-but/but.frcmod",
        "BUT = loadmol2 ../paprika/data/cb6-but/but.mol2",
        "model = loadpdb ../paprika/data/cb6-but/cb6-but.pdb",
    ]
    sys.build()
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == waters


def test_solvation_bind3p(clean_files):
    logger.debug("Solvating with 500 Bind3P waters in a truncated octahedron...")
    sys = TLeap()
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.target_waters = 2000
    sys.output_prefix = "solvate"
    sys.pbc_type = PBCBox.cubic
    sys.set_water_model("tip3p", model_type="bind3p")
    sys.template_lines = [
        "source leaprc.gaff",
        "loadamberparams ../paprika/data/cb6-but/cb6.frcmod",
        "CB6 = loadmol2 ../paprika/data/cb6-but/cb6.mol2",
        "BUT = loadmol2 ../paprika/data/cb6-but/but.mol2",
        "model = loadpdb ../paprika/data/cb6-but/cb6-but.pdb",
    ]
    sys.build()

    prmtop = os.path.join("./tmp/solvate.prmtop")
    inpcrd = os.path.join("./tmp/solvate.rst7")
    structure = pmd.load_file(prmtop, inpcrd, structure=True)
    water_index = [
        atom.index
        for atom in structure.topology.atoms()
        if atom.residue.name == "WAT" and atom.residue.index == 10
    ]
    system = structure.createSystem()
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    ][0]
    oxygen = nonbonded.getParticleParameters(water_index[0])
    assert (
        pytest.approx(oxygen[1].value_in_unit(openmm_unit.angstrom) - 3.1319, abs=1e-3)
        == 0.0
    )
    assert (
        pytest.approx(
            oxygen[2].value_in_unit(openmm_unit.kilocalorie_per_mole) - 0.1818, abs=1e-3
        )
        == 0.0
    )


@pytest.mark.slow
def test_solvation_spatial_size(clean_files):
    """Test that we can solvate CB6-BUT with a buffer size in Angstroms."""
    random_int = np.random.randint(10, 20)
    random_size = random_int * np.random.random_sample(1) + random_int
    logger.debug("Trying buffer size of {} A...".format(random_size[0]))
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.buffer_value = float(random_size[0])
    sys.output_prefix = "solvate"
    sys.pbc_type = PBCBox.cubic
    sys.build()
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == sys.target_waters


@pytest.mark.slow
def test_solvation_potassium_control(clean_files):
    """Test there is no potassium by default. A negative control."""
    waters = np.random.randint(1000, 10000)
    logger.debug("Trying {} waters with potassium...".format(waters))
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.counter_cation = "K+"
    sys.build()
    potassium = sp.check_output(
        ["grep -oh 'K+' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(potassium) == 0


@pytest.mark.slow
def test_solvation_with_additional_ions(clean_files):
    """Test that we can solvate CB6-BUT with additional ions."""
    waters = np.random.randint(1000, 10000)
    cations = ["LI", "Na+", "K+", "RB", "CS"]
    anions = ["F", "Cl-", "BR", "IOD"]
    n_cations = np.random.randint(1, 10)
    n_anions = np.random.randint(1, 10)
    random_cation = random.choice(cations)
    random_anion = random.choice(anions)
    logger.debug("Trying {} waters with additional ions...".format(waters))
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.neutralize = False
    sys.add_ions = [random_cation, n_cations, random_anion, n_anions]
    sys.build()
    # These should come in the RESIDUE_LABEL region of the prmtop and be before all the water.
    cation_number = sp.check_output(
        [
            "grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | "
            + "grep -oh '{} ' | wc -w".format(random_cation)
        ],
        shell=True,
    )
    anion_number = sp.check_output(
        [
            "grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | "
            + "grep -oh '{} ' | wc -w".format(random_anion)
        ],
        shell=True,
    )
    logger.debug("Expecting...")
    logger.debug("cation = {}\tn_cations={}".format(random_cation, n_cations))
    logger.debug("anion  = {}\t n_anions={}".format(random_anion, n_anions))
    logger.debug("Found...")
    logger.debug("             n_cations={}".format(cation_number))
    logger.debug("              n_anions={}".format(anion_number))

    assert int(cation_number) == n_cations and int(anion_number) == n_anions


def test_solvation_by_M_and_m(clean_files):
    """Test that we can solvate CB6-BUT through molarity and molality."""
    logger.debug("Trying 10 A buffer with 150 mM NaCl...")
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.buffer_value = 10.0
    sys.output_prefix = "solvate"
    sys.neutralize = False
    sys.pbc_type = PBCBox.rectangular
    sys.add_ions = ["NA", "0.150M", "CL", "0.150M", "K", "0.100m", "BR", "0.100m"]
    sys.build(clean_files=False)

    # Molarity Check
    obs_num_na = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | " + "grep -oh 'NA ' | wc -w"],
        shell=True,
    )
    obs_num_cl = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | " + "grep -oh 'CL ' | wc -w"],
        shell=True,
    )

    volume = sys.get_volume()
    volume_in_liters = volume * ANGSTROM_CUBED_TO_LITERS
    calc_num_na = np.ceil((6.022 * 10**23) * 0.150 * volume_in_liters)
    calc_num_cl = np.ceil((6.022 * 10**23) * 0.150 * volume_in_liters)
    assert int(obs_num_na) == int(calc_num_na)
    assert int(obs_num_cl) == int(calc_num_cl)

    # Molality Check
    obs_num_k = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | " + "grep -oh 'K ' | wc -w"],
        shell=True,
    )
    obs_num_br = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./tmp/solvate.prmtop | " + "grep -oh 'BR ' | wc -w"],
        shell=True,
    )
    calc_num_waters = sys.count_residues()["WAT"]
    calc_num_k = np.ceil(0.100 * calc_num_waters * 0.018)
    calc_num_br = np.ceil(0.100 * calc_num_waters * 0.018)
    assert int(obs_num_k) == int(calc_num_k)
    assert int(obs_num_br) == int(calc_num_br)


@pytest.mark.slow
def test_alignment_workflow(clean_files):
    """Test that we can solvate CB6-BUT after alignment."""
    cb6 = pmd.load_file(
        os.path.join(
            os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
        )
    )
    zalign(cb6, ":CB6", ":BUT", save=True, filename="./tmp/tmp.pdb")
    waters = np.random.randint(1000, 10000)
    sys = TLeap()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = "tmp.pdb"
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.build()
    logger.debug("Trying {} waters after alignment...".format(waters))
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == waters


def test_hydrogen_mass_repartitioning(clean_files):
    """Test that hydrogen mass is repartitioned."""
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")
    sys = TLeap()
    but_frcmod = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.frcmod")
    )
    but_mol2 = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.mol2")
    )

    sys.template_lines = [
        "source leaprc.gaff",
        f"loadamberparams {but_frcmod}",
        f"BUT = loadmol2 {but_mol2}",
        f"model = loadmol2 {but_mol2}",
    ]
    sys.output_path = temporary_directory
    sys.output_prefix = "but"
    sys.pbc_type = None
    sys.neutralize = False
    sys.build()

    but = pmd.load_file(
        os.path.join(temporary_directory, sys.output_prefix + ".prmtop")
    )
    assert np.allclose(but["@H="].atoms[0].mass, 1.008)

    sys.repartition_hydrogen_mass()
    but = pmd.load_file(
        os.path.join(temporary_directory, sys.output_prefix + ".prmtop")
    )
    assert np.allclose(but["@H="].atoms[0].mass, 3.024)


def test_multiple_pdb_files(clean_files):
    """
    Test that multiple `loadpdb` lines are carried through.
    Reference: https://github.com/GilsonLabUCSD/pAPRika/issues/141
    """
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")
    sys = TLeap()
    cb6_frcmod = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6.frcmod")
    )
    cb6_mol2 = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6.mol2")
    )
    but_frcmod = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.frcmod")
    )
    but_mol2 = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.mol2")
    )

    sys.template_lines = [
        "source leaprc.gaff",
        f"loadamberparams {cb6_frcmod}",
        f"CB6 = loadmol2 {cb6_mol2}",
        f"loadamberparams {but_frcmod}",
        f"BUT = loadmol2 {but_mol2}",
        "a = loadpdb cb6-but-notcentered.pdb",
        "b = loadpdb cb6-but-minimized.pdb",
        "model = combine { a b }",
    ]
    sys.output_path = temporary_directory
    sys.output_prefix = "multi"
    sys.pbc_type = None
    sys.neutralize = False
    sys.build(clean_files=False)

    with open(f"{temporary_directory}/multi.tleap.in", "r") as f:
        lines = f.read()

    assert "a = loadpdb cb6-but-notcentered.pdb" in lines
    assert "b = loadpdb cb6-but-minimized.pdb" in lines
    assert "combine" in lines
    assert "savepdb model multi.pdb" in lines


def test_conversions(clean_files):
    """Test that conversion methods work in Tleap module."""
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    but_frcmod = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.frcmod")
    )
    but_mol2 = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/but.mol2")
    )

    sys = TLeap()
    sys.template_lines = [
        "source leaprc.gaff",
        f"loadamberparams {but_frcmod}",
        f"BUT = loadmol2 {but_mol2}",
        f"model = loadmol2 {but_mol2}",
    ]
    sys.output_path = temporary_directory
    sys.output_prefix = "but"
    sys.pbc_type = None
    sys.neutralize = False
    sys.build()

    from paprika.build.system.utils import ConversionToolkit
    from paprika.utils import is_file_and_not_empty

    # Gromacs ParmEd test
    sys.convert_to_gromacs(overwrite=True)
    top_file = os.path.join(sys.output_path, sys.output_prefix + ".top")
    gro_file = os.path.join(sys.output_path, sys.output_prefix + ".gro")
    assert is_file_and_not_empty(top_file) is True
    assert is_file_and_not_empty(gro_file) is True

    # Gromacs InterMol test
    sys.convert_to_gromacs(toolkit=ConversionToolkit.InterMol)
    assert is_file_and_not_empty(top_file) is True
    assert is_file_and_not_empty(gro_file) is True
    top_file = os.path.join(sys.output_path, sys.output_prefix + "_from_amber.top")
    gro_file = os.path.join(sys.output_path, sys.output_prefix + "_from_amber.gro")
    assert is_file_and_not_empty(top_file) is True
    assert is_file_and_not_empty(gro_file) is True

    # CHARMM ParmEd test
    sys.convert_to_charmm()
    psf_file = os.path.join(sys.output_path, sys.output_prefix + ".psf")
    crd_file = os.path.join(sys.output_path, sys.output_prefix + ".crd")

    assert is_file_and_not_empty(psf_file) is True
    assert is_file_and_not_empty(crd_file) is True

    # CHARMM Intermol test
    sys.convert_to_charmm(toolkit=ConversionToolkit.InterMol)
    rtf_file = os.path.join(sys.output_path, sys.output_prefix + ".rtf")
    prm_file = os.path.join(sys.output_path, sys.output_prefix + ".prm")

    assert is_file_and_not_empty(psf_file) is True
    assert is_file_and_not_empty(crd_file) is True
    assert is_file_and_not_empty(rtf_file) is True
    assert is_file_and_not_empty(prm_file) is True

    # LAMMPS test
    sys.convert_to_lammps()
    input_file = os.path.join(sys.output_path, sys.output_prefix + ".input")
    lmp_file = os.path.join(sys.output_path, sys.output_prefix + ".lmp")
    assert is_file_and_not_empty(input_file) is True
    assert is_file_and_not_empty(lmp_file) is True

    # DESMOND test
    sys.convert_to_desmond()
    cms_file = os.path.join(sys.output_path, sys.output_prefix + ".cms")
    assert is_file_and_not_empty(cms_file) is True


def test_ion_solvation(clean_files):
    """Test that tleap module is solvating an ion properly."""
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")

    ions = ["Na+", "Cl-"]
    for ion in ions:
        # Write PDB file
        with open(os.path.join(temporary_directory, "ion.pdb"), "w") as f:
            f.writelines(
                "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n"
            )
            f.writelines(
                f"ATOM      1 {ion:3s}  {ion:3s}     1       0.000   0.000   0.000  1.00  0.00\n"
            )
            f.writelines("END")

        # Write MOL2 file
        with open(os.path.join(temporary_directory, "ion.mol2"), "w") as f:
            charge = 1.0 if "+" in ion else -1.0
            gaff_type = ion.split("+" if "+" in ion else "-")[0].lower()
            f.writelines("@<TRIPOS>MOLECULE\n")
            f.writelines(f"{ion}\n")
            f.writelines("1     0     1     0     0\n")
            f.writelines("SMALL\n")
            f.writelines("No Charge or Current Charge\n\n\n")
            f.writelines("@<TRIPOS>ATOM\n")
            f.writelines(
                f"      1 {ion:3s}          0.0000     0.0000     0.0000 {gaff_type}         1 {ion}       {charge:.6f}\n"
            )
            f.writelines("@<TRIPOS>BOND\n")
            f.writelines("@<TRIPOS>SUBSTRUCTURE\n")
            f.writelines(
                f"     1 {ion:3s}         1 TEMP              0 ****  ****    0 ROOT"
            )

        # Run tleap
        ion_sys = TLeap()
        ion_sys.output_path = temporary_directory
        ion_sys.target_waters = 1000
        ion_sys.pbc_type = PBCBox.cubic
        ion_sys.neutralize = False
        ion_sys.template_lines = [
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
            f"{ion} = loadmol2 ion.mol2",
            "model = loadpdb ion.pdb",
        ]
        ion_sys.build()

        # Make sure Amber files are generated
        assert (
            is_file_and_not_empty(os.path.join(temporary_directory, "build.prmtop"))
            is True
        )
        assert (
            is_file_and_not_empty(os.path.join(temporary_directory, "build.rst7"))
            is True
        )

        # Check number of atoms in files
        structure = pmd.load_file(
            os.path.join(temporary_directory, "build.prmtop"),
            os.path.join(temporary_directory, "build.rst7"),
        )
        assert len(structure.atoms) == 3001
