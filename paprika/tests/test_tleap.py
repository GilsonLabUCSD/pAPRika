"""
Tests tleap tools.
"""

import random as random
import shutil

import pytest

from paprika.align import *
from paprika.dummy import *
from paprika.tleap import *


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
    """ Test that we can solvate CB6-BUT using default settings. """
    waters = np.random.randint(100, 10000)
    log.debug("Trying {} waters with default settings...".format(waters))
    sys = System()
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


@pytest.mark.parametrize("shape", ["octahedral", "cubic"])
def test_solvation_shapes(shape, clean_files):
    """ Test that we can solvate CB6-BUT with a truncated octahedron. """
    waters = np.random.randint(1000, 10000)
    log.debug("Trying {} waters in a truncated octahedron...".format(waters))
    sys = System()
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


@pytest.mark.slow
def test_solvation_spatial_size(clean_files):
    """ Test that we can solvate CB6-BUT with an buffer size in Angstroms. """
    random_int = np.random.randint(10, 20)
    random_size = random_int * np.random.random_sample(1) + random_int
    log.debug("Trying buffer size of {} A...".format(random_size[0]))
    sys = System()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/cb6-but.pdb"
    )
    sys.buffer_value = float(random_size[0])
    sys.output_prefix = "solvate"
    sys.pbc_type = "cubic"
    sys.build()
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == sys.target_waters


@pytest.mark.slow
def test_solvation_potassium_control(clean_files):
    """ Test there is no potassium by default. A negative control. """
    waters = np.random.randint(1000, 10000)
    log.debug("Trying {} waters with potassium...".format(waters))
    sys = System()
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
    """ Test that we can solvate CB6-BUT with additional ions. """
    waters = np.random.randint(1000, 10000)
    cations = ["LI", "Na+", "K+", "RB", "CS"]
    anions = ["F", "Cl-", "BR", "IOD"]
    n_cations = np.random.randint(1, 10)
    n_anions = np.random.randint(1, 10)
    random_cation = random.choice(cations)
    random_anion = random.choice(anions)
    log.debug("Trying {} waters with additional ions...".format(waters))
    sys = System()
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
    log.debug("Expecting...")
    log.debug("cation = {}\tn_cations={}".format(random_cation, n_cations))
    log.debug("anion  = {}\t n_anions={}".format(random_anion, n_anions))
    log.debug("Found...")
    log.debug("             n_cations={}".format(cation_number))
    log.debug("              n_anions={}".format(anion_number))

    assert int(cation_number) == n_cations and int(anion_number) == n_anions


def test_solvation_by_M_and_m(clean_files):
    """ Test that we can solvate CB6-BUT through molarity and molality. """
    log.debug("Trying 10 A buffer with 150 mM NaCl...")
    sys = System()
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
    sys.pbc_type = "rectangular"
    sys.add_ions = ["NA", "0.150M", "CL", "0.150M", "K", "0.100m", "BR", "0.100m"]
    sys.build()

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
    calc_num_na = np.ceil((6.022 * 10 ** 23) * (0.150) * volume_in_liters)
    calc_num_cl = np.ceil((6.022 * 10 ** 23) * (0.150) * volume_in_liters)
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
    """ Test that we can solvate CB6-BUT after alignment. """
    cb6 = pmd.load_file(
        os.path.join(
            os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
        )
    )
    zalign(cb6, ":CB6", ":BUT", save=True, filename="./tmp/tmp.pdb")
    waters = np.random.randint(1000, 10000)
    sys = System()
    sys.template_file = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/tleap_solvate.in"
    )
    sys.output_path = "tmp"
    sys.loadpdb_file = "tmp.pdb"
    sys.target_waters = waters
    sys.output_prefix = "solvate"
    sys.build()
    log.debug("Trying {} waters after alignment...".format(waters))
    grepped_waters = sp.check_output(
        ["grep -oh 'WAT' ./tmp/solvate.prmtop | wc -w"], shell=True
    )
    assert int(grepped_waters) == waters


def test_add_dummy(clean_files):
    """ Test that dummy atoms get added correctly """
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")
    host_guest = pmd.load_file(
        os.path.join(
            os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
        ),
        structure=True,
    )
    host_guest = zalign(host_guest, ":BUT@C", ":BUT@C3", save=False)
    host_guest = add_dummy(host_guest, residue_name="DM1", z=-11.000, y=2.000, x=-1.500)

    host_guest.write_pdb(
        os.path.join(temporary_directory, "cb6-but-dum.pdb"), renumber=False
    )
    with open(os.path.join(temporary_directory, "cb6-but-dum.pdb"), "r") as f:
        lines = f.readlines()
        test_line1 = lines[123].rstrip()
        test_line2 = lines[124].rstrip()
    ref_line1 = "TER     123      BUT     2"
    ref_line2 = (
        "HETATM  123 DUM  DM1     3      -1.500   2.000 -11.000  0.00  0.00          PB"
    )
    assert ref_line1 == test_line1
    assert ref_line2 == test_line2

    write_dummy_frcmod(path=temporary_directory)
    write_dummy_mol2(path=temporary_directory, filename="dm1.mol2", residue_name="DM1")
    sys = System()
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
        "loadamberparams dummy.frcmod",
        "DM1 = loadmol2 dm1.mol2",
        "model = loadpdb cb6-but-dum.pdb",
    ]
    sys.output_path = temporary_directory
    sys.output_prefix = "cb6-but-dum"
    sys.pbc_type = None
    sys.neutralize = False
    sys.build()
    with open(
            os.path.join(os.path.dirname(__file__), "../data/cb6-but/REF_cb6-but-dum.rst7"),
            "r",
    ) as f:
        contents = f.read()
        reference = [float(i) for i in contents.split()[2:]]

    with open(os.path.join(temporary_directory, "cb6-but-dum.rst7"), "r") as f:
        contents = f.read()
        new = [float(i) for i in contents.split()[2:]]
    assert np.allclose(reference, new)


def test_hydrogen_mass_repartitioning(clean_files):
    """ Test that hydrogen mass is repartitioned. """
    temporary_directory = os.path.join(os.path.dirname(__file__), "tmp")
    sys = System()
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

    but = pmd.load_file(os.path.join(temporary_directory, sys.output_prefix + ".prmtop"))
    assert np.allclose(but["@H="].atoms[0].mass, 1.008)

    sys.repartition_hydrogen_mass()
    but = pmd.load_file(os.path.join(temporary_directory, sys.output_prefix + ".prmtop"))
    assert np.allclose(but["@H="].atoms[0].mass, 3.024)
