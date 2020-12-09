"""
Tests the dummy atoms utilities.
"""

import os
import shutil

import numpy as np
import parmed as pmd
import pytest

from paprika.build.align import zalign
from paprika.build.dummy import (
    add_dummy,
    extract_dummy_atoms,
    write_dummy_frcmod,
    write_dummy_mol2,
)
from paprika.build.system import TLeap


@pytest.fixture
def clean_files(directory=os.path.join(os.path.dirname(__file__), "tmp")):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


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


def test_extract_dummy():
    cb6 = pmd.load_file(
        os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb")
    )
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    aligned_cb6 = add_dummy(aligned_cb6, residue_name="DM1", z=-6.0)
    aligned_cb6 = add_dummy(aligned_cb6, residue_name="DM2", z=-9.0)
    aligned_cb6 = add_dummy(aligned_cb6, residue_name="DM3", z=-11.2, y=2.2)

    dm_index = []
    for atom in aligned_cb6.topology.atoms():
        if atom.residue.name in ["DM1", "DM2", "DM3"]:
            dm_index.append(atom.index)

    assert len(dm_index) != 0

    dummy_dict = extract_dummy_atoms(aligned_cb6, serial=False)

    assert all(dummy_dict["DM1"]["pos"] == [0.0, 0.0, -6.0])
    assert all(dummy_dict["DM2"]["pos"] == [0.0, 0.0, -9.0])
    assert all(dummy_dict["DM3"]["pos"] == [0.0, 2.2, -11.2])

    assert dummy_dict["DM1"]["mass"] == 208.0
    assert dummy_dict["DM2"]["mass"] == 208.0
    assert dummy_dict["DM3"]["mass"] == 208.0

    assert dummy_dict["DM1"]["idx"] == dm_index[0]
    assert dummy_dict["DM2"]["idx"] == dm_index[1]
    assert dummy_dict["DM3"]["idx"] == dm_index[2]
