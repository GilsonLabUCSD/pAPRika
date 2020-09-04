"""
Tests the dummy atoms utilities.
"""

import pytest
import os
import parmed as pmd
from paprika.align import zalign
import paprika.dummy as dummy


def test_extract_dummy():
    cb6 = pmd.load_file(os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.pdb"))
    aligned_cb6 = zalign(cb6, ":CB6", ":BUT")
    aligned_cb6 = dummy.add_dummy(aligned_cb6, residue_name="DM1", z=-6.0)
    aligned_cb6 = dummy.add_dummy(aligned_cb6, residue_name="DM2", z=-9.0)
    aligned_cb6 = dummy.add_dummy(aligned_cb6, residue_name="DM3", z=-11.2, y=2.2)

    dm_index = []
    for atom in aligned_cb6.topology.atoms():
        if atom.residue.name in ["DM1", "DM2", "DM3"]:
            dm_index.append(atom.index)

    assert (len(dm_index) != 0)

    dummy_dict = dummy.extract_dummy_atoms(aligned_cb6, serial=False)

    assert (all(dummy_dict["DM1"]["pos"] == [0.0, 0.0, -6.0]))
    assert (all(dummy_dict["DM2"]["pos"] == [0.0, 0.0, -9.0]))
    assert (all(dummy_dict["DM3"]["pos"] == [0.0, 2.2, -11.2]))

    assert (dummy_dict["DM1"]["mass"] == 208.0)
    assert (dummy_dict["DM2"]["mass"] == 208.0)
    assert (dummy_dict["DM3"]["mass"] == 208.0)

    assert (dummy_dict["DM1"]["idx"] == dm_index[0])
    assert (dummy_dict["DM2"]["idx"] == dm_index[1])
    assert (dummy_dict["DM3"]["idx"] == dm_index[2])
