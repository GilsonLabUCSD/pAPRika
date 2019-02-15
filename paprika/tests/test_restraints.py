"""
Tests the restraints utilities.
"""

import pytest

from paprika.restraints import *


def test_DAT_restraint():
    # Method 1
    log.info("### Testing restraint 1, Method 1")
    rest1 = DAT_restraint()
    rest1.amber_index = True
    rest1.continuous_apr = False
    rest1.auto_apr = False
    rest1.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest1.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest1.mask2 = ":BUT@C3"
    rest1.attach["target"] = 3.0
    rest1.attach["num_windows"] = 4
    rest1.attach["fc_initial"] = 0.0
    rest1.attach["fc_final"] = 3.0
    rest1.pull["fc"] = rest1.attach["fc_final"]
    rest1.pull["num_windows"] = 4
    rest1.pull["target_initial"] = rest1.attach["target"]
    rest1.pull["target_final"] = 6.0
    rest1.release["target"] = rest1.pull["target_final"]
    rest1.release["num_windows"] = rest1.attach["num_windows"]
    rest1.release["fc_initial"] = rest1.attach["fc_initial"]
    rest1.release["fc_final"] = rest1.attach["fc_final"]
    rest1.initialize()
    assert rest1.index1 == [13, 31, 49, 67, 85, 103]
    assert rest1.index2 == [119]
    assert rest1.index3 == None
    assert rest1.index4 == None
    assert np.allclose(
        rest1.phase["attach"]["force_constants"], np.array([0.0, 1.0, 2.0, 3.0])
    )
    assert np.allclose(rest1.phase["attach"]["targets"], np.array([3.0, 3.0, 3.0, 3.0]))
    assert np.allclose(
        rest1.phase["pull"]["force_constants"], np.array([3.0, 3.0, 3.0, 3.0])
    )
    assert np.allclose(rest1.phase["pull"]["targets"], np.array([3.0, 4.0, 5.0, 6.0]))
    assert np.allclose(
        rest1.phase["release"]["force_constants"], np.array([0.0, 1.0, 2.0, 3.0])
    )
    assert np.allclose(
        rest1.phase["release"]["targets"], np.array([6.0, 6.0, 6.0, 6.0])
    )
    window_list = create_window_list([rest1])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 1a
    log.info("### Testing restraint 2, Method 1a")
    rest2 = DAT_restraint()
    rest2.amber_index = True
    rest2.continuous_apr = False
    rest2.auto_apr = False
    rest2.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest2.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest2.mask2 = ":BUT@C3"
    rest2.mask3 = ":BUT@C"
    rest2.attach["target"] = 180.0
    rest2.attach["num_windows"] = 4
    rest2.attach["fc_final"] = 75.0
    rest2.pull["fc"] = rest2.attach["fc_final"]
    rest2.pull["num_windows"] = 4
    rest2.pull["target_final"] = 180.0
    rest2.release["target"] = rest2.pull["target_final"]
    rest2.release["num_windows"] = rest2.attach["num_windows"]
    rest2.release["fc_final"] = rest2.attach["fc_final"]
    rest2.initialize()
    assert rest2.index1 == [13, 31, 49, 67, 85, 103]
    assert rest2.index2 == [119]
    assert rest2.index3 == [109]
    assert rest2.index4 == None
    assert np.allclose(
        rest2.phase["attach"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(
        rest2.phase["attach"]["targets"], np.array([180.0, 180.0, 180.0, 180.0])
    )
    assert np.allclose(
        rest2.phase["pull"]["force_constants"], np.array([75.0, 75.0, 75.0, 75.0])
    )
    assert np.allclose(
        rest2.phase["pull"]["targets"], np.array([0.0, 60.0, 120.0, 180.0])
    )
    assert np.allclose(
        rest2.phase["release"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(
        rest2.phase["release"]["targets"], np.array([180.0, 180.0, 180.0, 180.0])
    )
    window_list = create_window_list([rest2])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 2 (Note auto_apr = True)
    log.info("### Testing restraint 3, Method 2")
    rest3 = DAT_restraint()
    rest3.amber_index = True
    rest3.continuous_apr = False
    rest3.auto_apr = True
    rest3.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest3.mask1 = ":CB6@O2"
    rest3.mask2 = ":CB6@O"
    rest3.mask3 = ":BUT@C3"
    rest3.mask4 = ":BUT@C"
    rest3.attach["target"] = 90.0
    rest3.attach["fc_increment"] = 25.0
    rest3.attach["fc_initial"] = 0.0
    rest3.attach["fc_final"] = 75.0
    rest3.pull["target_increment"] = 1.0
    rest3.pull["target_final"] = 93.0
    rest3.release["fc_final"] = 75.0
    rest3.initialize()
    assert rest3.index1 == [31]
    assert rest3.index2 == [13]
    assert rest3.index3 == [119]
    assert rest3.index4 == [109]
    assert np.allclose(
        rest3.phase["attach"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(
        rest3.phase["attach"]["targets"], np.array([90.0, 90.0, 90.0, 90.0])
    )
    assert np.allclose(
        rest3.phase["pull"]["force_constants"], np.array([75.0, 75.0, 75.0, 75.0])
    )
    assert np.allclose(
        rest3.phase["pull"]["targets"], np.array([90.0, 91.0, 92.0, 93.0])
    )
    assert np.allclose(
        rest3.phase["release"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(
        rest3.phase["release"]["targets"], np.array([93.0, 93.0, 93.0, 93.0])
    )
    window_list = create_window_list([rest3])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 2a
    log.info("### Testing restraint 4, Method 2a")
    rest4 = DAT_restraint()
    rest4.amber_index = True
    rest4.continuous_apr = False
    rest4.auto_apr = False
    rest4.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest4.mask1 = ":CB6@O2"
    rest4.mask2 = ":CB6@O"
    rest4.mask3 = ":BUT@C3"
    rest4.mask4 = ":BUT@C"
    rest4.attach["target"] = 0.0
    rest4.attach["fc_increment"] = 25.0
    rest4.attach["fc_final"] = 75.0
    rest4.pull["fc"] = 75.0
    rest4.pull["target_increment"] = 1.0
    rest4.pull["target_final"] = 3.0
    rest4.release["target"] = 3.0
    rest4.release["fc_increment"] = 25.0
    rest4.release["fc_final"] = 75.0
    rest4.initialize()
    assert rest4.index1 == [31]
    assert rest4.index2 == [13]
    assert rest4.index3 == [119]
    assert rest4.index4 == [109]
    assert np.allclose(
        rest4.phase["attach"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(rest4.phase["attach"]["targets"], np.array([0.0, 0.0, 0.0, 0.0]))
    assert np.allclose(
        rest4.phase["pull"]["force_constants"], np.array([75.0, 75.0, 75.0, 75.0])
    )
    assert np.allclose(rest4.phase["pull"]["targets"], np.array([0.0, 1.0, 2.0, 3.0]))
    assert np.allclose(
        rest4.phase["release"]["force_constants"], np.array([0.0, 25.0, 50.0, 75.0])
    )
    assert np.allclose(
        rest4.phase["release"]["targets"], np.array([3.0, 3.0, 3.0, 3.0])
    )
    window_list = create_window_list([rest4])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "p003",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 3
    log.info("### Testing restraint 5, Method 3")
    rest5 = DAT_restraint()
    rest5.amber_index = True
    rest5.continuous_apr = False
    rest5.auto_apr = False
    rest5.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest5.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest5.mask2 = ":BUT@C*"
    rest5.attach["target"] = 0.0
    rest5.attach["fraction_list"] = [0.0, 0.2, 0.5, 1.0]
    rest5.attach["fc_final"] = 5.0
    rest5.pull["fc"] = rest5.attach["fc_final"]
    rest5.pull["fraction_list"] = [0.0, 0.5, 1.0]
    rest5.pull["target_final"] = 1.0
    rest5.release["target"] = rest5.pull["target_final"]
    rest5.release["fraction_list"] = [0.0, 0.3, 0.6, 1.0]
    rest5.release["fc_final"] = rest5.attach["fc_final"]
    rest5.initialize()
    assert rest5.index1 == [13, 31, 49, 67, 85, 103]
    assert rest5.index2 == [109, 113, 115, 119]
    assert rest5.index3 == None
    assert rest5.index4 == None
    assert np.allclose(
        rest5.phase["attach"]["force_constants"], np.array([0.0, 1.0, 2.5, 5.0])
    )
    assert np.allclose(rest5.phase["attach"]["targets"], np.array([0.0, 0.0, 0.0, 0.0]))
    assert np.allclose(
        rest5.phase["pull"]["force_constants"], np.array([5.0, 5.0, 5.0])
    )
    assert np.allclose(rest5.phase["pull"]["targets"], np.array([0.0, 0.5, 1.0]))
    assert np.allclose(
        rest5.phase["release"]["force_constants"], np.array([0.0, 1.5, 3.0, 5.0])
    )
    assert np.allclose(
        rest5.phase["release"]["targets"], np.array([1.0, 1.0, 1.0, 1.0])
    )
    window_list = create_window_list([rest5])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "p000",
        "p001",
        "p002",
        "r000",
        "r001",
        "r002",
        "r003",
    ]

    # Method 4
    log.info("### Testing restraint 6, Method 4")
    rest6 = DAT_restraint()
    rest6.amber_index = True
    rest6.continuous_apr = False
    rest6.auto_apr = False
    rest6.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest6.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
    rest6.mask2 = ":BUT@C*"
    rest6.attach["target"] = 0.0
    rest6.attach["fraction_increment"] = 0.25
    rest6.attach["fc_final"] = 5.0
    rest6.pull["fc"] = rest6.attach["fc_final"]
    rest6.pull["fraction_increment"] = 0.5
    rest6.pull["target_final"] = 1.0
    rest6.release["target"] = rest6.pull["target_final"]
    rest6.release["fraction_increment"] = 0.33
    rest6.release["fc_final"] = rest6.attach["fc_final"]
    rest6.initialize()
    assert rest6.index1 == [13, 31, 49, 67, 85, 103]
    assert rest6.index2 == [109, 113, 115, 119]
    assert rest6.index3 == None
    assert rest6.index4 == None
    assert np.allclose(
        rest6.phase["attach"]["force_constants"], np.array([0.0, 1.25, 2.5, 3.75, 5.0])
    )
    assert np.allclose(
        rest6.phase["attach"]["targets"], np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    )
    assert np.allclose(
        rest6.phase["pull"]["force_constants"], np.array([5.0, 5.0, 5.0])
    )
    assert np.allclose(rest6.phase["pull"]["targets"], np.array([0.0, 0.5, 1.0]))
    ### Note, the 6.6 in the following test is wrong ... needs to get fixed.
    assert np.allclose(
        rest6.phase["release"]["force_constants"], np.array([0.0, 1.65, 3.3, 4.95, 6.6])
    )
    assert np.allclose(
        rest6.phase["release"]["targets"], np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    )
    window_list = create_window_list([rest6])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "a003",
        "a004",
        "p000",
        "p001",
        "p002",
        "r000",
        "r001",
        "r002",
        "r003",
        "r004",
    ]

    # Method 5 (Note continuous_apr = True)
    log.info("### Testing restraint 7, Method 5")
    rest7 = DAT_restraint()
    rest7.amber_index = True
    rest7.continuous_apr = True
    rest7.auto_apr = False
    rest7.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest7.mask1 = ":1@O,O1,:BUT@H1"
    rest7.mask2 = ":CB6@N"
    rest7.attach["target"] = 0.0
    rest7.attach["fc_list"] = [0.0, 0.5, 1.0, 2.0]
    rest7.pull["fc"] = 2.0
    rest7.pull["target_list"] = [0.0, 0.5, 1.0, 1.5]
    rest7.release["target"] = 1.5
    rest7.release["fc_list"] = [0.0, 0.66, 1.2, 2.0]
    rest7.initialize()
    assert rest7.index1 == [13, 14, 111]
    assert rest7.index2 == [3]
    assert rest7.index3 == None
    assert rest7.index4 == None
    assert np.allclose(
        rest7.phase["attach"]["force_constants"], np.array([0.0, 0.5, 1.0, 2.0])
    )
    assert np.allclose(rest7.phase["attach"]["targets"], np.array([0.0, 0.0, 0.0, 0.0]))
    assert np.allclose(
        rest7.phase["pull"]["force_constants"], np.array([2.0, 2.0, 2.0, 2.0])
    )
    assert np.allclose(rest7.phase["pull"]["targets"], np.array([0.0, 0.5, 1.0, 1.5]))
    assert np.allclose(
        rest7.phase["release"]["force_constants"], np.array([0.0, 0.66, 1.2, 2.0])
    )
    assert np.allclose(
        rest7.phase["release"]["targets"], np.array([1.5, 1.5, 1.5, 1.5])
    )
    window_list = create_window_list([rest7])
    assert window_list == [
        "a000",
        "a001",
        "a002",
        "p000",
        "p001",
        "p002",
        "p003",
        "r001",
        "r002",
        "r003",
    ]

    # Just Attach
    log.info("### Testing restraint 8, just attach")
    rest8 = DAT_restraint()
    rest8.amber_index = True
    rest8.continuous_apr = False
    rest8.auto_apr = False
    rest8.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest8.mask1 = ":CB6@O"
    rest8.mask2 = ":BUT@C3"
    rest8.attach["target"] = 0.0
    rest8.attach["num_windows"] = 4
    rest8.attach["fc_initial"] = 0.0
    rest8.attach["fc_final"] = 3.0
    rest8.initialize()
    assert rest8.index1 == [13]
    assert rest8.index2 == [119]
    assert rest8.index3 == None
    assert rest8.index4 == None
    assert np.allclose(
        rest8.phase["attach"]["force_constants"], np.array([0.0, 1.0, 2.0, 3.0])
    )
    assert np.allclose(rest8.phase["attach"]["targets"], np.array([0.0, 0.0, 0.0, 0.0]))
    assert rest8.phase["pull"]["force_constants"] == None
    assert rest8.phase["pull"]["targets"] == None
    assert rest8.phase["release"]["force_constants"] == None
    assert rest8.phase["release"]["targets"] == None
    window_list = create_window_list([rest8])
    assert window_list == ["a000", "a001", "a002", "a003"]

    # Just Pull
    log.info("### Testing restraint 9, just pull")
    rest9 = DAT_restraint()
    rest9.amber_index = True
    rest9.continuous_apr = False
    rest9.auto_apr = False
    rest9.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest9.mask1 = ":CB6@O"
    rest9.mask2 = ":BUT@C3"
    rest9.pull["fc"] = 3.0
    rest9.pull["num_windows"] = 4
    rest9.pull["target_initial"] = 0.0
    rest9.pull["target_final"] = 3.0
    rest9.initialize()
    assert rest9.index1 == [13]
    assert rest9.index2 == [119]
    assert rest9.index3 == None
    assert rest9.index4 == None
    assert rest9.phase["attach"]["force_constants"] == None
    assert rest9.phase["attach"]["targets"] == None
    assert np.allclose(
        rest9.phase["pull"]["force_constants"], np.array([3.0, 3.0, 3.0, 3.0])
    )
    assert np.allclose(rest9.phase["pull"]["targets"], np.array([0.0, 1.0, 2.0, 3.0]))
    assert rest9.phase["release"]["force_constants"] == None
    assert rest9.phase["release"]["targets"] == None
    window_list = create_window_list([rest9])
    assert window_list == ["p000", "p001", "p002", "p003"]

    # Just Release
    log.info("### Testing restraint 10, just release")
    rest10 = DAT_restraint()
    rest10.amber_index = True
    rest10.continuous_apr = False
    rest10.auto_apr = False
    rest10.topology = os.path.join(os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb")
    rest10.mask1 = ":CB6@O"
    rest10.mask2 = ":BUT@C3"
    rest10.release["target"] = 0.0
    rest10.release["num_windows"] = 3
    rest10.release["fc_initial"] = 0.0
    rest10.release["fc_final"] = 2.0
    rest10.initialize()
    assert rest10.index1 == [13]
    assert rest10.index2 == [119]
    assert rest10.index3 == None
    assert rest10.index4 == None
    assert rest10.phase["attach"]["force_constants"] == None
    assert rest10.phase["attach"]["targets"] == None
    assert rest10.phase["pull"]["force_constants"] == None
    assert rest10.phase["pull"]["targets"] == None
    assert np.allclose(
        rest10.phase["release"]["force_constants"], np.array([0.0, 1.0, 2.0])
    )
    assert np.allclose(rest10.phase["release"]["targets"], np.array([0.0, 0.0, 0.0]))
    window_list = create_window_list([rest10])
    assert window_list == ["r000", "r001", "r002"]

    # Test inconsistent continuous_apr:
    with pytest.raises(Exception) as e_info:
        window_list = create_window_list([rest7, rest8])

    # Test inconsistent windows:
    with pytest.raises(Exception) as e_info:
        window_list = create_window_list([rest1, rest10])


def test_static_DAT_restraint():
    structure = pmd.load_file(os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.prmtop"),
                              os.path.join(os.path.dirname(__file__), "../data/cb6-but/vac.rst7"))
    r = static_DAT_restraint(
        restraint_mask_list=[":BUT@C3", ":CB6@O"],
        num_window_list=[10, 10, 10],
        ref_structure=structure,
        force_constant=5.0,
        amber_index=True,
    )
