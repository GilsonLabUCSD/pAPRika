"""
Test that we can save and load restraints as JSON.
"""

import unittest
import os

import paprika
from paprika.restraints import *
from paprika.restraints_json import *


class TestJSON(unittest.TestCase):
    def test_save_and_load_single_restraint(self):
        """ Test we can save a simple restraint """
        rest1 = DAT_restraint()
        rest1.amber_index = True
        rest1.continuous_apr = False
        rest1.auto_apr = False
        rest1.topology = "./cb6-but/cb6-but-notcentered.pdb"
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
        save_restraints([rest1], "rest1.json")
        restraints = load_restraints("rest1.json")
        os.remove("rest1.json")
        assert rest1 == restraints[0]

    def test_save_and_load_list_restraint(self):
        """ Test we can save and load a list of restraints """
        rest1 = DAT_restraint()
        rest1.amber_index = True
        rest1.continuous_apr = False
        rest1.auto_apr = False
        rest1.topology = "./cb6-but/cb6-but-notcentered.pdb"
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

        rest2 = DAT_restraint()
        rest2.amber_index = True
        rest2.continuous_apr = False
        rest2.auto_apr = False
        rest2.topology = "./cb6-but/cb6-but-notcentered.pdb"
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

        save_restraints([rest1, rest2], "restraints.json")
        restraints = load_restraints("restraints.json")
        os.remove("restraints.json")
        assert rest1 == restraints[0]
        assert rest2 == restraints[1]


if __name__ == "__main__":
    log.debug("{}".format(paprika.__version__))
    unittest.main()
