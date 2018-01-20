"""
Tests the solvation of the system using `tleap`.
"""

import unittest
import warnings
import numpy as np
import logging as log
import subprocess as sp
import random as random
import parmed as pmd
import paprika
from paprika.openmm_simulate import *


class TestOpenMM(unittest.TestCase):
    def test_minimization_finishes(self):
        """ Test that we can minimize CB6-BUT with OpenMM. """

        simulation = OpenMM_GB_simulation()
        simulation.topology = '../test/cb6-but/vac.topo'
        simulation.min['platform'] = 'CPU'
        simulation.min['coordinates'] = '../test/cb6-but/vac.crds'
        result = simulation.minimize(save=False)
        state = result.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        self.assertAlmostEqual(energy, -827.9, places=1)


if __name__ == '__main__':
    log.debug('{}'.format(paprika.__version__))
    unittest.main()
