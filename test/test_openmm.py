"""
Tests basic OpenMM simulations.
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
from paprika.restraints import *
from paprika.utils import HAS_OPENMM


@unittest.skipUnless(HAS_OPENMM, 'Cannot run OpenMM tests without OpenMM installed.')
class TestOpenMM(unittest.TestCase):
    def test_minimization_finishes(self):
        """ Test that we can minimize CB6-BUT with OpenMM. """

        simulation = OpenMM_GB_simulation()
        simulation.topology = '../test/cb6-but/vac.topo'
        simulation.min['platform'] = 'CPU'
        simulation.min['coordinates'] = '../test/cb6-but/vac.crds'
        result, system = simulation.minimize(save=False)
        state = result.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        self.assertAlmostEqual(energy, -827.9, places=1)

    def test_soft_minimization(self):
        """ Test that we can minimize CB6-BUT with OpenMM, turning on interactions slowly. """

        sim = OpenMM_GB_simulation()
        sim.topology = '../test/cb6-but/vac.topo'
        sim.min['platform'] = 'CPU'
        sim.min['coordinates'] = '../test/cb6-but/vac.crds'

        simulation, system = sim.setup_system(sim.min)

        result = sim.turn_on_interactions_slowly(simulation, system)
        state = result.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        self.assertAlmostEqual(energy, -827.9, places=1)

        # How do we get all CustomExternalForces?

    # def test_openmm_restraint(self):
    #     """ Test that we can impose restraints with OpenMM. """
    #     simulation = OpenMM_GB_simulation()
    #     simulation.topology = '../test/cb6-but/vac.topo'
    #     simulation.md['platform'] = 'CPU'
    #     simulation.md['coordinates'] = '../test/cb6-but/vac.crds'
    #     simulation.md['steps'] = 100

    #     restraint = DAT_restraint()
    #     restraint.mask1 = ':BUT'
    #     # Need a way to get the system here...

    #  restrraint could/should be a list

    #     # setup_openmm_restraints(system?, restraint, phase='a', window=0)

    #     result, system = simulation.run_md(save=False)
    #     state = result.context.getState(getEnergy=True)
    #     energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    #     print(energy)
    #     # self.assertAlmostEqual(energy, -827.9, places=1)


if __name__ == '__main__':
    log.debug('{}'.format(paprika.__version__))
    unittest.main()