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


@unittest.skipUnless(HAS_OPENMM,
                     'Cannot run OpenMM tests without OpenMM installed.')
class TestOpenMM(unittest.TestCase):
    def test_minimization_finishes(self):
        """ Test that we can minimize CB6-BUT with OpenMM. """

        my_simulation = OpenMM_GB_simulation()
        my_simulation.topology = '../test/cb6-but/vac.topo'
        my_simulation.min['platform'] = 'CPU'
        my_simulation.min['coordinates'] = '../test/cb6-but/vac.crds'

        simulation, _ = my_simulation.setup_system(my_simulation.min)

        result = my_simulation.minimize(simulation, save=False)
        state = result.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        self.assertAlmostEqual(energy, -827.9, places=1)

    # def test_soft_minimization(self):
    #     """ Test that we can minimize CB6-BUT with OpenMM, turning on interactions slowly. """

    #     my_simulation = OpenMM_GB_simulation()
    #     my_simulation.topology = '../test/cb6-but/vac.topo'
    #     my_simulation.min['platform'] = 'CPU'
    #     my_simulation.min['coordinates'] = '../test/cb6-but/vac.crds'

    #     simulation, system = my_simulation.setup_system(my_simulation.min)

    #     result = my_simulation.turn_on_interactions_slowly(simulation, system)
    #     state = result.context.getState(getEnergy=True)
    #     energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    #     self.assertAlmostEqual(energy, -827.9, places=1)

    def test_openmm_single_restraint(self):
        """ Test that we can impose restraints with OpenMM. """
        my_simulation = OpenMM_GB_simulation()
        my_simulation.topology = '../test/cb6-but/vac.topo'
        my_simulation.md['platform'] = 'Reference'
        my_simulation.md['coordinates'] = '../test/cb6-but/vac.crds'
        my_simulation.md['steps'] = 1000

        simulation, system = my_simulation.setup_system(
            my_simulation.md, seed=42)

        restraint = DAT_restraint()
        restraint.structure_file = my_simulation.topology
        restraint.mask1 = ':BUT'
        restraint.mask2 = ':CB6'
        restraint.attach['target'] = 3.0
        restraint.attach['num_windows'] = 4
        restraint.attach['fc_initial'] = 0.0
        restraint.attach['fc_final'] = 3.0
        restraint.auto_apr = False
        restraint.initialize()

        my_simulation.add_openmm_restraints(
            system, [restraint], phase='attach', window=3)

        result = my_simulation.run_md(simulation, seed=42, save=False)
        state = result.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        self.assertAlmostEqual(energy, -750.0, places=1)


if __name__ == '__main__':
    log.debug('{}'.format(paprika.__version__))
    unittest.main()
