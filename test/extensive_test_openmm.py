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
from paprika.utils import HAS_OPENMM, decompose_openmm_energy


def extensive_test_openmm_single_restraint():
    """ Test that we can impose restraints with OpenMM. """
    my_simulation = OpenMM_GB_simulation()
    my_simulation.topology = '../test/cb6-but/vac.topo'
    my_simulation.md['platform'] = 'CPU'
    my_simulation.md['coordinates'] = '../test/cb6-but/vac.crds'
    my_simulation.md['steps'] = 10000
    my_simulation.md['output'] = 'tmp.nc'

    system = my_simulation.setup_system(my_simulation.md, seed=42)

    restraint = DAT_restraint()
    restraint.structure_file = my_simulation.topology
    restraint.mask1 = ':BUT'
    restraint.mask2 = ':CB6'
    restraint.attach['target'] = 10.0
    restraint.attach['num_windows'] = 4
    restraint.attach['fc_initial'] = 0.0
    restraint.attach['fc_final'] = 3.0
    restraint.auto_apr = False
    restraint.initialize()

    # log.debug(system.getForces())
    my_simulation.add_openmm_restraints(
        system, [restraint], phase='attach', window=3)
    simulation = my_simulation.setup_simulation(system, my_simulation.md)

    result = my_simulation.run_md(simulation, seed=42, save=True)
    state = result.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    # Needs to pass a structure...
    # print(decompose_openmm_energy(my_simulation.topology, my_simulation.context))


extensive_test_openmm_single_restraint()
