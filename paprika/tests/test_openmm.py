"""
Tests basic OpenMM simulations.
"""

import os

import pytest
import paprika

# Skip these tests if you cannot import paprika.openmm.
pytest.importorskip("paprika.openmm")

from paprika.openmm import *
from paprika.restraints import *
from paprika.tests import addons


@addons.using_openmm
@pytest.mark.xfail
def test_minimization_finishes():
    """ Test that we can minimize CB6-BUT with OpenMM. """

    my_simulation = OpenMM_GB_simulation()
    my_simulation.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.topo"
    )
    my_simulation.min["platform"] = "Reference"
    my_simulation.min["coordinates"] = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.crds"
    )

    system = my_simulation.setup_system(my_simulation.min)
    simulation = my_simulation.setup_simulation(system, my_simulation.min)

    result = my_simulation.minimize(simulation, save=False)
    state = result.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    np.testing.assert_almost_equal(energy, -821.306, decimal=2)


@addons.using_openmm
@pytest.mark.xfail
def test_soft_minimization():
    """ Test that we can minimize CB6-BUT with OpenMM, turning on interactions slowly. """

    my_simulation = OpenMM_GB_simulation()
    my_simulation.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.topo"
    )
    my_simulation.min["platform"] = "Reference"
    my_simulation.min["coordinates"] = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.crds"
    )
    my_simulation.min["max_iterations"] = 100

    system = my_simulation.setup_system(my_simulation.min)
    simulation = my_simulation.setup_simulation(system, my_simulation.min)

    result = my_simulation.turn_on_interactions_slowly(simulation, system)
    state = result.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    np.testing.assert_almost_equal(energy, -821.306, decimal=2)


@addons.using_openmm
@pytest.mark.xfail
def test_openmm_single_restraint():
    """ Test that we can impose restraints with OpenMM. """
    my_simulation = OpenMM_GB_simulation()
    my_simulation.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.topo"
    )
    my_simulation.md["platform"] = "Reference"
    my_simulation.md["coordinates"] = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.crds"
    )
    my_simulation.md["steps"] = 1000

    system = my_simulation.setup_system(my_simulation.md)

    restraint = DAT_restraint()
    restraint.topology = my_simulation.topology
    restraint.mask1 = ":BUT"
    restraint.mask2 = ":CB6"
    restraint.attach["target"] = 3.0
    restraint.attach["num_windows"] = 4
    restraint.attach["fc_initial"] = 0.0
    restraint.attach["fc_final"] = 3.0
    restraint.auto_apr = False
    restraint.initialize()

    my_simulation.add_openmm_restraints(system, [restraint], phase="attach", window=3)

    simulation = my_simulation.setup_simulation(system, my_simulation.md)

    result = my_simulation.run_md(simulation, save=False)
    state = result.context.getState(getEnergy=True)
    state.getPotentialEnergy() / unit.kilocalories_per_mole
    # np.testing.assert_almost_equal(energy, -709.6, decimal=1)


@addons.using_openmm
@pytest.mark.xfail
def test_openmm_two_restraints():
    """ Test that we can impose restraints with OpenMM. """
    my_simulation = OpenMM_GB_simulation()
    my_simulation.topology = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.topo"
    )
    my_simulation.md["platform"] = "Reference"
    my_simulation.md["coordinates"] = os.path.join(
        os.path.dirname(__file__), "../data/cb6-but/vac.crds"
    )
    my_simulation.md["steps"] = 1000

    system = my_simulation.setup_system(my_simulation.md)

    restraint_1 = DAT_restraint()
    restraint_1.topology = my_simulation.topology
    restraint_1.mask1 = ":BUT"
    restraint_1.mask2 = ":CB6"
    restraint_1.attach["target"] = 3.0
    restraint_1.attach["num_windows"] = 4
    restraint_1.attach["fc_initial"] = 0.0
    restraint_1.attach["fc_final"] = 3.0
    restraint_1.auto_apr = False
    restraint_1.initialize()

    restraint_2 = DAT_restraint()
    restraint_2.topology = my_simulation.topology
    restraint_2.mask1 = ":CB6@O"
    restraint_2.mask2 = ":CB6@O6"
    restraint_2.attach["target"] = 8.0
    restraint_2.attach["num_windows"] = 4
    restraint_2.attach["fc_initial"] = 0.0
    restraint_2.attach["fc_final"] = 30.0
    restraint_2.auto_apr = False
    restraint_2.initialize()

    my_simulation.add_openmm_restraints(
        system, [restraint_1, restraint_2], phase="attach", window=3
    )

    simulation = my_simulation.setup_simulation(system, my_simulation.md)

    result = my_simulation.run_md(simulation, save=True)
    state = result.context.getState(getEnergy=True)
    state.getPotentialEnergy() / unit.kilocalories_per_mole

    utils.decompose_openmm_energy(result)

    # np.testing.assert_almost_equal(energy, -709.6, decimal=1)
    # np.testing.assert_almost_equal(energies['total'], -709.6, decimal=1)
    # np.testing.assert_almost_equal(energies['restraint'], 1.6, decimal=1)
    # np.testing.assert_almost_equal(energies['non-restraint'], -711.1, decimal=1)

    for filename in ["md.nc", "md.csv"]:
        if os.path.isfile("./" + filename):
            os.remove("./" + filename)
