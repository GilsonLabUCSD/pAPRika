import logging as log

import numpy as np
from paprika.restraints import *

import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from mdtraj.reporters import NetCDFReporter


class OpenMM_GB_simulation():
    """Setup and run a GB simulation in OpenMM."""

    def __init__(self):

        self.path = './'
        self.coordinates = None
        self.topology = None
        self.phase = None
        self.window = None

        self.min = dict(
            coordinates=self.coordinates,
            prefix='minimize',
            forcefield=None,
            platform='CUDA',
            devices='1',
            precision='mixed',
            max_iterations=5000,
            reporter_frequency=1000,
            tolerance=1,
            solvent=app.HCT,
            salt=0.1 * unit.mole / unit.liter,
            nonbonded_method=app.NoCutoff,
            nonbonded_cutoff=None,
            constraints=app.HBonds,
            temperature=300 * unit.kelvin,
            friction=1.0 / unit.picoseconds,
            timestep=2.0 * unit.femtoseconds,
            soft=False)
        self.min['output'] = self.path + self.min['prefix'] + '.pdb'

        self.md = dict(
            coordinates=self.coordinates,
            minimized_coordinates=None,
            integrator=None,
            prefix='md',
            forcefield=None,
            platform='CUDA',
            devices='1',
            precision='mixed',
            reporter_frequency=1000,
            solvent=app.HCT,
            salt=0.1 * unit.mole / unit.liter,
            nonbonded_method=app.NoCutoff,
            nonbonded_cutoff=None,
            constraints=app.HBonds,
            temperature=300 * unit.kelvin,
            friction=1.0 / unit.picoseconds,
            timestep=2.0 * unit.femtoseconds,
            steps=10000)
        self.md['output'] = self.path + self.md['prefix'] + '.nc'
        self.md['data'] = self.path + self.md['prefix'] + '.csv'

    def setup_system(self, settings, seed=None):
        """
        Provide a way to create an OpenMM system object with minimization or MD settings.

        Parameters
        ----------
        settings : dict
            A dictionary containing simulation settings.

        Returns
        -------
        simulation : simtk.openmm.app.Simulation
            The simulation object.
        system : simtk.openmm.System
            The system object.
        """

        prmtop = app.AmberPrmtopFile(self.topology)
        app.AmberInpcrdFile(settings['coordinates'])
        self.integrator = mm.LangevinIntegrator(settings['temperature'],
                                                settings['friction'],
                                                settings['timestep'])

        if seed is not None:
            self.integrator.setRandomNumberSeed(seed)
        try:
            log.debug('Integrator random number seed: {}'.format(
                self.integrator.getRandomNumberSeed()))
        except AttributeError:
            pass
        system = prmtop.createSystem(
            nonbondedMethod=settings['nonbonded_method'],
            implicitSolvent=settings['solvent'],
            implicitSolventSaltConc=settings['salt'],
            constraints=settings['constraints'])

        # simulation = app.Simulation(prmtop.topology, system, integrator,
        #                             platform, prop)
        # simulation.context.setPositions(inpcrd.positions)
        return system

    def setup_simulation(self, system, settings):
        inpcrd = app.AmberInpcrdFile(settings['coordinates'])
        prmtop = app.AmberPrmtopFile(self.topology)
        platform, prop = self.setup_platform(settings)
        simulation = app.Simulation(prmtop.topology, system, self.integrator,
                                    platform, prop)
        simulation.context.setPositions(inpcrd.positions)
        return simulation

    def setup_platform(self, settings):
        """
        Set whether CPU or GPU should be used for OpenMM simulations, split off for readability.

        Parameters
        ----------
        settings : dict
            A dictionary containing simulation settings.

        Returns
        -------
        platform : simtk.openmm.mm.Platform
            The platform object.
        properties : dict
            The simulation settings for GPUs.

        """
        platform = mm.Platform.getPlatformByName(settings['platform'])
        log.debug('Platform: {}'.format(settings['platform']))
        if settings['platform'] == 'CUDA':
            properties = dict(
                CudaPrecision=settings['precision'],
                CudaDeviceIndex=settings['devices'])
            log.debug(properties)
        else:
            properties = None

        return platform, properties

    def turn_on_interactions_slowly(self, simulation, system):
        """
        Provide an interface to slowly turn on the nonbonded parameters.

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
            The simulation object.
        system : simtk.openmm.System
            The system object.

        Returns
        -------
        simulation : simtk.openmm.app.Simulation
            The simulation (after minimization).
        """
        # Phase 1: minimize with nonbonded interactions disabled.
        # This is the first 40% of the maximum iterations.
        log.debug('Minimization phase 1 for {} steps.'.format(
            int(0.4 * self.min['max_iterations'])))
        simulation.minimizeEnergy(
            maxIterations=int(0.4 * self.min['max_iterations']),
            tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole)
        # Phase 2: slowly turn on the nonbonded interactions.
        # This is the next 40% of the maximum iterations.
        # This increases the nonbonded interactions linearly, which is not
        # the same as using `IINC` in AMBER.
        log.debug('Minimization phase 2 for {} steps.'.format(
            int(0.4 * self.min['max_iterations'])))
        forces = {
            force.__class__.__name__: force
            for force in system.getForces()
        }
        nb_force = forces['NonbondedForce']
        for scale in np.linspace(0, 1.0,
                                 int(0.4 * self.min['max_iterations'])):
            for particle in range(nb_force.getNumParticles()):
                [charge, sigma,
                 epsilon] = nb_force.getParticleParameters(particle)
                nb_force.setParticleParameters(particle, charge * scale, sigma,
                                               epsilon * scale)
            simulation.minimizeEnergy(
                maxIterations=1,
                tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole)
        # Phase 3: minimize with nonbonded interactions at full strength.
        # This is the last 20% of the maximum iterations.
        log.debug('Minimization phase 3 for {} steps.'.format(
            int(0.2 * self.min['max_iterations'])))
        simulation.minimizeEnergy(
            maxIterations=int(0.2 * self.min['max_iterations']),
            tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole)
        return simulation

    def add_openmm_restraints(self, system, restraints, phase, window):
        """
        Loop through the instances of `DAT_restraint`, after an OpenMM `system` has been created and
        call the function to apply the restraint for a given phase and window of the calculation.

        Parameters
        ----------
        system : simtk.openmm.System
            The system object.
        restraints : list
            A list of restraints to be applied.
        phase : str
            The phase of the simulation, used with the window, to determine the restraint
            target value and strength.
        window : int
            The window of the simulation, used with the phase, to determine the restraint target
            value and strength.

        Returns
        -------
        system: simtk.openmm.System
            The same system object, with restraints applied.

        """
        for i, restraint in enumerate(restraints):
            # log.debug(system.getForces())
            system = setup_openmm_restraints(system, restraint, phase, window)
            # log.debug(system.getForces())

        return system

    def minimize(self, simulation, save=True):
        """
        Minimize with OpenMM.
        """

        log.info('Running OpenMM minimization...')

        if self.min['soft']:
            simulation = self.turn_on_interactions_slowly(system, simulation)
        else:
            simulation.minimizeEnergy(
                maxIterations=self.min['max_iterations'],
                tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole)

        if save:
            self.md['minimized_coordinates'] = simulation.context.getState(
                getPositions=True).getPositions()
            app.PDBFile.writeFile(simulation.topology,
                                  self.md['minimized_coordinates'],
                                  open(self.min['output'], 'w'))
        return simulation

    def run_md(self, simulation, seed=None, save=True):
        """
        Run MD with OpenMM.
        """

        if self.md['minimized_coordinates']:
            simulation.context.setPositions(self.md['minimized_coordinates'])
        if seed is not None:
            log.debug('Velocity random number seed: {}'.format(seed))
            simulation.context.setVelocitiesToTemperature(
                self.md['temperature'] * unit.kelvin, seed)

        else:
            simulation.context.setVelocitiesToTemperature(
                self.md['temperature'] * unit.kelvin)

        if save:
            reporter = NetCDFReporter(
                self.md['output'], self.md['reporter_frequency'], cell=False)
            simulation.reporters.append(reporter)
            simulation.reporters.append(
                app.StateDataReporter(
                    self.md['data'],
                    self.md['reporter_frequency'],
                    step=True,
                    potentialEnergy=True,
                    temperature=True,
                    density=True))

        log.info('Running OpenMM MD...')
        simulation.step(self.md['steps'])
        log.info('MD completed.')

        if save:
            reporter.close()
        return simulation
