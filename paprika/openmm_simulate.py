import logging as log
import subprocess as sp
from collections import OrderedDict

import numpy as np
import parmed as pmd
from paprika.restraints import *

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit
    from mdtraj.reporters import NetCDFReporter
    log.debug('OpenMM support: Yes')
except:
    log.debug('OpenMM support: No')


class OpenMM_GB_simulation():
    """Setup and run a GB simulation in AMBER or OpenMM.

    Parameters
    ----------
    phase : phase of the calculation used to specify the value of the restraints, should be
    "attach", "pull" or "release".
    window : window of the calculation used to specify the value of the restraints, should be
    numeric.
    """

    def __init__(self):

        self.path = './'
        self.coordinates = None
        self.topology = None
        self.phase = None
        self.window = None

        self.min = dict(coordinates=self.coordinates,
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

        self.md = dict(coordinates=self.coordinates,
                       minimized_coordinates=None,
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

    def turn_on_interactions_slowly(self, system, simulation):
            # Phase 1: minimize with nonbonded interactions disabled.
            # This is the first 40% of the maximum iterations.
            log.debug('Minimization phase 1 for {} steps.'.format(int(0.4 * self.min[
                'max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.4 * self.openmm['min']['max_iterations']),
                tolerance=self.openmm['min']['tolerance'] *
                unit.kilojoule / unit.mole
            )
            # Phase 2: slowly turn on the nonbonded interactions.
            # This is the next 40% of the maximum iterations.
            # This increases the nonbonded interactions linearly, which is not
            # the same as using `IINC` in AMBER.
            log.debug('Minimization phase 2 for {} steps.'.format(int(0.4 * self.min['max_iterations'])))
            for scale in np.linspace(0, 1.0, int(0.4 * self.min['max_iterations'])):
                log.debug('Scaling NB interactions to {0:0.4f} / 1.0'.format(scale))
                for particle in range(system.getNumParticles()):
                    [charge, sigma, epsilon] = mm.NonbondedForce.getParticleParameters(particle)
                    mm.NonbondedForce.setParticleParameters(particle, charge * scale, sigma,
                                                            epsilon * scale)
                simulation.minimizeEnergy(
                    maxIterations=1,
                    tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole
                )
            # Phase 3: minimize with nonbonded interactions at full strength.
            # This is the last 20% of the maximum iterations.
            log.debug('Minimization phase 3 for {} steps.'.format(int(0.2 * self.min[
                'max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.2 * self.min['max_iterations']),
                tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole
            )
            return simulation

    def add_openmm_restraints(self, system):
        """
        Loop through the instances of `DAT_restraint`, after an OpenMM `system` has been created and
        call the function to apply the restraint for a given phase and window of the calculation.
        """
        for i, restraint in enumerate(DAT_restraint.restraint_list):
            log.debug('Setting up restraint number {} in phase {} and window {}...'.format(i,
                                                                                           self.phase,
                                                                                           self.window))
            print(system)
            system = setup_openmm_restraints(system, restraint, self.phase, self.window)
            print(system)
            return system

    def minimize(self):
        """
        Run MD with OpenMM.
        """
        prmtop = pmd.load_file(self.topology, self.min['coordinates'])
        # I'm not sure why we need an integrator for minimization!
        integrator = mm.LangevinIntegrator(
            self.min['temperature'],
            self.min['friction'],
            self.min['timestep']
        )
        platform = mm.Platform.getPlatformByName(self.min['platform'])
        if self.min['platform'] == 'CUDA':
            prop = dict(CudaPrecision=self.min['precision'],
                        CudaDeviceIndex=self.min['devices'])
        else:
            prop = None

        if self.min['forcefield'] is not None:
            log.warning('We haven\'t tested running OpenMM with an external force field yet.')
            forcefield = app.ForceField(self.min['forcefield'])
            log.warning('Probably need to load a separate topology here...')
            system = forcefield.createSystem(
                nonbondedMethod=self.min['nonbonded_method'],
                implicitSolvent=self.min['solvent'],
                implicitSolventSaltConc=self.min['salt'],
                constraints=self.min['constraints']
               )
        else:
            system = prmtop.createSystem(
                nonbondedMethod=self.min['nonbonded_method'],
                implicitSolvent=self.min['solvent'],
                implicitSolventSaltConc=self.min['salt'],
                constraints=self.min['constraints']
                )
        simulation = app.Simulation(prmtop.topology, system, integrator, platform, prop)
        simulation.context.setPositions(prmtop.positions)
        system = self.add_openmm_restraints(system)
        log.info('Running OpenMM minimization...')

        if self.min['soft']:
            simulation = self.turn_on_interactions_slowly(system, simulation)
        else:
            simulation.minimizeEnergy(
                maxIterations=self.min['max_iterations'],
                tolerance=self.min['tolerance'] * unit.kilojoule / unit.mole
            )

        self.md['minimized_coordinates'] = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, self.md['minimized_coordinates'],
                              open(self.min['output'], 'w'))
        log.info('Minimization completed.')

    def run_md(self):
        """
        Run MD with OpenMM.
        """
        prmtop = pmd.load_file(self.topology, self.md['coordinates'])
        # I'm not sure why we need an integrator for minimization!
        integrator = mm.LangevinIntegrator(
            self.md['temperature'],
            self.md['friction'],
            self.md['timestep']
        )
        platform = mm.Platform.getPlatformByName(self.md['platform'])
        if self.md['platform'] == 'CUDA':
            prop = dict(CudaPrecision=self.min['precision'],
                        CudaDeviceIndex=self.min['devices'])
        else:
            prop = None

        if self.min['forcefield'] is not None:
            log.warning('We haven\'t tested running OpenMM with an external force field yet.')
            forcefield = app.ForceField(self.md['forcefield'])
            log.warning('Probably need to load a separate topology here...')
            system = forcefield.createSystem(
                nonbondedMethod=self.md['nonbonded_method'],
                implicitSolvent=self.md['solvent'],
                implicitSolventSaltConc=self.md['salt'],
                constraints=self.md['constraints']
            )
        else:
            system = prmtop.createSystem(
                nonbondedMethod=self.md['nonbonded_method'],
                implicitSolvent=self.md['solvent'],
                implicitSolventSaltConc=self.md['salt'],
                constraints=self.md['constraints']
                )
        simulation = app.Simulation(prmtop.topology, system, integrator, platform, prop)
        system = self.add_openmm_restraints(system)

        if self.md['minimized_coordinates']:
            simulation.context.setPositions(self.md['minimized_coordinates'])
        else:
            simulation.context.setPositions(prmtop.positions)
        simulation.context.setVelocitiesToTemperature(self.md['temperature'] * unit.kelvin)
        reporter = NetCDFReporter(self.md['output'], self.md['reporter_frequency'])
        simulation.reporters.append(reporter)
        simulation.reporters.append(app.StateDataReporter(self.md['data'],
                                                          self.md['reporter_frequency'],
                                                          step=True,
                                                          potentialEnergy=True,
                                                          temperature=True,
                                                          density=True))

        log.info('Running OpenMM MD...')
        simulation.step(self.md['steps'])
        log.info('MD completed.')
        reporter.close()