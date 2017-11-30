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
    log.debug('OpenMM support: Yes')
except:
    log.debug('OpenMM support: No')


def _amber_write_input_file(filename, dictionary, title='Input.'):
    """
    Process a dictionary to generate a pmemd input file.
    """

    log.debug('Writing {}'.format(filename))
    with open(filename, 'w') as f:
        f.write("{}\n".format(title))
        for namelist in ['cntrl', 'ewald']:
            if dictionary[namelist] is not None:
                f.write(" &{}\n".format(namelist))
                for key, val in dictionary[namelist].items():
                    if val is not None:
                        f.write("  {:15s} {:s},\n".format(
                            key + ' =', str(val)))
                f.write(" /\n")
        if dictionary['cntrl']['nmropt'] == 1:
            if dictionary['wt'] is not None:
                for line in dictionary['wt']:
                    f.write(" " + line + "\n")
            f.write(" &wt type = 'END', /\n")
            if dictionary['DISANG'] is not None:
                f.write("DISANG = {}\n".format(dictionary['DISANG']))
                f.write("LISTOUT = POUT\n\n")
        if dictionary['GROUP'] is not None:
            f.write("{:s}".format(dictionary['GROUP']))


def _apply_openmm_restraints(system):
         """
        Loop through the instances of `DAT_restraint`, after an OpenMM `system` has been created and
        call the function to apply the restraint for a given phase and window of the calculation.
        Parameters
        ----------
        restraint: a `DAT_restraint` object
        system: OpenMM System Class
        """
        for i, restraint in enumerate(DAT_restraint.get_instances()):
            log.debug('Setting up restraint number {} in phase {} and window {}...'.format(
                i, self.phase, self.window))
            # self.setup_openmm_restraints(system, restraint, self.phase, self.window)

class pme_simulation():
    """
    All PME/PBC simulation methods.
    """

    def __init__(self):
        self.todo = 'todo'


class gb_simulation(phase=None, window=None):
    """Setup and run a GB simulation in AMBER or OpenMM.

    Parameters
    ----------
    phase: phase of the calculation used to specify the value of the restraints, should be "attach", "pull" or "release".
    window: window of the calculation used to specify the value of the restraints, should be numeric.
    """

    def __init__(self):

        # Setup simulation directory and files
        self.path = './'
        self.amber_executable = 'pmemd'
        self.toplogy = self.path + 'prmtop'
        self.restraints = self.path + 'restraints.in'

        self.amber = {}
        self.amber['min'] = {}
        self.amber['min']['prefix'] = 'minimize'
        self.amber['min']['input'] = self.path + self.amber['min']['prefix'] + '.in'
        self.amber['min']['output'] = self.path + self.amber['min']['prefix'] + '.out'
        self.amber['min']['restart'] = self.path + self.amber['min']['prefix'] + '.rst'
        self.amber['min']['info'] = self.path + self.amber['min']['prefix'] + '.inf'
        self.amber['min']['coords'] = self.path + self.amber['min']['prefix'] + '.inpcrd'

        self.amber['min']['cntrl'] = OrderedDict()
        # Minimize
        self.amber['min']['cntrl']['imin'] = 1
        # Coordinates will be read
        self.amber['min']['cntrl']['ntx'] = 1
        # Not a restart
        self.amber['min']['cntrl']['irest'] = 0
        # Maximum number of cycles
        self.amber['min']['cntrl']['maxcyc'] = 5000
        # Number of steps until switching from SD to CG
        self.amber['min']['cntrl']['ncyc'] = 1000
        # Number of steps to output human-readable energy
        self.amber['min']['cntrl']['ntpr'] = 100
        # Output format
        self.amber['min']['cntrl']['ntxo'] = 1
        # Whether all interactions are calculated
        self.amber['min']['cntrl']['ntf'] = 1
        # Whether SHAKE is performed
        self.amber['min']['cntrl']['ntc'] = 1
        # Nonbonded cutoff
        self.amber['min']['cntrl']['cut'] = 999.0
        # Flag for generalized Born
        self.amber['min']['cntrl']['igb'] = 1
        self.amber['min']['cntrl']['ntr'] = None
        self.amber['min']['cntrl']['restraint_wt'] = None
        self.amber['min']['cntrl']['restraintmask'] = None
        self.amber['min']['cntrl']['nmropt'] = 1
        self.amber['min']['cntrl']['pencut'] = -1

        self.amber['min']['ewald'] = None
        self.amber['min']['wt'] = None    # or []
        self.amber['min']['DISANG'] = self.restraints
        self.amber['min']['GROUP'] = None    # or []

        self.amber['md'] = {}
        self.amber['md']['prefix'] = 'md'
        self.amber['md']['input'] = self.path + self.amber['md']['prefix'] + '.in'
        self.amber['md']['output'] = self.path + self.amber['md']['prefix'] + '.out'
        self.amber['md']['restart'] = self.path + self.amber['md']['prefix'] + '.rst'
        self.amber['md']['info'] = self.path + self.amber['md']['prefix'] + '.inf'
        self.amber['md']['coords'] = self.path + self.amber['md']['prefix'] + '.inpcrd'
        self.amber['md']['ref_coords'] = None

        self.amber['md']['cntrl'] = OrderedDict()
        self.amber['md']['cntrl']['imin'] = 0
        self.amber['md']['cntrl']['ntx'] = 1
        self.amber['md']['cntrl']['irest'] = 0
        self.amber['md']['cntrl']['dt'] = 0.002
        self.amber['md']['cntrl']['nstlim'] = 5000
        self.amber['md']['cntrl']['ntpr'] = 500
        self.amber['md']['cntrl']['ntwe'] = 0
        self.amber['md']['cntrl']['ntwr'] = 5000
        self.amber['md']['cntrl']['ntxo'] = 1
        self.amber['md']['cntrl']['ntwx'] = 500
        self.amber['md']['cntrl']['ioutfm'] = 1
        self.amber['md']['cntrl']['ntf'] = 2
        self.amber['md']['cntrl']['ntc'] = 2
        self.amber['md']['cntrl']['igb'] = 1
        self.amber['md']['cntrl']['cut'] = 999.0
        self.amber['md']['cntrl']['ntt'] = 3
        self.amber['md']['cntrl']['gamma_ln'] = 1.0
        self.amber['md']['cntrl']['ig'] = -1
        self.amber['md']['cntrl']['temp0'] = 298.15
        self.amber['md']['cntrl']['ntr'] = None
        self.amber['md']['cntrl']['restraint_wt'] = None
        self.amber['md']['cntrl']['restraintmask'] = None
        self.amber['md']['cntrl']['nmropt'] = 1
        self.amber['md']['cntrl']['pencut'] = -1

        self.amber['md']['ewald'] = None
        self.amber['md']['wt'] = None
        self.amber['md']['DISANG'] = self.restraints
        self.amber['md']['GROUP'] = None

        self.openmm = {}
        self.openmm['forcfield'] = None
        self.openmm['min'] = {}
        self.openmm['min']['coords'] = self.amber['min']['coords']
        self.openmm['min']['prefix'] = 'minimize'
        self.openmm['min']['output'] = self.path + self.openmm['min']['prefix'] + '.pdb'
        self.openmm['min']['positions'] = None
        self.openmm['min']['platform'] = 'CUDA'  # CUDA, OpenCL, CPU, Reference
        self.openmm['min']['precision'] = 'mixed'
        self.openmm['min']['max_iterations'] = 5000
        self.openmm['min']['reporter_frequency'] = 1000
        self.openmm['min']['tolerance'] = 1
        # Hawkins-Cramer-Truhlar (igb = 1 in AMBER)
        self.openmm['min']['solvent'] = app.HCT
        self.openmm['min']['salt_conc'] = 0.1 * unit.mole / unit.liter
        self.openmm['min']['nonbonded_method'] = app.NoCutoff
        self.openmm['min']['nonbonded_cutoff'] = None
        self.openmm['min']['constraints'] = app.HBonds
        self.openmm['min']['temperature'] = 300 * unit.kelvin
        self.openmm['min']['friction'] = 1.0 / unit.picoseconds
        self.openmm['min']['timestep'] = 2.0 * unit.femtoseconds

        self.openmm['md'] = {}
        self.openmm['md']['coords'] = self.amber['md']['coords']
        self.openmm['md']['prefix'] = 'md'
        self.openmm['md']['output'] = self.path + self.openmm['md']['prefix'] + '.nc'
        self.openmm['md']['data'] = self.path + self.openmm['md']['prefix'] + '.csv'

        self.openmm['md']['platform'] = 'CUDA'
        self.openmm['md']['precision'] = 'mixed'
        self.openmm['md']['reporter_frequency'] = 1000
        self.openmm['md']['solvent'] = app.HCT
        self.openmm['md']['salt_conc'] = 0.1 * unit.mole / unit.liter
        self.openmm['md']['nonbonded_method'] = app.NoCutoff
        self.openmm['md']['nonbonded_cutoff'] = None
        self.openmm['md']['constraints'] = app.HBonds
        self.openmm['md']['temperature'] = 300 * unit.kelvin
        self.openmm['md']['friction'] = 1.0 / unit.picoseconds
        self.openmm['md']['timestep'] = 2.0 * unit.femtoseconds
        self.openmm['md']['steps'] = 10000

    def minimize_with_amber(self, soft=False):
        """
        Minimize the system.

        If soft=True, slowly turn on non-bonded interactions during minimization
        so that the restraints get enforced first.

        """

        if soft:
            self.amber['min']['wt'] = [
                "&wt type = 'NB', istep1=0, istep2=2000, value1 = 0.0, value2=0.0, IINC=50, /",
                "&wt type = 'NB', istep1=2000, istep2=4000, value1 = 0.0, value2=1.0, IINC=50, /"]
        _amber_write_input_file(
            self.amber['min']['input'], self.amber['min'], title='GB Minimization.')

        log.info('Running GB Minimization at {}'.format(self.path))
        sp.call("{6} -O -p {0} -c {1} -ref {1} -i {2}\
                 -o {3} -r {4} -inf {5}"
                .format(self.toplogy, self.amber['min']['coords'], self.amber['min']['input'],
                        self.amber['min']['output'], self.amber['min']['restart'], self.amber['min']['info'], self.amber_executable), cwd=self.path, shell=True)
        log.debug('TODO: Does the above command need `shell=True` ...?')
        log.debug('TODO: Catch errors here...')

    def run_md_with_amber(self):
        """
        Run MD with AMBER.
        """

        _amber_write_input_file(self.amber['md']['input'], self.amber['md'], title='MD.')
        log.info('Running AMBER MD at {}'.format(self.path))

        if self.amber['md']['ref_coords']:
            log.warn('Check the specificiation of reference coordinates is correct.')
            execution_string = "{6} -O -p {0} -c {1} -ref {7} -i {2}\
                     -o {3} -r {4} -inf {5}".format(self.toplogy, self.amber['md']['coords'], self.amber['md']['input'],
                                                    self.amber['md']['output'], self.amber['md'][
                                                        'restart'], self.amber['md']['restart'],
                                                    self.amber_executable, self.amber['md']['ref_coords'])

            log.debug(execution_string)
            sp.call(execution_string, cwd=self.path)

        else:
            execution_string = "{6} -O -p {0} -c {1} -i {2}\
                     -o {3} -r {4} -inf {5}".format(self.toplogy, self.amber['md']['coords'], self.amber['md']['input'],
                                                    self.amber['md']['output'], self.amber['md'][
                                                        'restart'], self.amber['md']['restart'],
                                                    self.amber_executable)

            log.debug(execution_string)
            sp.call(execution_string, cwd=self.path)
        log.info('Minimization completed...')

    def minimize_with_openmm(self, soft=False):
        """
        Run MD with OpenMM.
        (The `soft` option can probably be folded into a class variable to tidy this up.)
        """
        prmtop = pmd.load_file(self.toplogy, self.openmm['min']['coords'])
        # I'm not sure if I need an integrator for minimization!
        integrator = mm.LangevinIntegrator(
            self.openmm['min']['temperature'],
            self.openmm['min']['friction'],
            self.openmm['min']['timestep']
        )
        platform = mm.Platform.getPlatformByName(
            self.openmm['min']['platform'])
        if self.openmm['min']['platform'] == 'CUDA':
            prop = dict(CudaPrecision=self.openmm['min']['precision'])
        else:
            prop = None

        if self.openmm['forcfield'] is not None:
            log.warn(
                'We haven\'t tested running OpenMM with an external force field yet.')
            forcefield = app.ForceField(self.opennmm['forcefield'])
            log.warn('Probably need to load a separate topology here...')
            system = forcefield.createSystem(
                nonbondedMethod=self.openmm['min']['nonbonded_method'],
                implicitSolvent=self.openmm['min']['solvent'],
                implicitSolventSaltConc=self.openmm['min']['salt_conc'],
                constraints=self.openmm['min']['constraints'],
                platform=platform,
                prop=prop)
        else:
            system = prmtop.createSystem(
                nonbondedMethod=self.openmm['min']['nonbonded_method'],
                implicitSolvent=self.openmm['min']['solvent'],
                implicitSolventSaltConc=self.openmm['min']['salt_conc'],
                constraints=self.openmm['min']['constraints'],
                platform=platform,
                prop=prop)
        simulation = app.Simulation(prmtop.topology, system)
        simulation.context.setPositions(prmtop.positions)
        self.add_openmm_restraints(system)
        log.info('Running OpenMM minimization...')

        if soft:
            # Phase 1: minimize with nonbonded interactions disabled.
            # This is the first 40% of the maximum iterations.
            log.debug('Minimization phase 1 for {} steps.'
                      .format(int(0.4 * self.openmm['min']['max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.4 * self.openmm['min']['max_iterations']),
                tolerance=self.openmm['min']['tolerance'] *
                unit.kilojoule / unit.mole
            )
            # Phase 2: slowly turn on the nonbonded interactions.
            # This is the next 40% of the maximum iterations.
            # This increases the nonbonded interactions linearly, which is not
            # the same as using `IINC` in AMBER.
            log.debug('Minimization phase 2 for {} steps.'
                      .format(int(0.4 * self.openmm['min']['max_iterations'])))
            for scale in np.linspace(0, 1.0, int(0.4 * self.openmm['min']['max_iterations'])):
                log.debug(
                    'Scaling NB interactions to {0:0.4f} / 1.0'.format(scale))
                for particle in range(system.getNumParticles()):
                    [charge, sigma, epsilon] = mm.NonbondedForce.getParticleParameters(
                        particle)
                    mm.NonbondedForce.setParticleParameters(
                        particle, charge * scale, sigma, epsilon * scale)
                simulation.minimizeEnergy(
                    maxIterations=1,
                    tolerance=self.openmm['min']['tolerance'] *
                    unit.kilojoule / unit.mole
                )
            # Phase 3: minimize with nonbonded interactions at full strength.
            # This is the last 20% of the maximum iterations.
            log.debug('Minimization phase 3 for {} steps.'
                      .format(int(0.2 * self.openmm['min']['max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.2 * self.openmm['min']['max_iterations']),
                tolerance=self.openmm['min']['tolerance'] *
                unit.kilojoule / unit.mole
            )
            state = simulation.context.getState(
                getEnergy=True, getPositions=True)
        else:
            simulation.minimizeEnergy(
                maxIterations=self.openmm['min']['max_iterations'],
                tolerance=self.openmm['min']['tolerance'] *
                unit.kilojoule / unit.mole
            )
        self.openmm['min']['positions'] = simulation.context.getState(
            getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, self.openmm['min']['positions'], open(
            self.openmm['min']['output'], 'w'))
        log.info('Minimization completed...')

    def run_md_with_openmm(self):
        """
        Run MD with OpenMM.
        """
        prmtop = pmd.load_file(self.toplogy, self.openmm['md']['coords'])
        # I'm not sure if I need an integrator for minimization!
        integrator = mm.LangevinIntegrator(
            self.openmm['md']['temperature'],
            self.openmm['md']['friction'],
            self.openmm['md']['timestep']
        )
        platform = mm.Platform.getPlatformByName(self.openmm['md']['platform'])
        if self.openmm['md']['platform'] == 'CUDA':
            prop = dict(CudaPrecision=self.openmm['md']['precision'])
        else:
            prop = None
        if self.openmm['forcfield'] is not None:
            log.warn(
                'We haven\'t tested running OpenMM with an external force field yet.')
            forcefield = app.ForceField(self.openmm['forcefield'])
            log.warn('Probably need to load a separate topology here...')
            system = forcefield.createSystem(
                nonbondedMethod=self.openmm['md']['nonbonded_method'],
                implicitSolvent=self.openmm['md']['solvent'],
                implicitSolventSaltConc=self.openmm['md']['salt_conc'],
                constraints=self.openmm['md']['constraints'],
                platform=platform,
                prop=prop)
        else:
            system = prmtop.createSystem(
                nonbondedMethod=self.openmm['md']['nonbonded_method'],
                implicitSolvent=self.openmm['md']['solvent'],
                implicitSolventSaltConc=self.openmm['md']['salt_conc'],
                constraints=self.openmm['md']['constraints'],
                platform=platform,
                prop=prop)
        simulation = app.Simulation(prmtop.topology, system)
        self.add_openmm_restraints(system)

        if self.openmm['min']['positions']:
            simulation.context.setPositions(self.openmm['min']['positions'])
        else:
            simulation.context.setPositions(prmtop.positions)
        simulation.context.setVelocitiesToTemperature(
            self.openmm['md']['temperature'] * unit.kelvin)
        reporter = mm.NetCDFReporter(
            self.openmm['md']['output'], self.openmm['md']['reporter_frequency'])
        simulation.reporters.append(reporter)
        simulation.reporters.append(app.StateDataReporter(self.openmm['md']['data'],
                                                          self.openmm['md']['reporter_frequency'],
                                                          step=True,
                                                          potentialEnergy=True,
                                                          temperature=True,
                                                          density=True))

        log.info('Running OpenMM MD...')
        simulation.step(self.openmm['md']['steps'])
        log.info('MD completed...')
        reporter.close()

    
