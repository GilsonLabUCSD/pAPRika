import logging as log
import os as os
import subprocess as sp
from collections import OrderedDict

import numpy as np
import parmed as pmd

try:
    # OpenMM Imports
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit
    log.debug('OpenMM support: Yes')
except:
    log.debug('OpenMM support: No')


def _write_input_file(filename, input_dict, title='Input.'):
    """
    Process a dictionary to generate a pmemd input file.
    """

    log.debug('Writing {}'.format(filename))
    with open(filename, 'w') as f:
        f.write("{}\n".format(title))
        for namelist in ['cntrl','ewald']:
            if input_dict[namelist] is not None:
                f.write(" &{}\n".format(namelist))
                for key,val in input_dict[namelist].items():
                    if val is not None:
                        f.write("  {:15s} {:s},\n".format(key+' =',str(val)))
                f.write(" /\n")
        if input_dict['cntrl']['nmropt'] == 1:          ### Need to check if this exists first!!
            if input_dict['wt'] is not None:
                for line in  input_dict['wt']:
                    f.write(" "+line+"\n")
            f.write(" &wt type = 'END', /\n")
            if input_dict['DISANG'] is not None:
                f.write("DISANG = {}\n".format(input_dict['DISANG']))
                f.write("LISTOUT = POUT\n\n")
        if input_dict['GROUP'] is not None:
            f.write("{:s}".format(input_dict['GROUP']))



class pme_simulation():
    """
    All PME/PBC simulation methods.
    """

    def __init__(self):
        self.todo = 'todo'


class gb_simulation():
    """
    All GB simulation methods.
    """

    def __init__(self):

        # Setup simulation directory and files
        self.path = '.'
        self.toplogy = 'prmtop'
        self.coordinates = 'inpcrd'
        self.restraints = 'restraints.in'
        self.pdb = None

        # See AMBER manual for definitions
        # I'm not really pleased with the format here, but I want to keep
        # the order of the keys, so I have to go line by line for the
        # OrderedDict.  If I make the dict all at once, it won't be in order.
        # This format should allow flexible editing of the input files
        
        self.amber = {}
        # Minimization input file settings
        self.amber['minimization'] = {}
        self.amber['minimization']['cntrl'] =                     OrderedDict()
        # Minimize
        self.amber['minimization']['cntrl']['imin'] =             1 
        # Coordinates will be read
        self.amber['minimization']['cntrl']['ntx'] =              1
        # Not a restart
        self.amber['minimization']['cntrl']['irest'] =            0
        # Maximum number of cycles
        self.amber['minimization']['cntrl']['maxcyc'] =           5000
        # Number of steps until switching from SD to CG
        self.amber['minimization']['cntrl']['ncyc'] =             1000
        # Number of steps to output human-readable energy
        self.amber['minimization']['cntrl']['ntpr'] =             100
        # Output format
        self.amber['minimization']['cntrl']['ntxo'] =             1
        # Whether all interactions are calculated
        self.amber['minimization']['cntrl']['ntf'] =              1
        # Whether SHAKE is performed
        self.amber['minimization']['cntrl']['ntc'] =              1
        # Nonbonded cutoff
        self.amber['minimization']['cntrl']['cut'] =              999.0
        # Flag for generalized Born
        self.amber['minimization']['cntrl']['igb'] =              1
        self.amber['minimization']['cntrl']['ntr'] =              None
        self.amber['minimization']['cntrl']['restraint_wt'] =     None
        self.amber['minimization']['cntrl']['restraintmask'] =    None
        self.amber['minimization']['cntrl']['nmropt'] =           1
        self.amber['minimization']['cntrl']['pencut'] =           -1
        self.amber['minimization']['ewald'] =                     None    # or OrderedDict() (Add the other namelists?)
        self.amber['minimization']['wt'] =                        None    # or []
        self.amber['minimization']['DISANG'] =                    self.restraints
        self.amber['minimization']['GROUP'] =                     None    # or []

        # Production input file settings
        self.amber['md'] = {}
        self.amber['md']['cntrl'] =                       OrderedDict()
        self.amber['md']['cntrl']['imin'] =               0
        self.amber['md']['cntrl']['ntx'] =                1
        self.amber['md']['cntrl']['irest'] =              0
        self.amber['md']['cntrl']['dt'] =                 0.002
        self.amber['md']['cntrl']['nstlim'] =             5000
        self.amber['md']['cntrl']['ntpr'] =               500
        self.amber['md']['cntrl']['ntwe'] =               0
        self.amber['md']['cntrl']['ntwr'] =               5000
        self.amber['md']['cntrl']['ntxo'] =               1
        self.amber['md']['cntrl']['ntwx'] =               500
        self.amber['md']['cntrl']['ioutfm'] =             1
        self.amber['md']['cntrl']['ntf'] =                2
        self.amber['md']['cntrl']['ntc'] =                2
        self.amber['md']['cntrl']['igb'] =                1
        self.amber['md']['cntrl']['cut'] =                999.0
        self.amber['md']['cntrl']['ntt'] =                3
        self.amber['md']['cntrl']['gamma_ln'] =           1.0
        self.amber['md']['cntrl']['ig'] =                 -1
        self.amber['md']['cntrl']['temp0'] =              298.15
        self.amber['md']['cntrl']['ntr'] =                None
        self.amber['md']['cntrl']['restraint_wt'] =       None
        self.amber['md']['cntrl']['restraintmask'] =      None
        self.amber['md']['cntrl']['nmropt'] =             1
        self.amber['md']['cntrl']['pencut'] =             -1
        self.amber['md']['ewald'] =                       None    # or OrderedDict() (Add the other namelists?)
        self.amber['md']['wt'] =                          None
        self.amber['md']['DISANG'] =                      self.restraints
        self.amber['md']['GROUP'] =                       None

        self.openmm                                       = {}
        self.openmm['forcfield']                          = None
        self.openmm['minimization']                       = {}
        # Define the platform to use: CUDA, OpenCL, CPU, or Reference.
        self.openmm['minimization']['platform']           = 'CUDA'
        self.openmm['minimization']['max_iterations']     = 5000
        self.openmm['minimization']['reporter_frequency'] = 1000
        self.openmm['minimization']['tolerance']          = 1
        self.openmm['minimization']['solvent']            = app.HCT  # Hawkins-Cramer-Truhlar (igb = 1 in AMBER)
        self.openmm['minimization']['nonbonded_method']   = app.NoCutoff
        self.openmm['minimization']['nonbonded_cutoff']   = None
        self.openmm['minimization']['constraints']        = app.HBonds

        #self.print_time = 1.0 #ps.  Figure out a way to update dicts if this is reset by user. 
        # offer way for user to just provide the input file?

    def minimize(self, soft=False):
        """
        Minimize the system.

            If soft=True, slowly turn on non-bonded interactions during minimization
            so that the restraints get enforced first.

        """

        if soft:
            self.amber['minimization']['wt'] = [ "&wt type = 'NB', istep1=0, istep2=2000, value1 = 0.0, value2=0.0, IINC=50, /",
                                                "&wt type = 'NB', istep1=2000, istep2=4000, value1 = 0.0, value2=1.0, IINC=50, /" ]

        _write_input_file(self.path+'/gb_minimize.in', self.amber['minimization'], title='GB Minimization.')

        ### ToDo: Use nice error catching method dave figured out.
        log.debug('Running GB Minimization at {}'.format(self.path))
        sp.call("pmemd -O -p {0} -c {1} -ref {1} -i gb_minimize.in\
                 -o mdout.gb_minimize -r rst.gb_minimize -inf mdinfo.gb_minimize".format(self.prmtop,
                 self.coordinates), cwd=self.path, shell=True)

    def setup_openmm_restraints(self, system, restraint, phase, window):
        """
        Add particle restraints with OpenMM.
        This should probably go into `restraints.py`.
        """
        if restraints is None:
            log.error('Restraints requested but unable to setup restraints in OpenMM.')
            sys.exit(1)
        # http://docs.openmm.org/7.1.0/api-c++/generated/OpenMM.CustomExternalForce.html
        # It's possible we might need to use `periodicdistance`.
        positional_restraint = mm.CustomExternalForce('k * ((x - x_0)**2 + (y - y_0)**2 + ' \
                                                      '(z - z_0)**2)')
        bond_restraint       = mm.CustomBondForce('k * (r - r_0)**2')
        angle_restraint      = mm.CustomAngleForce('k * (theta - theta_0)**2')
        dihedral_restraint   = mm.CustomTorsionForce('k * (theta - theta_0**2')

        positional_restraint.addPerParticleParameter('k')
        positional_restraint.addPerParticleParameter('x_0')
        positional_restraint.addPerParticleParameter('y_0')
        positional_restraint.addPerParticleParameter('z_0')
        
        bond_restraint.addPerBondParameter('k')
        bond_restraint.addPerBondParameter('r_0')

        angle_restraint.addPerAngleParameter('k')
        angle_restraint.addPerAngleParameter('theta_0')

        dihedral_restraint.addPerTorsionParameter('k')
        dihedral_restraint.addPerTorsionParameter('theta_0')

        # Maybe we can move this section somewhere else and create the restraints before the 
        # OpenMM system is defined.
        system.addForce(positional_restraint)
        system.addForce(bond_restraint)
        system.addForce(angle_restraint)
        system.addForce(dihedral_restraint)

        # Have to figure out if group restraint...

        if restraint.mask1 is group:

        if restraint.mask1 is not None and restraint.mask2 is not None and \
        restraint.mask3 is None and restraint.mask4 is None:
            bond_restraint = True
        if restraint.mask1 is not None and restraint.mask2 is not None and \
            restraint.mask3 is not None and restraint.mask4 is None:
            angle_restraint = True
        if restraint.mask1 is not None and restraint.mask2 is not None and \
            restraint.mask3 is not None and restraint.mask4 is not None:
            torsion_restraint = True
        
        if bond_restraint:
            bond_restraint.addBond(restraint.index1, restraint.index2, 
                                   [restraint.phase[phase]['force_constants'][window],
                                   [restraint.phase[phase]['targets'][window]]])
        

    def minimize_openmm(self, soft=False):
        """
        Run MD with OpenMM.
        (The `soft` option can probably be folded into a class variable to tidy this up.)
        """
        prmtop = app.AmberPrmtopFile(self.toplogy)
        inpcrd = app.AmberInpcrdFile(self.coordinates)
        platform = mm.Platform.getPlatformByName(self.openmm['minimization']['platform'])
        prop = dict(CudaPrecision='mixed')


        # Call setup_openmm_restraints?



        if self.openmm['forcfield'] is not None:
            # In this case, it may be easier to not go through AMBER file types, but for now
            # using an alternative file still uses the AMBER file types as intermediaries.
            log.warn('We haven\'t tested running OpenMM with an external force field yet.')
            forcefield = app.ForceField(self.opennmm['forcefield'])
            log.warn('Probably need to load a separate topology here...')
            system = forcefield.createSystem(
                implicitSolvent=self.openmm['minimization']['solvent'],
                nonbondedMethod=self.openmm['minimization']['nonbonded_method'],
                constraints=self.openmm['minimization']['constraints'],
                platform=platform,
                prop=prop)
        else:    
            system = prmtop.createSystem(
                implicitSolvent=self.openmm['minimization']['solvent'],
                nonbondedMethod=self.openmm['minimization']['nonbonded_method'],
                constraints=self.openmm['minimization']['constraints'],
                platform=platform,
                prop=prop)
        simulation = app.Simulation(prmtop.topology, system)
        simulation.context.setPositions(inpcrd)

        if soft:
            # Phase 1: minimize with nonbonded interactions disabled.
            # This is the first 40% of the maximum iterations.
            log.debug('Minimization phase 1 for {} steps.'
                .format(int(0.4 * self.openmm['minimization']['max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.4 * self.openmm['minimization']['max_iterations']),
                tolerance=self.openmm['minimization']['tolerance']*unit.kilojoule/unit.mole
                )
            # Phase 2: slowly turn on the nonbonded interactions.
            # This is the next 40% of the maximum iterations.
            # This increases the nonbonded interactions linearly, which is not
            # the same as using `IINC` in AMBER.
            log.debug('Minimization phase 2 for {} steps.'
                .format(int(0.4 * self.openmm['minimization']['max_iterations'])))
            for scale in np.linspace(0, 1.0, int(0.4 * self.openmm['minimization']['max_iterations'])):
                log.debug('Scaling NB interactions to {0:0.4f} / 1.0'.format(scale))
                for particle in range(system.getNumParticles()):
                    [charge, sigma, epsilon] = mm.NonbondedForce.getParticleParameters(particle)
                    mm.NonbondedForce.setParticleParameters(particle, charge * scale, sigma, epsilon * scale)
                simulation.minimizeEnergy(
                    maxIterations=1,
                    tolerance=self.openmm['minimization']['tolerance']*unit.kilojoule/unit.mole
                )
            # Phase 3: minimize with nonbonded interactions at full strength.
            # This is the last 20% of the maximum iterations.
            log.debug('Minimization phase 3 for {} steps.'
                .format(int(0.2 * self.openmm['minimization']['max_iterations'])))
            simulation.minimizeEnergy(
                maxIterations=int(0.2 * self.openmm['minimization']['max_iterations']),
                tolerance=self.openmm['minimization']['tolerance']*unit.kilojoule/unit.mole
                )
            state = simulation.context.getState(getEnergy=True, getPositions=True)
        else:
            simulation.minimizeEnergy(
                maxIterations=self.openmm['minimization']['max_iterations'],
                tolerance=self.openmm['minimization']['tolerance']*unit.kilojoule/unit.mole
                )



    def run_md(self, input_file, input_crds, output_suffix, ref_crds=None):
        """
        Run MD.
        """

        _write_input_file(self.path+'/'+input_file, self.amber['md'], title='MD.')

        if ref_crds is None:
            ref_crds = input_crds

        ### ToDo: Use nice error catching method dave figured out.
        ### ToDo: Is this the best way to do file naming?
        log.debug('Running MD at {}'.format(self.path))
        sp.call("pmemd -O -p {0} -c {1} -ref {2} -i {3} -o mdout.{4} -r rst.{4} -x traj.{4} -e mden.{4} -inf mdinfo.{4}\
                ".format(self.toplogy,self.coordinates,ref_crds,input_file,output_suffix), cwd=self.path, shell=True)


    def run_md_openmm(self):
        """
        Run MD with OpenMM.
        """
        if restraints:
            # Figure out if there are restraints.
            pass

        else:
            openmm = pmd.load_file(self.toplogy, self.coordinates)
            system = openmm.createSystem(
                nonbondedMethod=app.PME,
                constraints=app.HBonds
            )
            integrator = mm.Langevin

    def amber_to_pdb(self):
        
        # To run with the SMIRNOFF force field, we must use PDB as an intermediary format.
        # This will print CONECT records between H1 and H2, if using TIP3P. Possibly others.
        self.pdb = '.'.join(self.toplogy.split('.')[0:-1]) + '.pdb'
        pdb_input = '.'.join(self.toplogy.split('.')[0:-1]) + '.in'
        log.info('Converting AMBER coordinates and topology to PDB format.')
        log.debug('Calling `cpptraj`...')
        with open(pdb_input, 'w') as file:
            file.write('parm {}\n'.format(self.toplogy))
            file.write('trajin {}\n'.format(self.coordinates))
            file.write('trajout {} conect\n'.format(self.pdb))
        sp.check_call(['cpptraj', '-i', pdb_input])