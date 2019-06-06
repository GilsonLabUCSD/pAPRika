"""
This class is going to contain the simulation setup...
"""

import pkg_resources
import shutil
import textwrap
import parmed as pmd
import os as os
import subprocess as sp
import simtk.openmm as openmm
import simtk.unit as unit

from pathlib import Path
from paprika import align
from paprika.dummy import add_dummy
from paprika.tleap import System
from paprika.restraints import static_DAT_restraint, DAT_restraint
from paprika.restraints.restraints import create_window_list
from paprika.restraints.read_yaml import read_yaml
from paprika.io import save_restraints

import logging

logger = logging.getLogger(__name__)


def _get_installed_benchmarks():
    _installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        _installed_benchmarks[entry_point.name] = entry_point.load()
    return _installed_benchmarks

def read_openmm_system_from_xml(filename):
    with open(filename, "rb") as file:
        return openmm.XmlSerializer.deserialize(file.read().decode())

class Setup(object):
    """
    The Setup class provides a wrapper function around the preparation of the host-guest system and the application of restraints.
    """

    def __init__(self, host, guest, backend="openmm"):
        self.host = host
        self.guest = guest
        self.backend = backend

        self.directory = Path("benchmarks").joinpath(self.host).joinpath(self.guest)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.window_list = []

        if self.backend == "amber":
            raise NotImplementedError

        installed_benchmarks = get_benchmarks()
        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks)
        self.benchmark_path = host_yaml.parent
        self.host_yaml = read_yaml(host_yaml)
        self.guest_yaml = read_yaml(guest_yaml)

        self.build_desolvated_windows()
        # Read in XML from Simon
        self.read_solvated_windows()
        self.dummy_atom_indices = self.add_dummy_atoms()
        # Initialize the restraints
        self.static_restraints, \
        self.conformational_restraints, \
        self.wall_restraints, \
        self.guest_restraints = self.initialize_restraints()
        for window in self.window_list:
            self.initialize_calculation(window)
        # self.build_bound()
        # self.static_restraints, \
        # self.conformational_restraints, \
        # self.wall_restraints, \
        # self.guest_restraints = self.initialize_restraints()
        # self.build_windows()
        # self.openmm_wrapper()

    def parse_yaml(self, installed_benchmarks):
        """
        Read the YAML recipe for the host and guest.

        Returns
        -------

        """
        try:
            host_yaml = installed_benchmarks["host_guest_pairs"][self.host]["yaml"]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for host: {self.host}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        try:
            guest_yaml = installed_benchmarks["host_guest_pairs"][self.host][self.guest]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        return host_yaml, guest_yaml

    def build_bound(self):
        """
        This is going to change to accommodate SMIRNOFF99Frosst.

        Returns
        -------
        """
        system = System()
        system.output_path = self.directory
        system.output_prefix = f"{self.host}-{self.guest}-unaligned"
        system.pbc_type = None
        system.neutralize = False
        system.template_lines = [
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
            f"loadamberparams {self.benchmark_path.joinpath(self.host_yaml['frcmod'])}",
            f"loadamberparams {self.benchmark_path.joinpath('dummy.frcmod')}",
            f"{self.host.upper()} = loadmol2 {self.benchmark_path.joinpath(self.host_yaml['structure'])}",
            f"{self.guest.upper()} = loadmol2 {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['structure'])}",
            f"DM1 = loadmol2 {self.benchmark_path.joinpath('dm1.mol2')}",
            f"DM2 = loadmol2 {self.benchmark_path.joinpath('dm2.mol2')}",
            f"DM3 = loadmol2 {self.benchmark_path.joinpath('dm3.mol2')}",
            f"model = loadpdb {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['complex'])}",
        ]
        system.build()

        source_prmtop = self.directory.joinpath(f"{self.host}-{self.guest}-unaligned.prmtop")
        source_inpcrd = self.directory.joinpath(f"{self.host}-{self.guest}-unaligned.rst7")
        structure = pmd.load_file(str(source_prmtop), str(source_inpcrd),
                                  structure=True)
        guest_angle_restraint = self.guest_yaml["restraints"][-1]["restraint"]["atoms"].split()                          
        aligned_structure = align.zalign(structure, guest_angle_restraint[1], guest_angle_restraint[2])
        destination_prmtop = self.directory.joinpath(f"{self.host}-{self.guest}.prmtop")
        destination_inpcrd = self.directory.joinpath(f"{self.host}-{self.guest}.rst7")
        destination_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")
        aligned_structure.save(str(destination_prmtop), overwrite=True)
        aligned_structure.save(str(destination_inpcrd), overwrite=True)
        aligned_structure.save(str(destination_pdb), overwrite=True)

    def align(self, coordinates, topology):
        structure = pmd.load_file(str(coordinates), structure=True)
        intermediate_pdb = self.directory.joinpath(f"tmp.pdb")
        destination_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

        guest_angle_restraint_mask = self.guest_yaml["restraints"][-1]["restraint"]["atoms"].split()
        aligned_structure = align.zalign(structure, guest_angle_restraint_mask[1], guest_angle_restraint_mask[2])
        aligned_structure.save(str(intermediate_pdb), overwrite=True)
        p = create_pdb_with_conect(intermediate_pdb,
                               amber_prmtop=topology,
                               output_pdb=destination_pdb)
        if p:
            os.remove(intermediate_pdb)


    def _create_dummy_restraint(self, initial_structure):

        windows = [
            self.host_yaml["calculation"]["windows"]["attach"],
            self.host_yaml["calculation"]["windows"]["pull"],
            None,
        ]
        restraint = self.guest_yaml["restraints"][0]

        guest_restraint = DAT_restraint()
        guest_restraint.auto_apr = True
        guest_restraint.continuous_apr = True
        guest_restraint.amber_index = False if self.backend == "openmm" else True
        guest_restraint.topology = str(initial_structure)
        guest_restraint.mask1 = "@1"
        guest_restraint.mask2 = "@2"

        guest_restraint.attach["target"] = restraint["restraint"]["attach"][
            "target"
        ]
        guest_restraint.attach["fc_final"] = restraint["restraint"]["attach"][
            "force_constant"
        ]
        guest_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
            "lambda"
        ]["attach"]

        guest_restraint.pull["target_final"] = self.host_yaml["calculation"][
            "target"
        ]["pull"]
        guest_restraint.pull["num_windows"] = windows[1]

        guest_restraint.initialize()

        return guest_restraint

    def build_desolvated_windows(self):
        initial_structure = self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['complex'])
        initial_topology = self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['prmtop'])
        self.align(coordinates=initial_structure,
                   topology=initial_topology)
        logger.debug("Setting up dummy restraint to build window list.")
        _dummy_restraint = self._create_dummy_restraint(initial_structure = str(initial_structure))
        self.window_list = create_window_list([_dummy_restraint])

        for window in self.window_list:
            logger.debug(f"Translating guest in window {window}...")
            self.directory.joinpath("windows").joinpath(window).mkdir(parents=True, exist_ok=True)
            self.translate(window, restraint = _dummy_restraint)


    def read_solvated_windows(self):
        for window in self.window_list:
            logger.debug(f"Reading solvated OpenMM System...")
            try:
                system = read_openmm_system_from_xml(self.directory.joinpath("windows").joinpath(window).joinpath(
                "system.xml"))

            except:
                logger.warning(f"Missing system.xml in {window}")


    def add_dummy_atoms(self):
        # First add dummy atoms to structure
        guest_angle_restraint_mask = self.guest_yaml["restraints"][-1]["restraint"]["atoms"].split()
        for window in self.window_list:
            logger.debug(f"Adding dummy atoms to structures...")
            try:
                structure = pmd.load_file(str(self.directory.joinpath("windows").joinpath(window).joinpath(
                    "output.pdb")))
                offset_coordinates = structure[guest_angle_restraint_mask[1]].coordinates

                structure = add_dummy(structure, x = offset_coordinates[0][0],
                          y = offset_coordinates[0][1],
                          z = offset_coordinates[0][2] - 6.0,
                                      residue_name="DM1")
                structure = add_dummy(structure, x = offset_coordinates[0][0],
                          y = offset_coordinates[0][1],
                          z = offset_coordinates[0][2] - 9.0,
                                      residue_name="DM2")
                structure = add_dummy(structure, x = offset_coordinates[0][0],
                        y = offset_coordinates[0][1] + 2.2,
                        z = offset_coordinates[0][2] - 11.2,
                                      residue_name="DM3")
                structure.save(str(self.directory.joinpath("windows").joinpath(window).joinpath(
                    "dummy.pdb")), overwrite=True)
            except:
                logger.warning(f"Missing output.pdb in {window}")

            # Add dummy atoms to System
            try:
                system = read_openmm_system_from_xml(self.directory.joinpath("windows").joinpath(window).joinpath(
                    "system.xml"))
                dummy_atom_indices = []
                dummy_atom_indices.append([system.addParticle(mass = 207) for _ in range(3)])

                for force_index in range(system.getNumForces()):
                    force = system.getForce(force_index)
                    if not isinstance(force, openmm.NonbondedForce):
                        continue
                    force.addParticle(0.0, 1.0, 0.0)
                    force.addParticle(0.0, 1.0, 0.0)
                    force.addParticle(0.0, 1.0, 0.0)

                for atom in structure.atoms:
                    if atom.name == "DUM":

                        positional_restraint = openmm.CustomExternalForce('k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
                        positional_restraint.addPerParticleParameter('k')
                        positional_restraint.addPerParticleParameter('x0')
                        positional_restraint.addPerParticleParameter('y0')
                        positional_restraint.addPerParticleParameter('z0')
                        # I haven't found a way to get this to use ParmEd's unit library here.
                        # ParmEd correctly reports `atom.positions` as units of Ångstroms.
                        # But then we can't access atom indices.
                        # Using `atom.xx` works for coordinates, but is unitless.

                        k = 50.0 * unit.kilocalories_per_mole / unit.angstroms**2
                        x0 = 0.1 * atom.xx * unit.nanometers
                        y0 = 0.1 * atom.xy * unit.nanometers
                        z0 = 0.1 * atom.xz * unit.nanometers
                        positional_restraint.addParticle(atom.idx, [k, x0, y0, z0])
                        system.addForce(positional_restraint)

                system_xml = openmm.XmlSerializer.serialize(system)
                with open(self.directory.joinpath("windows").joinpath(window).joinpath("dummy.xml"), "w") as file:
                    file.write(system_xml)
            except:
                logger.warning(f"Missing system.xml")

        return dummy_atom_indices

    def initialize_restraints(self):

        windows = [
            self.host_yaml["calculation"]["windows"]["attach"],
            self.host_yaml["calculation"]["windows"]["pull"],
            None,
        ]
        # structure = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")
        structure = self.directory.joinpath("windows").joinpath("a000").joinpath("dummy.pdb")

        static_restraints = []
        for restraint in self.host_yaml["calculation"]["restraints"]["static"]:
            static = static_DAT_restraint(
                restraint_mask_list=restraint["restraint"]["atoms"].split(),
                num_window_list=windows,
                ref_structure=str(structure),
                force_constant=restraint["restraint"]["force_constant"],
                amber_index=False if self.backend == "openmm" else True,
            )
            static_restraints.append(static)

        conformational_restraints = []
        if self.host_yaml["calculation"]["restraints"]["conformational"]:
            for conformational in self.host_yaml["calculation"]["restraints"][
                "conformational"
            ]:
                raise NotImplementedError
        else:
            logger.debug("Skipping conformational restraints...")

        wall_restraints = []
        if self.host_yaml["calculation"]["restraints"]["wall"]:
            for wall in self.host_yaml["calculation"]["restraints"]["wall"]:
                raise NotImplementedError
        else:
            logger.debug("Skipping wall restraints...")

        guest_restraints = []
        for restraint in self.guest_yaml["restraints"]:
            mask = restraint["restraint"]["atoms"].split()

            guest_restraint = DAT_restraint()
            guest_restraint.auto_apr = True
            guest_restraint.continuous_apr = True
            guest_restraint.amber_index = False if self.backend == "openmm" else True
            guest_restraint.topology = str(structure)
            guest_restraint.mask1 = mask[0]
            guest_restraint.mask2 = mask[1]
            guest_restraint.mask3 = mask[2] if len(mask) > 2 else None
            guest_restraint.mask4 = mask[3] if len(mask) > 3 else None

            guest_restraint.attach["target"] = restraint["restraint"]["attach"][
                "target"
            ]
            guest_restraint.attach["fc_final"] = restraint["restraint"]["attach"][
                "force_constant"
            ]
            guest_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
                "lambda"
            ]["attach"]

            guest_restraint.pull["target_final"] = self.host_yaml["calculation"][
                "target"
            ]["pull"]
            guest_restraint.pull["num_windows"] = windows[1]

            guest_restraint.initialize()
            guest_restraints.append(guest_restraint)

        return static_restraints,        conformational_restraints,        wall_restraints,         guest_restraints

    def build_windows(self):
        save_restraints(self.static_restraints + self.conformational_restraints + self.wall_restraints + self.guest_restraints, self.directory.joinpath("restraints.json"))
        window_list = create_window_list(self.guest_restraints)
        
        for window in window_list:
            self.directory.joinpath("windows").joinpath(window).mkdir(parents=True, exist_ok=True)
            self.translate(window)
            self.solvate(window)

            if self.backend == "amber":
                # Write the restraint file in each window.
                raise NotImplementedError
            # from paprika.restraints.amber_restraints import amber_restraint_line
            # for window in window_list:
            #     with open(self.directory.joinpath("windows").joinpath(window).joinpath("disang.rest"), "a") as file:
            #         for restraint in self.static_restraints + self.conformational_restraints + self.wall_restraints + self.guest_restraints:
            #             string = amber_restraint_line(restraint, window)
            #             if string is not None:
            #                 file.write(string)

            else:
                self.initialize_calculation(window)
    
    def translate(self, window, restraint):
        window_path = self.directory.joinpath("windows").joinpath(window)
        if window[0] == "a":
            # Copy the initial structure.
            source_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")
            shutil.copy(source_pdb, window_path)
        elif window[0] == "p":
            # Translate the guest.
            source_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

            structure = pmd.load_file(str(source_pdb), structure = True)
            target_difference = restraint.phase['pull']['targets'][int(window[\
                1:])] - restraint.pull['target_initial']
            for atom in structure.atoms:
                if atom.residue.name == self.guest.upper():
                    atom.xz += target_difference

            intermediate_pdb = window_path.joinpath(f"tmp.pdb")
            destination_pdb = window_path.joinpath(f"{self.host}-{self.guest}.pdb")
            structure.save(str(intermediate_pdb), overwrite=True)

            p = create_pdb_with_conect(intermediate_pdb,
                                   amber_prmtop=self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['prmtop']),
                                   output_pdb=destination_pdb)
            if p:
                os.remove(intermediate_pdb)



        elif window[0] == "r":
            # Copy the final pull window.
            source_pdb = self.directory.joinpath("windows").joinpath(f"p{self.host_yaml['calculation']['windows']['pull']:03d}").joinpath(f"{self.host}-{self.guest}.pdb")
            shutil.copy(source_pdb, window_path)

    def solvate(self, window):
        window_path = self.directory.joinpath("windows").joinpath(window)
        source_prmtop = window_path.joinpath(f"{self.host}-{self.guest}.prmtop")
        source_inpcrd = window_path.joinpath(f"{self.host}-{self.guest}.rst7")
        pdb = window_path.joinpath(f"{self.host}-{self.guest}.pdb")
        
        if not Path.exists(pdb):
            structure = pmd.load_file(str(source_prmtop), str(source_inpcrd), structure = True)
            structure.save(str(pdb), overwrite=True)
        if Path.exists(window_path.joinpath(f"{self.host}-{self.guest}-sol.prmtop")):
            logger.debug("Found solvated structure. Skipping.")
            return
    
        system = System()
        system.output_path = window_path
        system.output_prefix = f"{self.host}-{self.guest}-sol"
        
        system.target_waters = self.host_yaml["calculation"]["waters"]
        system.neutralize = True
        # system.add_ions = ["Na+", 6, "Cl-", 6]
        # This should be read from the YAML file.
        system.template_lines = [
        "source leaprc.gaff",
        "source leaprc.water.tip3p",
        f"loadamberparams {self.benchmark_path.joinpath(self.host_yaml['frcmod'])}",
        f"loadamberparams {self.benchmark_path.joinpath('dummy.frcmod')}",
        f"{self.host.upper()} = loadmol2 {self.benchmark_path.joinpath(self.host_yaml['structure'])}",
        f"{self.guest.upper()} = loadmol2 {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['structure'])}",
        f"DM1 = loadmol2 {self.benchmark_path.joinpath('dm1.mol2')}",
        f"DM2 = loadmol2 {self.benchmark_path.joinpath('dm2.mol2')}",
        f"DM3 = loadmol2 {self.benchmark_path.joinpath('dm3.mol2')}",
        f"model = loadpdb {self.host}-{self.guest}.pdb",
        ]
        system.build()

    def initialize_calculation(self, window): 
        if self.backend == "amber":
            # Write simulation input files in each directory
            raise NotImplementedError

        window_path = self.directory.joinpath("windows").joinpath(window)
        system = read_openmm_system_from_xml(window_path.joinpath("dummy.xml"))

        for restraint in self.static_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=10)
        for restraint in self.conformational_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=11)
        for restraint in self.guest_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=12)
        for restraint in self.wall_restraints:
            # Custom force here...
            raise NotImplementedError

        system_xml=openmm.XmlSerializer.serialize(system)

        with open(window_path.joinpath("restraints.xml"), "w") as file:
            file.write(system_xml)

def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks

def apply_openmm_restraints(system, restraint, window, ForceGroup=None):
    if window[0] == "a":
        phase = "attach"
    elif window[0] == "p":
        phase = "pull"
    elif window[0] == "r":
        phase = "release"
    window_number = int(window[1:])

    if restraint.mask2 and not restraint.mask3:
        if not restraint.group1 and not restraint.group2:
            bond_restraint = openmm.CustomBondForce('k * (r - r_0)^2')
            bond_restraint.addPerBondParameter('k')
            bond_restraint.addPerBondParameter('r_0')

            r_0 = restraint.phase[phase]["targets"][window_number] * unit.angstroms
            k = restraint.phase[phase]["force_constants"][window_number] * unit.kilocalories_per_mole / \
                 unit.angstrom**2
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])
            system.addForce(bond_restraint)
        else:
            # Probably needs openmm.CustomCentroidBondForce (?)
            raise NotImplementedError
        if ForceGroup:
            bond_restraint.setForceGroup(ForceGroup)

    elif restraint.mask3 and not restraint.mask4:
        if not restraint.group1 and not restraint.group2 and not restraint.group3:
            angle_restraint = openmm.CustomAngleForce('k * (theta - theta_0)^2')
            angle_restraint.addPerAngleParameter('k')
            angle_restraint.addPerAngleParameter('theta_0')

            theta_0 = restraint.phase[phase]["targets"][window_number] * unit.degrees
            k = restraint.phase[phase]["force_constants"][window_number] * unit.kilocalories_per_mole / unit.radian**2
            angle_restraint.addAngle(restraint.index1[0], restraint.index2[0], restraint.index3[0],
            [k, theta_0])
            system.addForce(angle_restraint)
        else:
            # Probably needs openmm.CustomCentroidAngleForce (?)
            raise NotImplementedError
        if ForceGroup:
            angle_restraint.setForceGroup(ForceGroup)

    elif restraint.mask4:
        if not restraint.group1 and not restraint.group2 and not restraint.group3 and not restraint.group4:
            dihedral_restraint = openmm.CustomTorsionForce('k * (theta - theta_0)^2')
            dihedral_restraint.addPerTorsionParameter('k')
            dihedral_restraint.addPerTorsionParameter('theta_0')

            theta_0 = restraint.phase[phase]["targets"][window_number] * unit.degrees
            k = restraint.phase[phase]["force_constants"][window_number] * unit.kilocalories_per_mole / unit.radian**2
            dihedral_restraint.addTorsion(restraint.index1[0], restraint.index2[0], restraint.index3[0],
                                          restraint.index4[0],
            [k, theta_0])
            system.addForce(dihedral_restraint)
        else:
            # Probably needs openmm.CustomCentroidTorsionForce (?)
            raise NotImplementedError
        if ForceGroup:
            dihedral_restraint.setForceGroup(ForceGroup)
    return system



def create_pdb_with_conect(input_pdb, amber_prmtop, output_pdb):
    """
    Create a PDB file containing CONECT records.
    `cpptraj` must be in your PATH.
    Parameters
    ----------
    input_pdb : str
        Existing structure from e.g., Mobley's Benchmark Sets repository
    amber_prmtop : str
        AMBER (or other) parameters for the residues in the solvated PDB file
    output_pdb : str
        Output PDB file name
    path : str
        Directory for input and output files
    """
    logger.info(f'Creating {output_pdb} with CONECT records...')
    cpptraj = \
        f'''
    parm {amber_prmtop.resolve()}
    trajin {input_pdb.resolve()}
    trajout {output_pdb.resolve()} conect
    '''
    path = output_pdb.resolve().parent

    cpptraj_input = output_pdb.resolve().with_suffix(".in")
    cpptraj_output = output_pdb.resolve().with_suffix(".out")
    logger.debug(f"Writing {cpptraj_input}")

    with open(cpptraj_input, 'w') as file:
        file.write(cpptraj)
    with open(cpptraj_output, 'w') as file:
        p = sp.Popen(
            ['cpptraj', '-i', cpptraj_input],
            cwd=path,
            stdout=file,
            stderr=file)
        output, error = p.communicate()
    if p.returncode == 0:
        logger.debug('PDB file written by cpptraj.')
        return True
    elif p.returncode == 1:
        logger.error('Error returned by cpptraj.')
        logger.error(f'Output: {output}')
        logger.error(f'Error: {error}')
        p = sp.Popen(['cat', cpptraj_output], cwd=path, stdout=sp.PIPE)
        for line in p.stdout:
            logger.error(line.decode("utf-8").strip(), )
    else:
        logger.error(f'Output: {output}')
        logger.error(f'Error: {error}')