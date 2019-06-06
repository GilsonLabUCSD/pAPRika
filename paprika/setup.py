"""
This class is going to contain the simulation setup...
"""

import pkg_resources
import shutil
import parmed as pmd
import os as os
import simtk.openmm as openmm
import simtk.unit as unit


from pathlib import Path
from paprika import align
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

    def __init__(self, host, guest, backend="openmm", directory_path="benchmarks"):
        self.host = host
        self.guest = guest
        self.backend = backend
        self.directory = Path(directory_path).joinpath(self.host).joinpath(self.guest)
        self.desolvated_window_paths = []
        self.window_list = []

        if self.backend == "amber":
            raise NotImplementedError

        self.directory.mkdir(parents=True, exist_ok=True)
        installed_benchmarks = get_benchmarks()
        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks)
        self.benchmark_path = host_yaml.parent
        self.host_yaml = read_yaml(host_yaml)
        self.guest_yaml = read_yaml(guest_yaml)

        # Here, we build desolvated windows and pass the files to the Property Estimator.
        # These files are stored in `self.desolvated_window_paths`.
        self.build_desolvated_windows()
        # Now, we read in the solvated windows from the Property Estimator
        self.add_dummy_atoms()
        self.static_restraints, self.conformational_restraints, self.wall_restraints, self.guest_restraints = (
            self.initialize_restraints()
        )
        for window in self.window_list:
            self.initialize_calculation(window)

    def parse_yaml(self, installed_benchmarks):
        """
        Read the YAML recipe for the host and guest.

        Returns
        -------

        """
        try:
            host_yaml = installed_benchmarks["host_guest_systems"][self.host]["yaml"]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for host: {self.host}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        try:
            guest_yaml = installed_benchmarks["host_guest_systems"][self.host][self.guest]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError

        return host_yaml, guest_yaml

    def align(self, input_pdb):
        structure = pmd.load_file(str(input_pdb), structure=True)
        intermediate_pdb = self.directory.joinpath(f"tmp.pdb")
        destination_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

        guest_angle_restraint_mask = self.guest_yaml["restraints"][-1]["restraint"][
            "atoms"
        ].split()
        aligned_structure = align.zalign(
            structure, guest_angle_restraint_mask[1], guest_angle_restraint_mask[2]
        )
        aligned_structure.save(str(intermediate_pdb), overwrite=True)

        # Save aligned PDB file with CONECT records.
        positions_pdb = openmm.app.PDBFile(str(intermediate_pdb))
        topology_pdb = openmm.app.PDBFile(str(input_pdb))

        positions = positions_pdb.positions
        topology = topology_pdb.topology
        with open(destination_pdb, "w") as file:
            openmm.app.PDBFile.writeFile(topology, positions, file)
        os.remove(intermediate_pdb)

    def build_desolvated_windows(self):
        initial_structure = self.benchmark_path.joinpath(self.guest).joinpath(
            self.guest_yaml["complex"]
        )
        self.align(input_pdb=initial_structure)
        logger.debug("Setting up dummy restraint to build window list.")
        _dummy_restraint = self._create_dummy_restraint(
            initial_structure=str(initial_structure)
        )
        self.window_list = create_window_list([_dummy_restraint])

        for window in self.window_list:
            logger.debug(f"Translating guest in window {window}...")
            self.directory.joinpath("windows").joinpath(window).mkdir(
                parents=True, exist_ok=True
            )
            self.translate(window, topology_pdb=initial_structure, restraint=_dummy_restraint)

            window_pdb_file_name = f"{self.host}-{self.guest}.pdb"

            self.desolvated_window_paths.append(
                str(
                    self.directory.joinpath("windows")
                    .joinpath(window)
                    .joinpath(window_pdb_file_name)
                )
            )

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

        guest_restraint.attach["target"] = restraint["restraint"]["attach"]["target"]
        guest_restraint.attach["fc_final"] = restraint["restraint"]["attach"][
            "force_constant"
        ]
        guest_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
            "lambda"
        ]["attach"]

        guest_restraint.pull["target_final"] = self.host_yaml["calculation"]["target"][
            "pull"
        ]
        guest_restraint.pull["num_windows"] = windows[1]

        guest_restraint.initialize()

        return guest_restraint

    def translate(self, window, topology_pdb, restraint):
        window_path = self.directory.joinpath("windows").joinpath(window)
        if window[0] == "a":
            # Copy the initial structure.
            source_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")
            shutil.copy(source_pdb, window_path)
        elif window[0] == "p":
            # Translate the guest.
            source_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

            structure = pmd.load_file(str(source_pdb), structure=True)
            target_difference = (
                restraint.phase["pull"]["targets"][int(window[1:])]
                - restraint.pull["target_initial"]
            )
            for atom in structure.atoms:
                if atom.residue.name == self.guest.upper():
                    atom.xz += target_difference

            intermediate_pdb = window_path.joinpath(f"tmp.pdb")
            destination_pdb = window_path.joinpath(f"{self.host}-{self.guest}.pdb")
            structure.save(str(intermediate_pdb), overwrite=True)


            input_pdb = openmm.app.PDBFile(str(intermediate_pdb))
            topology_pdb = openmm.app.PDBFile(str(topology_pdb))
            positions = input_pdb.positions
            topology = topology_pdb.topology
            with open(destination_pdb, "w") as file:
                openmm.app.PDBFile.writeFile(topology, positions, file)
            os.remove(intermediate_pdb)

        elif window[0] == "r":
            # Copy the final pull window.
            source_pdb = (
                self.directory.joinpath("windows")
                .joinpath(f"p{self.host_yaml['calculation']['windows']['pull']:03d}")
                .joinpath(f"{self.host}-{self.guest}.pdb")
            )
            shutil.copy(source_pdb, window_path)

    def add_dummy_atoms(
        self,
        input_pdb="output.pdb",
        input_xml="system.xml",
        output_pdb="output.pdb",
        output_xml="output.xml",
    ):
        # First add dummy atoms to structure
        guest_angle_restraint_mask = self.guest_yaml["restraints"][-1]["restraint"][
            "atoms"
        ].split()
        logger.debug(f"Adding dummy atoms to {input_pdb}")
        try:
            structure = pmd.load_file(input_pdb)
            offset_coordinates = structure[guest_angle_restraint_mask[1]].coordinates
            input_pdb_file = openmm.app.PDBFile(input_pdb)

            positions = input_pdb_file.positions
            positions.append(
                openmm.Vec3(
                    offset_coordinates[0][0],
                    offset_coordinates[0][1],
                    offset_coordinates[0][2] - 6.0,
                )
                * unit.angstrom
            )
            positions.append(
                openmm.Vec3(
                    offset_coordinates[0][0],
                    offset_coordinates[0][1],
                    offset_coordinates[0][2] - 9.0,
                )
                * unit.angstrom
            )
            positions.append(
                openmm.Vec3(
                    offset_coordinates[0][0],
                    offset_coordinates[0][1] + 2.2,
                    offset_coordinates[0][2] - 11.2,
                )
                * unit.angstrom
            )

            topology = input_pdb_file.topology
            for dummy_index in range(1, 4):
                dummy_chain = topology.addChain(None)
                dummy_residue = topology.addResidue(f"DM{dummy_index}", dummy_chain)
                topology.addAtom(f"DUM", None, dummy_residue)

            with open(output_pdb, "w") as file:
                openmm.app.PDBFile.writeFile(topology, positions, file)

        except:
            logger.warning(f"Missing {input_pdb}")

        # Add dummy atoms to System
        try:
            system = read_openmm_system_from_xml(input_xml)
            [system.addParticle(mass=207) for _ in range(3)]

            for force_index in range(system.getNumForces()):
                force = system.getForce(force_index)
                if not isinstance(force, openmm.NonbondedForce):
                    continue
                force.addParticle(0.0, 1.0, 0.0)
                force.addParticle(0.0, 1.0, 0.0)
                force.addParticle(0.0, 1.0, 0.0)

            for atom in structure.atoms:
                if atom.name == "DUM":

                    positional_restraint = openmm.CustomExternalForce(
                        "k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
                    )
                    positional_restraint.addPerParticleParameter("k")
                    positional_restraint.addPerParticleParameter("x0")
                    positional_restraint.addPerParticleParameter("y0")
                    positional_restraint.addPerParticleParameter("z0")
                    # I haven't found a way to get this to use ParmEd's unit library here.
                    # ParmEd correctly reports `atom.positions` as units of Ångstroms.
                    # But then we can't access atom indices.
                    # Using `atom.xx` works for coordinates, but is unitless.

                    k = 50.0 * unit.kilocalories_per_mole / unit.angstroms ** 2
                    x0 = 0.1 * atom.xx * unit.nanometers
                    y0 = 0.1 * atom.xy * unit.nanometers
                    z0 = 0.1 * atom.xz * unit.nanometers
                    positional_restraint.addParticle(atom.idx, [k, x0, y0, z0])
                    system.addForce(positional_restraint)

            system_xml = openmm.XmlSerializer.serialize(system)
            with open(output_xml, "w") as file:
                file.write(system_xml)
        except:
            logger.warning(f"Missing {input_xml}")

    def initialize_restraints(self, PDB_filename="output.pdb"):

        windows = [
            self.host_yaml["calculation"]["windows"]["attach"],
            self.host_yaml["calculation"]["windows"]["pull"],
            None,
        ]
        structure = (
            self.directory.joinpath("windows").joinpath("a000").joinpath(PDB_filename)
        )

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

        return (
            static_restraints,
            conformational_restraints,
            wall_restraints,
            guest_restraints,
        )

    def initialize_calculation(self, window, input_xml="system.xml", output_xml="system.xml"):
        if self.backend == "amber":
            # Write simulation input files in each directory
            raise NotImplementedError
        try:
            system = read_openmm_system_from_xml(input_xml)
        except:
            logger.warning(f"Cannot read XML from {input_xml}")

        for restraint in self.static_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=10)
        for restraint in self.conformational_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=11)
        for restraint in self.guest_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=12)
        for restraint in self.wall_restraints:
            # Custom force here...
            raise NotImplementedError

        system_xml = openmm.XmlSerializer.serialize(system)

        with open(output_xml, "w") as file:
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
            bond_restraint = openmm.CustomBondForce("k * (r - r_0)^2")
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")

            r_0 = restraint.phase[phase]["targets"][window_number] * unit.angstroms
            k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.angstrom ** 2
            )
            bond_restraint.addBond(restraint.index1[0], restraint.index2[0], [k, r_0])
            system.addForce(bond_restraint)
        else:
            bond_restraint = openmm.CustomCentroidBondForce(
                2, "k * (distance(g1, g2) - r_0)^2"
            )
            bond_restraint.addPerBondParameter("k")
            bond_restraint.addPerBondParameter("r_0")
            r_0 = restraint.phase[phase]["targets"][window_number] * unit.angstroms
            k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.angstrom ** 2
            )
            g1 = bond_restraint.addGroup(restraint.index1)
            g2 = bond_restraint.addGroup(restraint.index2)
            bond_restraint.addBond([g1, g2], [k, r_0])
            system.addForce(bond_restraint)

        if ForceGroup:
            bond_restraint.setForceGroup(ForceGroup)

    elif restraint.mask3 and not restraint.mask4:
        if not restraint.group1 and not restraint.group2 and not restraint.group3:
            angle_restraint = openmm.CustomAngleForce("k * (theta - theta_0)^2")
            angle_restraint.addPerAngleParameter("k")
            angle_restraint.addPerAngleParameter("theta_0")

            theta_0 = restraint.phase[phase]["targets"][window_number] * unit.degrees
            k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.radian ** 2
            )
            angle_restraint.addAngle(
                restraint.index1[0],
                restraint.index2[0],
                restraint.index3[0],
                [k, theta_0],
            )
            system.addForce(angle_restraint)
        else:
            # Probably needs openmm.CustomCentroidAngleForce (?)
            raise NotImplementedError
        if ForceGroup:
            angle_restraint.setForceGroup(ForceGroup)

    elif restraint.mask4:
        if (
            not restraint.group1
            and not restraint.group2
            and not restraint.group3
            and not restraint.group4
        ):
            dihedral_restraint = openmm.CustomTorsionForce("k * (theta - theta_0)^2")
            dihedral_restraint.addPerTorsionParameter("k")
            dihedral_restraint.addPerTorsionParameter("theta_0")

            theta_0 = restraint.phase[phase]["targets"][window_number] * unit.degrees
            k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.radian ** 2
            )
            dihedral_restraint.addTorsion(
                restraint.index1[0],
                restraint.index2[0],
                restraint.index3[0],
                restraint.index4[0],
                [k, theta_0],
            )
            system.addForce(dihedral_restraint)
        else:
            # Probably needs openmm.CustomCentroidTorsionForce (?)
            raise NotImplementedError
        if ForceGroup:
            dihedral_restraint.setForceGroup(ForceGroup)
    return system
