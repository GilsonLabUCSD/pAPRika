"""
This class contains a simulation setup wrapper for use with the OpenFF Evaluator.
"""

import logging
import os
import shutil
import subprocess as sp
from pathlib import Path

import numpy as np
import parmed as pmd
import pkg_resources
import pytraj as pt
import simtk.openmm as openmm
import simtk.unit as unit

from paprika import align
from paprika.restraints import static_DAT_restraint, DAT_restraint
from paprika.restraints.read_yaml import read_yaml
from paprika.restraints.restraints import create_window_list

logger = logging.getLogger(__name__)
_PI_ = np.pi

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

    def __init__(self, host, guest=None,
                 backend="openmm", directory_path="benchmarks",
                 additional_benchmarks=None, generate_gaff_files=False, gaff_version="gaff2",
                 guest_orientation=None, build=True):
        self.host = host
        self.guest = guest if guest is not None else "release"
        self.backend = backend
        self.directory = Path(directory_path).joinpath(self.host).joinpath(f"{self.guest}-{guest_orientation}" if
                                                                           guest_orientation is not None else
                                                                           f"{self.guest}")
        self.desolvated_window_paths = []
        self.window_list = []

        if self.backend == "amber":
            # Generate `frcmod` and dummy atom files.
            raise NotImplementedError

        self.directory.mkdir(parents=True, exist_ok=True)
        installed_benchmarks = get_benchmarks()
        if additional_benchmarks is not None:
            installed_benchmarks.update(additional_benchmarks)

        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks, guest_orientation)

        self.benchmark_path = host_yaml.parent
        self.host_yaml = read_yaml(host_yaml)
        if guest:
            self.guest_yaml = read_yaml(guest_yaml["yaml"])

        if build:
            # Here, we build desolvated windows and pass the files to the OpenFF Evaluator.
            # These files are stored in `self.desolvated_window_paths`.
            self.build_desolvated_windows(guest_orientation)
            if generate_gaff_files:
                generate_gaff(mol2_file=self.benchmark_path.joinpath(self.host_yaml["structure"]),
                              residue_name=self.host_yaml["resname"],
                              output_name=self.host,
                              directory_path=self.directory,
                              gaff=gaff_version)
                if guest:
                    generate_gaff(mol2_file=self.benchmark_path.joinpath(
                        self.guest).joinpath(self.guest_yaml["structure"]),
                                  output_name=self.guest,
                                  residue_name=self.guest_yaml["name"],
                                  directory_path=self.directory,
                                  gaff=gaff_version)
        if not build:
            self.populate_window_list(input_pdb=os.path.join(self.directory, f"{self.host}-{self.guest}.pdb" if self.guest is not None
            else f"{self.host}.pdb"))

    def parse_yaml(self, installed_benchmarks, guest_orientation):
        """
        Read the YAML recipe for the host and guest.

        Returns
        -------

        """
        try:
            if guest_orientation:

                host_yaml = installed_benchmarks["host_guest_systems"][self.host]["yaml"][guest_orientation]
            else:
                host_yaml = installed_benchmarks["host_guest_systems"][self.host]["yaml"]["p"]

        except KeyError:
            logger.error(f"Cannot find YAML recipe for host: {self.host}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        try:
            guest_yaml = installed_benchmarks["host_guest_systems"][self.host][self.guest]
        except KeyError:
            if self.guest == "release":
                guest_yaml = None
            else:
                logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
                logger.debug(installed_benchmarks)
                raise FileNotFoundError

        return host_yaml, guest_yaml

    def align(self, input_pdb):
        structure = pmd.load_file(str(input_pdb), structure=True)
        intermediate_pdb = self.directory.joinpath(f"tmp.pdb")
        destination_pdb = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

        if not self.guest == "release":
            # Align the host-guest complex so the first guest atom is at (0, 0, 0) and the second guest atom lies
            # along the positive z-axis.
            guest_angle_restraint_mask = self.guest_yaml["restraints"]["guest"][-1]["restraint"][
                "atoms"
            ].split()
            aligned_structure = align.zalign(
                structure, guest_angle_restraint_mask[1], guest_angle_restraint_mask[2]
            )
            aligned_structure.save(str(intermediate_pdb), overwrite=True)

        else:
            # Create a PDB file just for the host.

            host = pmd.load_file(str(input_pdb), structure=True)
            host_coordinates = host[f":{self.host_yaml['resname'].upper()}"].coordinates
            # Cheap way to get the center of geometry
            offset_coordinates = pmd.geometry.center_of_mass(host_coordinates,
                                                             masses=np.ones(len(host_coordinates)))

            # Find the principal components, take the two largest, and find the vector orthogonal to that
            # (should be cross-product right hand rule, I think). Use that vector to align with the z-axis.
            # This may not generalize to non-radially-symmetric host molecules.

            aligned_coords = np.empty_like(structure.coordinates)
            for atom in range(len(structure.atoms)):
                aligned_coords[atom] = structure.coordinates[atom] - offset_coordinates
            structure.coordinates = aligned_coords

            inertia_tensor = np.dot(structure.coordinates.transpose(), structure.coordinates)
            eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
            order = np.argsort(eigenvalues)
            axis_3, axis_2, axis_1 = eigenvectors[:, order].transpose()

            dummy_axis = np.cross(axis_1, axis_2)

            self._add_dummy_to_PDB(input_pdb=input_pdb,
                                   output_pdb=intermediate_pdb,
                                   offset_coordinates=offset_coordinates,
                                   dummy_atom_tuples=[(0, 0, 0),
                                                      (dummy_axis[0], dummy_axis[1], dummy_axis[2])])
            structure = pmd.load_file(str(intermediate_pdb), structure=True)

            for atom in structure.atoms:
                atom.mass = 1.0

            aligned_structure = align.zalign(
                structure, ":DM1", ":DM2"
            )
            aligned_structure["!:DM1&!:DM2"].save(str(intermediate_pdb),
                                                  overwrite=True)


        # Save aligned PDB file with CONECT records.
        positions_pdb = openmm.app.PDBFile(str(intermediate_pdb))
        topology_pdb = openmm.app.PDBFile(str(input_pdb))

        positions = positions_pdb.positions
        topology = topology_pdb.topology
        with open(destination_pdb, "w") as file:
            openmm.app.PDBFile.writeFile(topology, positions, file)
        os.remove(intermediate_pdb)

    def populate_window_list(self, input_pdb):
        logger.debug("Setting up dummy restraint to build window list.")
        _dummy_restraint = self._create_dummy_restraint(
            initial_structure=str(input_pdb),
        )
        self.window_list = create_window_list([_dummy_restraint])
        return _dummy_restraint

    def build_desolvated_windows(self, guest_orientation):
        if self.guest != "release":
            if not guest_orientation:
                initial_structure = self.benchmark_path.joinpath(self.guest).joinpath(
                    self.guest_yaml["complex"]
                )
            else:
                base_name = Path(self.guest_yaml["complex"]).stem
                orientation_structure = base_name + f"-{guest_orientation}.pdb"
                initial_structure = self.benchmark_path.joinpath(self.guest).joinpath(
                    orientation_structure
                )

        else:
            initial_structure = self.directory.joinpath(self.benchmark_path.joinpath(self.host_yaml["structure"]))
            host = pt.iterload(str(initial_structure), str(initial_structure))
            host.save(str(self.directory.joinpath(f"{self.host}.pdb")), overwrite=True, options='conect')
            initial_structure = str(self.directory.joinpath(f"{self.host}.pdb"))

        self.align(input_pdb=initial_structure)

        _dummy_restraint = self.populate_window_list(input_pdb=initial_structure)

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

        if self.guest != "release":
            windows = [
                self.host_yaml["calculation"]["windows"]["attach"],
                self.host_yaml["calculation"]["windows"]["pull"],
                None,
            ]
        else:
            windows = [
                None,
                None,
                self.host_yaml["calculation"]["windows"]["release"]
            ]


        guest_restraint = DAT_restraint()
        guest_restraint.auto_apr = True
        guest_restraint.continuous_apr = True
        guest_restraint.amber_index = False if self.backend == "openmm" else True
        guest_restraint.topology = str(initial_structure)
        guest_restraint.mask1 = "@1"
        guest_restraint.mask2 = "@2"

        if self.guest != "release":
            restraint = self.guest_yaml["restraints"]["guest"][0]
            guest_restraint.attach["target"] = restraint["restraint"]["attach"][
                "target"
            ]
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
        else:
            # Remember, the purpose of this *fake* restraint is *only* to figure out how many windows to make,
            # so we can use the OpenFF Evaluator to solvate the structures for us. To figure out how many winodws
            # we need, just setting the lambda values should be sufficient.
            guest_restraint.auto_apr = False
            guest_restraint.continuous_apr = False

            guest_restraint.release["target"] = 1.0
            guest_restraint.release["fc_final"] = 1.0
            guest_restraint.release["fraction_list"] = self.host_yaml["calculation"][
                "lambda"
            ]["release"]

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
            try:
                # Copy the final pull window, if it exists
                source_pdb = (
                    self.directory.joinpath("windows")
                    .joinpath(f"p{self.host_yaml['calculation']['windows']['pull']:03d}")
                    .joinpath(f"{self.host}-{self.guest}.pdb")
                )
                shutil.copy(source_pdb, window_path)
            except FileNotFoundError:
                # Copy the initial structure, assuming we are doing a standalone release calculation.
                shutil.copy(self.directory.joinpath(f"{self.host}-{self.guest}.pdb"),
                window_path)

    def _add_dummy_to_PDB(self, input_pdb, output_pdb, offset_coordinates,
                          dummy_atom_tuples):
        input_pdb_file = openmm.app.PDBFile(input_pdb)

        positions = input_pdb_file.positions

        # When we pass in a guest, we have multiple coordinates and the function expects to address the first guest
        # atom coordinates.
        # When we pass in the center of mass of the host, we'll only have one set of coordinates.
        if len(np.shape(offset_coordinates)) < 2:
            offset_coordinates = [offset_coordinates, ]

        for index, dummy_atom_tuple in enumerate(dummy_atom_tuples):

            positions.append(
                openmm.Vec3(
                    offset_coordinates[0][0] + dummy_atom_tuple[0],
                    offset_coordinates[0][1] + dummy_atom_tuple[1],
                    offset_coordinates[0][2] + dummy_atom_tuple[2],
                )
                * unit.angstrom
            )

        topology = input_pdb_file.topology
        for dummy_index in range(len(dummy_atom_tuples)):
            dummy_chain = topology.addChain(None)
            dummy_residue = topology.addResidue(f"DM{dummy_index + 1}", dummy_chain)
            topology.addAtom(f"DUM", None, dummy_residue)

        with open(output_pdb, "w") as file:
            openmm.app.PDBFile.writeFile(topology, positions, file)

    def _add_dummy_to_System(self, system, dummy_atom_tuples):
        [system.addParticle(mass=207) for _ in range(len(dummy_atom_tuples))]

        for force_index in range(system.getNumForces()):
            force = system.getForce(force_index)
            if not isinstance(force, openmm.NonbondedForce):
                continue
            force.addParticle(0.0, 1.0, 0.0)
            force.addParticle(0.0, 1.0, 0.0)
            force.addParticle(0.0, 1.0, 0.0)

        return system

    def add_dummy_atoms(
        self,
        reference_pdb="reference.pdb",
        solvated_pdb="output.pdb",
        solvated_xml="system.xml",
        dummy_pdb="output.pdb",
        dummy_xml="output.xml",
    ):

        reference_structure = pmd.load_file(reference_pdb, structure=True)

        # Determine the offset coordinates for the new dummy atoms.
        if self.guest == "release":

            host_coordinates = reference_structure[f":{self.host_yaml['resname'].upper()}"].coordinates
            # Cheap way to get the center of geometry
            offset_coordinates = pmd.geometry.center_of_mass(host_coordinates,
                                                             masses=np.ones(len(host_coordinates)))

        else:
            guest_angle_restraint_mask = self.guest_yaml["restraints"]["guest"][-1]["restraint"][
                "atoms"
            ].split()

            offset_coordinates = reference_structure[f':{self.guest_yaml["name"].upper()} | :{self.host_yaml["resname"].upper()}']\
                [guest_angle_restraint_mask[1]].coordinates

        # First add dummy atoms to structure
        logger.debug(f"Adding dummy atoms to {solvated_pdb}")
        try:

            self._add_dummy_to_PDB(solvated_pdb, dummy_pdb, offset_coordinates,
                                   dummy_atom_tuples=[(0, 0, -6.0),
                                                      (0, 0, -9.0),
                                                      (0, 2.2, -11.2)])
        except FileNotFoundError:
            logger.warning(f"Missing {solvated_pdb}")


        self._wrap(dummy_pdb)

        # Add dummy atoms to System
        if solvated_xml is not None:

            try:

                system = read_openmm_system_from_xml(solvated_xml)
                system = self._add_dummy_to_System(system, dummy_atom_tuples=[(0, 0, -6.0),
                                                                              (0, 0, -9.0),
                                                                              (0, 2.2, -11.2)])
                system_xml = openmm.XmlSerializer.serialize(system)
                with open(dummy_xml, "w") as file:
                    file.write(system_xml)

            except FileNotFoundError:
                logger.warning(f"Missing {solvated_xml}")

    @staticmethod
    def _wrap(file, mask=":DM3"):
        logging.info(f"Re-wrapping {file} to avoid pulling near periodic boundaries.")
        structure = pmd.load_file(file, structure=True)

        anchor = structure[mask]
        anchor_z = anchor.atoms[0].xz

        for atom in structure.atoms:
            atom.xz -= anchor_z - 2.0

        structure.save(file, overwrite=True)


    def initialize_restraints(self, structure="output.pdb"):

        if self.guest != "release":

            windows = [
                self.host_yaml["calculation"]["windows"]["attach"],
                self.host_yaml["calculation"]["windows"]["pull"],
                None,
            ]
        else:
            windows = [None,
                       None,
                       self.host_yaml["calculation"]["windows"]["release"]
            ]

        static_restraints = []
        for restraint in self.host_yaml["restraints"]["static"]:

            static = static_DAT_restraint(
                restraint_mask_list=restraint["restraint"]["atoms"].split(),
                num_window_list=windows,
                ref_structure=str(structure),
                force_constant=restraint["restraint"]["force_constant"],
                amber_index=False if self.backend == "openmm" else True,
            )
            static_restraints.append(static)

        conformational_restraints = []
        if self.host_yaml["restraints"]["conformational"]:

            for conformational in self.host_yaml["restraints"][
                "conformational"
            ]:

                mask = conformational["restraint"]["atoms"].split()

                conformational_restraint = DAT_restraint()
                conformational_restraint.auto_apr = True
                conformational_restraint.continuous_apr = True
                conformational_restraint.amber_index = False if self.backend == "openmm" else True
                conformational_restraint.topology = str(structure)
                conformational_restraint.mask1 = mask[0]
                conformational_restraint.mask2 = mask[1]
                conformational_restraint.mask3 = mask[2] if len(mask) > 2 else None
                conformational_restraint.mask4 = mask[3] if len(mask) > 3 else None

                if self.guest != "release":

                    conformational_restraint.attach["target"] = conformational["restraint"][
                        "target"
                    ]
                    conformational_restraint.attach["fc_final"] = conformational["restraint"][
                        "force_constant"
                    ]
                    conformational_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
                        "lambda"
                    ]["attach"]

                    conformational_restraint.pull["target_final"] = conformational["restraint"][
                        "target"
                    ]
                    conformational_restraint.pull["num_windows"] = windows[1]

                else:

                    conformational_restraint.auto_apr = False
                    conformational_restraint.continuous_apr = False

                    conformational_restraint.release["target"] = conformational["restraint"][
                        "target"
                    ]
                    conformational_restraint.release["fc_final"] = conformational["restraint"][
                        "force_constant"
                    ]
                    conformational_restraint.release["fraction_list"] = self.host_yaml["calculation"][
                        "lambda"
                    ]["release"]


                conformational_restraint.initialize()
                conformational_restraints.append(conformational_restraint)
        else:
            logger.debug("Skipping conformational restraints...")

        symmetry_restraints = []
        if self.guest != "release" and "symmetry_correction" in self.guest_yaml:
            for symmetry in self.guest_yaml["symmetry_correction"]["restraints"]:
                symmetry_restraint = DAT_restraint()
                symmetry_restraint.auto_apr = True
                symmetry_restraint.continuous_apr = True
                symmetry_restraint.amber_index = False if self.backend == "openmm" else True
                symmetry_restraint.topology = str(structure)
                symmetry_restraint.mask1 = symmetry["atoms"].split()[0]
                symmetry_restraint.mask2 = symmetry["atoms"].split()[1]
                symmetry_restraint.mask3 = symmetry["atoms"].split()[2]

                symmetry_restraint.attach["fc_final"] = symmetry["force_constant"]
                symmetry_restraint.attach["fraction_list"] = [1.0] * len(self.host_yaml["calculation"][
                        "lambda"
                    ]["attach"])
                # This target should be overridden by the custom values.
                symmetry_restraint.attach["target"] = 999.99
                symmetry_restraint.custom_restraint_values["r2"] = 91
                symmetry_restraint.custom_restraint_values["r3"] = 91
                # 0 force constant between 91 degrees and 180 degrees.
                symmetry_restraint.custom_restraint_values["rk3"] = 0.0
                symmetry_restraint.initialize()

                symmetry_restraints.append(symmetry_restraint)

        else:
            logger.debug("Skipping symmetry restraints...")

        wall_restraints = []
        if self.guest != "release" and "wall_restraints" in self.guest_yaml['restraints']:
            for wall in self.guest_yaml["restraints"]["wall_restraints"]:
                wall_restraint = DAT_restraint()
                wall_restraint.auto_apr = True
                wall_restraint.continuous_apr = True
                wall_restraint.amber_index = False if self.backend == "openmm" else True
                wall_restraint.topology = str(structure)
                wall_restraint.mask1 = wall["restraint"]["atoms"].split()[0]
                wall_restraint.mask2 = wall["restraint"]["atoms"].split()[1]

                wall_restraint.attach["fc_final"] = wall["restraint"]["force_constant"]
                wall_restraint.attach["fraction_list"] = [1.0] * len(self.host_yaml["calculation"][
                                                                             "lambda"
                                                                         ]["attach"])
                wall_restraint.attach["target"] = wall["restraint"]["target"]
                # Minimum distance is 0 Angstrom
                wall_restraint.custom_restraint_values["r1"] = 0
                wall_restraint.custom_restraint_values["r2"] = 0
                # Harmonic force constant beyond target distance.
                wall_restraint.custom_restraint_values["rk2"] = wall["restraint"]["force_constant"]
                wall_restraint.custom_restraint_values["rk3"] = wall["restraint"]["force_constant"]
                wall_restraint.initialize()

                wall_restraints.append(wall_restraint)

        else:
            logger.debug("Skipping wall restraints...")

        guest_restraints = []
        for restraint in [] if not hasattr(self, 'guest_yaml') else self.guest_yaml["restraints"]["guest"]:
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

            guest_restraint.pull["target_final"] = restraint["restraint"]["pull"][
                "target"
            ]
            guest_restraint.pull["num_windows"] = windows[1]

            guest_restraint.initialize()
            guest_restraints.append(guest_restraint)

        return (
            static_restraints,
            conformational_restraints,
            symmetry_restraints,
            wall_restraints,
            guest_restraints,
        )

    def initialize_calculation(self, window, structure_path="output.pdb",
                               input_xml="system.xml", output_xml="system.xml"):
        if self.backend == "amber":
            # Write simulation input files in each directory
            raise NotImplementedError
        try:
            system = read_openmm_system_from_xml(input_xml)
        except FileNotFoundError:
            logger.warning(f"Cannot read XML from {input_xml}")

        # Apply the positional restraints.
        structure = pmd.load_file(structure_path, structure=True)

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
                positional_restraint.setForceGroup(15)

        for restraint in self.static_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=10)
        for restraint in self.conformational_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=11)
        for restraint in self.guest_restraints:
            system = apply_openmm_restraints(system, restraint, window, ForceGroup=12)
        for restraint in self.symmetry_restraints:
            system = apply_openmm_restraints(system, restraint, window, flat_bottom=True, ForceGroup=13)
        for restraint in self.wall_restraints:
            system = apply_openmm_restraints(system, restraint, window, flat_bottom=True, ForceGroup=14)

        system_xml = openmm.XmlSerializer.serialize(system)

        with open(output_xml, "w") as file:
            file.write(system_xml)

def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks

def apply_openmm_restraints(system, restraint, window, flat_bottom=False, ForceGroup=None):
    if window[0] == "a":
        phase = "attach"
    elif window[0] == "p":
        phase = "pull"
    elif window[0] == "r":
        phase = "release"
    window_number = int(window[1:])

    if flat_bottom and phase == "attach" and restraint.mask3:
        flat_bottom_force = openmm.CustomAngleForce('step(-(theta - theta_0)) * k * (theta - theta_0)^2')
        # If theta is greater than theta_0, then the argument to step is negative, which means the force is off.
        flat_bottom_force.addPerAngleParameter("k")
        flat_bottom_force.addPerAngleParameter("theta_0")

        theta_0 = 91.0 * unit.degrees
        k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.radian ** 2
        )
        flat_bottom_force.addAngle(
            restraint.index1[0],
            restraint.index2[0],
            restraint.index3[0],
            [k, theta_0],
        )
        system.addForce(flat_bottom_force)
        if ForceGroup:
            flat_bottom_force.setForceGroup(ForceGroup)

        return system
    elif flat_bottom and phase == "attach" and not restraint.mask3:
        flat_bottom_force = openmm.CustomBondForce('step((r - r_0)) * k * (r - r_0)^2')
        # If x is greater than x_0, then the argument to step is positive, which means the force is on.
        flat_bottom_force.addPerBondParameter("k")
        flat_bottom_force.addPerBondParameter("r_0")

        r_0 = restraint.phase[phase]["targets"][window_number] * unit.angstrom
        k = (
                restraint.phase[phase]["force_constants"][window_number]
                * unit.kilocalories_per_mole
                / unit.radian ** 2
        )
        flat_bottom_force.addBond(
            restraint.index1[0],
            restraint.index2[0],
            [k, r_0],
        )
        system.addForce(flat_bottom_force)
        if ForceGroup:
            flat_bottom_force.setForceGroup(ForceGroup)

        return system

    elif flat_bottom and phase == "pull":
        return system
    elif flat_bottom and phase == "release":
        return system

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
            dihedral_restraint = openmm.CustomTorsionForce(f"k * min(min(abs(theta - theta_0), abs(theta - theta_0 + 2 * {_PI_})), abs(theta - theta_0 - 2 * {_PI_}))^2")
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

def generate_gaff(mol2_file, residue_name, output_name=None, need_gaff_atom_types=True, generate_frcmod=True,
                  directory_path="benchmarks", gaff="gaff2"):

    if output_name is None:
        output_name = mol2_file.stem

    if need_gaff_atom_types:
        _generate_gaff_atom_types(mol2_file=mol2_file,
                                  residue_name=residue_name,
                                  output_name=output_name,
                                  gaff=gaff,
                                  directory_path=directory_path)
        logging.debug("Checking to see if we have a multi-residue MOL2 file that should be converted "
                      "to single-residue...")
        structure = pmd.load_file(os.path.join(directory_path, f"{output_name}.{gaff}.mol2"), structure=True)
        if len(structure.residues) > 1:
            structure[":1"].save("tmp.mol2")
            if os.path.exists("tmp.mol2"):
                os.rename("tmp.mol2", os.path.join(directory_path, f"{output_name}.{gaff}.mol2"))
                logging.debug("Saved single-residue MOL2 file for `tleap`.")
            else:
                raise RuntimeError("Unable to convert multi-residue MOL2 file to single-residue for `tleap`.")

        if generate_frcmod:
            _generate_frcmod(mol2_file=f'{output_name}.{gaff}.mol2',
                             gaff=gaff,
                             output_name=output_name,
                             directory_path=directory_path)
        else:
            raise NotImplementedError()

def _generate_gaff_atom_types(mol2_file, residue_name, output_name, gaff="gaff2", directory_path="benchmarks"):
    
    p = sp.Popen(["antechamber", "-i", str(mol2_file), "-fi", "mol2",
              "-o", f"{output_name}.{gaff}.mol2", "-fo", "mol2",
              "-rn", f"{residue_name.upper()}",
              "-at", f"{gaff}",
              "-an", "no",
              "-dr", "no",
              "-pf", "yes"], cwd=directory_path)
    p.communicate()

    files = ["ANTECHAMBER_AC.AC", "ANTECHAMBER_AC.AC0",
             "ANTECHAMBER_BOND_TYPE.AC", "ANTECHAMBER_BOND_TYPE.AC0",
             "ATOMTYPE.INF"]
    files = [directory_path.joinpath(i) for i in files]
    for file in files:
        if file.exists():
            logger.debug(f"Removing temporary file: {file}")
            file.unlink()
            
    if not os.path.exists(f"{output_name}.{gaff}.mol2"):
        # Try with the newer (AmberTools 19) version of `antechamber` which doesn't have the `-dr` flag
        p = sp.Popen(["antechamber", "-i", str(mol2_file), "-fi", "mol2",
                      "-o", f"{output_name}.{gaff}.mol2", "-fo", "mol2",
                      "-rn", f"{residue_name.upper()}",
                      "-at", f"{gaff}",
                      "-an", "no",
                      "-pf", "yes"], cwd=directory_path)
        p.communicate()

        files = ["ANTECHAMBER_AC.AC", "ANTECHAMBER_AC.AC0",
                 "ANTECHAMBER_BOND_TYPE.AC", "ANTECHAMBER_BOND_TYPE.AC0",
                 "ATOMTYPE.INF"]
        files = [directory_path.joinpath(i) for i in files]
        for file in files:
            if file.exists():
                logger.debug(f"Removing temporary file: {file}")
                file.unlink()


def _generate_frcmod(mol2_file, gaff, output_name, directory_path="benchmarks"):
    sp.Popen(["parmchk2", "-i", str(mol2_file), "-f", "mol2",
              "-o", f"{output_name}.{gaff}.frcmod",
              "-s", f"{gaff}"
              ], cwd=directory_path)
