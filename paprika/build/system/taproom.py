import json
import os
import shutil
from typing import Any, Dict, List

import mdtraj
import numpy
import openmm
import openmm.app as app
import parmed
from joblib import Parallel, delayed
from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from tqdm.auto import tqdm

from paprika.evaluator import Setup
from paprika.io import PaprikaEncoder, save_restraints
from paprika.restraints import DAT_restraint, create_window_list, parse_window
from paprika.restraints.openmm import apply_dat_restraint, apply_positional_restraints
from paprika.taproom import get_benchmarks, read_yaml_schema
from paprika.utils import index_from_mask


class BuildTaproomAPR:
    """A class to generate APR files from Taproom database in pAPRika.

    TODO: Implement an option to create files with implicit solvent.
    TODO: Implement an option to use GAFF force field.

    Parameters
    ----------
    host_code: str
        The 3-letter Taproom code for the host molecule.
    guest_code: str
        The 3-letter Taproom code for the guest molecule.
    n_water: int
        Number of Water molecules.
    build_folder: str
        Temporary folder to save intermediate files.
    working_folder: str
        The main folder to write the APR structure files.

    Example
    -------
    >>> from paprika.build.system import BuildTaproomAPR
    >>> from openff.toolkit import ForceField
    >>>
    >>> # Select OpenFF force field version
    >>> force_field = ForceField("openff-2.0.0.offxml")
    >>>
    >>> # Initiate system object
    >>> system = BuildTaproomAPR(host_code="bcd", guest_code="hex", n_water=3000, force_field=force_field)
    >>>
    >>> # Build APR files
    >>> system.build_system()
    >>>
    >>> # We can extend the default `r_final` specified in Taproom
    >>> from openff.units import unit
    >>> system.extend_pull_distance(extend_by=6*unit.angstrom)
    >>> system.build_system()
    >>>
    >>> # Creating these files can take 10-20 minutes on one core. We can speed things up with more CPU cores
    >>> system.build_system(n_cpus=4)
    """

    def __init__(
        self,
        host_code: str,
        guest_code: str,
        n_water: int,
        force_field: ForceField,
        build_folder: str = "build_files",
        working_folder: str = "simulations",
    ):
        self._host_code = host_code
        self._guest_code = guest_code
        self._n_water = n_water
        self._force_field = force_field
        self._build_folder = build_folder
        self._working_folder = working_folder
        self._host_metadata = None
        self._guest_metadata = None
        self._orientations = None
        self._water_mol = None
        self._water_intrcg = None
        self._restraints = None
        self._n_cpus = 1

        self._initialize()

    def _initialize(self):
        # Create folder
        os.makedirs(self._build_folder, exist_ok=True)
        os.makedirs(self._working_folder, exist_ok=True)

        # Load Host-Guest system from Taproom
        taproom = get_benchmarks()
        self._host_metadata = taproom["host_guest_systems"][self._host_code]
        self._guest_metadata = taproom["host_guest_systems"][self._host_code][
            self._guest_code
        ]
        self._orientations = list(self._host_metadata["yaml"].keys())
        self._host_yaml_schema = read_yaml_schema(
            self._host_metadata["yaml"][self._orientations[0]]
        )
        self._guest_yaml_schema = read_yaml_schema(self._guest_metadata["yaml"])

        # Water Molecule
        self._water_mol = Molecule.from_smiles("O")
        self._water_intrcg = Interchange.from_smirnoff(
            force_field=self._force_field,
            topology=[self._water_mol] * self._n_water,
        )

        # Load restraints
        self._restraints = {
            "static_restraints": self._unnest_restraint_specs(
                self._host_yaml_schema["restraints"]["static"]
            ),
            "conformational_restraints": self._unnest_restraint_specs(
                self._host_yaml_schema["restraints"]["conformational"]
            ),
            "guest_restraints": self._unnest_restraint_specs(
                self._guest_yaml_schema["restraints"]["guest"]
            ),
            "wall_restraints": self._unnest_restraint_specs(
                self._guest_yaml_schema["restraints"]["wall_restraints"]
            ),
            "symmetry_restraints": self._unnest_restraint_specs(
                self._guest_yaml_schema["symmetry_correction"]["restraints"]
            ),
        }

    @staticmethod
    def _unnest_restraint_specs(
        restraint_specs: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """A helper method to un-nest restraint lists parsed from a taproom
        yaml file.

        Parameters
        ----------
        restraint_specs
            The restraint specs to un-nest.
        """
        return [
            value["restraint"]
            for value in restraint_specs
            if value["restraint"] is not None
        ]

    @staticmethod
    def _restraints_to_dict(restraints: List[DAT_restraint]):
        """Converts a list of ``paprika`` restraint objects to
        a list of JSON compatible dictionary representations
        """

        return [
            json.loads(json.dumps(restraint.__dict__, cls=PaprikaEncoder))
            for restraint in restraints
        ]

    def _solvate_and_add_dummy(
        self,
        complex_path: str,
        solvated_path: str,
        unique_molecules: List[Molecule],
        offset: float,
    ):
        """Solvate a PDB file with PackMol through OpenFF-Interchange.

        Parameters
        ----------
        complex_path: str
            The file path of the complex PDB.
        solvated_path: str
            The output file path for the solvated complex.
        unique_molecules: List[Molecule]
            List of unique molecules (to generate the OpenFF Topology)
        offset: float
            An offset to place the dummy atoms.

        Returns
        -------
        system_intrcg: openff.interchange.Interchange
            The solvated system as OpenFF-Interchange object.
        """
        rectangular_box = numpy.asarray(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
            ]
        )

        # 01 - Solvate structure
        pdbfile = app.PDBFile(complex_path)
        solute_topology = Topology.from_openmm(
            pdbfile.topology, unique_molecules=unique_molecules
        )
        solute_intrcg = Interchange.from_smirnoff(
            force_field=self._force_field,
            topology=solute_topology,
            charge_from_molecules=unique_molecules,
        )
        solvated_topology = pack_box(
            molecules=[self._water_mol],
            number_of_copies=[self._n_water],
            solute=solute_intrcg.topology,
            box_shape=rectangular_box,
            mass_density=0.95 * unit.grams / unit.milliliters,
            center_solute="ORIGIN",
        )
        solute_intrcg.box = solvated_topology.box_vectors
        self._water_intrcg.box = solvated_topology.box_vectors
        solvated_topology.to_file(solvated_path)

        # 02 - Add Dummy Atoms to PDB
        input_structure = parmed.load_file(solvated_path, structure=True)
        Setup.add_dummy_atoms_to_structure(
            input_structure,
            dummy_atom_offsets=[
                numpy.array([0, 0, -offset]),
                numpy.array([0, 0, -3.0 - offset]),
                numpy.array([0, 2.2, -5.2 - offset]),
            ],
            offset_coordinates=numpy.zeros(3),
        )

        # 03 - Shift structure to avoid issues with PBC
        input_structure.coordinates += numpy.array(
            [
                input_structure.box[0] * 0.5,
                input_structure.box[1] * 0.5,
                -input_structure.coordinates[-1, 2] + 5.0,
            ]
        )

        # 04 - Write PDB for solvated system
        with open(solvated_path, "w") as f:
            app.PDBFile.writeFile(
                input_structure.topology,
                input_structure.positions,
                f,
                keepIds=True,
            )

        # 05 - Combine interchange objects
        system_intrcg = solute_intrcg + self._water_intrcg
        system_intrcg.box = solvated_topology.box_vectors

        return system_intrcg

    @staticmethod
    def _create_system_and_add_dummy(
        system_intrcg: Interchange, system_output_path: str
    ):
        """Convert `Interchange` object to OpenMM System and add dummy atoms.

        Parameters
        ----------
        system_intrcg: Interchange
            The Interchange object to convert.
        system_output_path: str
            The file path to save the OpenMM System to XML file.
        """
        openmm_system = system_intrcg.to_openmm()

        for _ in range(3):
            openmm_system.addParticle(mass=207)

        for i, force in enumerate(openmm_system.getForces()):
            if isinstance(force, openmm.NonbondedForce):
                force.addParticle(0.0, 1.0, 0.0)
                force.addParticle(0.0, 1.0, 0.0)
                force.addParticle(0.0, 1.0, 0.0)

        with open(system_output_path, "w") as f:
            f.write(openmm.XmlSerializer.serialize(openmm_system))

    def _build_pull_structures(
        self,
        i,
        orient,
        complex_path,
        guest_atom_indices,
        guest_orientation_mask,
        pull_distance,
        offset,
        n_windows,
        unique_molecules,
    ):
        """Function that translates guest molecules from host that is to be wrapped in `delayed` for parallelism."""

        folder = f"{self._working_folder}/pull-{orient}/p{i:03}"
        os.makedirs(folder, exist_ok=True)

        # 01 - Prepare complex structure
        structure = Setup.prepare_complex_structure(
            complex_path,
            guest_atom_indices,
            guest_orientation_mask,
            pull_distance=pull_distance,
            pull_window_index=i,
            n_pull_windows=n_windows["pull"],
        )
        complex_prepared_path = (
            f"{folder}/{self._host_code}-{self._guest_code}-{orient}.pdb"
        )
        with open(complex_prepared_path, "w") as f:
            app.PDBFile.writeFile(
                structure.topology,
                structure.positions,
                f,
                keepIds=True,
            )

        # 02 - Solvate structure
        complex_solvated_path = f"{folder}/restrained.pdb"
        host_guest_system_intrcg = self._solvate_and_add_dummy(
            complex_prepared_path,
            complex_solvated_path,
            unique_molecules=unique_molecules,
            offset=offset,
        )

        if i == 0:
            # 03 - Create Host-Guest OpenMM System with Dummy Atoms
            system_output_path = f"{self._build_folder}/{self._host_code}-{self._guest_code}-dum-solv.xml"
            self._create_system_and_add_dummy(
                host_guest_system_intrcg, system_output_path
            )

        # 04 - Clean up
        os.remove(complex_prepared_path)

    def _build_apr_structures(self):
        """Build and prepare the APR structures and windows."""
        host_resname = self._host_yaml_schema["resname"]
        n_windows = self._host_yaml_schema["calculation"]["windows"]

        # Create OpenFF Molecule instances of molecules
        guest_mol = Molecule.from_file(
            str(
                self._host_metadata["path"]
                .joinpath(self._guest_yaml_schema["name"])
                .joinpath(self._guest_yaml_schema["structure"]["sdf"])
            )
        )
        host_mol = Molecule.from_file(
            str(
                self._host_metadata["path"].joinpath(
                    self._host_yaml_schema["structure"]["sdf"]
                )
            )
        )

        # --------------------------------------------------------------------- #
        # Prepare Host-Guest Complex
        # --------------------------------------------------------------------- #
        print("Generating files for the `pull` phase.")
        for orient in self._orientations:
            # 01 - Load complex structure
            complex_path = str(
                self._host_metadata["path"]
                .joinpath(self._guest_yaml_schema["name"])
                .joinpath(self._guest_yaml_schema["complex"])
            ).replace(".pdb", f"-{orient}.pdb")
            structure = parmed.load_file(complex_path, structure=True)

            # 02 - Get Guest indices and mask
            guest_atom_indices = index_from_mask(
                structure, f":{self._guest_yaml_schema['name'].upper()}"
            )
            G1 = self._guest_yaml_schema["aliases"][3]["G1"]
            G2 = self._guest_yaml_schema["aliases"][4]["G2"]
            guest_orientation_mask = f"{G1} {G2}"

            # 03 - Initial `r` values
            r_initial = self._guest_yaml_schema["restraints"]["guest"][0]["restraint"][
                "attach"
            ]["target"]
            r_final = self._guest_yaml_schema["restraints"]["guest"][0]["restraint"][
                "pull"
            ]["target"]
            pull_distance = (r_final - r_initial).m_as(unit.angstrom)
            offset = self._guest_yaml_schema["restraints"]["guest"][0]["restraint"][
                "attach"
            ]["target"].m_as(unit.angstrom)

            # --------------------------------------------------------------------- #
            # Prepare `pull` windows
            # --------------------------------------------------------------------- #
            Parallel(n_jobs=self._n_cpus)(
                delayed(self._build_pull_structures)(
                    i,
                    orient,
                    complex_path,
                    guest_atom_indices,
                    guest_orientation_mask,
                    pull_distance,
                    offset,
                    n_windows,
                    [host_mol, guest_mol],
                )
                for i in tqdm(range(n_windows["pull"]))
            )

            # --------------------------------------------------------------------- #
            # Prepare `attach` windows - Copy PDB from p000
            # --------------------------------------------------------------------- #
            print("Generating files for the `attach` phase.")
            for i in tqdm(range(n_windows["attach"])):
                folder = f"{self._working_folder}/attach-{orient}/a{i:03}"
                os.makedirs(folder, exist_ok=True)
                shutil.copy(
                    f"{self._working_folder}/pull-{orient}/p000/restrained.pdb",
                    f"{folder}/restrained.pdb",
                )

        # --------------------------------------------------------------------- #
        # Prepare Host-only Structure
        # --------------------------------------------------------------------- #
        print("Generating files for the `release` phase.")
        # 01 - Remove Guest molecule from complex structure
        complex_solvate_path = f"{self._working_folder}/pull-p/p000/restrained.pdb"
        complex_structure = parmed.load_file(complex_solvate_path, structure=True)
        host_atom_indices = index_from_mask(complex_structure, mask=f":{host_resname}")

        mdtraj_trajectory = mdtraj.load_pdb(complex_solvate_path)
        host_trajectory = mdtraj_trajectory.atom_slice(host_atom_indices)
        host_trajectory.save(f"{self._build_folder}/host_input.pdb")

        # 02 - Align host molecule
        host_structure = Setup.prepare_host_structure(
            f"{self._build_folder}/host_input.pdb"
        )
        output_coordinate_path = f"{self._build_folder}/host_input_aligned.pdb"
        with open(output_coordinate_path, "w") as file:
            app.PDBFile.writeFile(
                host_structure.topology, host_structure.positions, file, True
            )

        # 03 - Solvate host molecule and add dummy atoms
        host_solvated_path = f"{self._build_folder}/{self._host_code}-dum-solv.pdb"
        host_system_intrcg = self._solvate_and_add_dummy(
            output_coordinate_path,
            host_solvated_path,
            unique_molecules=[host_mol],
            offset=offset,
        )

        # 04 - Create Host-only OpenMM System with Dummy Atoms
        system_output_path = f"{self._build_folder}/{self._host_code}-dum-solv.xml"
        self._create_system_and_add_dummy(host_system_intrcg, system_output_path)

        # --------------------------------------------------------------------- #
        # Prepare `release` windows - Copy PDB
        # --------------------------------------------------------------------- #
        for i in tqdm(range(n_windows["release"])):
            folder = f"{self._working_folder}/release/r{i:03}"
            os.makedirs(folder, exist_ok=True)
            shutil.copy(
                host_solvated_path,
                f"{folder}/restrained.pdb",
            )

    def _apply_attach_restraints(self):
        """Apply restraints for the `attach` phase."""

        print("Applying restraints for `attach` phase...")
        attach_lambdas = self._host_yaml_schema["calculation"]["lambda"]["attach"]
        n_windows = self._host_yaml_schema["calculation"]["windows"]

        for orient in self._orientations:
            attach_folder = f"{self._working_folder}/attach-{orient}"
            complex_path = f"{self._working_folder}/pull-{orient}/p000/restrained.pdb"

            static_restraints = Setup.build_static_restraints(
                complex_path,
                n_attach_windows=n_windows["attach"],
                n_pull_windows=None,
                n_release_windows=None,
                restraint_schemas=self._restraints["static_restraints"],
            )
            conformational_restraints = Setup.build_conformational_restraints(
                complex_path,
                attach_lambdas=attach_lambdas,
                n_pull_windows=None,
                release_lambdas=None,
                restraint_schemas=self._restraints["conformational_restraints"],
            )
            guest_restraints = Setup.build_guest_restraints(
                complex_path,
                attach_lambdas=attach_lambdas,
                n_pull_windows=None,
                restraint_schemas=self._restraints["guest_restraints"],
            )
            symmetry_restraints = Setup.build_symmetry_restraints(
                complex_path,
                n_attach_windows=n_windows["attach"],
                restraint_schemas=self._restraints["symmetry_restraints"],
            )
            wall_restraints = Setup.build_wall_restraints(
                complex_path,
                n_attach_windows=n_windows["attach"],
                restraint_schemas=self._restraints["wall_restraints"],
            )

            symmetry_restraints = (
                [] if symmetry_restraints is None else symmetry_restraints
            )
            wall_restraints = [] if wall_restraints is None else wall_restraints
            guest_restraints = [] if guest_restraints is None else guest_restraints

            restraints_dictionary = {
                "static": self._restraints_to_dict(static_restraints),
                "conformational": self._restraints_to_dict(conformational_restraints),
                "symmetry": self._restraints_to_dict(symmetry_restraints),
                "wall": self._restraints_to_dict(wall_restraints),
                "guest": self._restraints_to_dict(guest_restraints),
            }

            with open(f"{attach_folder}/restraints.json", "w") as file:
                json.dump(restraints_dictionary, file)

            save_restraints(
                conformational_restraints + guest_restraints,
                filepath=f"{attach_folder}/apr_restraints.json",
            )

            # Apply restraints
            attach_windows = create_window_list(guest_restraints)
            Parallel(n_jobs=self._n_cpus)(
                delayed(self._apply_attach_to_system)(
                    window,
                    complex_path,
                    attach_folder,
                    static_restraints,
                    conformational_restraints,
                    guest_restraints,
                    symmetry_restraints,
                    wall_restraints,
                )
                for window in tqdm(attach_windows)
            )

    def _apply_pull_restraints(self):
        """Apply restraints for the `pull` phase."""

        print("Applying restraints for `pull` phase...")
        attach_lambdas = self._host_yaml_schema["calculation"]["lambda"]["attach"]
        n_windows = self._host_yaml_schema["calculation"]["windows"]

        for orient in self._orientations:
            pull_folder = f"{self._working_folder}/pull-{orient}"
            complex_path = f"{self._working_folder}/pull-{orient}/p000/restrained.pdb"

            static_restraints = Setup.build_static_restraints(
                complex_path,
                n_attach_windows=n_windows["attach"],
                n_pull_windows=n_windows["pull"],
                n_release_windows=None,
                restraint_schemas=self._restraints["static_restraints"],
            )
            conformational_restraints = Setup.build_conformational_restraints(
                complex_path,
                attach_lambdas=attach_lambdas,
                n_pull_windows=n_windows["pull"],
                release_lambdas=None,
                restraint_schemas=self._restraints["conformational_restraints"],
            )
            guest_restraints = Setup.build_guest_restraints(
                complex_path,
                attach_lambdas=attach_lambdas,
                n_pull_windows=n_windows["pull"],
                restraint_schemas=self._restraints["guest_restraints"],
            )
            guest_restraints = [] if guest_restraints is None else guest_restraints

            # Remove the `attach` phases from the restraints as these restraints are
            # only being used for the pull phase.
            for restraint in (
                static_restraints + conformational_restraints + guest_restraints
            ):
                for key in restraint.phase["attach"]:
                    restraint.phase["attach"][key] = None

            restraints_dictionary = {
                "static": self._restraints_to_dict(static_restraints),
                "conformational": self._restraints_to_dict(conformational_restraints),
                "symmetry": None,
                "wall": None,
                "guest": self._restraints_to_dict(guest_restraints),
            }

            with open(f"{pull_folder}/restraints.json", "w") as file:
                json.dump(restraints_dictionary, file)

            save_restraints(
                conformational_restraints + guest_restraints,
                filepath=f"{pull_folder}/apr_restraints.json",
            )

            # Apply restraints
            pull_windows = create_window_list(guest_restraints)
            Parallel(n_jobs=self._n_cpus)(
                delayed(self._apply_pull_to_system)(
                    window,
                    complex_path,
                    pull_folder,
                    static_restraints,
                    conformational_restraints,
                    guest_restraints,
                )
                for window in tqdm(pull_windows)
            )

    def _apply_release_restraints(self):
        """Apply restraints for the `release` phase."""

        print("Applying restraints for `release` phase...")
        release_lambdas = self._host_yaml_schema["calculation"]["lambda"]["release"]
        n_windows = self._host_yaml_schema["calculation"]["windows"]

        release_folder = f"{self._working_folder}/release"
        host_solvated_path = f"{release_folder}/r000/restrained.pdb"

        static_restraints = Setup.build_static_restraints(
            host_solvated_path,
            n_attach_windows=None,
            n_pull_windows=None,
            n_release_windows=n_windows["release"],
            restraint_schemas=self._restraints["static_restraints"],
        )
        conformational_restraints = Setup.build_conformational_restraints(
            host_solvated_path,
            attach_lambdas=None,
            n_pull_windows=None,
            release_lambdas=release_lambdas,
            restraint_schemas=self._restraints["conformational_restraints"],
        )

        restraints_dictionary = {
            "static": self._restraints_to_dict(static_restraints),
            "conformational": self._restraints_to_dict(conformational_restraints),
            "symmetry": None,
            "wall": None,
            "guest": None,
        }

        with open(f"{release_folder}/restraints.json", "w") as file:
            json.dump(restraints_dictionary, file)

        save_restraints(
            conformational_restraints,
            filepath=f"{release_folder}/apr_restraints.json",
        )

        # Apply restraints
        system_path = f"{self._build_folder}/{self._host_code}-dum-solv.xml"
        release_windows = create_window_list(conformational_restraints)
        Parallel(n_jobs=self._n_cpus)(
            delayed(self._apply_release_to_system)(
                window,
                host_solvated_path,
                system_path,
                release_folder,
                static_restraints,
                conformational_restraints,
            )
            for window in tqdm(release_windows)
        )

    @staticmethod
    def _apply_attach_to_system(
        window,
        complex_path,
        attach_folder,
        static_restraints,
        conformational_restraints,
        guest_restraints,
        symmetry_restraints,
        wall_restraints,
    ):
        """Function that applies `attach` restraints that is to be wrapped in `delayed` for parallelism."""

        window_number, phase = parse_window(window)
        folder = f"{attach_folder}/{window}"
        os.makedirs(folder, exist_ok=True)

        system_path = f"{folder}/restrained.xml"
        with open(system_path, "r") as file:
            system = openmm.XmlSerializer.deserialize(file.read())

        for restraint in static_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=10)
        for restraint in conformational_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=11)
        for restraint in guest_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=12)
        for restraint in symmetry_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=13)
        for restraint in wall_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=14)

        apply_positional_restraints(complex_path, system, force_group=15)

        new_system_path = f"{folder}/restrained.xml"
        with open(new_system_path, "w") as file:
            file.write(openmm.XmlSerializer.serialize(system))

    @staticmethod
    def _apply_pull_to_system(
        window,
        complex_path,
        pull_folder,
        static_restraints,
        conformational_restraints,
        guest_restraints,
    ):
        """Function that applies `pull` restraints that is to be wrapped in `delayed` for parallelism."""

        window_number, phase = parse_window(window)
        folder = f"{pull_folder}/{window}"
        os.makedirs(folder, exist_ok=True)

        system_path = f"{folder}/restrained.xml"
        with open(system_path, "r") as file:
            system = openmm.XmlSerializer.deserialize(file.read())

        for restraint in static_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=10)
        for restraint in conformational_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=11)
        for restraint in guest_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=12)

        apply_positional_restraints(complex_path, system, force_group=15)

        new_system_path = f"{folder}/restrained.xml"
        with open(new_system_path, "w") as file:
            file.write(openmm.XmlSerializer.serialize(system))

    @staticmethod
    def _apply_release_to_system(
        window,
        host_solvated_path,
        system_path,
        release_folder,
        static_restraints,
        conformational_restraints,
    ):
        """Function that applies `release` restraints that is to be wrapped in `delayed` for parallelism."""

        window_number, phase = parse_window(window)
        folder = f"{release_folder}/{window}"
        os.makedirs(folder, exist_ok=True)

        with open(system_path, "r") as file:
            system = openmm.XmlSerializer.deserialize(file.read())

        for restraint in static_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=10)
        for restraint in conformational_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=11)

        apply_positional_restraints(host_solvated_path, system, force_group=15)

        new_system_path = f"{folder}/restrained.xml"
        with open(new_system_path, "w") as file:
            file.write(openmm.XmlSerializer.serialize(system))

    def extend_pull_distance(self, extend_by: unit.Quantity):
        """Extend the pull distance further than what is in the original Taproom Metadata.

        Parameters
        ----------
        extend_by: openff.unit.units.Quantity
            The extended distance to pull the guest molecule to.
        """

        # Determine dr between windows
        pull_distance = (
            self._restraints["guest_restraints"][0]["pull"]["target"]
            - self._restraints["guest_restraints"][0]["attach"]["target"]
        ).m_as(unit.angstrom)
        n_pull_windows = self._host_yaml_schema["calculation"]["windows"]["pull"]
        dr = pull_distance / n_pull_windows

        # Update distance
        self._restraints["guest_restraints"][0]["pull"]["target"] = extend_by
        new_pull_distance = (
            self._restraints["guest_restraints"][0]["pull"]["target"]
            - self._restraints["guest_restraints"][0]["attach"]["target"]
        ).m_as(unit.angstrom)
        updated_n_windows = int(new_pull_distance / dr)

        self._host_yaml_schema["calculation"]["windows"]["pull"] = updated_n_windows

    def build_system(self, n_cpus: int = 1, clean_files: bool = False):
        """Build the APR files in the order:

        (1) Generate Structures and add dummy atoms
        (2) Apply `attach` restraints
        (3) Apply `pull` restraints
        (4) Apply `release` restraints

        Parameters
        ----------
        n_cpus: int
            Number of CPUs to spread the workload.
        clean_files: bool
            Option to delete temporary files.
        """

        self._n_cpus = n_cpus

        self._build_apr_structures()
        self._apply_attach_restraints()
        self._apply_pull_restraints()
        self._apply_release_restraints()

        if clean_files:
            shutil.rmtree(self._build_folder)
