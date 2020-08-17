"""
This class contains a simulation setup wrapper for use with the OpenFF Evaluator.
"""

import logging
import os
import subprocess
from typing import Any, Dict, List, Optional

import numpy as np
import parmed as pmd
import pkg_resources

from paprika import align
from paprika.restraints import DAT_restraint, static_DAT_restraint

logger = logging.getLogger(__name__)
_PI_ = np.pi


def _get_installed_benchmarks():
    _installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        _installed_benchmarks[entry_point.name] = entry_point.load()
    return _installed_benchmarks


class Setup(object):
    """
    The Setup class provides a wrapper function around the preparation of the host-guest
    system and the application of restraints.
    """

    @classmethod
    def prepare_host_structure(
        cls, coordinate_path: str, host_atom_indices: Optional[List[int]] = None
    ) -> pmd.Structure:
        """Attempts to align the cavity of the host along the z-axis ready
        for the APR calculation.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which contains the host molecule.
        host_atom_indices
            The indices of the host atoms in the coordinate file if the file. This may
            be used if the coordinate file contains more than a single molecule.
        """
        import parmed.geometry

        # noinspection PyTypeChecker
        structure = pmd.load_file(coordinate_path, structure=True)

        # Extract the host from the full structure.
        if not host_atom_indices:
            host_atom_indices = range(len(structure.atoms))

        host_structure = structure[
            "@" + ",".join(map(lambda x: str(x + 1), host_atom_indices))
        ]

        # noinspection PyTypeChecker
        center_of_mass: np.ndarray = pmd.geometry.center_of_mass(
            host_structure.coordinates, masses=np.ones(len(host_structure.coordinates))
        )

        # Remove the COM from the host coordinates to make alignment easier.
        structure.coordinates -= center_of_mass
        host_structure.coordinates -= center_of_mass

        # Find the principal components of the host, take the two largest, and find
        # the vector orthogonal to that. Use that vector to align with the z-axis.
        # This may not generalize to non-radially-symmetric host molecules.
        inertia_tensor = np.dot(
            host_structure.coordinates.transpose(), host_structure.coordinates
        )

        eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
        order = np.argsort(eigenvalues)

        _, axis_2, axis_1 = eigenvectors[:, order].transpose()

        cavity_axis = np.cross(axis_1, axis_2)

        # Add dummy atoms which will be used to align the structure.
        cls.add_dummy_atoms_to_structure(structure, [np.array([0, 0, 0]), cavity_axis])

        # Give atoms uniform mass so that the align code uses the center
        # of geometry rather than the center of mass.
        for atom in structure.atoms:
            atom.mass = 1.0

        aligned_structure = align.zalign(structure, ":DM1", ":DM2")

        # Return a copy of the original structure where things like the masses
        # have not been changed and dummy atoms not added.

        # noinspection PyTypeChecker
        structure: pmd.Structure = pmd.load_file(coordinate_path, structure=True)
        structure.coordinates = aligned_structure["!:DM1&!:DM2"].coordinates

        return structure

    @classmethod
    def prepare_complex_structure(
        cls,
        coordinate_path: str,
        guest_atom_indices: List[int],
        guest_orientation_mask: str,
        pull_distance: float,
        pull_window_index: int,
        n_pull_windows: int,
    ) -> pmd.Structure:
        """Prepares the coordinates of both the host and the guest molecule in a
        given APR window.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which contains the guest molecule
            bound to the host.
        guest_atom_indices
        guest_orientation_mask
        pull_distance
        pull_window_index
        n_pull_windows
        """

        # Align the host-guest complex so the first guest atom is at (0, 0, 0) and the
        # second guest atom lies along the positive z-axis.
        # noinspection PyTypeChecker
        structure: pmd.Structure = pmd.load_file(coordinate_path, structure=True)

        (
            guest_orientation_mask_0,
            guest_orientation_mask_1,
        ) = guest_orientation_mask.split(" ")

        aligned_structure = align.zalign(
            structure, guest_orientation_mask_0, guest_orientation_mask_1
        )

        # n_pull_windows = self.host_yaml["calculation"]["windows"]["pull"]

        # target_initial = self.guest_yaml["restraints"]["guest"][0]["attach"]["target"]
        # target_final = self.host_yaml["calculation"]["target"]["pull"]

        target_distance = np.linspace(0.0, pull_distance, n_pull_windows)[
            pull_window_index
        ]
        target_difference = target_distance

        for guest_index in guest_atom_indices:
            aligned_structure.atoms[guest_index].xz += target_difference

        return aligned_structure

    @staticmethod
    def add_dummy_atoms_to_structure(
        structure: pmd.Structure,
        dummy_atom_offsets: List[np.ndarray],
        offset_coordinates: Optional[np.ndarray] = None,
    ):

        if offset_coordinates is None:
            offset_coordinates = np.zeros(3)

        full_coordinates = np.vstack(
            [
                structure.coordinates,
                *[
                    offset_coordinates + dummy_atom_offset
                    for dummy_atom_offset in dummy_atom_offsets
                ],
            ]
        )

        for index in range(len(dummy_atom_offsets)):
            structure.add_atom(pmd.Atom(name="DUM"), f"DM{index + 1}", 1)

        structure.positions = full_coordinates

    # def __init__(self, host, guest=None,
    #              backend="openmm", directory_path="benchmarks",
    #              additional_benchmarks=None, generate_gaff_files=False,
    #              gaff_version="gaff2",
    #              guest_orientation=None, build=True):
    #     self.host = host
    #     self.guest = guest if guest is not None else "release"
    #     self.backend = backend
    #     self.directory = Path(directory_path).joinpath(self.host).joinpath(
    #     f"{self.guest}-{guest_orientation}" if
    #     guest_orientation is not None else
    #     f"{self.guest}")
    #     self.desolvated_window_paths = []
    #     self.window_list = []
    #
    #     if self.backend == "amber":
    #         # Generate `frcmod` and dummy atom files.
    #         raise NotImplementedError
    #
    #     self.directory.mkdir(parents=True, exist_ok=True)
    #     installed_benchmarks = get_benchmarks()
    #     if additional_benchmarks is not None:
    #         installed_benchmarks.update(additional_benchmarks)
    #
    #     host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks,
    #     guest_orientation)
    #
    #     self.benchmark_path = host_yaml.parent
    #     self.host_yaml = read_yaml(host_yaml)
    #     if guest:
    #         self.guest_yaml = read_yaml(guest_yaml["yaml"])
    #
    #     if build:
    #         # Here, we build desolvated windows and pass the files to the OpenFF
    #         Evaluator.
    #         # These files are stored in `self.desolvated_window_paths`.
    #         self.build_desolvated_windows(guest_orientation)
    #         if generate_gaff_files:
    #             generate_gaff(mol2_file=self.benchmark_path.joinpath(
    #             self.host_yaml["structure"]),
    #                           residue_name=self.host_yaml["resname"],
    #                           output_name=self.host,
    #                           directory_path=self.directory,
    #                           gaff=gaff_version)
    #             if guest:
    #                 generate_gaff(mol2_file=self.benchmark_path.joinpath(
    #                     self.guest).joinpath(self.guest_yaml["structure"]),
    #                               output_name=self.guest,
    #                               residue_name=self.guest_yaml["name"],
    #                               directory_path=self.directory,
    #                               gaff=gaff_version)
    #     if not build:
    #         self.populate_window_list(input_pdb=os.path.join(self.directory,
    #         f"{self.host}-{self.guest}.pdb" if self.guest is not None
    #         else f"{self.host}.pdb"))

    # def parse_yaml(self, installed_benchmarks, guest_orientation):
    #     """
    #     Read the YAML recipe for the host and guest.
    #
    #     Returns
    #     -------
    #
    #     """
    #     try:
    #         if guest_orientation:
    #
    #             host_yaml = installed_benchmarks["host_guest_systems"][self.host]
    #             ["yaml"][guest_orientation]
    #         else:
    #             host_yaml = installed_benchmarks["host_guest_systems"][self.host]
    #             ["yaml"]["p"]
    #
    #     except KeyError:
    #         logger.error(f"Cannot find YAML recipe for host: {self.host}")
    #         logger.debug(installed_benchmarks)
    #         raise FileNotFoundError
    #     try:
    #         guest_yaml = installed_benchmarks["host_guest_systems"][self.host]
    #         [self.guest]
    #     except KeyError:
    #         if self.guest == "release":
    #             guest_yaml = None
    #         else:
    #             logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
    #             logger.debug(installed_benchmarks)
    #             raise FileNotFoundError
    #
    #     return host_yaml, guest_yaml

    @classmethod
    def build_static_restraints(
        cls,
        coordinate_path: str,
        n_attach_windows: Optional[int],
        n_pull_windows: Optional[int],
        n_release_windows: Optional[int],
        restraint_schemas: List[Dict[str, Any]],
        use_amber_indices: bool = False,
    ) -> List[DAT_restraint]:
        """

        Parameters
        ----------
        coordinate_path
        n_attach_windows
        n_pull_windows
        n_release_windows
        restraint_schemas

            This should be a dictionary of the form:

                * atoms: List[str]
                * force_constant: float

        use_amber_indices

        Returns
        -------

        """

        if n_pull_windows is not None:
            assert n_attach_windows is not None

        # restraint_schemas = self.host_yaml["restraints"]["static"]

        static_restraints: List[DAT_restraint] = []

        n_windows = [n_attach_windows, n_pull_windows, n_release_windows]

        for restraint_schema in restraint_schemas:

            static = static_DAT_restraint(
                restraint_mask_list=restraint_schema["atoms"].split(),
                num_window_list=n_windows,
                ref_structure=coordinate_path,
                force_constant=restraint_schema["force_constant"],
                amber_index=use_amber_indices,
                continuous_apr=False,
            )
            static_restraints.append(static)

        return static_restraints

    @classmethod
    def build_conformational_restraints(
        cls,
        coordinate_path: str,
        attach_lambdas: Optional[List[float]],
        n_pull_windows: Optional[int],
        release_lambdas: Optional[List[float]],
        restraint_schemas: List[Dict[str, Any]],
        use_amber_indices: bool = False,
    ) -> List[DAT_restraint]:
        """

        Parameters
        ----------
        coordinate_path
        attach_lambdas
        n_pull_windows
        release_lambdas
        restraint_schemas

            This should be a dictionary of the form:

                * atoms: List[str]
                * force_constant: float
                * target: float

        use_amber_indices

        Returns
        -------

        """

        if n_pull_windows is not None:
            assert attach_lambdas is not None

        # restraint_schemas = self.host_yaml["restraints"]["conformational"]

        # attach_lambdas = self.host_yaml["calculation"]["lambda"]["attach"]
        # release_lambdas = self.host_yaml["calculation"]["lambda"]["release"]

        restraints = []

        for restraint_schema in restraint_schemas:

            mask = restraint_schema["atoms"].split()

            restraint = DAT_restraint()
            restraint.amber_index = use_amber_indices
            restraint.topology = coordinate_path
            restraint.mask1 = mask[0]
            restraint.mask2 = mask[1]
            restraint.mask3 = mask[2] if len(mask) > 2 else None
            restraint.mask4 = mask[3] if len(mask) > 3 else None
            restraint.auto_apr = True
            restraint.continuous_apr = False

            if attach_lambdas:

                restraint.attach["target"] = restraint_schema["target"]
                restraint.attach["fc_final"] = restraint_schema["force_constant"]
                restraint.attach["fraction_list"] = attach_lambdas

            if n_pull_windows:

                restraint.pull["target_final"] = restraint_schema["target"]
                restraint.pull["num_windows"] = n_pull_windows

            if release_lambdas:

                restraint.auto_apr = False
                restraint.release["target"] = restraint_schema["target"]
                restraint.release["fc_final"] = restraint_schema["force_constant"]
                restraint.release["fraction_list"] = release_lambdas

            restraint.initialize()
            restraints.append(restraint)

        return restraints

    @classmethod
    def build_symmetry_restraints(
        cls,
        coordinate_path: str,
        n_attach_windows: int,
        restraint_schemas: List[Dict[str, Any]],
        use_amber_indices: bool = False,
    ) -> List[DAT_restraint]:
        """

        Parameters
        ----------
        coordinate_path
        n_attach_windows
        restraint_schemas

            This should be a dictionary of the form:

                * atoms: List[str]
                * force_constant: float

        use_amber_indices

        Returns
        -------

        """

        # restraint_schemas = self.guest_yaml["symmetry_correction"]["restraints"]

        restraints = []

        for restraint_schema in restraint_schemas:

            restraint = DAT_restraint()
            restraint.auto_apr = True
            restraint.continuous_apr = False
            restraint.amber_index = use_amber_indices
            restraint.topology = coordinate_path
            restraint.mask1 = restraint_schema["atoms"].split()[0]
            restraint.mask2 = restraint_schema["atoms"].split()[1]
            restraint.mask3 = restraint_schema["atoms"].split()[2]

            restraint.attach["fc_final"] = restraint_schema["force_constant"]
            restraint.attach["fraction_list"] = [1.0] * n_attach_windows

            # This target should be overridden by the custom values.
            restraint.attach["target"] = 999.99
            restraint.custom_restraint_values["r2"] = 91
            restraint.custom_restraint_values["r3"] = 91

            # 0 force constant between 91 degrees and 180 degrees.
            restraint.custom_restraint_values["rk3"] = 0.0
            restraint.initialize()

            restraints.append(restraint)

        return restraints

    @classmethod
    def build_wall_restraints(
        cls,
        coordinate_path: str,
        n_attach_windows: int,
        restraint_schemas: List[Dict[str, Any]],
        use_amber_indices: bool = False,
    ) -> List[DAT_restraint]:
        """

        Parameters
        ----------
        coordinate_path
        n_attach_windows
        restraint_schemas

            This should be a dictionary of the form:

                * atoms: List[str]
                * force_constant: float
                * target: float

        use_amber_indices

        Returns
        -------

        """

        # restraint_schemas = self.guest_yaml["restraints"]["wall_restraints"]

        restraints = []

        for restraint_schema in restraint_schemas:

            restraint = DAT_restraint()
            restraint.auto_apr = True
            restraint.continuous_apr = False
            restraint.amber_index = use_amber_indices
            restraint.topology = coordinate_path
            restraint.mask1 = restraint_schema["atoms"].split()[0]
            restraint.mask2 = restraint_schema["atoms"].split()[1]

            restraint.attach["fc_final"] = restraint_schema["force_constant"]
            restraint.attach["fraction_list"] = [1.0] * n_attach_windows
            restraint.attach["target"] = restraint_schema["target"]

            # Minimum distance is 0 Angstrom
            restraint.custom_restraint_values["r1"] = 0
            restraint.custom_restraint_values["r2"] = 0

            # Harmonic force constant beyond target distance.
            restraint.custom_restraint_values["rk2"] = restraint_schema[
                "force_constant"
            ]
            restraint.custom_restraint_values["rk3"] = restraint_schema[
                "force_constant"
            ]

            restraint.initialize()

            restraints.append(restraint)

        return restraints

    @classmethod
    def build_guest_restraints(
        cls,
        coordinate_path: str,
        attach_lambdas: List[float],
        n_pull_windows: Optional[int],
        restraint_schemas: List[Dict[str, Any]],
        use_amber_indices: bool = False,
    ) -> List[DAT_restraint]:
        """

        Parameters
        ----------
        coordinate_path
        attach_lambdas
        n_pull_windows
        restraint_schemas

            This should be a dictionary of the form:

                * atoms: List[str]
                * force_constant: float
                * target: float

        use_amber_indices

        Returns
        -------

        """

        # restraint_schemas = self.guest_yaml["restraints"]["guest"]
        # attach_lambdas = self.host_yaml["calculation"]["lambda"]["attach"]

        restraints = []

        for restraint_schema in restraint_schemas:

            mask = restraint_schema["atoms"].split()

            guest_restraint = DAT_restraint()
            guest_restraint.auto_apr = True
            guest_restraint.continuous_apr = False
            guest_restraint.amber_index = use_amber_indices
            guest_restraint.topology = coordinate_path
            guest_restraint.mask1 = mask[0]
            guest_restraint.mask2 = mask[1]
            guest_restraint.mask3 = mask[2] if len(mask) > 2 else None
            guest_restraint.mask4 = mask[3] if len(mask) > 3 else None

            guest_restraint.attach["target"] = restraint_schema["attach"]["target"]
            guest_restraint.attach["fc_final"] = restraint_schema["attach"][
                "force_constant"
            ]
            guest_restraint.attach["fraction_list"] = attach_lambdas

            if n_pull_windows:

                guest_restraint.pull["target_final"] = restraint_schema["pull"][
                    "target"
                ]
                guest_restraint.pull["num_windows"] = n_pull_windows

            guest_restraint.initialize()
            restraints.append(guest_restraint)

        return restraints

    @classmethod
    def apply_openmm_positional_restraints(
        cls, coordinate_path, system, force_group: int = 15
    ):
        """

        Parameters
        ----------
        coordinate_path
        system
        force_group

        Returns
        -------

        """

        from simtk import openmm, unit

        # noinspection PyTypeChecker
        structure: pmd.Structure = pmd.load_file(coordinate_path, structure=True)

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
                # But then we can't access atom indices. Using `atom.xx` works for
                # coordinates, but is unitless.
                k = 50.0 * unit.kilocalories_per_mole / unit.angstroms ** 2

                x0 = 0.1 * atom.xx * unit.nanometers
                y0 = 0.1 * atom.xy * unit.nanometers
                z0 = 0.1 * atom.xz * unit.nanometers

                positional_restraint.addParticle(atom.idx, [k, x0, y0, z0])
                system.addForce(positional_restraint)
                positional_restraint.setForceGroup(force_group)


def apply_openmm_restraints(
    system, restraint, phase, window_number, flat_bottom=False, force_group=None
):

    from simtk import openmm, unit

    assert phase in {"attach", "pull", "release"}

    if flat_bottom and phase == "attach" and restraint.mask3:
        flat_bottom_force = openmm.CustomAngleForce(
            "step(-(theta - theta_0)) * k * (theta - theta_0)^2"
        )
        # If theta is greater than theta_0, then the argument to step is negative,
        # which means the force is off.
        flat_bottom_force.addPerAngleParameter("k")
        flat_bottom_force.addPerAngleParameter("theta_0")

        theta_0 = 91.0 * unit.degrees
        k = (
            restraint.phase[phase]["force_constants"][window_number]
            * unit.kilocalories_per_mole
            / unit.radian ** 2
        )
        flat_bottom_force.addAngle(
            restraint.index1[0], restraint.index2[0], restraint.index3[0], [k, theta_0],
        )
        system.addForce(flat_bottom_force)
        if force_group:
            flat_bottom_force.setForceGroup(force_group)

        return

    elif flat_bottom and phase == "attach" and not restraint.mask3:
        flat_bottom_force = openmm.CustomBondForce("step((r - r_0)) * k * (r - r_0)^2")
        # If x is greater than x_0, then the argument to step is positive, which means
        # the force is on.
        flat_bottom_force.addPerBondParameter("k")
        flat_bottom_force.addPerBondParameter("r_0")

        r_0 = restraint.phase[phase]["targets"][window_number] * unit.angstrom
        k = (
            restraint.phase[phase]["force_constants"][window_number]
            * unit.kilocalories_per_mole
            / unit.radian ** 2
        )
        flat_bottom_force.addBond(
            restraint.index1[0], restraint.index2[0], [k, r_0],
        )
        system.addForce(flat_bottom_force)
        if force_group:
            flat_bottom_force.setForceGroup(force_group)

        return

    elif flat_bottom and phase == "pull":
        return
    elif flat_bottom and phase == "release":
        return

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

        if force_group:
            bond_restraint.setForceGroup(force_group)

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
        if force_group:
            angle_restraint.setForceGroup(force_group)

    elif restraint.mask4:
        if (
            not restraint.group1
            and not restraint.group2
            and not restraint.group3
            and not restraint.group4
        ):
            dihedral_restraint = openmm.CustomTorsionForce(
                f"k * min(min(abs(theta - theta_0), abs(theta - theta_0 + 2 * "
                f"{_PI_})), abs(theta - theta_0 - 2 * {_PI_}))^2"
            )
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
        if force_group:
            dihedral_restraint.setForceGroup(force_group)


def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks


def generate_gaff(
    mol2_file,
    residue_name,
    output_name=None,
    need_gaff_atom_types=True,
    generate_frcmod=True,
    directory_path="benchmarks",
    gaff="gaff2",
):

    if output_name is None:
        output_name = mol2_file.stem

    if need_gaff_atom_types:
        _generate_gaff_atom_types(
            mol2_file=mol2_file,
            residue_name=residue_name,
            output_name=output_name,
            gaff=gaff,
            directory_path=directory_path,
        )
        logging.debug(
            "Checking to see if we have a multi-residue MOL2 file that should be converted "
            "to single-residue..."
        )
        structure = pmd.load_file(
            os.path.join(directory_path, f"{output_name}.{gaff}.mol2"), structure=True
        )
        if len(structure.residues) > 1:
            structure[":1"].save("tmp.mol2")
            if os.path.exists("tmp.mol2"):
                os.rename(
                    "tmp.mol2",
                    os.path.join(directory_path, f"{output_name}.{gaff}.mol2"),
                )
                logging.debug("Saved single-residue MOL2 file for `tleap`.")
            else:
                raise RuntimeError(
                    "Unable to convert multi-residue MOL2 file to single-residue for `tleap`."
                )

        if generate_frcmod:
            _generate_frcmod(
                mol2_file=f"{output_name}.{gaff}.mol2",
                gaff=gaff,
                output_name=output_name,
                directory_path=directory_path,
            )
        else:
            raise NotImplementedError()


def _generate_gaff_atom_types(
    mol2_file, residue_name, output_name, gaff="gaff2", directory_path="benchmarks"
):

    p = subprocess.Popen(
        [
            "antechamber",
            "-i",
            str(mol2_file),
            "-fi",
            "mol2",
            "-o",
            f"{output_name}.{gaff}.mol2",
            "-fo",
            "mol2",
            "-rn",
            f"{residue_name.upper()}",
            "-at",
            f"{gaff}",
            "-an",
            "no",
            "-dr",
            "no",
            "-pf",
            "yes",
        ],
        cwd=directory_path,
    )
    p.communicate()

    files = [
        "ANTECHAMBER_AC.AC",
        "ANTECHAMBER_AC.AC0",
        "ANTECHAMBER_BOND_TYPE.AC",
        "ANTECHAMBER_BOND_TYPE.AC0",
        "ATOMTYPE.INF",
    ]
    files = [directory_path.joinpath(i) for i in files]
    for file in files:
        if file.exists():
            logger.debug(f"Removing temporary file: {file}")
            file.unlink()

    if not os.path.exists(f"{output_name}.{gaff}.mol2"):
        # Try with the newer (AmberTools 19) version of `antechamber` which doesn't have the `-dr` flag
        p = subprocess.Popen(
            [
                "antechamber",
                "-i",
                str(mol2_file),
                "-fi",
                "mol2",
                "-o",
                f"{output_name}.{gaff}.mol2",
                "-fo",
                "mol2",
                "-rn",
                f"{residue_name.upper()}",
                "-at",
                f"{gaff}",
                "-an",
                "no",
                "-pf",
                "yes",
            ],
            cwd=directory_path,
        )
        p.communicate()

        files = [
            "ANTECHAMBER_AC.AC",
            "ANTECHAMBER_AC.AC0",
            "ANTECHAMBER_BOND_TYPE.AC",
            "ANTECHAMBER_BOND_TYPE.AC0",
            "ATOMTYPE.INF",
        ]
        files = [directory_path.joinpath(i) for i in files]
        for file in files:
            if file.exists():
                logger.debug(f"Removing temporary file: {file}")
                file.unlink()


def _generate_frcmod(mol2_file, gaff, output_name, directory_path="benchmarks"):
    subprocess.Popen(
        [
            "parmchk2",
            "-i",
            str(mol2_file),
            "-f",
            "mol2",
            "-o",
            f"{output_name}.{gaff}.frcmod",
            "-s",
            f"{gaff}",
        ],
        cwd=directory_path,
    )
