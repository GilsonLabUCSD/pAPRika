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


class Setup:
    """
    The Setup class provides a wrapper function around the preparation of the host-guest
    system and the application of restraints.
    """

    @classmethod
    def prepare_host_structure(
        cls, coordinate_path: str, host_atom_indices: Optional[List[int]] = None
    ) -> pmd.Structure:
        """Prepares the coordinates of a host molecule ready for the release phase of
        an APR calculation.
        This currently involves aligning the cavity of the host along the z-axis, and
        positioning the host so that its center of geometry at (0, 0, 0).

        Notes
        -----
        The hosts cavity axis is for now determined as the vector orthogonal to the
        two largest principal components of the host.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which contains the host molecule.
        host_atom_indices
            The (0-based) indices of the host atoms in the coordinate file if the
            file. This may be used if the coordinate file contains more than a single
            molecule.

        Returns
        -------
            A ParmEd structure which contains only the aligned and centered host.
            If ``host_atom_indices`` are provided, the structure will contain only
            the referenced host atoms.
        """

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
        """Prepares the coordinates of a host molecule ready for the pull (+ attach)
        phase of an APR calculation.

        This currently involves aligning the complex so that the guest molecule (or
        rather, the guest atoms specified by ``guest_orientation_mask``) is aligned
        with the z-axis, and the first atom specified by ``guest_orientation_mask`` is
        positioned at (0, 0, 0).

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which contains the guest molecule
            bound to the host.
        guest_atom_indices
            The (0-based) indices of the atoms in the coordinate file which
            correspond to the guest molecule.
        guest_orientation_mask
            The string mask which describes which guest atoms will be restrained to
            keep the molecule aligned to the z-axis and at a specific distance from
            the host. This should be of the form 'X Y' where X Y are ParmEd selectors
            for the two guest atoms to restrain relative to the dummy atoms.
        pull_distance
            The total distance that the guest will be pulled along the z-axis during
            the pull phase in units of Angstroms.
        pull_window_index
            The index of the window to prepare coordinates for. This determines the
            distance to place the guest from the host in the returned structure. An
            index of zero corresponds to the guest in its initial position, while an
            ``index = n_pull_windows - `` corresponds to the guest positioned
            a distance of ``pull_distance`` away from the host.
        n_pull_windows
            The total number of pull windows being used in the calculation. This
            will determine the distance to move the guest at each window.

        Returns
        -------
            A ParmEd structure which contains the aligned and positioned host and
            guest molecules.
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
        """A convenience method to add a number of dummy atoms to an existing
        ParmEd structure, and to position those atoms at a specified set of positions.

        Parameters
        ----------
        structure
            The structure to add the dummy atoms to.
        dummy_atom_offsets
            The list of positions (defined by a 3-d numpy array) of the dummy atoms
            to add.
        offset_coordinates
            An optional amount to offset each of the dummy atom positions by with
            shape=(3,).
        """

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
        """A method to convert a set of static restraints defined by their 'schemas'
        into corresponding ``DAT_restraint``objects.

        Each 'schema' should be a dictionary with:

            * an ``atoms`` entry with a value of the atom selection make which specifies
              which atoms the restraint will apply to
            * a ``force_constant`` entry which specifies the force constant of the
              restraint.

        These 'schemas` map directly to the 'restraints -> static -> restraint'
        dictionaries specified in the `taproom` host YAML files.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which the restraints will be applied to.
            This should contain either the host or the complex, the dummy atoms and
            and solvent.
        n_attach_windows
            The total number of attach windows being used in the APR calculation.
        n_pull_windows
            The total number of pull windows being used in the APR calculation.
        n_release_windows
            The total number of release windows being used in the APR calculation.
        restraint_schemas
            The list of dictionaries which provide the settings to use for each
            static restraint to add.
        use_amber_indices
            Whether to use amber based (i.e. starting from 1) restraint indices or
            OpenMM based (i.e. starting from 0) indices.

        Returns
        -------
            The constructed static restraint objects.
        """

        if n_pull_windows is not None:
            assert n_attach_windows is not None

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
        """A method to convert a set of conformational restraints defined by their
        'schemas' into corresponding ``DAT_restraint``objects.

        Each 'schema' should be a dictionary with:

            * an ``atoms`` entry with a value of the atom selection make which specifies
              which atoms the restraint will apply to
            * a ``force_constant`` entry which specifies the force constant of the
              restraint.
            * a ``target`` entry which specifies the target value of the restraint.

        These 'schemas` map directly to the 'restraints -> conformational -> restraint'
        dictionaries specified in the ``taproom`` host YAML files.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which the restraints will be applied to.
            This should contain either the host or the complex, the dummy atoms and
            and solvent.
        attach_lambdas
            The values 'lambda' being used during the attach phase of the APR
            calculation.
        n_pull_windows
            The total number of pull windows being used in the APR calculation.
        release_lambdas
            The values 'lambda' being used during the release phase of the APR
            calculation.
        restraint_schemas
            The list of dictionaries which provide the settings to use for each
            conformational restraint to add.
        use_amber_indices
            Whether to use amber based (i.e. starting from 1) restraint indices or
            OpenMM based (i.e. starting from 0) indices.

        Returns
        -------
            The constructed conformational restraint objects.
        """

        if n_pull_windows is not None:
            assert attach_lambdas is not None

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
        """A method to convert a set of symmetry restraints defined by their 'schemas'
        into corresponding ``DAT_restraint``objects.

        Each 'schema' should be a dictionary with:

            * an ``atoms`` entry with a value of the atom selection make which specifies
              which atoms the restraint will apply to
            * a ``force_constant`` entry which specifies the force constant of the
              restraint.

        These 'schemas` map directly to the 'restraints -> symmetry_correction
        -> restraint' dictionaries specified in the `taproom` guest YAML files.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which the restraints will be applied to.
            This should contain either the host or the complex, the dummy atoms and
            and solvent.
        n_attach_windows
            The total number of attach windows being used in the APR calculation.
        restraint_schemas
            The list of dictionaries which provide the settings to use for each
            symmetry restraint to add.
        use_amber_indices
            Whether to use amber based (i.e. starting from 1) restraint indices or
            OpenMM based (i.e. starting from 0) indices.

        Returns
        -------
            The constructed symmetry restraint objects.
        """

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
        """A method to convert a set of wall restraints defined by their 'schemas'
        into corresponding ``DAT_restraint``objects.

        Each 'schema' should be a dictionary with:

            * an ``atoms`` entry with a value of the atom selection make which specifies
              which atoms the restraint will apply to
            * a ``force_constant`` entry which specifies the force constant of the
              restraint.
            * a ``target`` entry which specifies the target value of the restraint.

        These 'schemas` map directly to the 'restraints -> wall_restraints -> restraint'
        dictionaries specified in the `taproom` guest YAML files.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which the restraints will be applied to.
            This should contain either the host or the complex, the dummy atoms and
            and solvent.
        n_attach_windows
            The total number of attach windows being used in the APR calculation.
        restraint_schemas
            The list of dictionaries which provide the settings to use for each
            wall restraint to add.
        use_amber_indices
            Whether to use amber based (i.e. starting from 1) restraint indices or
            OpenMM based (i.e. starting from 0) indices.

        Returns
        -------
            The constructed wall restraint objects.
        """

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
        """A method to convert a set of guest restraints defined by their 'schemas'
        into corresponding ``DAT_restraint``objects.

        Each 'schema' should be a dictionary with:
            * an ``atoms`` entry with a value of the atom selection make which specifies
              which atoms the restraint will apply to

        and additionally a nested ``attach`` and ``pull`` dictionary with

            * a ``force_constant`` entry which specifies the force constant of the
              restraint.
            * a ``target`` entry which specifies the target value of the restraint.

        These 'schemas` map directly to the 'restraints -> guest -> restraint'
        dictionaries specified in the `taproom` guest YAML files.

        Parameters
        ----------
        coordinate_path
            The path to the coordinate file which the restraints will be applied to.
            This should contain either the host or the complex, the dummy atoms and
            and solvent.
        attach_lambdas
            The values 'lambda' being used during the attach phase of the APR
            calculation.
        n_pull_windows
            The total number of pull windows being used in the APR calculation.
        restraint_schemas
            The list of dictionaries which provide the settings to use for each
            wall restraint to add.
        use_amber_indices
            Whether to use amber based (i.e. starting from 1) restraint indices or
            OpenMM based (i.e. starting from 0) indices.

        Returns
        -------
            The constructed wall restraint objects.
        """

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


def get_benchmarks():
    """
    Determine the installed `taproom` benchmarks.
    """
    installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        installed_benchmarks[entry_point.name] = entry_point.load()

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
