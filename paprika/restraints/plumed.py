import logging
import os

import numpy as np
from paprika.dummy import extract_dummy_atoms
from paprika.restraints.utils import get_bias_potential_type, parse_window
from paprika.utils import get_key, return_parmed_structure
from parmed.structure import Structure as ParmedStructureClass

logger = logging.getLogger(__name__)

_PI_ = np.pi


class Plumed:
    """
    This class converts restraints generated in pAPRika DAT_restraints into Plumed restraints.

    Example:
    -------
        >>> plumed = Plumed()
        >>> plumed.file_name = 'plumed.dat'
        >>> plumed.path = './windows'
        >>> plumed.window_list = window_list
        >>> plumed.restraint_list = restraint_list
        >>> plumed.dump_to_file()

    plumed.dat:
        UNITS LENGTH=A ENERGY=kcal/mol TIME=ns
        # Collective variables
        c1 : DISTANCE ATOMS=175,150 NOPBC
        c2 : ANGLE    ATOMS=176,175,150 NOPBC
        c3 : ANGLE    ATOMS=175,150,165 NOPBC
        # Bias potential
        RESTRAINT   ARG=c7  AT=  6.000 KAPPA=  10.00
        RESTRAINT   ARG=c8  AT=  3.142 KAPPA= 200.00
        RESTRAINT   ARG=c9  AT=  3.142 KAPPA= 200.00
    """

    @property
    def file_name(self):
        """
        str: The Plumed file name where the restraints will be written to.
        """
        return self._file_name

    @file_name.setter
    def file_name(self, value: str):
        self._file_name = value

    @property
    def window_list(self):
        """
        list: The list of APR windows
        """
        return self._window_list

    @window_list.setter
    def window_list(self, value: list):
        self._window_list = value

    @property
    def restraint_list(self):
        """
        list: The list of restraints to be converted to Plumed.
        """
        return self._restraint_list

    @restraint_list.setter
    def restraint_list(self, value: list):
        self._restraint_list = value

    @property
    def path(self):
        """
        os.PathLike: The parent directory that contains the simulation windows.
        """
        return self._path

    @path.setter
    def path(self, value: str):
        self._path = value

    @property
    def uses_legacy_k(self):
        """
        bool: Option to specify whether the force constant parsed into DAT_restraint()
        is Amber-style or Gromacs/NAMD-style. Amber-style force constants have their value
        multiplied by a factor of 1/2 whereas Gromacs/NAMD-style does not. Plumed follows
        the Gromacs/NAMD-style for the force constant and the equations below demonstrates
        this point.

            * Amber:  U = K x (x - x0)²    * Plumed: U = 1/2 x k x (x - x0)²
            --> K(Amber) = 1/2 k(Plumed)

        i.e. "uses_legacy_k" is set to True (default) the force constants will be multiplied by 2.
        """
        return self._uses_legacy_k

    @uses_legacy_k.setter
    def uses_legacy_k(self, value: bool):
        self._uses_legacy_k = value

    @property
    def units(self):
        """
        dict: Dictionary of units for Plumed as strings. The dictionary requires the key values
        of 'energy', 'length, 'time'.
        """
        return self._units

    @units.setter
    def units(self, value: dict):
        self._units = value

    def __init__(self):
        self._file_name = "plumed.dat"
        self._restraint_list = None
        self._window_list = None
        self._path = "./"
        self._uses_legacy_k = True
        self.k_factor = 1.0
        self._units = None
        self.header_line = None
        self.group_index = None
        self.group_atoms = None

    def _initialize(self):
        # Set factor for spring constant
        if self.uses_legacy_k:
            self.k_factor = 2.0

        # Check units
        # NOTE: We have to resort to strings until (if) we migrate all quantities to Pint/SimTK units
        if self.units is None:
            self.units = {"energy": "kcal/mol", "length": "A", "time": "ns"}
        _check_plumed_units(self.units)

        # header line
        self.header_line = f"UNITS LENGTH={self.units['length']} ENERGY={self.units['energy']} TIME={self.units['time']}"

    def dump_to_file(self):

        self._initialize()

        # Loop over APR windows
        for windows in self.window_list:

            window, phase = parse_window(windows)

            # Check if file exist and write header line
            with open(os.path.join(self.path, windows, self.file_name), "w") as file:
                file.write(self.header_line + "\n")

            cv_index = 1
            cv_list = []
            bias_list = []

            self.group_index = 1
            self.group_atoms = {}

            # Parse each restraint in the list
            for restraint in self.restraint_list:
                # Skip restraint if the target or force constant is not defined.
                # Example: wall restraints only used during the attach phase.
                try:
                    target = restraint.phase[phase]["targets"][window]
                    force_constant = (
                        restraint.phase[phase]["force_constants"][window]
                        * self.k_factor
                    )
                except TypeError:
                    continue

                # Convert list to comma-separated string
                atom_index = self._get_atom_indices(restraint)
                atom_string = ",".join(map(str, atom_index))

                # Determine collective variable type
                colvar_type = "distance"
                if len(atom_index) == 3:
                    colvar_type = "angle"
                    target *= _PI_ / 180.0
                elif len(atom_index) == 4:
                    colvar_type = "torsion"
                    target *= _PI_ / 180.0

                # Determine bias type for this restraint
                bias_type = get_bias_potential_type(restraint, phase, window)

                # Append string to lists
                cv_list.append(
                    f"c{cv_index}: {colvar_type.upper()} ATOMS={atom_string} NOPBC\n"
                )
                bias_list.append(
                    f"{bias_type.upper()} ARG=c{cv_index} AT={target:.4f} KAPPA={force_constant:.2f}\n"
                )

                # Increment colvar index
                cv_index += 1

            # Write collective variables to file
            self._write_colvar_to_file(windows, cv_list, bias_list)

    def _write_colvar_to_file(self, window, cv_list, bias_list):
        with open(os.path.join(self.path, window, self.file_name), "a") as file:
            if len(self.group_atoms) != 0:
                file.write("# Centroid groups\n")
                for key, value in self.group_atoms.items():
                    file.write(f"{key}: COM ATOMS={value}\n")

            file.write("# Collective variables\n")
            for line in cv_list:
                file.write(line)

            file.write("# Bias potentials\n")
            for line in bias_list:
                file.write(line)

    def _get_atom_indices(self, restraint):
        # Check atom index setting
        index_shift = 0
        if not restraint.amber_index:
            index_shift = 1
            logger.debug("Atom indices starts from 0 --> shifting indices by 1.")

        # Collect DAT atom indices
        atom_index = []
        for i in range(4):
            ii = i + 1

            if not eval(f"restraint.index{ii}"):
                if ii in [1, 2]:
                    raise Exception("There must be at least two atoms in a restraint.")
            elif not eval(f"restraint.group{ii}"):
                atom_index.append(eval(f"restraint.index{ii}")[0] + index_shift)
            else:
                igr1 = ""
                for index in eval(f"restraint.index{ii}"):
                    igr1 += "{},".format(index + index_shift)

                if not get_key(self.group_atoms, igr1):
                    self.group_atoms[f"g{self.group_index}"] = igr1
                    self.group_index += 1

                atom_index.append(get_key(self.group_atoms, igr1)[0])

        return atom_index

    def add_dummy_atoms_to_file(self, structure):
        # Extract dummy atoms
        for windows in self.window_list:
            if isinstance(structure, str):
                structure = return_parmed_structure(structure)
            elif isinstance(structure, ParmedStructureClass):
                pass
            else:
                raise Exception(
                    "add_dummy_atoms_to_file does not support the type associated with structure: "
                    + type(structure)
                )

            dummy_atoms = extract_dummy_atoms(structure, serial=True)

            # Write dummy atom info to plumed file
            plumed_file = os.path.join(self.path, windows, self.file_name)
            if os.path.isfile(plumed_file):
                with open(plumed_file, "a") as file:
                    _write_dummy_to_file(file, dummy_atoms)
            else:
                raise Exception(f"ERROR: '{plumed_file}' file does not exists!")


def _check_plumed_units(units):
    """
    Checks the specified units and makes sure that its supported.
    """
    if units["energy"] not in ["kj/mol", "kcal/mol"]:
        raise Exception(
            f"Specified unit for energy ({units['energy']}) is not supported."
        )

    if units["length"] not in ["nm", "A"]:
        raise Exception(
            f"Specified unit for length ({units['length']}) is not supported."
        )

    if units["time"] not in ["ps", "fs", "ns"]:
        raise Exception(f"Specified unit for time ({units['time']}) is not supported.")


def _write_dummy_to_file(file, dummy_atoms, kpos=100.0):
    """
    Append to the "plumed.dat" file the dummy atoms colvar definition and position restraints

    Parameters
    ----------
    file : class '_io.TextIOWrapper'
        The file object handle to save the plumed file.
    dummy_atoms : dict
        Dictionary containing information about the dummy atoms.
    kpos : float
        Spring constant used to restrain dummy atoms (kcal/mol/A^2).

    Output
    ------
    dm1: POSITION ATOM = 170 NOPBC
    dm2: POSITION ATOM = 171 NOPBC
    dm3: POSITION ATOM = 172 NOPBC
    RESTRAINT...
    ARG = dm1.x, dm1.y, dm1.z, dm2.x, dm2.y, dm2.z, dm3.x, dm3.y, dm3.z
    AT = 19.68, 20.3, 26.9, 19.68, 20.3, 23.9, 19.68, 22.5, 21.7
    KAPPA = 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0
    ...
    RESTRAINT

    """

    file.write("# Dummy Atoms\n")
    file.write(f"dm1: POSITION ATOM={dummy_atoms['DM1']['idx']} NOPBC\n")
    file.write(f"dm2: POSITION ATOM={dummy_atoms['DM2']['idx']} NOPBC\n")
    file.write(f"dm3: POSITION ATOM={dummy_atoms['DM3']['idx']} NOPBC\n")

    arg = "dm1.x,dm1.y,dm1.z," "dm2.x,dm2.y,dm2.z," "dm3.x,dm3.y,dm3.z,"

    at = (
        f"{dummy_atoms['DM1']['pos'][0]:0.3f},"
        f"{dummy_atoms['DM1']['pos'][1]:0.3f},"
        f"{dummy_atoms['DM1']['pos'][2]:0.3f},"
    )
    at += (
        f"{dummy_atoms['DM2']['pos'][0]:0.3f},"
        f"{dummy_atoms['DM2']['pos'][1]:0.3f},"
        f"{dummy_atoms['DM2']['pos'][2]:0.3f},"
    )
    at += (
        f"{dummy_atoms['DM3']['pos'][0]:0.3f},"
        f"{dummy_atoms['DM3']['pos'][1]:0.3f},"
        f"{dummy_atoms['DM3']['pos'][2]:0.3f},"
    )

    kappa = (
        f"{kpos:0.1f},{kpos:0.1f},{kpos:0.1f},"
        f"{kpos:0.1f},{kpos:0.1f},{kpos:0.1f},"
        f"{kpos:0.1f},{kpos:0.1f},{kpos:0.1f},"
    )

    file.write(f"RESTRAINT ...\n")
    file.write(f"ARG={arg}\n")
    file.write(f"AT={at}\n")
    file.write(f"KAPPA={kappa}\n")
    file.write(f"LABEL=dummy\n")
    file.write(f"... RESTRAINT\n")
