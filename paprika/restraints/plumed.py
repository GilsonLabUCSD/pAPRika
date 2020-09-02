import os
import logging
import numpy as np
import parmed as pmd

from paprika.utils import get_key, return_parmed_structure
from paprika.restraints.utils import get_bias_potential_type, parse_window
from parmed.structure import Structure as ParmedStructureClass

logger = logging.getLogger(__name__)

_PI_ = np.pi


class Plumed(object):
    """
    This class is converts restraints generated in pAPRika DAT_restraints into Plumed restraints.

    Example:
    -------
        >>> plumed = Plumed()
        >>> plumed.file_name = 'plumed.dat'
        >>> plumed.path = './windows'
        >>> plumed.window_list = window_list
        >>> plumed.restraint_list = restraint_list
        >>> plumed.dump_to_file()


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
    def legacy_k(self):
        """
        bool: Option to specify whether the force constant parsed into DAT_restraint()
        is halved or not (i.e. value is already multiplied by a factor of 1/2).
        """
        return self._legacy_k

    @legacy_k.setter
    def legacy_k(self, value: bool):
        self._legacy_k = value

    @property
    def energy_unit(self):
        """
        str: Units for energies in Plumed.
        """
        return self._energy_unit

    @energy_unit.setter
    def energy_unit(self, value: bool):
        self._energy_unit = value

    @property
    def length_unit(self):
        """
        str: Units for length in Plumed.
        """
        return self._length_unit

    @length_unit.setter
    def length_unit(self, value: bool):
        self._length_unit = value

    @property
    def time_unit(self):
        """
        str: Units for time in Plumed.
        """
        return self._time_unit

    @time_unit.setter
    def time_unit(self, value: bool):
        self._time_unit = value

    def __init__(self):
        self._file_name = "plumed.dat"
        self._restraint_list = None
        self._window_list = None
        self._path = './'
        self._legacy_k = True
        self.k_factor = 1.0
        self._energy_unit = "kcal/mol"
        self._length_unit = "A"
        self._time_unit = "ns"
        self.header_line = None

    def _initialize(self):
        # Set factor for spring constant
        if self.legacy_k:
            self.k_factor = 2.0

        # Check units
        # NOTE: We have to resort to strings until we migrate all quantities to Pint/SimTK units
        _check_plumed_units(self.energy_unit, self.length_unit, self.time_unit)

        # header line
        self.header_line = f"UNITS LENGTH={self.length_unit} ENERGY={self.energy_unit} TIME={self.time_unit}"

    def dump_to_file(self):

        self._initialize()

        # Loop over APR windows
        for windows in self.window_list:

            window, phase = parse_window(windows)

            # Check if file exist and write header line
            try:
                with open(os.path.join(self.path, windows, self.file_name), "w") as file:
                    file.write(self.header_line + "\n")

            except IOError:
                raise Exception("Cannot write header to file. Please check the specified path exist.")

            print(f"Converting pAPRika restraint in {windows} to Plumed")

            cv_index = 1

            center_list = []
            colvar_list = []
            bias_list = []

            group_index = 1
            group_atoms = {}

            # Parse each restraint in the list
            for restraint in self.restraint_list:
                # Skip restraint if the target or force constant is not defined.
                # Example: wall restraints only used during the attach phase.
                try:
                    target = restraint.phase[phase]["targets"][window]
                    force_constant = restraint.phase[phase]["force_constants"][window] * self.k_factor
                except TypeError:
                    continue

                # Check atom index setting
                index_shift = 0
                if not restraint.amber_index:
                    index_shift = 1
                    logger.info("Atom indices starts from 0 --> shifting indices by 1.")

                # Collect DAT atom indices
                atom_index = []

                if not restraint.index1:
                    raise Exception("There must be at least two atoms in a restraint.")
                elif not restraint.group1:
                    atom_index.append(restraint.index1[0] + index_shift)
                else:
                    igr1 = ""
                    for index in restraint.index1:
                        igr1 += "{},".format(index + index_shift)

                    if not get_key(group_atoms, igr1):
                        group_atoms[f"g{group_index}"] = igr1
                        group_index += 1

                    atom_index.append(get_key(group_atoms, igr1)[0])

                if not restraint.index2:
                    raise Exception("There must be at least two atoms in a restraint.")
                elif not restraint.group2:
                    atom_index.append(restraint.index2[0] + index_shift)
                else:
                    igr2 = ""
                    for index in restraint.index2:
                        igr2 += "{},".format(index + index_shift)

                    if not get_key(group_atoms, igr2):
                        group_atoms[f"g{group_index}"] = igr2
                        group_index += 1

                    atom_index.append(get_key(group_atoms, igr2)[0])

                if not restraint.index3:
                    pass
                elif not restraint.group3:
                    atom_index.append(restraint.index3[0] + index_shift)
                else:
                    igr3 = ""
                    for index in restraint.index3:
                        igr3 += "{},".format(index + index_shift)

                    if not get_key(group_atoms, igr3):
                        group_atoms[f"g{group_index}"] = igr3
                        group_index += 1

                    atom_index.append(get_key(group_atoms, igr3)[0])

                if not restraint.index4:
                    pass
                elif not restraint.group4:
                    atom_index.append(restraint.index4[0] + index_shift)
                else:
                    igr4 = ""
                    for index in restraint.index4:
                        igr4 += "{},".format(index + index_shift)

                    if not get_key(group_atoms, igr4):
                        group_atoms[f"g{group_index}"] = igr4
                        group_index += 1

                    atom_index.append(get_key(group_atoms, igr4)[0])

                # Convert list to comma-separated string
                atom_string = ','.join(map(str, atom_index))

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
                if bias_type is "harmonic":
                    bias_type = "restraint"

                # Append string to lists
                colvar_list.append(f"c{cv_index:<2}: {colvar_type.upper():<8} ATOMS={atom_string} NOPBC\n")
                bias_list.append(f"{bias_type.upper():<11} ARG=c{cv_index:<2} AT={target:>7.3f} KAPPA={force_constant:>7.2f}\n")

                # Increment colvar index
                cv_index += 1

            # Write info to file
            with open(os.path.join(self.path, windows, self.file_name), "a") as file:
                if len(group_atoms) != 0:
                    file.write("# Center groups\n")
                    for key, value in group_atoms.items():
                        file.write(f"{key:<3}: COM ATOMS={value}\n")

                file.write("# Collective variables\n")
                for line in colvar_list:
                    file.write(line)

                file.write("# Bias potentials\n")
                for line in bias_list:
                    file.write(line)

    def add_dummy_atoms_to_file(self, structure):
        from paprika.dummy import extract_dummy_atoms
        from paprika.utils import return_parmed_structure
        
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

            # Write dummy atom info to 'plumed.dat'
            plumed_file = os.path.join(self.path, windows, self.file_name)
            if os.path.isfile(plumed_file):
                with open(plumed_file, "a") as file:
                    _write_dummy_to_file(file, dummy_atoms)
            else:
                raise Exception(
                    f"ERROR: '{plumed_file}' file does not exists!"
                )


def _check_plumed_units(energy, length, time):
    """
    Checks the specified units and makes sure that its supported
    """
    if energy not in ['kj/mol', 'kcal/mol']:
        raise Exception(f"Specified unit for energy ({energy}) is not supported.")

    if length not in ['nm', 'A']:
        raise Exception(f"Specified unit for length ({length}) is not supported.")

    if time not in ['ps', 'fs', 'ns']:
        raise Exception(f"Specified unit for time ({time}) is not supported.")


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
