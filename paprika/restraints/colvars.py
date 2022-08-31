import logging
import os

import numpy as np
from openff.units import unit as openff_unit

from paprika.restraints.plumed import Plumed
from paprika.restraints.utils import get_bias_potential_type, parse_window
from paprika.utils import check_unit, get_key

logger = logging.getLogger(__name__)

_PI_ = np.pi


class Colvars(Plumed):
    """
    This class converts restraints generated with :class:`paprika.restraints.DAT_restraints` into collective variables
    (`Colvars`) restraints that is available as a plugin in `NAMD` and `LAMMPS`.

    .. note ::
        The ``Colvars`` module is described in the reference below and the source code is available on Github
        https://github.com/Colvars/colvars

        `Fiorin, G., Klein, M. L. & Hénin, J. Using collective variables to drive molecular dynamics simulations.
        Mol. Phys. 111, 3345–3362 (2013).`

    Examples
    --------
        >>> colvars = Colvars()
        >>> colvars.file_name = 'colvars.tcl'
        >>> colvars.path = './windows'
        >>> colvars.window_list = window_list
        >>> colvars.restraint_list = restraint_file
        >>> colvars.dump_to_file()

    The commands above will write the restraints to ``windows/*/colvars.tcl`` and contains the `Colvars`-style
    restraints

    .. code-block::

        ColvarsTrajFrequency 500
        ColvarsRestartFrequency 50000
        # Collective variables
        colvar {
          name c1
          distance {
            forceNoPBC yes
            group1 { atomNumbers 123 }
            group2 { atomNumbers 1 }
          }
        }
        # Bias potentials
        harmonic {
          colvars c1
          centers 6.0
          forceConstant 10.0000
        }

    The positional restraints on dummy atoms, however, is not added automatically. This restraints on the dummy
    atoms can be added to ``windows/*/colvars.tcl`` using the code below.

        >>> for window in window_list:
        >>>     structure = pmd.load_file("topology.prmtop", "coordinates.rst7")
        >>>     colvars.add_dummy_atoms_to_file(structure, window)

    This appends the file with the following

    .. code-block::

        # Dummy atom position restraints
        colvar {
          name dummyAtoms
          cartesian {
            atoms { atomNumbers 123 124 125 }
          }
        }
        harmonic {
          colvars dummyAtoms
          centers ( 0.0, 0.0, -6.0, 0.0, 0.0, -9.0, 0.0, 2.2, -11.2 )
          forceConstant 100.00
        }
    """

    @property
    def output_freq(self) -> int:
        """int: The frequency at which the `colvars` will be printed to ``*.colvars.traj`` file."""
        return self._output_freq

    @output_freq.setter
    def output_freq(self, value: int):
        self._output_freq = value

    def __init__(self):

        super().__init__()

        self._file_name = "colvars.dat"
        self._output_freq = 500
        self._colvars_factor = {}

    def _initialize(self):
        # Set factor for spring constant
        if self.uses_legacy_k:
            self.k_factor = 2.0

        # header line
        self.header_line = (
            f"ColvarsTrajFrequency {self.output_freq}\n"
            f"ColvarsRestartFrequency {self.output_freq*100} "
        )

    def dump_to_file(self):
        """
        Write the `Colvars`-style restraints to file.
        """

        self._initialize()

        # Loop over APR windows
        for window in self.window_list:

            window_number, phase = parse_window(window)

            # Check if file exist and write header line
            with open(os.path.join(self.path, window, self.file_name), "w") as file:
                file.write(self.header_line + "\n")

            cv_index = 1
            cv_dict = {}
            cv_lines = []
            bias_lines = []

            # Parse each restraint in the list
            for restraint in self.restraint_list:
                # Skip restraint if the target or force constant is not defined.
                # Example: wall restraints only used during the attach phase.
                try:
                    target = restraint.phase[phase]["targets"][window_number]
                    force_constant = (
                        restraint.phase[phase]["force_constants"][window_number]
                        * self.k_factor
                    )
                except TypeError:
                    continue

                # Get atom indices in space separated string
                atom_index = self._get_atom_indices(restraint)
                atom_string = " ".join(map(str, atom_index))

                # Determine bias type for half-harmonic potentials
                bias_type, restraint_values = get_bias_potential_type(
                    restraint, phase, window_number, return_values=True
                )
                if bias_type == "upper_walls":
                    target = restraint_values["r3"]
                    force_constant = restraint_values["rk3"] * self.k_factor
                elif bias_type == "lower_walls":
                    target = restraint_values["r2"]
                    force_constant = restraint_values["rk2"] * self.k_factor

                # Convert units to the correct type for COLVAR module
                energy_units = openff_unit.kcal / openff_unit.mole
                if restraint.restraint_type == "distance":
                    target = target.to(openff_unit.angstrom)
                    force_constant = force_constant.to(
                        energy_units / openff_unit.angstrom**2
                    )
                elif (
                    restraint.restraint_type == "angle"
                    or restraint.restraint_type == "torsion"
                ):
                    target = target.to(openff_unit.degrees)
                    force_constant = force_constant.to(
                        energy_units / openff_unit.degrees**2
                    )

                # Append cv to list
                # The code below prevents duplicate cv definition.
                # While not necessary, it makes the plumed file cleaner.
                if not get_key(cv_dict, atom_string):
                    cv_key = f"c{cv_index}"
                    cv_dict[cv_key] = atom_string

                    cv_template_lines = [
                        "colvar {",
                        f"  name {cv_key}",
                        f"  {'dihedral' if restraint.restraint_type == 'torsion' else restraint.restraint_type} {{",
                        "    forceNoPBC yes",
                        f"    group1 {{ atomNumbers {atom_index[0]} }}",
                        f"    group2 {{ atomNumbers {atom_index[1]} }}",
                    ]

                    if restraint.restraint_type == "angle":
                        cv_template_lines += [
                            f"    group3 {{ atomNumbers {atom_index[2]} }}"
                        ]

                    if restraint.restraint_type == "torsion":
                        cv_template_lines += [
                            f"    group3 {{ atomNumbers {atom_index[2]} }}",
                            f"    group4 {{ atomNumbers {atom_index[3]} }}",
                        ]

                    cv_template_lines += ["  }", "}"]

                    cv_lines.append(cv_template_lines)

                    bias_lines.append(
                        self._get_bias_block(bias_type, cv_key, target, force_constant)
                    )

                else:
                    cv_key = get_key(cv_dict, atom_string)[0]

                    bias_lines.append(
                        self._get_bias_block(bias_type, cv_key, target, force_constant)
                    )

                # Increment cv index
                cv_index += 1

            # Write collective variables to file
            self._write_colvar_to_file(window, cv_lines, bias_lines)

    def _write_colvar_to_file(self, window, cv_list, bias_list):
        with open(os.path.join(self.path, window, self.file_name), "a") as file:
            file.write("# Collective variables\n")
            for colvar in cv_list:
                for line in colvar:
                    file.write(line + "\n")

            file.write("# Bias potentials\n")
            for bias in bias_list:
                for line in bias:
                    file.write(line + "\n")

    @staticmethod
    def _get_bias_block(bias_type, cv_key, target, force_constant):
        bias_template_lines = []

        if bias_type == "restraint":
            bias_template_lines = [
                "harmonic {",
                f"  colvars {cv_key}",
                f"  centers {target.magnitude:.4f}",
                f"  forceConstant {force_constant.magnitude:.4f}",
                "}",
            ]
        elif bias_type == "upper_walls":
            bias_template_lines = [
                "harmonicWalls {",
                f"  colvars {cv_key}",
                f"  upperWalls {target.magnitude:.4f}",
                f"  upperWallConstant {force_constant.magnitude:.4f}",
                "}",
            ]
        elif bias_type == "lower_walls":
            bias_template_lines = [
                "harmonicWalls {",
                f"  colvars {cv_key}",
                f"  lowerWalls {target.magnitude:.4f}",
                f"  lowerWallConstant {force_constant.magnitude:.4f}",
                "}",
            ]

        return bias_template_lines

    @staticmethod
    def _get_atom_indices(restraint):
        # Check atom index setting
        index_shift = 0
        if not restraint.amber_index:
            index_shift = 1
            logger.debug("Atom indices starts from 0 --> shifting indices by 1.")

        # Collect DAT atom indices
        atom_index = []

        if not restraint.group1:
            atom_index.append("{}".format(restraint.index1[0] + index_shift))
        else:
            igr1 = ""
            for index in restraint.index1:
                igr1 += "{} ".format(index + index_shift)

            atom_index.append(igr1)

        if not restraint.group2:
            atom_index.append("{}".format(restraint.index2[0] + index_shift))
        else:
            igr2 = ""
            for index in restraint.index2:
                igr2 += "{} ".format(index + index_shift)

            atom_index.append(igr2)

        if restraint.index3 and not restraint.group3:
            atom_index.append("{}".format(restraint.index3[0] + index_shift))
        elif restraint.group3:
            igr3 = ""
            for index in restraint.index3:
                igr3 += "{} ".format(index + index_shift)

            atom_index.append(igr3)

        if restraint.index4 and not restraint.group4:
            atom_index.append("{}".format(restraint.index4[0] + index_shift))
        elif restraint.group4:
            igr4 = ""
            for index in restraint.index4:
                igr4 += "{} ".format(index + index_shift)

            atom_index.append(igr4)

        return atom_index

    @staticmethod
    def _write_dummy_to_file(file, dummy_atoms, kpos=100.0):
        """
        Append to the "colvars.dat" file the dummy atoms colvar definition and position restraints

        Parameters
        ----------
        file : class '_io.TextIOWrapper'
            The file object handle to save the plumed file.
        dummy_atoms : dict
            Dictionary containing information about the dummy atoms.
        kpos : float or openff.unit.Quantity
            Spring constant used to restrain dummy atoms (default for float: kcal/mol/A^2).

        Examples
        --------
        .. code :: tcl

            colvar {
                name dummy
                cartesian {
                    atoms { atomNumbers 1 2 3 }
                }
            }
            harmonic {
                colvars dummy
                centers ( 0.0, 0.0, -6.0, 0.0, 0.0, -9.0, 0.0, 2.2, -11.2)
                forceConstant 100.0
            }
        """
        # Check k units
        kpos = check_unit(
            kpos,
            base_unit=openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2,
        )

        # Get dummy atom indices
        dummy_indices = [dummy_atoms[key]["idx"] for key in dummy_atoms.keys()]
        dummy_index_string = " ".join(map(str, dummy_indices))

        cv_template_lines = [
            "colvar {",
            "  name dummyAtoms",
            "  cartesian {",
            f"    atoms {{ atomNumbers {dummy_index_string} }}",
            "  }",
            "}",
        ]

        # Get dummy atom positions
        dummy_position = []
        for dummy in dummy_atoms.keys():
            dummy_position += list(dummy_atoms[dummy]["pos"])
        dummy_position_string = ", ".join(map(str, dummy_position))

        bias_template_lines = [
            "harmonic {",
            "  colvars dummyAtoms",
            f"  centers ( {dummy_position_string} )",
            f"  forceConstant {kpos.magnitude:.2f}",
            "}",
        ]

        file.write("# Dummy atom position restraints\n")

        # Write colvar to file
        for line in cv_template_lines:
            file.write(line + "\n")

        # Write bias potential to file
        for line in bias_template_lines:
            file.write(line + "\n")
