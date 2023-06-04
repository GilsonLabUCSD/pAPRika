import logging
import os
import re
import subprocess

import numpy
import parmed

from paprika.build.system.utils import (
    ANGSTROM_CUBED_TO_LITERS,
    N_A,
    ConversionToolkit,
    PBCBox,
)

logger = logging.getLogger(__name__)

# TODO: refactor TLeap class, implement a build class for PSFGEN, PackMol and add TopoTools/VMD support


class TLeap(object):
    """
    Class for building AMBER prmtop/rst7 files with ``TLeap``.

    Parameters
    ----------
    add_ions : list, default=None
        A list of additional ions to be added to the system, specified by the residue name and
        an indication of the amount: an **integer number**, the **molarity** (M), or the **molality** (m).
        For example, the following list shows valid examples: add_ions = ['MG', 2, 'Cl-', 4, 'K', '0.150M',
        'BR', '0.150M', 'NA', '0.050m', 'F', '0.050m'].
    add_ion_residues : list, default=None
        A property formatted list of ions to add to the system.  The format is the same as ``add_ions``
        except that only integer amounts are allowed.  This is set by ``set_additional_ions``.
    buffer_value : float, default=12.0
        The desired solvation buffer value. This will add a layer of water around your solute with the
        minimum distance to the periodic box edge defined by ``buffer_value``. If ``target_waters`` is set,
        it will override this value.
    buffer_value_history : list, default=[0]
        The history of ``buffer_value`` adjustments stored in a list.
    counter_cation : str, default='Na+'
        If ``neutralize`` is set to ``True``, this cation will be used.
    counter_anion : str, default='Cl-'
        If ``neutralize`` is set to ``True``, this anion will be used.
    cycles_since_last_exp_change : int, default=0
        A count of the number of ``buffer_value`` adjustments since the last exponent change. Should
        always start at 0.
    exponent : int, default=`
        The initial value used for dynamically adjusting the buffer value.
    loadpdb_file : str, default=None
        If specified, this PDB file will be loaded via ``loadpdb`` instead of whatever file is specified
        in the template. This file should be present in ``output_path``, or if you specify it with a path,
        the path must be relative to ``output_path``.
    manual_switch_thresh : int, default=``target_waters``**(1./3.)
        The threshold difference between waters actually added and target_waters for which we can
        safely switch to manual removal of waters. If too large, then there will be big air pockets
        in our box. If too small, it will take forever to converge. If the user does not set this
        value, it will set proportionately to the target_waters.
    max_cycles : int, default=50
        The maximum number of ``buffer_value`` adjustment cycles.
    min_exponent_limit : int, default=-5
        The minimum limit that we let exponent get to. If it gets too small, it has no effect on
        the number of waters added.
    neutralize : bool, default=True
        Whether to neutralize the system with counterions.
    output_path : str, default='./'
        Directory path where ``TLeap`` input files will be written and executed. Associated input
        files (``frcmod``, ``mol2``, ``pdb``, etc) need to present in that path.
    output_prefix : str, default='build'
        Append a prefix to files created by this module.
    pbc_type : :class:`PBCBox`, default=PBCBox.cubic
        Type of solvation (`cubic`, `rectangular`, `octahedral`, or `None`).
    solvent_molecular_mass : float, default=0.018 (for water)
        `kg/mol` for the solvent. Used for computing the molality.
    target_waters : int, default=None
        The desired number of waters to solvate your system with. If specified, this will override any
        ``buffer_value`` settings.
    template_file : str, default=None
        Name of template file to read boilerplate ``TLeap`` input (e.g., `source leaprc`... lines). Any
        frcmod/mol2/pdb files which are loaded by the template must be present in ``output_path``.
    template_lines : list, default=None
        The list of ``TLeap`` commands which are needed to build the system. Either the ``template_file``
        or the ``template_lines`` needs to be set prior to running ``build()``. Files called for loading
        should be present in output_path (see above).
    unit : str, default='model'
        The ``TLeap`` unit name.
    waters_added_history : list, default=[0]
        The history of number of waters added which correlate to the values in ``buffer_value_history``.
    model : dict, default={"library": "tip3p", "frcmod": "tip3p", "water_box": "TIP3PBOX"}
        The type of water model to solvate the system. Use ``TLeap.getAmberWater()`` to specify the water model.
    waters_to_remove : list, default=None
        A list of the water residues to be manually removed after solvation.
    write_save_lines : bool, default=True
        Whether or not to do ``saveamberparm`` and ``savepdb`` when ``tleap`` is executed. During optimization
        loops it speeds up the process to set as ``False``, but it must be returned to ``True`` for the final
        execution.

    Example
    -------
        >>> from paprika.build.system import TLeap
        >>> from paprika.build.system.utils import PBCBox
        >>>
        >>> # Initialize TLeap object
        >>> system = TLeap()
        >>> system.output_path = "windows"
        >>> system.output_prefix = "host-guest-sol"
        >>> system.target_waters = 2000
        >>> system.pbc_type = PBCBox.cubic
        >>> system.set_water_model("tip3p", model_type="force-balance")
        >>> system.neutralize = True
        >>>
        >>> # template_lines lists all structure files etc.
        >>> system.template_lines = [
        >>>     "source leaprc.gaff",
        >>>     "loadamberparams host.frcmod",
        >>>     "loadamberparams guest.frcmod",
        >>>     "loadamberparams dummy.frcmod",
        >>>     "GST = loadmol2 guest.mol2",
        >>>     "HST = loadmol2 host.mol2",
        >>>     "DM1 = loadmol2 dm1.mol2",
        >>>     "DM2 = loadmol2 dm2.mol2",
        >>>     "DM3 = loadmol2 dm3.mol2",
        >>>     "model = loadpdb host-guest.pdb",
        >>> ]
        >>>
        >>> # Run TLeap
        >>> system.build()
        >>>
        >>> # Repartition masses of Hydrogen atoms
        >>> system.repartition_hydrogen_mass()
        >>>
        >>> # Convert AMBER files to GROMACS
        >>> system.convert_to_gromacs()

    """

    @property
    def add_ions(self):
        """list: A list of ions to add to the system."""
        return self._add_ions

    @add_ions.setter
    def add_ions(self, value):
        self._add_ions = value

    @property
    def add_ion_residues(self):
        """list: A property formatted list of ions to add to the system."""
        return self._add_ion_residues

    @add_ion_residues.setter
    def add_ion_residues(self, value):
        self._add_ion_residues = value

    @property
    def buffer_value(self):
        """float: The desired solvation buffer value."""
        return self._buffer_value

    @buffer_value.setter
    def buffer_value(self, value):
        self._buffer_value = value

    @property
    def buffer_value_history(self):
        """list: A list of the history of ``buffer_value`` adjustments."""
        return self._buffer_value_history

    @buffer_value_history.setter
    def buffer_value_history(self, value):
        self._buffer_value_history = value

    @property
    def counter_cation(self):
        """str: AMBER-style residue name for the cation."""
        return self._counter_cation

    @counter_cation.setter
    def counter_cation(self, value):
        self._counter_cation = value

    @property
    def counter_anion(self):
        """str: AMBER-style residue name for the anion."""
        return self._counter_anion

    @counter_anion.setter
    def counter_anion(self, value):
        self._counter_anion = value

    @property
    def cycles_since_last_exp_change(self):
        """int: Number of ``buffer_value`` adjustments since the last exponent change."""
        return self._cycles_since_last_exp_change

    @cycles_since_last_exp_change.setter
    def cycles_since_last_exp_change(self, value):
        self._cycles_since_last_exp_change = value

    @property
    def exponent(self):
        """int: The initial value used for adjusting the ``buffer_value``."""
        return self._exponent

    @exponent.setter
    def exponent(self, value):
        self._exponent = value

    @property
    def loadpdb_file(self):
        """os.PathLike: The PDB file to load using the ``TLeap`` command ``loadpdb``."""
        return self._loadpdb_file

    @loadpdb_file.setter
    def loadpdb_file(self, value):
        self._loadpdb_file = value

    @property
    def manual_switch_thresh(self):
        """int: The threshold difference between waters actually added and target_waters
        for which we can safely switch to manual removal of waters."""
        return self._manual_switch_thresh

    @manual_switch_thresh.setter
    def manual_switch_thresh(self, value):
        self._manual_switch_thresh = value

    @property
    def max_cycles(self):
        """int: The maximum number of ``buffer_value`` adjustment cycles."""
        return self._max_cycles

    @max_cycles.setter
    def max_cycles(self, value):
        self._max_cycles = value

    @property
    def min_exponent_limit(self):
        """int: The minimum limit that we let exponent get to."""
        return self._min_exponent_limit

    @min_exponent_limit.setter
    def min_exponent_limit(self, value):
        self._min_exponent_limit = value

    @property
    def neutralize(self):
        """str: AMBER-style residue name for the anion."""
        return self._neutralize

    @neutralize.setter
    def neutralize(self, value):
        self._neutralize = value

    @property
    def output_path(self):
        """os.PathLike: Directory path where ``TLeap`` input files will be written and executed.
        Associated input files (``frcmod``, ``mol2``, ``pdb``, etc) need to be present in that path.
        """
        return self._output_path

    @output_path.setter
    def output_path(self, value):
        self._output_path = value

    @property
    def output_prefix(self):
        """str: The name for the output files."""
        return self._output_prefix

    @output_prefix.setter
    def output_prefix(self, value):
        self._output_prefix = value

    @property
    def pbc_type(self):
        """:class:`PBCBox`: The periodic box type for the solvated system."""
        return self._pbc_type

    @pbc_type.setter
    def pbc_type(self, value):
        self._pbc_type = value

    @property
    def solvent_molecular_mass(self):
        """float: The molecular mass in kg/mol for the solvent molecule."""
        return self._solvent_molecular_mass

    @solvent_molecular_mass.setter
    def solvent_molecular_mass(self, value):
        self._solvent_molecular_mass = value

    @property
    def target_waters(self):
        """int: Number of water molecules to add to the system."""
        return self._target_waters

    @target_waters.setter
    def target_waters(self, value):
        self._target_waters = value

    @property
    def template_file(self):
        """os.PathLike: Specify a ``TLeap`` template file."""
        return self._template_file

    @template_file.setter
    def template_file(self, value):
        self._template_file = value

    @property
    def template_lines(self):
        """list: List of user-based commands as strings for ``TLeap``."""
        return self._template_lines

    @template_lines.setter
    def template_lines(self, value):
        self._template_lines = value

    @property
    def unit(self):
        """str: The `unit` name for the structure in ``TLeap``."""
        return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value

    @property
    def waters_added_history(self):
        """list: The history of number of waters added."""
        return self._waters_added_history

    @waters_added_history.setter
    def waters_added_history(self, value):
        self._waters_added_history = value

    @property
    def water_model(self):
        """dict: The water model for ``TLeap`` to use stored as a dictionary."""
        return self._water_model

    @water_model.setter
    def water_model(self, value):
        self._water_model = value

    @property
    def waters_to_remove(self):
        """list: A list of water molecules to remove after solvation."""
        return self._waters_to_remove

    @waters_to_remove.setter
    def waters_to_remove(self, value):
        self._waters_to_remove = value

    @property
    def write_save_lines(self):
        """bool: Whether to run ``saveamberparm`` and ``savepdb`` when ``TLeap`` is executed."""
        return self._write_save_lines

    @write_save_lines.setter
    def write_save_lines(self, value):
        self._write_save_lines = value

    def __init__(self):

        # User Settings: Defaults
        self._template_file = None
        self._template_lines = None
        self._loadpdb_file = None
        self._pbc_type = PBCBox.cubic
        self._buffer_value = 12.0
        self._target_waters = None
        self._water_model = {
            "library": "tip3p",
            "frcmod": "tip3p",
            "water_box": "TIP3PBOX",
        }
        self._neutralize = True
        self._counter_cation = "Na+"
        self._counter_anion = "Cl-"
        self._add_ions = None
        self._output_path = "../"
        self._output_prefix = "build"

        # Advanced Settings: Defaults
        self._unit = "model"
        self._exponent = 1
        self._min_exponent_limit = -5
        self._cycles_since_last_exp_change = 0
        self._max_cycles = 50
        self._manual_switch_thresh = None
        self._waters_to_remove = None
        self._add_ion_residues = None
        self._solvent_molecular_mass = 0.018
        self._buffer_value_history = [0]
        self._waters_added_history = [0]
        self._write_save_lines = True

    def build(self, clean_files: bool = True):
        """
        Build the ``TLeap`` system.

        Parameters
        ----------
        clean_files : bool, optional
            Whether to delete log files after completion.
        """

        logger.debug("Running tleap.build() in {}".format(self.output_path))

        # Check input
        if self.template_file and self.template_lines:
            raise Exception("template_file and template_lines cannot both be specified")
        elif self.template_file:
            with open(self.template_file, "r") as f:
                self.template_lines = f.read().splitlines()
        elif self.template_lines:
            for i, line in enumerate(self.template_lines):
                self.template_lines[i] = line.rstrip()
        else:
            raise Exception(
                "Either template_file or template_lines needs to be specified"
            )

        # Filter out any interfering lines
        self.filter_template()

        # Either just write/run tleap, or do solvate
        if self.pbc_type is None:
            self.write_input()
            self.run()
        else:
            self.solvate()

        # Cleanup
        if clean_files:
            os.remove(os.path.join(self.output_path, self.output_prefix + ".tleap.in"))
            os.remove(os.path.join(self.output_path, "leap.log"))

            if self.water_model["frcmod"] == "bind3p":
                os.remove(os.path.join(self.output_path, "frcmod.bind3p"))

    def filter_template(self):
        """
        Filter out any ``template_lines`` that may interfere with solvation.
        """

        filtered_lines = []

        number_of_pdb_files = len(re.findall("loadpdb", "".join(self.template_lines)))
        for line in self.template_lines:
            # Find loadpdb line, replace pdb file if necessary, set unit name
            if re.search("loadpdb", line):
                words = line.rstrip().replace("=", " ").split()
                self.unit = words[0]
                filtered_lines.append(
                    "{} = loadpdb {}\n".format(
                        self.unit,
                        self.loadpdb_file
                        if self.loadpdb_file and number_of_pdb_files == 1
                        else words[2],
                    )
                )
            elif re.search("combine", line):
                words = line.rstrip().replace("=", " ").split()
                self.unit = words[0]
                logger.debug(
                    f"Found `combine` keyword and reassigning `self.unit` to {self.unit}..."
                )
                filtered_lines.append(line)
            # Remove any included solvation and ionization commands if pbc_type
            # is not None
            elif self.pbc_type is not None:
                if not re.search(
                    r"^\s*(addions|addions2|addionsrand|desc|quit|solvate|save)",
                    line,
                    re.IGNORECASE,
                ):
                    filtered_lines.append(line)
            else:
                filtered_lines.append(line)

        self.template_lines = filtered_lines

        # Add the water leaprc library to template_lines
        if self.pbc_type and self.water_model:
            # We need to source the water model after all other modules are loaded
            # because some modules will overwrite the water model definition.
            # Example case: source leaprc.protein.ff14SB contains tip3p parameters so
            # we need to load this first before the water frcmod file.
            insert_index = 0
            for idx, line in enumerate(self.template_lines):
                if "source" in line:
                    insert_index = idx + 1
            self.template_lines.insert(
                insert_index, f"loadamberparams frcmod.{self.water_model['frcmod']}"
            )
            self.template_lines.insert(
                insert_index, f"source leaprc.water.{self.water_model['library']}"
            )

    def write_input(self):
        """
        Write a ``TLeap`` input file based on ``template_lines`` and other things we have set.
        """

        file_path = os.path.join(self.output_path, self.output_prefix + ".tleap.in")
        if not os.path.exists(os.path.dirname(file_path)):
            try:
                os.makedirs(os.path.dirname(file_path))
            except OSError:
                raise

        # If Bind3P, write frcmod.bind3p file
        if self.water_model["frcmod"] == "bind3p":
            frcmod_path = os.path.join(self.output_path, "frcmod.bind3p")
            with open(frcmod_path, "w") as f:
                f.write(
                    """\
This is the additional/replacement parameter set for Bind3P water
MASS
OW    16.0
HW     1.008   0.000

BOND
OW-HW  553.0    0.9572      Bind3P water
HW-HW  553.0    1.5136      Bind3P water

ANGLE

DIHE

NONBON
  OW          1.7577  0.1818             Bind3P water model
  HW          0.0000  0.0000             Bind3P water model
"""
                )

        with open(file_path, "w") as f:
            for line in self.template_lines:
                f.write(line + "\n")

            if self.pbc_type == PBCBox.cubic:
                f.write(
                    "solvatebox {} {} {} iso\n".format(
                        self.unit, self.water_model["water_box"], self.buffer_value
                    )
                )
            elif self.pbc_type == PBCBox.rectangular:
                f.write(
                    "solvatebox {} {} {{10.0 10.0 {}}}\n".format(
                        self.unit, self.water_model["water_box"], self.buffer_value
                    )
                )
            elif self.pbc_type == PBCBox.octahedral:
                f.write(
                    "solvateoct {} {} {} iso\n".format(
                        self.unit, self.water_model["water_box"], self.buffer_value
                    )
                )
            elif self.pbc_type is None:
                f.write("# Skipping solvation ...\n")
            else:
                raise Exception(
                    "Incorrect pbc_type value provided: "
                    + str(self.pbc_type)
                    + ". Only `cubic`, `rectangular`, `octahedral`, and None are valid"
                )
            if self.neutralize:
                f.write("addionsrand {} {} 0\n".format(self.unit, self.counter_cation))
                f.write("addionsrand {} {} 0\n".format(self.unit, self.counter_anion))
            # Additional ions should be specified as a list, with residue name and number of ions in pairs, like ['NA',
            # 5] for five additional sodium ions. By this point, if the user specified a molality or molarity,
            # it should already have been converted into a number.
            if self.add_ion_residues:
                for residue, amount in zip(
                    self.add_ion_residues[0::2], self.add_ion_residues[1::2]
                ):
                    f.write("addionsrand {} {} {}\n".format(self.unit, residue, amount))
            if self.waters_to_remove:
                for water_number in self.waters_to_remove:
                    f.write(
                        "remove {} {}.{}\n".format(self.unit, self.unit, water_number)
                    )

            # Note, the execution of tleap is assumed to take place in the
            # same directory as all the associated input files, so we won't
            # put directory paths on the saveamberparm or savepdb commands.
            if self.output_prefix and self.write_save_lines:
                f.write("savepdb {} {}.pdb\n".format(self.unit, self.output_prefix))
                f.write(
                    "saveamberparm {} {}.prmtop {}.rst7\n".format(
                        self.unit, self.output_prefix, self.output_prefix
                    )
                )
            else:
                pass
            f.write("desc {}\n".format(self.unit))
            f.write("quit\n")

    def run(self):
        """
        Execute ``TLeap``.

        Returns
        -------
        output : list
            The tleap stdout returned as a list.
        """

        self.check_for_leap_log()
        file_name = self.output_prefix + ".tleap.in"

        output = subprocess.Popen(
            ["tleap", "-s ", "-f ", file_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=self.output_path,
        )
        output = output.stdout.read().decode().splitlines()

        self.grep_leap_log()
        return output

    def grep_leap_log(self):
        """
        Check for a few keywords in the ``TLeap`` output.
        """
        try:
            with open(self.output_path + "leap.log", "r") as file:
                for line in file.readlines():
                    if re.search(
                        "ERROR|WARNING|Warning|duplicate|FATAL|Could|Fatal|Error", line
                    ):
                        logger.warning(
                            "It appears there was a problem with solvation: check `leap.log`..."
                        )
        except BaseException:
            return

    def check_for_leap_log(self, log_file="leap.log"):
        """
        Check if `leap.log` exists in ``output_path``, and if so, delete so the current run doesn't append.

        Parameters
        ----------
        log_file : str, optional
            Name of the tleap logfile. Default: leap.log
        """
        log_file_path = os.path.join(self.output_path, log_file)
        try:
            os.remove(log_file_path)
            logger.debug("Deleted existing leap logfile: " + log_file_path)
        except OSError:
            pass

    def solvate(self):
        """
        Solvate a structure with an exact number of waters or buffer size.
        """

        # If buffer_value is set but not target_waters, figure out what
        # target_waters value corresponds to buffer_value. We do this because
        # the solvate() algorithm is focused on getting target_waters correct.
        if self.target_waters is None:
            self.target_waters = self.count_waters()
        # If target_waters is set, we assume that overrides buffer_value. We'll
        # set buffer_value = 2.0 because that seems to lead to smooth convergence
        # of the solvate() routine. Setting this to 2.0 prevents issues when
        # solvating a single ion.
        else:
            self.buffer_value = 2.0

        # If the user has not set manual_switch_thresh, we will set it proportionately
        # to the target_waters. This will control when we can start manually deleting
        # waters rather than adjusting the buffer_value.
        if self.manual_switch_thresh is None:
            self.manual_switch_thresh = int(numpy.ceil(self.target_waters ** (1.0 / 3.0)))
            if self.manual_switch_thresh < 12:
                self.manual_switch_thresh = 12
            logger.debug(
                "manual_switch_thresh is set to: {:.0f}".format(
                    self.manual_switch_thresh
                )
            )

        if self.add_ions:
            self.set_additional_ions()

        # Speed up the initial optimization loop by not writing prmtops and
        # pdbs
        self.write_save_lines = False

        # First, a coarse adjustment...
        # This will run for 50 iterations or until we (a) have more waters than the target and (b) are within ~12 waters
        # of the target (that can be manually removed).
        cycle = 0
        while cycle < self.max_cycles:
            # Start by not manually removing any water, just adjusting the buffer value to get near the target number of
            # waters...

            # Find out how many waters for *this* buffer value...
            waters = self.count_waters()
            self.waters_added_history.append(waters)
            # noinspection PyTypeChecker
            self.buffer_value_history.append(self.buffer_value)
            logger.debug(
                "Cycle {:02.0f} {:.0f} {:10.7f} {:6.0f} ({:6.0f})".format(
                    cycle, self.exponent, self.buffer_value, waters, self.target_waters
                )
            )

            # Possible location of a switch to adjust the buffer_value by polynomial
            # fit approach.

            # If we've nailed it, break!
            if waters == self.target_waters:
                # Run one more time and save files
                self.final_solvation_run()
                return
            # If we are close, go to fine adjustment...
            elif (
                waters > self.target_waters
                and (waters - self.target_waters) < self.manual_switch_thresh
            ):
                self.remove_waters_manually()
                # Run once more and save files
                self.final_solvation_run()
                return
            # Otherwise, try to keep adjusting the number of waters...
            else:
                self.adjust_buffer_value()
                # Now that we're close, let's re-evaluate how many ions to add, in case the volume has changed a lot.
                # (This could be slow and run less frequently...)
                if self.add_ions and cycle % 10 == 0:
                    self.set_additional_ions()
                cycle += 1

        if cycle >= self.max_cycles and waters > self.target_waters:
            logger.debug(
                "The added waters ({}) didn't reach the manual_switch_thresh ({}) with max_cycles ({}), "
                "but we'll try manual removal anyway.".format(
                    waters, self.target_waters, self.max_cycles
                )
            )
            self.remove_waters_manually()
            self.final_solvation_run()

        if cycle >= self.max_cycles and waters < self.target_waters:
            raise Exception(
                "Automatic adjustment of the buffer value resulted in fewer waters \
                added than targeted by `buffer_water`. Try increasing manual_switch_thresh"
            )
        else:
            raise Exception(
                "Automatic adjustment of the buffer value was unable to converge on \
                a solution to within the specified manual_switch_thresh: {:.0f}".format(
                    self.manual_switch_thresh
                )
            )

    def final_solvation_run(self):
        """
        Run the solvation with ``write_save_lines=True``.
        """

        self.write_save_lines = True

        # We'll put this in a loop in case addionsrand does something different than the
        # last round, which indicated success. We only check waters, because addionsrand
        # can replace water molecules sometimes, which messes up the count.
        for i in range(50):
            residues = self.count_residues(print_results=True)
            if residues["WAT"] == self.target_waters:
                break
            if i == 49:
                raise Exception(
                    "Unable to add the correct waters at 50 cycles during final_solvation_run()"
                )
            logger.info(
                "The final solvation step added the wrong number of waters. Repeating ..."
            )

    def count_waters(self):
        """
        Quickly check the number of waters added for a given buffer size.

        Returns
        -------
        waters : int
            Number of water molecules.
        """
        waters = self.count_residues()["WAT"]
        return waters

    def count_residues(self, print_results=False):
        """
        Run and parse ``TLeap`` output and return a dictionary of residues in the structure.

        Parameters
        ----------
        print_results: bool, optional
            Whether to print residues to log file.

        Returns
        -------
        residues : dict
            Dictionary of added residues and their number
        """
        for attempt in range(10):
            self.write_input()
            output = self.run()
            # Return a dictionary of {'RES' : number of RES}
            residues = {}
            for line in output:
                # Is this line a residue from `desc` command?
                match = re.search("^R<(.*) ", line)
                if match:
                    residue_name = match.group(1)
                    # If this residue is not in the dictionary, initialize and
                    # set the count to 1.
                    if residue_name not in residues:
                        residues[residue_name] = 1
                    # If this residue is in the dictionary, increment the count
                    # each time we find an instance.
                    elif residue_name in residues:
                        residues[residue_name] += 1
            if residues:
                break
            if attempt == 9:
                raise Exception(
                    "tleap was unable to successfully create the system after 10 attempts."
                    + " Investigate the leap.log for errors."
                )

        # logger.debug(residues)
        if print_results:
            for key, value in sorted(residues.items()):
                logger.info("{:10s} {:10.0f}".format(key, value))

        return residues

    def set_additional_ions(self):
        """
        Determine whether additional ions (instead of or in addition to neutralization) are requested...

        **Sets:**

        ``self.add_ion_residues`` : list
            A processed list of ions and their amounts that can be passed to ``TLeap``.

        """
        if not self.add_ions:
            return None
        if len(self.add_ions) < 2:
            raise Exception("No amount specified for additional ions.")
        if len(self.add_ions) % 2 == 1:
            raise Exception(
                "The 'add_ions' list requires an even number of elements. "
                "Make sure there is a residue mask followed by a value for "
                "each ion to be added (or molarity ending in 'M' or molality ending in 'm')."
            )
        self.add_ion_residues = []
        for ion, amount in zip(self.add_ions[0::2], self.add_ions[1::2]):
            self.add_ion_residues.append(ion)
            if isinstance(amount, int):
                self.add_ion_residues.append(amount)
            elif isinstance(amount, str) and amount[-1] == "m":
                # User specifies molality...
                # number to add = (molality) x (number waters) x (kg/mol
                # solvent)
                number_to_add = int(
                    numpy.ceil(
                        float(amount[:-1])
                        * self.target_waters
                        * self.solvent_molecular_mass
                    )
                )
                self.add_ion_residues.append(number_to_add)
            elif isinstance(amount, str) and amount[-1] == "M":
                # User specifies molarity...
                volume = self.get_volume()
                if volume is None:
                    raise Exception(
                        "The volume of the system could not be found and thus "
                        "the correct ion count could not be determined."
                    )
                number_of_atoms = float(amount[:-1]) * N_A
                liters = volume * ANGSTROM_CUBED_TO_LITERS
                number_to_add = int(numpy.ceil(number_of_atoms * liters))
                self.add_ion_residues.append(number_to_add)
            else:
                raise Exception("Unanticipated error calculating how many ions to add.")

    def get_volume(self):
        """
        Run and parse ``TLeap`` output and return the volume of the structure.

        Returns
        -------
        volume : float
            The volume of the structure in cubic angstroms.

        """
        output = self.run()
        # Return the total simulation volume
        for line in output:
            line = line.strip()
            if "Volume" in line:
                match = re.search("Volume(.*)", line)
                volume = float(match.group(1)[1:-4])
                return volume
        logger.warning("Could not determine total simulation volume.")
        return None

    def remove_waters_manually(self):
        """
        Remove a few water molecules manually with ``TLeap`` to exactly match a desired number of waters.
        """

        cycle = 0
        max_cycles = 100
        waters = self.waters_added_history[-1]
        while waters > self.target_waters:
            # Retrieve excess water residues
            water_surplus = waters - self.target_waters
            water_residues = self.list_waters()

            # THIS IS HACKY. But if we've gone > 5 cycles, we're probably
            # in a ping-pong convergence problem. So we'll try to solve it
            # by increasing the water_surplus value, manually, and increase
            # it as we go through more cycles.
            if cycle > 5:
                additional_water = int(float(cycle) / 5.0)
                water_surplus += additional_water
                logger.debug(
                    "Detected trouble with manually removing water. Increasing the number"
                    "of surplus waters by {}".format(additional_water)
                )

            self.waters_to_remove = water_residues[-1 * water_surplus :]
            logger.debug("Manually removing waters... {}".format(self.waters_to_remove))

            # Get counts for all residues
            residues = self.count_residues()

            # Check if we reached target
            waters = residues["WAT"]
            if waters == self.target_waters:
                return
            cycle += 1
            if cycle > max_cycles:
                raise Exception(
                    "Solvation failed due to an unanticipated problem with water removal."
                )

    def list_waters(self):
        """
        Run and parse ``TLeap`` output and return the a list of water residues.

        Returns
        -------
        water_residues : list
            A list of the water residues in the structure.

        """
        output = self.run()

        # Return a list of residue numbers for the waters
        water_residues = []
        for line in output:
            # Is this line a water?
            match = re.search("^R<WAT (.*)>", line)
            if match:
                water_residues.append(match.group(1))
        return water_residues

    def adjust_buffer_value(self):
        """
        Determine whether to increase or decrease the buffer thickness to match a desired number of waters.

        **Sets:**

        ``self.buffer_value`` : float
            A new buffer size to try.
        ``self.exponent`` : int
            Adjusts the order of magnitude of buffer value changes.
        ``self.cycles_since_last_exp_change`` : int
            Resets this value to 0 when exponent value is changed to help with the logic gates.

        """

        # If the number of waters was less than the target and is now greater than the target, make the buffer smaller
        # smaller
        if (
            self.waters_added_history[-2]
            < self.target_waters
            < self.waters_added_history[-1]
        ):
            # If its been more than one round since last exponent change,
            # change exponent
            if self.cycles_since_last_exp_change > 1:
                logger.debug("Adjustment loop 1a")
                self.exponent -= 1
                self.buffer_value = self.buffer_value_history[-1] + -5 * (
                    10 ** self.exponent
                )
                self.cycles_since_last_exp_change = 0
            else:
                logger.debug("Adjustment loop 1b")
                self.buffer_value = self.buffer_value_history[-1] + -1 * (
                    10 ** self.exponent
                )
                self.cycles_since_last_exp_change += 1
        # If the number of waters was greater than the target and is now less
        # than the target, make the buffer bigger
        elif (
            self.waters_added_history[-2]
            > self.target_waters
            > self.waters_added_history[-1]
        ):
            # If its been more than one round since last exponent change,
            # change exponent
            if self.cycles_since_last_exp_change > 1:
                logger.debug("Adjustment loop 2a")
                self.exponent -= 1
                self.buffer_value = self.buffer_value_history[-1] + 5 * (
                    10 ** self.exponent
                )
                self.cycles_since_last_exp_change = 0
            else:
                logger.debug("Adjustment loop 2b")
                self.buffer_value = self.buffer_value_history[-1] + 1 * (
                    10 ** self.exponent
                )
                self.cycles_since_last_exp_change += 1
        # If the last two rounds of solvation have too many waters, make the
        # buffer smaller...
        elif (
            self.waters_added_history[-2] > self.target_waters
            and self.waters_added_history[-1] > self.target_waters
        ):
            logger.debug("Adjustment loop 3")
            self.buffer_value = self.buffer_value_history[-1] + -1 * (
                10 ** self.exponent
            )
            self.cycles_since_last_exp_change += 1
        # If the last two rounds of solvation had too few waters, make the
        # buffer bigger...
        elif (
            self.waters_added_history[-2] < self.target_waters
            and self.waters_added_history[-1] < self.target_waters
        ):
            logger.debug("Adjustment loop 4")
            self.buffer_value = self.buffer_value_history[-1] + 1 * (
                10 ** self.exponent
            )
            self.cycles_since_last_exp_change += 1
        else:
            raise Exception(
                "The buffer_values search died due to an unanticipated set of variable values"
            )

        if self.exponent <= self.min_exponent_limit:
            raise Exception(
                "Automatic adjustment of the buffer value failed to get near enough to target_waters "
                "before the exponent value was below {:.0f}. Try increasing manual_switch_thresh.".format(
                    self.exponent
                )
            )

    def repartition_hydrogen_mass(self, options=None):
        """
        Repartitions the masses of Hydrogen atoms in the system by a factor 3 and
        overwrites the `.prmtop` file: `output_path/output_prefix.prmtop`

        Parameters
        ----------
        options : str, optional
            Optional keyword(s) for the repartitioning and following the
            :class:`parmed.tools.actions` documentation the usage is '[<mass>] [dowater]'.

        """

        prmtop = os.path.join(self.output_path, self.output_prefix + ".prmtop")

        structure = parmed.load_file(prmtop, structure=True)

        # noinspection PyTypeChecker
        parmed.tools.actions.HMassRepartition(structure, arg_list=options).execute()

        structure.save(prmtop, overwrite=True)

    def convert_to_gromacs(
        self,
        overwrite=False,
        output_path=None,
        output_prefix=None,
        toolkit=ConversionToolkit.InterMol,
    ):
        """
        Convert AMBER topology and coordinate files to GROMACS format.

        Parameters
        ----------
        overwrite: bool, optional, default=False
            Option to overwrite GROMACS ``.top`` and ``.gro`` files if they already
            exists in the folder.
        output_path: str, optional, default=None
            Alternate directory path where the AMBER files are located. Default is the
            `path` parsed to the :class:`TLeap` object.
        output_prefix: str, optional, default=None
            Alternate file name prefix for the Amber files. Default is the `prefix` parsed
            to the :class:`TLeap` object.
        toolkit: :class:`ConversionToolkit`, default=ConversionToolkit.InterMol
            Option to choose the toolkit for converting the AMBER files, ParmEd or InterMol
        """

        if output_path is None:
            output_path = self.output_path

        if output_prefix is None:
            output_prefix = self.output_prefix

        file_name = os.path.join(output_path, output_prefix)

        # Check if Amber Topology file(s) exist
        topology = numpy.array([f"{file_name}.{ext}" for ext in ["prmtop", "parm7"]])
        check_topology = [os.path.isfile(file) for file in topology]
        if not any(check_topology):
            raise FileNotFoundError("Cannot find any AMBER topology file.")

        # Get the first topology file in the list of file(s) that exists
        prmtop = topology[check_topology][0]

        # Check if Amber Coordinate file(s) exist
        if toolkit == ConversionToolkit.ParmEd:
            coordinates = numpy.array(
                [
                    f"{file_name}.{ext}"
                    for ext in ["rst7", "inpcrd", "rst", "crd", "pdb"]
                ]
            )
        else:  # InterMol does not support PDB for loading Amber files
            coordinates = numpy.array(
                [f"{file_name}.{ext}" for ext in ["rst7", "inpcrd", "rst", "crd"]]
            )
        check_coordinates = [os.path.isfile(file) for file in coordinates]
        if not any(check_coordinates):
            raise FileNotFoundError("Cannot find any AMBER coordinates file.")

        # Get the first coordinate file in the list of file(s) that exists
        inpcrd = coordinates[check_coordinates][0]

        # Gromacs file names
        top_file = f"{file_name}.top"
        gro_file = f"{file_name}.gro"

        # Convert with ParmEd
        if toolkit == ConversionToolkit.ParmEd:
            # Load Amber files
            structure = parmed.load_file(prmtop, inpcrd, structure=True)

            if overwrite:
                structure.save(top_file, format="gromacs", overwrite=True)
                structure.save(gro_file, format="gro", overwrite=True)
            else:
                if not os.path.isfile(top_file):
                    structure.save(top_file, format="gromacs")
                else:
                    logger.info(f"Topology file {top_file} exists, skipping writing file.")

                if not os.path.isfile(gro_file):
                    structure.save(gro_file, format="gro")
                else:
                    logger.info(
                        f"Coordinates file {gro_file} exists, skipping writing file."
                    )

        # Convert with InterMol
        elif toolkit == ConversionToolkit.InterMol:
            from intermol.convert import _load_amber, _save_gromacs

            # Load Amber files
            system, prefix, prmtop_in, crd_in, amb_structure = _load_amber(
                [prmtop, inpcrd]
            )

            # Save Gromacs files
            output_status = dict()
            _save_gromacs(system, file_name, output_status)

            if output_status["gromacs"] != "Converted":
                raise Exception(
                    f"Converting AMBER to GROMACS with Intermol unsuccessful: {output_status['gromacs']}"
                )

    def convert_to_charmm(
        self,
        overwrite=False,
        output_path=None,
        output_prefix=None,
        toolkit=ConversionToolkit.InterMol,
    ):
        """
        Convert AMBER topology and coordinate files to CHARMM format.

        Parameters
        ----------
        overwrite: bool, optional, default=False
            Option to overwrite CHARMM ``.psf`` and ``.pdb`` files if they already
            exists in the folder.
        output_path: str, optional, default=None
            Alternate directory path where the AMBER files are located. Default is the
            `path` parsed to the :class:`TLeap` object.
        output_prefix: str, optional, default=None
            Alternate file name prefix for the Amber files. Default is the `prefix` parsed
            to the :class:`TLeap` object.
        toolkit: :class:`ConversionToolkit`, default=ConversionToolkit.InterMol
            Option to choose the toolkit for converting the AMBER files, ParmEd or InterMol
        """

        if output_path is None:
            output_path = self.output_path

        if output_prefix is None:
            output_prefix = self.output_prefix

        file_name = os.path.join(output_path, output_prefix)

        # Check if Amber Topology file(s) exist
        topology = numpy.array([f"{file_name}.{ext}" for ext in ["prmtop", "parm7"]])
        check_topology = [os.path.isfile(file) for file in topology]
        if not any(check_topology):
            raise FileNotFoundError("Cannot find any AMBER topology file.")

        # Get the first topology file in the list of file(s) that exists
        prmtop = topology[check_topology][0]

        # Check if Amber Coordinate file(s) exist
        if toolkit == ConversionToolkit.ParmEd:
            coordinates = numpy.array(
                [
                    f"{file_name}.{ext}"
                    for ext in ["rst7", "inpcrd", "rst", "crd", "pdb"]
                ]
            )
        else:  # InterMol does not support PDB for loading Amber files
            coordinates = numpy.array(
                [f"{file_name}.{ext}" for ext in ["rst7", "inpcrd", "rst", "crd"]]
            )
        check_coordinates = [os.path.isfile(file) for file in coordinates]
        if not any(check_coordinates):
            raise FileNotFoundError("Cannot find any AMBER coordinates file.")

        # Get the first coordinate file in the list of file(s) that exists
        inpcrd = coordinates[check_coordinates][0]

        # Charmm file names
        psf_file = f"{file_name}.psf"
        crd_file = f"{file_name}.crd"

        # Convert with ParmEd
        if toolkit == ConversionToolkit.ParmEd:
            # Load Amber files
            structure = parmed.load_file(prmtop, inpcrd, structure=True)

            if overwrite:
                structure.save(psf_file, format="psf", overwrite=True)
                structure.save(crd_file, format="charmmcrd", overwrite=True)
            else:
                if not os.path.isfile(psf_file):
                    structure.save(psf_file, format="psf")
                else:
                    logger.info(f"Topology file {psf_file} exists, skipping writing file.")

                if not os.path.isfile(crd_file):
                    structure.save(crd_file)
                else:
                    logger.info(
                        f"Coordinates file {crd_file} exists, skipping writing file."
                    )

        # Convert with InterMol
        elif toolkit == ConversionToolkit.InterMol:
            from intermol.convert import _load_amber, _save_charmm

            # Load Amber files
            system, prefix, prmtop_in, crd_in, amb_structure = _load_amber(
                [prmtop, inpcrd]
            )

            # Save Charmm files
            output_status = dict()
            _, _ = _save_charmm(amb_structure, file_name, output_status)

            if output_status["charmm"] != "Converted":
                raise Exception(
                    f"Converting AMBER to CHARMM with Intermol unsuccessful: {output_status['charmm']}"
                )

    def convert_to_lammps(
        self,
        pair_style=None,
        output_path=None,
        output_prefix=None,
    ):
        """
        Convert AMBER topology and coordinate files to LAMMPS format.

        Parameters
        ----------
        pair_style: str, optional, default=None
            LAMMPS option for ``pair_style``, if set to ``None`` (default) the setting will be for a periodic Ewald
            simulation.
        output_path: str, optional, default=None
            Alternate directory path where the AMBER files are located. Default is the
            `path` parsed to the :class:`TLeap` object.
        output_prefix: str, optional, default=None
            Alternate file name prefix for the Amber files. Default is the `prefix` parsed
            to the :class:`TLeap` object.
        """
        if output_path is None:
            output_path = self.output_path

        if output_prefix is None:
            output_prefix = self.output_prefix

        file_name = os.path.join(output_path, output_prefix)

        # Check if Amber Topology file(s) exist
        topology = numpy.array([f"{file_name}.{ext}" for ext in ["prmtop", "parm7"]])
        check_topology = [os.path.isfile(file) for file in topology]
        if not any(check_topology):
            raise FileNotFoundError("Cannot find any AMBER topology file.")

        # Get the first topology file in the list of file(s) that exists
        prmtop = topology[check_topology][0]

        # Check if Amber Coordinate file(s) exist
        coordinates = numpy.array(
            [f"{file_name}.{ext}" for ext in ["rst7", "inpcrd", "rst", "crd"]]
        )
        check_coordinates = [os.path.isfile(file) for file in coordinates]
        if not any(check_coordinates):
            raise FileNotFoundError("Cannot find any AMBER coordinates file.")

        # Get the first coordinate file in the list of file(s) that exists
        inpcrd = coordinates[check_coordinates][0]

        # Pair style settings
        args = dict()
        if pair_style is None:
            args[
                "lmp_settings"
            ] = "pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1e-8"
        else:
            args["lmp_settings"] = pair_style

        # Convert
        from intermol.convert import _load_amber, _save_lammps

        # Load Amber files
        system, prefix, prmtop_in, crd_in, amb_structure = _load_amber([prmtop, inpcrd])

        # Save Lammps files
        output_status = dict()
        _save_lammps(system, file_name, output_status, args)

        if output_status["lammps"] != "Converted":
            raise Exception(
                f"Converting AMBER to LAMMPS with Intermol unsuccessful: {output_status['lammps']}"
            )

    def convert_to_desmond(
        self,
        output_path=None,
        output_prefix=None,
    ):
        """
        Convert AMBER topology and coordinate files to DESMOND format.

        Parameters
        ----------
        output_path: str, optional, default=None
            Alternate directory path where the AMBER files are located. Default is the
            `path` parsed to the :class:`TLeap` object.
        output_prefix: str, optional, default=None
            Alternate file name prefix for the Amber files. Default is the `prefix` parsed
            to the :class:`TLeap` object.
        """
        if output_path is None:
            output_path = self.output_path

        if output_prefix is None:
            output_prefix = self.output_prefix

        file_name = os.path.join(output_path, output_prefix)

        # Check if Amber Topology file(s) exist
        topology = numpy.array([f"{file_name}.{ext}" for ext in ["prmtop", "parm7"]])
        check_topology = [os.path.isfile(file) for file in topology]
        if not any(check_topology):
            raise FileNotFoundError("Cannot find any AMBER topology file.")

        # Get the first topology file in the list of file(s) that exists
        prmtop = topology[check_topology][0]

        # Check if Amber Coordinate file(s) exist
        coordinates = numpy.array(
            [f"{file_name}.{ext}" for ext in ["rst7", "inpcrd", "rst", "crd"]]
        )
        check_coordinates = [os.path.isfile(file) for file in coordinates]
        if not any(check_coordinates):
            raise FileNotFoundError("Cannot find any AMBER coordinates file.")

        # Get the first coordinate file in the list of file(s) that exists
        inpcrd = coordinates[check_coordinates][0]

        # Convert
        from intermol.convert import _load_amber, _save_desmond

        # Load Amber files
        system, prefix, prmtop_in, crd_in, amb_structure = _load_amber([prmtop, inpcrd])

        # Save Desmond files
        output_status = dict()
        _save_desmond(system, file_name, output_status)

        if output_status["desmond"] != "Converted":
            raise Exception(
                f"Converting AMBER to DESMOND with Intermol unsuccessful: {output_status['desmond']}"
            )

    def set_water_model(self, model, model_type=None):
        """
        Specify the water model for ``TLeap`` to use. The water model is stored as a dictionary in
        the variable ``water_model``.

        Parameters
        ----------
        model: str
            The water model to use, models supported are ["spc", "opc", "tip3p", "tip4p"].
            Strings are case-insensitive.
        model_type: str
            The particular type of the water model, default is ``None`` and the key in parenthesis is the water box
            used in ``TLeap``. Strings are case-insensitive.
            * spc: None (SPCBOX), "flexible" (SPCFWBOX), "quantum" (QSPCFWBOX)
            * opc: None (OPCBOX), "three-point" (OPC3BOX)
            * tip3p: None (TIP3PBOX), "flexible" (TIP3PFBOX), "force-balance" (FB3BOX), "bind3p" (TIP3PBOX)
            * tip4p: None (TIP4PBOX), "ewald" (TIP4PEWBOX), "force-balance" (FB4BOX)
        """
        if model.lower() not in ["spc", "opc", "tip3p", "tip4p"]:
            raise KeyError(f"Water model {model} is not supported.")

        library = None
        frcmod = None
        water_box = None

        if model.lower() == "spc":
            library = "spce"
            frcmod = "spce"
            if model_type is None:
                water_box = "SPCBOX"
            elif model_type.lower() == "flexible":
                water_box = "SPCFWBOX"
            elif model_type.lower() == "quantum":
                water_box = "QSPCFWBOX"
            else:
                raise KeyError(
                    f"Water type {model_type} is not supported for SPC water model."
                )

        if model.lower() == "opc":
            library = "opc"
            if model_type is None:
                frcmod = "opc"
                water_box = "OPCBOX"
            elif model_type.lower() == "three-point":
                frcmod = "opc3"
                water_box = "OPC3BOX"
            else:
                raise KeyError(
                    f"Water type {model_type} is not supported for OPC water model."
                )

        if model.lower() == "tip3p":
            library = "tip3p"
            water_box = "TIP3PBOX"
            if model_type is None:
                frcmod = "tip3p"
            elif model_type.lower() == "bind3p":
                frcmod = "bind3p"
            elif model_type.lower() == "force-balance":
                frcmod = "tip3pfb"
                water_box = "FB3BOX"
            elif model_type.lower() == "flexible":
                frcmod = "tip3p"
                water_box = "TIP3PFBOX"
            else:
                raise KeyError(
                    f"Water type {model_type} is not supported for TIP3P water model."
                )

        if model.lower() == "tip4p":
            library = "tip4pew"
            if model_type is None:
                frcmod = "tip4p"
                water_box = "TIP4PBOX"
            elif model_type.lower() == "force-balance":
                frcmod = "tip4pfb"
                water_box = "FB4BOX"
            elif model_type.lower() == "ewald":
                frcmod = "tip4pew"
                water_box = "TIP4PEWBOX"
            else:
                raise KeyError(
                    f"Water type {model_type} is not supported for TIP4P water model."
                )

        self.water_model = {
            "library": library,
            "frcmod": frcmod,
            "water_box": water_box,
        }
