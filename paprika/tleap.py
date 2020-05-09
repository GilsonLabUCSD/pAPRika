import os as os
import re as re
import subprocess as sp
import logging as log
import numpy as np
import parmed as pmd


N_A = 6.0221409 * 10 ** 23
ANGSTROM_CUBED_TO_LITERS = 1 * 10 ** -27


class System(object):
    """
    Class for building AMBER prmtop/rst7 files with tleap.

    Parameters
    ----------
    template_file : str
        Name of template file to read boilerplate `tleap` input (e.g., `source leaprc`... lines). Any
        frcmod/mol2/pdb files which are loaded by the template must be present in output_path.
        Default: None
    template_lines : lst
        The list of tleap commands which are needed to build the system. Either the template_file or
        the template_lines needs to be set prior to running build(). Files called for loading should
        be present in output_path (see above). Default: None
    loadpdb_file : str
        If specified, this PDB file will be loaded via loadpdb instead of whatever file is specified
        in the template. This file should be present in output_path, or if you specify it with a path,
        the path must be relative to output_path. Default: None
    pbc_type : str
        Type of solvation (cubic, rectangular, octahedral, or None). Default: cubic
    buffer_value : float
        The desired solvation buffer value. This will add a layer of water around your solute with the
        minimum distance to the periodic box edge defined by buffer_value. If target_waters is set, it
        will override this value. Default: 12.0
    target_waters : int
        The desired number of waters to solvate your system with. If specified, this will override any
        buffer_value settings. Default: None
    water_box : str
        The type of water box, e.g., TIP3PBOX (see AMBER manual for acceptable values). Default: TIP3PBOX
    neutralize : bool
        Whether to neutralize the system with counterions. Default: True
    counter_cation : str
        If neutralize=True, this positive ion will be used. Default: Na+
    counter_anion : str
        If neuralize=True, this negative ion will be used. Default: Cl-
    add_ions : list
        A list of additional ions to be added to the system, specified by the residue name and
        an indication of the amount: an integer number, the molarity (M), or the molality (m).
        For example, the following list shows valid examples:
        add_ions = ['MG', 2, 'Cl-', 4, 'K', '0.150M', 'BR', '0.150M', 'NA', '0.050m', 'F', '0.050m']
        Default: None
    output_path : str
        Directory path where tleap input files will be written and executed. Associated input
        files (frcmod, mol2, pdb, etc) need to present in that path.
    output_prefix : str
        Append a prefix to files created by this module. Default: 'build'
    unit : str
        The tleap unit name. Default: 'model'
    exponent : int
        The initial value used for dynamically adjusting the buffer value. Default: 1
    min_exponent_limit : int
        The minimum limit that we let exponent get to. If it gets too small, it has no effect on
        the number of waters added. Default: -5
    cyc_since_last_exp_change : int
        A count of the number of buffer_value adjustments since the last exponent change. Should
        always start at 0. Default: 0
    max_cycles : int
        The maximum number of buffer_value adjustment cycles.  Default: 50
    manual_switch_thresh : int
        The threshold difference between waters actually added and target_waters for which we can
        safely switch to manual removal of waters. If too large, then there will be big air pockets
        in our box. If too small, it will take forever to converge. If the user does not set this
        value, it will set proportionately to the target_waters. Default: (target_waters)**(1./3.)
    waters_to_remove : list
        A list of the water residues to be manually removed after solvation. Default: None
    add_ion_residues : list
        A property formated list of ions to add to the system.  The format is the same as add_ions
        except that only integer amounts are allowed.  This is set by set_additional_ions. Default: None
    kg_per_mol_solvent : float
        kg/mol for the solvent. Used for equation computing molality. Default for water: 0.018
    buffer_val_history : list
        The history of buffer_value adjustments stored in a list. Default: [0]
    wat_added_history : list
        The history of number of waters added which correlate to the values in buffer_val_history.
        Default: [0]
    write_save_lines : bool
        Whether or not to do saveamberparm and savepdb when tleap is executed. During optimization
        loops it speeds up the process to set as False, but it must be returned to True for the final
        execution. Default: True

    """

    def __init__(self):

        # User Settings: Defaults
        self.template_file = None
        self.template_lines = None
        self.loadpdb_file = None
        self.pbc_type = "cubic"
        self.buffer_value = 12.0
        self.target_waters = None
        self.water_box = "TIP3PBOX"
        self.neutralize = True
        self.counter_cation = "Na+"
        self.counter_anion = "Cl-"
        self.add_ions = None
        self.output_path = "./"
        self.output_prefix = "build"

        # Advanced Settings: Defaults
        self.unit = "model"
        self.exponent = 1
        self.min_exponent_limit = -5
        self.cyc_since_last_exp_change = 0
        self.max_cycles = 50
        self.manual_switch_thresh = None
        self.waters_to_remove = None
        self.add_ion_residues = None
        self.kg_per_mol_solvent = 0.018
        self.buffer_val_history = [0]
        self.wat_added_history = [0]
        self.write_save_lines = True

    def build(self):
        """
        Build the tleap.System

        """

        log.debug("Running tleap.build() in {}".format(self.output_path))

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

    def filter_template(self):
        """
        Filter out any template_lines that may interfere with solvation.

        """

        filtered_lines = []
        for line in self.template_lines:
            # Find loadpdb line, replace pdb file if necessary, set unit name
            if re.search("loadpdb", line):
                words = line.rstrip().replace("=", " ").split()
                if self.loadpdb_file is None:
                    self.loadpdb_file = words[2]
                self.unit = words[0]
                filtered_lines.append(
                    "{} = loadpdb {}".format(self.unit, self.loadpdb_file)
                )
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

    def write_input(self):
        """
        Write a tleap input file based on template_lines and other things we have set.

        """

        file_path = os.path.join(self.output_path, self.output_prefix + ".tleap.in")
        if not os.path.exists(os.path.dirname(file_path)):
            try:
                os.makedirs(os.path.dirname(file_path))
            except OSError as e:
                raise

        with open(file_path, "w") as f:
            for line in self.template_lines:
                f.write(line + "\n")

            if self.pbc_type == "cubic":
                f.write(
                    "solvatebox {} {} {} iso\n".format(
                        self.unit, self.water_box, self.buffer_value
                    )
                )
            elif self.pbc_type == "rectangular":
                f.write(
                    "solvatebox {} {} {{10.0 10.0 {}}}\n".format(
                        self.unit, self.water_box, self.buffer_value
                    )
                )
            elif self.pbc_type == "octahedral":
                f.write(
                    "solvateoct {} {} {} iso\n".format(
                        self.unit, self.water_box, self.buffer_value
                    )
                )
            elif self.pbc_type is None:
                f.write("# Skipping solvation ...\n")
            else:
                raise Exception(
                    "Incorrect pbctype value provided: "
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
        Execute `tleap`.

        Returns
        -------
        output : list
            The tleap stdout returned as a list.

        """

        self.check_for_leap_log()
        file_name = self.output_prefix + ".tleap.in"

        output = sp.Popen(
            ["tleap", "-s ", "-f ", file_name],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=self.output_path,
        )
        output = output.stdout.read().decode().splitlines()

        self.grep_leap_log()
        return output

    def grep_leap_log(self):
        """
        Check for a few keywords in the `tleap` output.

        """
        try:
            with open(self.output_path + "leap.log", "r") as file:
                for line in file.readlines():
                    if re.search(
                        "ERROR|WARNING|Warning|duplicate|FATAL|Could|Fatal|Error", line
                    ):
                        log.warning(
                            "It appears there was a problem with solvation: check `leap.log`..."
                        )
        except BaseException:
            return

    def check_for_leap_log(self, log_file="leap.log"):
        """
        Check if `leap.log` exists in output_path, and if so, delete so the current run doesn't append.

        Parameters
        ----------
        log_file : str
            Name of the tleap logfile. Default: leap.log

        """
        log_file_path = os.path.join(self.output_path, log_file)
        try:
            os.remove(log_file_path)
            log.debug("Deleted existing leap logfile: " + log_file_path)
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
        # set buffer_value = 1.0 because that seems to lead to smooth convergence
        # of the solvate() routine.
        else:
            self.buffer_value = 1.0

        # If the user has not set manual_switch_thresh, we will set it proportionately
        # to the target_waters. This will control when we can start manually deleting
        # waters rather than adjusting the buffer_value.
        if self.manual_switch_thresh is None:
            self.manual_switch_thresh = int(np.ceil(self.target_waters ** (1. / 3.)))
            if self.manual_switch_thresh < 12:
                self.manual_switch_thresh = 12
            log.debug(
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
            self.wat_added_history.append(waters)
            self.buffer_val_history.append(self.buffer_value)
            log.debug(
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
            log.debug(
                "The added waters ({}) didn't reach the manual_switch_thresh ({}) with max_cycles ({}), but we'll try manual removal anyway.".format(
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
        Run the solvation with write_save_lines = True.

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
            log.info(
                "The final solvation step added the wrong number of waters. Repeating ..."
            )

    def count_waters(self):
        """
        Quickly check the number of waters added for a given buffer size.

        Returns
        -------
        waters : int

        """
        waters = self.count_residues()["WAT"]
        return waters

    def count_residues(self, print_results=False):
        """
        Run and parse `tleap` output and return a dictionary of residues in the structure.

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

        # log.debug(residues)
        if print_results:
            for key, value in sorted(residues.items()):
                log.info("{:10s} {:10.0f}".format(key, value))

        return residues

    def set_additional_ions(self):
        """
        Determine whether additional ions (instead of or in addition to neutralization) are requested...

        Sets
        -------
        self.add_ion_residues : list
            A processed list of ions and their amounts that can be passed to `tleap`

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
                    np.ceil(
                        float(amount[:-1])
                        * self.target_waters
                        * self.kg_per_mol_solvent
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
                number_to_add = int(np.ceil(number_of_atoms * liters))
                self.add_ion_residues.append(number_to_add)
            else:
                raise Exception("Unanticipated error calculating how many ions to add.")

    def get_volume(self):
        """
        Run and parse `tleap` output and return the volume of the structure.

        Returns
        -------
        volume : float
            The volume of the structure in cubic angstroms

        """
        output = self.run()
        # Return the total simulation volume
        for line in output:
            line = line.strip()
            if "Volume" in line:
                match = re.search("Volume(.*)", line)
                volume = float(match.group(1)[1:-4])
                return volume
        log.warning("Could not determine total simulation volume.")
        return None

    def remove_waters_manually(self):
        """
        Remove a few water molecules manually with `tleap` to exactly match a desired number of waters.

        """

        cycle = 0
        max_cycles = 100
        waters = self.wat_added_history[-1]
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
                log.debug(
                    "Detected trouble with manually removing water. Increasing the number of surplus waters by {}".format(
                        additional_water
                    )
                )

            self.waters_to_remove = water_residues[-1 * water_surplus :]
            log.debug("Manually removing waters... {}".format(self.waters_to_remove))

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
        Run and parse `tleap` output and return the a list of water residues.

        Returns
        -------
        water_residues : list
            A list of the water residues in the structure

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

        Sets
        -------
        self.buffer_value : float
            A new buffer size to try
        self.exponent : int
            Adjusts the order of magnitue of buffer value changes
        self.cyc_since_last_exp_change : int
            Resets this value to 0 when exponent value is changed to help with the logic gates.

        """

        # If the number of waters was less than the target and is now greater than the target, make the buffer smaller
        # smaller
        if (
            self.wat_added_history[-2] < self.target_waters
            and self.wat_added_history[-1] > self.target_waters
        ):
            # If its been more than one round since last exponent change,
            # change exponent
            if self.cyc_since_last_exp_change > 1:
                log.debug("Adjustment loop 1a")
                self.exponent -= 1
                self.buffer_value = self.buffer_val_history[-1] + -5 * (
                    10 ** self.exponent
                )
                self.cyc_since_last_exp_change = 0
            else:
                log.debug("Adjustment loop 1b")
                self.buffer_value = self.buffer_val_history[-1] + -1 * (
                    10 ** self.exponent
                )
                self.cyc_since_last_exp_change += 1
        # If the number of waters was greater than the target and is now less
        # than the target, make the buffer bigger
        elif (
            self.wat_added_history[-2] > self.target_waters
            and self.wat_added_history[-1] < self.target_waters
        ):
            # If its been more than one round since last exponent change,
            # change exponent
            if self.cyc_since_last_exp_change > 1:
                log.debug("Adjustment loop 2a")
                self.exponent -= 1
                self.buffer_value = self.buffer_val_history[-1] + 5 * (
                    10 ** self.exponent
                )
                self.cyc_since_last_exp_change = 0
            else:
                log.debug("Adjustment loop 2b")
                self.buffer_value = self.buffer_val_history[-1] + 1 * (
                    10 ** self.exponent
                )
                self.cyc_since_last_exp_change += 1
        # If the last two rounds of solvation have too many waters, make the
        # buffer smaller...
        elif (
            self.wat_added_history[-2] > self.target_waters
            and self.wat_added_history[-1] > self.target_waters
        ):
            log.debug("Adjustment loop 3")
            self.buffer_value = self.buffer_val_history[-1] + -1 * (10 ** self.exponent)
            self.cyc_since_last_exp_change += 1
        # If the last two rounds of solvation had too few waters, make the
        # buffer bigger...
        elif (
            self.wat_added_history[-2] < self.target_waters
            and self.wat_added_history[-1] < self.target_waters
        ):
            log.debug("Adjustment loop 4")
            self.buffer_value = self.buffer_val_history[-1] + 1 * (10 ** self.exponent)
            self.cyc_since_last_exp_change += 1
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
        overwrites the prmtop file: "self.output_path/self.output_prefix.prmtop"

        Parameters
        ----------
        options : str
            Optional keyword(s) for the repartitioning and following the 
            parmed.tools.actions documentation the usage is '[<mass>] [dowater]'.

        """

        prmtop = os.path.join(self.output_path, self.output_prefix + ".prmtop")

        structure = pmd.load_file(prmtop, structure=True)

        pmd.tools.actions.HMassRepartition(structure, arg_list=options).execute()

        structure.save(prmtop, overwrite=True)

