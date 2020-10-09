import abc
import logging
import os
import subprocess as sp
from collections import OrderedDict

import numpy as np
import parmed as pmd

from paprika.utils import get_dict_without_keys

from .base_class import BaseSimulation

logger = logging.getLogger(__name__)


class NAMD(BaseSimulation, abc.ABC):
    """
    A wrapper that can be used to set NAMD simulation parameters.

    .. todo ::
        Add wrapper support for ``Drude`` polarizable FF, steered-MD and alchemical calculations (perhaps as nested
        class).
    """

    @property
    def colvars_file(self):
        """os.PathLike: The name of the `Colvars`-style restraints file."""
        return self._colvars_file

    @colvars_file.setter
    def colvars_file(self, value):
        self._colvars_file = value

    @property
    def tcl_forces(self):
        """os.PathLike: The name of a TclForces scripts for applying a user defined force
        (e.g., can use this to generate a funnel potential)."""
        return self._tcl_forces

    @tcl_forces.setter
    def tcl_forces(self, value):
        self._tcl_forces = value

    @property
    def control(self):
        """dict: Dictionary for the output control of the MD simulation (frequency of energy, trajectory etc)."""
        return self._control

    @control.setter
    def control(self, value):
        self._control = value

    @property
    def nb_method(self):
        """dict:Dictionary for the non-bonded method options."""
        return self._nb_method

    @nb_method.setter
    def nb_method(self, value):
        self._nb_method = value

    @property
    def constraints(self):
        """dict: Dictionary for the bond constraint options."""
        return self._constraints

    @constraints.setter
    def constraints(self, value):
        self._constraints = value

    @property
    def thermostat(self):
        """dict: Dictionary for the thermostat options."""
        return self._thermostat

    @thermostat.setter
    def thermostat(self, value):
        self._thermostat = value

    @property
    def barostat(self):
        """dict: Dictionary for the barostat options."""
        return self._barostat

    @barostat.setter
    def barostat(self, value):
        self._barostat = value

    @property
    def temperature(self) -> float:
        """float: Desired temperature of the system (Kelvin). Default is 298.15 K."""
        return self._temperature

    @temperature.setter
    def temperature(self, value: float):
        self._temperature = value

    @property
    def pressure(self) -> float:
        """float: Desired pressure of the system (bar). Default is 1.01325 bar."""
        return self._pressure

    @pressure.setter
    def pressure(self, value: float):
        self._pressure = value

    @property
    def implicit(self) -> dict:
        """dict: Dictionary for implicit solvent options."""
        return self._implicit

    @implicit.setter
    def implicit(self, value: dict):
        self._implicit = value

    @property
    def charmm_parameters(self) -> list:
        """list: A list os.PathLike CHARMM parameters if running simulations with the CHARMM force field."""
        return self._charmm_parameters

    @charmm_parameters.setter
    def charmm_parameters(self, value: list):
        self._charmm_parameters = value

    @property
    def prefix(self):
        """
        str: The prefix for file names generated from this simulation.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, new_prefix):
        self._prefix = new_prefix
        self.input = new_prefix + ".in"
        self.logfile = new_prefix + ".log"

    @property
    def cell_basis_vectors(self):
        """
        dict: Dictionary for the cell basis vectors. If not defined the basis vectors and the
        origin is automatically set based on the information from the structure file.
        """
        return self._cell_basis_vectors

    @cell_basis_vectors.setter
    def cell_basis_vectors(self, value):
        self._cell_basis_vectors = value

    @property
    def previous(self) -> str:
        """
        str: Specify the prefix file name for the previous simulation. This is needed to
        get the coordinates, velocities and box information. Unlike AMBER and GROMACS,
        NAMD requires the restart files specified inside the configuration file instead
        of the command-line argument and the information is stored in different files.
        """
        return self._previous

    @previous.setter
    def previous(self, value: str):
        self._previous = value

    @property
    def custom_run_commands(self):
        """
        list: A list of NAMD-supported tcl commands if a custom run is preferred. One
        example of this usage is for a slow-heating MD run.

        .. code-block :: python

            simulation.custom_run_commands = [
                "for {set temp 0} {$temp <= 300} {incr temp 50} {",
                "langevinTemp $temp",
                "reinitvels $temp",
                "run 50000",
                "}",
            ]
        """
        return self._custom_run_commands

    @custom_run_commands.setter
    def custom_run_commands(self, value):
        self._custom_run_commands = value

    def __init__(self):

        super().__init__()

        # I/O
        self._colvars_file = None
        self._temperature = 298.15
        self._pressure = 1.01325

        # File names
        self.input = self._prefix + ".in"
        self.logfile = self._prefix + ".log"
        self._previous = None
        self._custom_run_commands = None
        self._tcl_forces = None
        self._charmm_parameters = []

        self._constraints = OrderedDict()
        self._constraints["rigidBonds"] = "all"
        self._constraints["usesettle"] = "on"

        self._nb_method = OrderedDict()
        self._nb_method["exclude"] = "scaled1-4"
        self._nb_method["1-4scaling"] = 0.833333
        self._nb_method["cutoff"] = 9.0
        self._nb_method["LJCorrection"] = "yes"
        self._nb_method["nonbondedFreq"] = 1
        self._nb_method["fullElectFrequency"] = 2
        self._nb_method["switching"] = "off"
        self._nb_method["wrapAll"] = "on"

        # Placeholder dict
        self._control = OrderedDict()
        self._implicit = OrderedDict()
        self._thermostat = OrderedDict()
        self._barostat = OrderedDict()
        self._cell_basis_vectors = OrderedDict()

    def _config_min(self):
        """
        Configure input settings for minimization.
        """
        self.control["outputname"] = self.prefix
        self.control["minimize"] = 5000

        if ".prmtop" in self.topology:
            self.control["readexclusions"] = "yes"
            self.control["scnb"] = 2.0
        elif ".psf" in self.topology:
            self.nb_method["1-4scaling"] = 1.0
            self.nb_method["cutoff"] = 12.0
            self.nb_method["switching"] = "on"
            self.nb_method["switchdist"] = 10.0
            self.nb_method["exclude"] = "scaled1-4"
            self.nb_method["LJCorrection"] = "off"

        self.control["restartFreq"] = 5000
        self.control["dcdFreq"] = 500
        self.control["xstFreq"] = 500
        self.control["outputEnergies"] = 500
        self.control["outputPressure"] = 500

        if "namd3" in self.executable:
            self.nb_method["CUDASOAintegrate"] = "off"

    def _config_md(self):
        """
        Configure input settings for MD.
        """
        self.constraints["timestep"] = 2.0

        if ".prmtop" in self.topology:
            self.control["readexclusions"] = "yes"
            self.control["scnb"] = 2.0
        elif ".psf" in self.topology:
            self.nb_method["1-4scaling"] = 1.0
            self.nb_method["cutoff"] = 12.0
            self.nb_method["switching"] = "on"
            self.nb_method["switchdist"] = 10.0
            self.nb_method["exclude"] = "scaled1-4"
            self.nb_method["LJCorrection"] = "no"

        self.control["outputname"] = self.prefix
        self.control["restartFreq"] = 5000
        self.control["dcdFreq"] = 500
        self.control["xstFreq"] = 500
        self.control["outputEnergies"] = 500
        self.control["outputPressure"] = 500
        self.control["run"] = 5000

        if "namd3" in self.executable:
            self.nb_method["CUDASOAintegrate"] = "on"
        else:
            self.nb_method["margin"] = 1.0
            self.constraints["stepspercycle"] = 20

        self.thermostat["langevin"] = "on"
        self.thermostat["langevinDamping"] = 1.0
        self.thermostat["langevinTemp"] = self.temperature
        self.thermostat["langevinHydrogen"] = "off"

    def config_gb_min(self):
        """
        Configure input settings for minimization with implicit solvent.
        """
        self.title = "Implicit solvent Minimization"

        self._config_min()

        self.implicit["GBIS"] = "on"
        self.implicit["solventDielectric"] = 78.5
        self.implicit["ionConcentration"] = 0.0

        self.nb_method["cutoff"] = 999.0
        self.nb_method["LJCorrection"] = "no"
        self.nb_method["nonbondedFreq"] = 2
        self.nb_method["fullElectFrequency"] = 4

        if ".psf" in self.topology:
            self.nb_method["switching"] = "off"
            del self.nb_method["switchdist"]

    def config_gb_md(self):
        """
        Configure input settings for MD with implicit solvent.
        """
        self.title = "Implicit solvent MD Simulation"

        self._config_md()

        self.implicit["GBIS"] = "on"
        self.implicit["solventDielectric"] = 78.5
        self.implicit["ionConcentration"] = 0.0

        self.constraints["timestep"] = 1.0

        self.nb_method["cutoff"] = 999.0
        self.nb_method["LJCorrection"] = "no"
        self.nb_method["nonbondedFreq"] = 2
        self.nb_method["fullElectFrequency"] = 4

        if ".psf" in self.topology:
            self.nb_method["switching"] = "off"
            del self.nb_method["switchdist"]

    def config_pbc_min(self, calculate_cell_vectors=True):
        """
        Configure input setting for an energy minimization run with periodic boundary conditions.

        Parameters
        ----------
        calculate_cell_vectors: bool, optional, default=True
            Calculate cell basis vectors based on `self.coordinates`. This is usually needed at the start of a
            simulation. When continuing a simulation NAMD reads in the ``*.xst`` file to get the basis vectors.

        """
        self.title = "PBC Minimization"

        self._config_min()

        self.nb_method["PME"] = "yes"
        self.nb_method["PMEGridSpacing"] = 1.0

        if calculate_cell_vectors:
            self._get_cell_basis_vectors()

    def config_pbc_md(
        self, ensemble="npt", barostat="langevin", calculate_cell_vectors=False
    ):
        """
        Configure input setting for a MD run with periodic boundary conditions.

        Parameters
        ----------
        ensemble: str, optional, default='npt'
            Configure MD setting to use ``nvt`` or ``npt``.
        barostat: str, optional, default='langevin'
            Option to choose one of two barostat implemented in NAMD: (1) `Langevin piston` or (2) `Berendsen`.
        calculate_cell_vectors: bool, optional, default=False
            Calculate cell basis vectors based on `self.coordinates`. This is usually needed at the start of a
            simulation. When continuing a simulation NAMD reads in the ``*.xsc`` file to get the cell basis vectors.

        """
        if ensemble.lower() not in ["nvt", "npt"]:
            raise Exception(f"Thermodynamic ensemble {ensemble} is not supported.")

        self.title = f"{ensemble.upper()} MD Simulation"

        self._config_md()
        self.nb_method["PME"] = "yes"
        self.nb_method["PMEGridSpacing"] = 1.0

        if ensemble.lower() == "npt":
            self.barostat["useGroupPressure"] = "yes"
            self.barostat["useFlexibleCell"] = "no"
            self.barostat["useConstantArea"] = "no"

            if barostat.lower() == "langevin":
                self.barostat["langevinPiston"] = "on"
                self.barostat["langevinPistonTarget"] = self.pressure
                self.barostat["langevinPistonPeriod"] = 200.0
                self.barostat["langevinPistonDecay"] = 100.0
                self.barostat["langevinPistonTemperature"] = self.temperature

            elif barostat.lower() == "berendsen":
                self.barostat["BerendsenPressure"] = "on"
                self.barostat["BerendsenPressureTarget"] = self.pressure
                self.barostat["BerendsenPressureCompressibility"] = 4.57e-5
                self.barostat["BerendsenPressureRelaxationTime"] = 100
                self.barostat["BerendsenPressureFreq"] = 10

        if calculate_cell_vectors:
            self._get_cell_basis_vectors()

    def _get_cell_basis_vectors(self):
        """
        Function to calculate the PBC cell basis vectors (needed when running a simulation for the first time).
        """
        structure = pmd.load_file(
            os.path.join(self.path, self.topology),
            os.path.join(self.path, self.coordinates),
            structure=True,
        )
        coordinates = structure.coordinates
        masses = np.ones(len(coordinates))
        center = pmd.geometry.center_of_mass(coordinates, masses)

        self.cell_basis_vectors["cellOrigin"] = list(center)

        # Check if box info exists
        if structure.box_vectors:
            import simtk.unit as unit

            self.cell_basis_vectors["cellBasisVector1"] = structure.box_vectors[
                0
            ].value_in_unit(unit.angstroms)
            self.cell_basis_vectors["cellBasisVector2"] = structure.box_vectors[
                1
            ].value_in_unit(unit.angstroms)
            self.cell_basis_vectors["cellBasisVector3"] = structure.box_vectors[
                2
            ].value_in_unit(unit.angstroms)

        else:
            self.cell_basis_vectors["cellBasisVector1"] = [
                structure.coordinates[:, 0].max() - structure.coordinates[:, 0].min(),
                0.0,
                0.0,
            ]
            self.cell_basis_vectors["cellBasisVector2"] = [
                0.0,
                structure.coordinates[:, 1].max() - structure.coordinates[:, 1].min(),
                0.0,
            ]
            self.cell_basis_vectors["cellBasisVector3"] = [
                0.0,
                0.0,
                structure.coordinates[:, 2].max() - structure.coordinates[:, 2].min(),
            ]

    @staticmethod
    def _write_dict_to_conf(f, dictionary):
        """
        Write dictionary to file, following NAMD format.

        Parameters
        ----------
        f : TextIO
            File where the dictionary should be written.
        dictionary : dict
            Dictionary of values.

        """
        for key, val in dictionary.items():
            if val is not None and not isinstance(val, list):
                f.write("{:30s} {:s}\n".format(key, str(val)))
            elif isinstance(val, list):
                f.write("{:30s} ".format(key))
                for i in val:
                    f.write("{:s} ".format(str(i)))
                f.write("\n")

    def _write_input_file(self):
        """
        Write the MD configuration to file.
        """

        logger.debug("Writing input file: {}".format(self.input))

        with open(os.path.join(self.path, self.input), "w") as conf:
            conf.write("## {:s}\n".format(self.title))

            # Topology and Coordinates files
            # AMBER PRMTOP/RST7
            if ".prmtop" in self.topology:
                conf.write("# AMBER files\n")
                conf.write("{:30s} {:s}\n".format("amber", "yes"))
                conf.write("{:30s} {:s}\n".format("parmfile", self.topology))

                if any([ext in self.coordinates for ext in [".rst7", ".inpcrd"]]):
                    conf.write("{:30s} {:s}\n".format("ambercoor", self.coordinates))
                elif ".pdb" in self.coordinates:
                    conf.write("{:30s} {:s}\n".format("coordinates", self.coordinates))
                else:
                    raise FileExistsError(
                        f"Coordinates file {self.coordinates} not does not exist."
                    )

                conf.write(
                    "{:30s} {:s}\n".format(
                        "readexclusions", self.control["readexclusions"]
                    )
                )
                conf.write("{:30s} {:s}\n".format("scnb", str(self.control["scnb"])))

            # GROMACS TOP/GRO
            elif ".top" in self.topology:
                conf.write("# GROMACS files\n")
                conf.write("{:30s} {:s}\n".format("gromacs", "on"))
                conf.write("{:30s} {:s}\n".format("grotopfile", self.topology))

                if ".gro" in self.coordinates:
                    conf.write("{:30s} {:s}\n".format("grocoorfile", self.coordinates))
                elif ".pdb" in self.coordinates:
                    conf.write("{:30s} {:s}\n".format("coordinates", self.coordinates))
                else:
                    raise FileExistsError(
                        f"Coordinates file {self.coordinates} not does not exist."
                    )

            # CHARMM PSF/PDB
            elif ".psf" in self.topology:
                conf.write("# CHARMM files\n")
                conf.write("{:30s} {:s}\n".format("paraTypeXplor", "on"))
                conf.write("{:30s} {:s}\n".format("structure", self.topology))

                if ".pdb" in self.coordinates:
                    conf.write("{:30s} {:s}\n".format("coordinates", self.coordinates))
                else:
                    raise FileExistsError(
                        f"Coordinates file {self.coordinates} not does not exist."
                    )

                if not self.charmm_parameters:
                    raise FileExistsError(
                        "CHARMM parameter file(s) not specified, please specify the parameters "
                        "with the variable 'namd.charmm_parameters'."
                    )

                conf.write("{:30s} {:s}\n".format("paraTypeCharmm ", "on"))
                for parameter in self.charmm_parameters:
                    conf.write("{:30s} {:s}\n".format("parameters", parameter))

            conf.write("# File I/O and run control\n")
            conf.write("{:30s} {:s}\n".format("outputname", self.control["outputname"]))

            # Load in previous files (.coor, .vel and/or .xsc)
            if self.previous is not None:
                # Coordinates
                bincoordinates = self.previous + ".restart.coor"
                if not os.path.isfile(os.path.join(self.path, bincoordinates)):
                    bincoordinates = self.previous + ".coor"

                conf.write("{:30s} {:s}\n".format("bincoordinates", bincoordinates))

                # Velocities
                binvelocities = self.previous + ".restart.vel"
                if not os.path.isfile(os.path.join(self.path, binvelocities)):
                    binvelocities = self.previous + ".vel"

                conf.write("{:30s} {:s}\n".format("binvelocities", binvelocities))

                # PBC box information
                if not self.implicit:
                    extendedsystem = self.previous + ".restart.xsc"
                    if not os.path.isfile(os.path.join(self.path, extendedsystem)):
                        extendedsystem = self.previous + ".xsc"

                    conf.write("{:30s} {:s}\n".format("extendedSystem", extendedsystem))

            else:
                conf.write("{:30s} {:s}\n".format("temperature", str(self.temperature)))

            conf.write("# Output frequency\n")
            self._write_dict_to_conf(
                conf,
                get_dict_without_keys(
                    self.control,
                    "minimize",
                    "run",
                    "outputname",
                    "readexclusions",
                    "scnb",
                ),
            )

            if self.implicit:
                conf.write("# Generalized Born Implicit Solvent\n")
                self._write_dict_to_conf(conf, self.implicit)

            if self.cell_basis_vectors:
                conf.write("# Cell basis vectors\n")
                self._write_dict_to_conf(conf, self.cell_basis_vectors)

            conf.write("# Nonbonded options\n")
            self._write_dict_to_conf(conf, self.nb_method)

            conf.write("# constraints\n")
            self._write_dict_to_conf(conf, self.constraints)

            conf.write("# Temperature coupling\n")
            if self.custom_run_commands:
                if any(["langevinTemp" in line for line in self.custom_run_commands]):
                    self._write_dict_to_conf(
                        conf, get_dict_without_keys(self.thermostat, "langevinTemp")
                    )
                else:
                    self._write_dict_to_conf(conf, self.thermostat)
            else:
                self._write_dict_to_conf(conf, self.thermostat)

            if self.barostat:
                conf.write("# Pressure coupling\n")
                self._write_dict_to_conf(conf, self.barostat)

            if self.colvars_file:
                conf.write("# Collective variables\n")
                conf.write("{:30s} {:s}\n".format("colvars", "on"))
                conf.write("{:30s} {:s}\n".format("colvarsConfig", self.colvars_file))

            if self.plumed_file:
                conf.write("# Plumed Restraints\n")
                conf.write("{:30s} {:s}\n".format("plumed", "on"))
                conf.write("{:30s} {:s}\n".format("plumedfile", self.plumed_file))

            if self.tcl_forces:
                conf.write("# TclForces\n")
                conf.write("{:30s} {:s}\n".format("tclForces", "on"))
                conf.write("{:30s} {:s}\n".format("tclForcesScript", self.tcl_forces))

            if not self.custom_run_commands:
                if "minimize" in self.control:
                    conf.write("# Minimization\n")
                    conf.write(
                        "{:30s} {:s}\n".format(
                            "minimize", str(self.control["minimize"])
                        )
                    )
                elif "run" in self.control:
                    conf.write("# Molecular dynamics\n")
                    conf.write("{:30s} {:s}\n".format("run", str(self.control["run"])))
            else:
                conf.write("# Custom Run\n")
                for line in self.custom_run_commands:
                    conf.write(line + "\n")

    def run(self, overwrite=True, fail_ok=False):
        """
        Method to run Molecular Dynamics simulation with NAMD.

        Parameters
        ----------
        overwrite: bool, optional, default=False
            Whether to overwrite simulation files.
        fail_ok: bool, optional, default=False
            Whether a failing simulation should stop execution of ``pAPRika``.

        """

        if overwrite or not self.check_complete():
            if "Minimization" in self.title:
                logger.info("Running Minimization at {}".format(self.path))
            elif "MD Simulation" in self.title:
                if self.thermostat and self.barostat:
                    logger.info("Running NPT MD at {}".format(self.path))
                elif not self.barostat:
                    logger.info("Running NVT MD at {}".format(self.path))

            # Check restraints file
            if self.colvars_file and self.plumed_file:
                raise Exception(
                    "Cannot use both Colvars-style and Plumed-style restraints at the same time."
                )

            # Set Plumed kernel library to path
            if self.plumed_file:
                self._set_plumed_kernel()

            # Write input file
            self._write_input_file()

            # Create executable list
            exec_list = self.executable.split()

            if "+p" not in self.executable:
                exec_list += ["+p", str(self.n_threads)]

            if "+devices" not in self.executable and self.gpu_devices is not None:
                exec_list += ["+devices", str(self.gpu_devices)]

            exec_list += [self.input]

            logger.debug("Exec line: " + " ".join(exec_list))

            # Execute
            namd_output = sp.Popen(
                exec_list,
                cwd=self.path,
                stdout=open(os.path.join(self.path, self.logfile), "w"),
                stderr=sp.PIPE,
                env=os.environ,
            )
            namd_stderr = namd_output.stderr.read().splitlines()

            if namd_stderr and any(
                ["ERROR:" in line.decode("utf-8").strip() for line in namd_stderr]
            ):
                logger.info("STDERR received from NAMD execution")
                for line in namd_stderr:
                    logger.error(line)

            # Check completion status
            if "minimize" in self.control and self.check_complete():
                logger.info("Minimization completed...")
            elif self.check_complete():
                logger.info("MD completed ...")
            else:
                logger.info(
                    "Simulation did not complete when executing the following ...."
                )
                logger.info(" ".join(exec_list))
                if not fail_ok:
                    raise Exception(
                        "Exiting due to failed simulation! Check logging info."
                    )
        else:
            logger.info(
                "Completed output detected ... Skipping. Use: run(overwrite=True) to overwrite"
            )

    def check_complete(self, alternate_file=None):
        """
        Check for the string "WallClock" in ``self.output`` file. If "WallClock" is found, then
        the simulation completed.

        Parameters
        ----------
        alternate_file : os.PathLike, optional, default=None
            If present, check for "WallClock" in this file rather than ``self.output``.
            Default: None

        Returns
        -------
        complete : bool
            True if "WallClock" is found in file. False, otherwise.
        """
        # Assume not completed
        complete = False

        if alternate_file:
            output_file = alternate_file
        else:
            output_file = os.path.join(self.path, self.logfile)

        if os.path.isfile(output_file):
            with open(output_file, "r") as f:
                strings = f.read()
                if "WallClock:" in strings:
                    complete = True

        if complete:
            logger.debug("{} has TIMINGS".format(output_file))
        else:
            logger.debug("{} does not have TIMINGS".format(output_file))

        return complete
