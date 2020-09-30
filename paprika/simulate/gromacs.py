import abc
import glob
import logging
import os
import subprocess as sp
from collections import OrderedDict

from .base_class import BaseSimulation

logger = logging.getLogger(__name__)


class GROMACS(BaseSimulation, abc.ABC):
    """
    A wrapper that can be used to set `GROMACS` simulation parameters.

    .. todo ::
        possibly modify this module to use the official python wrapper of `GROMACS`.

    """

    @property
    def index_file(self) -> str:
        """os.PathLike: `GROMACS` index file that specifies ``groups`` in the system. This is optional in a `GROMACS`
        simulation."""
        return self._index_file

    @index_file.setter
    def index_file(self, value: str):
        self._index_file = value

    @property
    def checkpoint(self) -> str:
        """os.PathLike: Checkpoint file (extension is .cpt) for starting a simulation from a previous state."""
        return self._checkpoint

    @checkpoint.setter
    def checkpoint(self, value: str):
        self._checkpoint = value

    @property
    def control(self):
        """dict: Dictionary for the output control of the MD simulation (frequency of energy, trajectory etc)."""
        return self._control

    @control.setter
    def control(self, value):
        self._control = value

    @property
    def nb_method(self):
        """dict: Dictionary for the nonbonded method options."""
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
    def tc_groups(self) -> list:
        """list: List of groups to apply thermostat based on the groups defined in the ``index_file``."""
        return self._tc_groups

    @tc_groups.setter
    def tc_groups(self, value: list):
        self._tc_groups = value

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
    def prefix(self):
        """str: The prefix for file names generated from this simulation."""
        return self._prefix

    @prefix.setter
    def prefix(self, new_prefix):
        self._prefix = new_prefix
        self.input = new_prefix + ".mdp"
        self.output = new_prefix + ".mdout"
        self.logfile = new_prefix + ".log"
        self.tpr = new_prefix + ".tpr"

    @property
    def custom_mdrun_command(self) -> str:
        """Custom commands for ``mdrun``. The default commands parsed to ``mdrun`` if all the variables are defined is

        .. code::

            gmx mdrun -deffnm ``prefix`` -nt ``n_threads`` -gpu_id ``gpu_devices`` -plumed ``plumed.dat``

        This is useful depending on how `GROMACS` was compiled, e.g. if `GROMACS` is compiled with the MPI library the
        you will need to use the command below:

        .. code::

            mpirun -np 6 gmx_mpi mdrun -deffnm ``prefix`` -ntomp 1 -gpu_id 0 -plumed ``plumed.dat``
        """
        return self._custom_mdrun_command

    @custom_mdrun_command.setter
    def custom_mdrun_command(self, value: str):
        self._custom_mdrun_command = value

    @property
    def grompp_maxwarn(self) -> int:
        """"int: Maximum number of warnings for `GROMPP` to ignore. default=1."""
        return self._grompp_maxwarn

    @grompp_maxwarn.setter
    def grompp_maxwarn(self, value: int):
        self._grompp_maxwarn = value

    def __init__(self):

        super().__init__()

        # I/O
        self._index_file = None
        self._plumed_file = None
        self._plumed_kernel_library = None
        self._custom_mdrun_command = None
        self._tc_groups = None
        self._temperature = 298.15
        self._pressure = 1.01325
        self._grompp_maxwarn = 1

        # File names
        self.input = self._prefix + ".mdp"
        self.output = self._prefix + ".mdout"
        self._checkpoint = None
        self.logfile = self._prefix + ".log"
        self.tpr = self._prefix + ".tpr"

        # Input file
        self._control = OrderedDict()
        self._control["nsteps"] = 5000
        self._control["nstxout"] = 500
        self._control["nstlog"] = 500
        self._control["nstenergy"] = 500
        self._control["nstcalcenergy"] = 500

        self._constraints = OrderedDict()
        self._constraints["constraint_algorithm"] = "lincs"
        self._constraints["constraints"] = "h-bonds"
        self._constraints["lincs_iter"] = 1
        self._constraints["lincs_order"] = 4

        self._nb_method = OrderedDict()
        self._nb_method["cutoff-scheme"] = "Verlet"
        self._nb_method["ns_type"] = "grid"
        self._nb_method["nstlist"] = 10
        self._nb_method["rlist"] = 0.9
        self._nb_method["rcoulomb"] = 0.9
        self._nb_method["rvdw"] = 0.9
        self._nb_method["coulombtype"] = "PME"
        self._nb_method["pme_order"] = 4
        self._nb_method["fourierspacing"] = 0.16
        self._nb_method["vdwtype"] = "Cut-off"
        self._nb_method["DispCorr"] = "EnerPres"
        self._nb_method["pbc"] = "xyz"

        # temperature dict
        self._thermostat = OrderedDict()
        self._barostat = OrderedDict()

    def _config_min(self):
        """
        Configure input settings for a minimization run.
        """
        self.control["continuation"] = "no"
        self.control["integrator"] = "steep"
        self.control["emtol"] = 10.0
        self.control["emstep"] = 0.01
        self.control["nsteps"] = 5000

    def config_pbc_min(self):
        """
        Configure input settings for minimization in periodic boundary conditions.
        """
        self._config_min()
        self.title = "PBC Minimization"
        self.control["continuation"] = "no"
        self.nb_method["nstlist"] = 1

    def _config_md(self):
        self.control["dt"] = 0.002
        self.control["integrator"] = "md"

        self.thermostat["tcoupl"] = "v-rescale"
        self.thermostat["tc-grps"] = "System"
        self.thermostat["ref_t"] = self.temperature
        self.thermostat["tau_t"] = 0.1
        if self.tc_groups:
            self.thermostat["tc-grps"] = self.tc_groups
            self.thermostat["tau_t"] = [0.1 for i in range(len(self.tc_groups))]
            self.thermostat["ref_t"] = [
                self.temperature for i in range(len(self.tc_groups))
            ]

    def config_pbc_md(self, ensemble="npt"):
        """
        Configure input setting for a MD run with periodic boundary conditions.

        Parameters
        ----------
        ensemble: str, optional, default='npt'
            Configure MD setting to use ``nvt`` or ``npt``.

        """
        if ensemble.lower() not in ["nvt", "npt"]:
            raise Exception(f"Thermodynamic ensemble {ensemble} is not supported.")

        self._config_md()
        self.title = f"{ensemble.upper()} MD Simulation"

        if ensemble.lower() == "nvt":
            self.barostat["pcoupl"] = "no"
            self.thermostat["gen_vel"] = "yes"
            self.thermostat["gen_temp"] = self.temperature
            self.thermostat["gen_seed"] = -1

        elif ensemble.lower() == "npt":
            self.thermostat["gen_vel"] = "no"
            self.barostat["pcoupl"] = "Berendsen"
            self.barostat["pcoupltype"] = "isotropic"
            self.barostat["tau_p"] = 5.0
            self.barostat["ref_p"] = self.pressure
            self.barostat["compressibility"] = 4.5e-5

    @staticmethod
    def _write_dict_to_mdp(f, dictionary):
        """
        Write dictionary to file, following GROMACS format.

        Parameters
        ----------
        f : TextIO
            File where the dictionary should be written.
        dictionary : dict
            Dictionary of values.

        """
        for key, val in dictionary.items():
            if val is not None and not isinstance(val, list):
                f.write("{:25s} {:s}\n".format(key, "= " + str(val)))
            elif isinstance(val, list):
                f.write("{:25s} {:s}".format(key, "= "))
                for i in val:
                    f.write("{:s} ".format(str(i)))
                f.write("\n")

    def _write_input_file(self):
        """
        Write the input file specification to file.
        """
        logger.debug("Writing {}".format(self.input))
        with open(os.path.join(self.path, self.input), "w") as mdp:
            mdp.write("{:25s} {:s}\n".format("title", "= " + self.title))

            mdp.write("; Run control\n")
            self._write_dict_to_mdp(mdp, self.control)

            mdp.write("; Nonbonded options\n")
            self._write_dict_to_mdp(mdp, self.nb_method)

            mdp.write("; Bond constraints\n")
            self._write_dict_to_mdp(mdp, self.constraints)

            if self.thermostat:
                mdp.write("; Temperature coupling\n")
                self._write_dict_to_mdp(mdp, self.thermostat)

            if self.barostat:
                mdp.write("; Pressure coupling\n")
                self._write_dict_to_mdp(mdp, self.barostat)

    def run(self, run_grompp=True, overwrite=False, fail_ok=False):
        """
        Method to run Molecular Dynamics simulation with `GROMACS`.

        Parameters
        ----------
        run_grompp: bool, optional, default=True
            Run `GROMPP` to generate ``.tpr`` file before running `MDRUN`
        overwrite: bool, optional, default=False
            Whether to overwrite simulation files.
        fail_ok: bool, optional, default=False
            Whether a failing simulation should stop execution of ``pAPRika``.

        """

        if overwrite or not self.check_complete():
            # Check the type of simulation: Minimization, NVT or NPT
            if self.control["integrator"] in ["steep", "cg", "l-bfgs"]:
                logger.info("Running Minimization at {}".format(self.path))
            elif self.control["integrator"] in [
                "md",
                "md-vv",
                "md-vv-avek",
                "sd",
                "bd",
            ]:
                if self.thermostat and self.barostat:
                    logger.info("Running NPT MD at {}".format(self.path))
                elif not self.barostat:
                    logger.info("Running NVT MD at {}".format(self.path))

            # Set Plumed kernel library to path
            self._set_plumed_kernel()

            # create executable list for GROMPP
            # gmx grompp -f npt.mdp -c coordinates.gro -p topology.top -t checkpoint.cpt -o npt.tpr -n index.ndx
            if run_grompp:

                # Clean previously generated files
                for file in glob.glob(os.path.join(self.path, f"{self.prefix}*")):
                    os.remove(file)

                # Write MDF input file
                self._write_input_file()

                # GROMPP list
                grompp_list = [self.executable, "grompp"]

                grompp_list += [
                    "-f",
                    self.input,
                    "-p",
                    self.topology,
                    "-c",
                    self.coordinates,
                    "-o",
                    self.tpr,
                    "-po",
                    self.output,
                    "-maxwarn",
                    str(self.grompp_maxwarn),
                ]
                if self.checkpoint:
                    grompp_list += ["-t", self.checkpoint]

                if self.index_file:
                    grompp_list += ["-n", self.index_file]

                # Run GROMPP
                grompp_output = sp.Popen(
                    grompp_list,
                    cwd=self.path,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    env=os.environ,
                )
                grompp_stdout = grompp_output.stdout.read().splitlines()
                grompp_stderr = grompp_output.stderr.read().splitlines()

                # Report any stdout/stderr which are output from execution
                if grompp_stdout:
                    logger.info("STDOUT received from GROMACS execution")
                    for line in grompp_stdout:
                        logger.info(line)

                # Not sure how to do this more efficiently/elegantly, "subprocess" seems to treat everything
                # Gromacs spits out from "grompp" as an error.
                if grompp_stderr and any(
                    ["Error" in line.decode("utf-8").strip() for line in grompp_stderr]
                ):
                    logger.info("STDERR received from GROMACS execution")
                    for line in grompp_stderr:
                        logger.error(line)

            # create executable list for MDRUN
            # gmx_mpi mdrun -v -deffnm npt -nt 6 -gpu_id 0 -plumed plumed.dat
            mdrun_list = []

            # Add any user specified command
            if self.custom_mdrun_command is not None:
                if self.executable not in self.custom_mdrun_command:
                    mdrun_list += [self.executable]

                if "mdrun" not in self.custom_mdrun_command:
                    mdrun_list += ["mdrun"]

                mdrun_list += self.custom_mdrun_command.split()

                # Output prefix
                if "-deffnm" not in self.custom_mdrun_command:
                    mdrun_list += ["-deffnm", self.prefix]

                # Add number of threads if not already specified in custom
                if not any(
                    [
                        cpu in self.custom_mdrun_command
                        for cpu in ["-nt", "-ntomp", "-ntmpi", "-ntomp_pme"]
                    ]
                ):
                    mdrun_list += [
                        "-ntomp" if "mpi" in self.executable else "-nt",
                        str(self.n_threads),
                    ]

                # Add gpu id if not already specified in custom
                if (
                    self.gpu_devices is not None
                    and "-gpu_id" not in self.custom_mdrun_command
                ):
                    mdrun_list += ["-gpu_id", str(self.gpu_devices)]

                # Add plumed file if not already specified in custom
                if self.plumed_file and "-plumed" not in self.custom_mdrun_command:
                    mdrun_list += ["-plumed", self.plumed_file]

            else:
                mdrun_list += [self.executable, "mdrun", "-deffnm", self.prefix]

                # Add number of threads
                mdrun_list += [
                    "-ntomp" if "mpi" in self.executable else "-nt",
                    str(self.n_threads),
                ]

                # Add gpu id
                if self.gpu_devices is not None:
                    mdrun_list += ["-gpu_id", str(self.gpu_devices)]

                # Add plumed file
                if self.plumed_file is not None:
                    mdrun_list += ["-plumed", self.plumed_file]

            # Run MDRUN
            mdrun_output = sp.Popen(
                mdrun_list,
                cwd=self.path,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                env=os.environ,
            )
            mdrun_out = mdrun_output.stdout.read().splitlines()
            mdrun_err = mdrun_output.stderr.read().splitlines()

            # Report any stdout/stderr which are output from execution
            if mdrun_out:
                logger.info("STDOUT received from MDRUN execution")
                for line in mdrun_out:
                    logger.info(line)

            # Same reasoning as before for "grompp".
            if mdrun_err and any(
                ["Error" in line.decode("utf-8").strip() for line in mdrun_err]
            ):
                logger.info("STDERR received from MDRUN execution")
                for line in mdrun_err:
                    logger.error(line)

            # Check completion status
            if (
                self.control["integrator"] in ["steep", "cg", "l-bfgs"]
                and self.check_complete()
            ):
                logger.info("Minimization completed...")
            elif self.check_complete():
                logger.info("Simulation completed...")
            else:
                logger.info(
                    "Simulation did not complete when executing the following ...."
                )
                logger.info(" ".join(mdrun_list))
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
        Check for the string "step N" in ``self.output`` file. If "step N" is found, then
        the simulation completed.

        Parameters
        ----------
        alternate_file : os.PathLike, optional, default=None
            If present, check for "step N" in this file rather than ``self.output``.
            Default: None

        Returns
        -------
        complete : bool
            True if "step N" is found in file. False, otherwise.

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
                if (
                    f" step {self.control['nsteps']} " in strings
                    or "Finished mdrun" in strings
                ):
                    complete = True

        if complete:
            logger.debug("{} has TIMINGS".format(output_file))
        else:
            logger.debug("{} does not have TIMINGS".format(output_file))

        return complete
