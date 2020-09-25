import abc
import logging
import os
import subprocess as sp
from collections import OrderedDict
from sys import platform

logger = logging.getLogger(__name__)


class BaseSimulation(object):
    """
    A base class for a Molecular Dynamics simulation wrapper.
    """

    @property
    def path(self):
        """os.PathLike: The path for the creation and execution of files for this protocol."""
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

    @property
    def executable(self):
        """str: The AMBER executable that will be used.

        .. note::
            This could be made safer by making an ``ENUM`` of ``sander``, ``pmemd`` and ``pmemd.cuda``.
        """
        return self._executable

    @executable.setter
    def executable(self, value):
        self._executable = value

    @property
    def n_threads(self) -> int:
        """int: Number of threads to use for the simulation."""
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value: int):
        self._n_threads = value

    @property
    def gpu_devices(self):
        """int or str: A wrapper around the environmental variable ``CUDA_VISIBLE_DEVICES``."""
        return self._gpu_devices

    @gpu_devices.setter
    def gpu_devices(self, value):
        self._gpu_devices = value

    @property
    def plumed_file(self) -> str:
        """os.PathLike: The name of the Plumed-style restraints file.

        .. note::
            When running simulations with a MD engine that supports Plumed and other restraint module(s),
            you can only use one or the other. If both are specified, an ``Exception`` will be thrown.
        """
        return self._plumed_file

    @plumed_file.setter
    def plumed_file(self, value: str):
        self._plumed_file = value

    @property
    def plumed_kernel_library(self) -> str:
        """os.PathLike: The file path of the Plumed kernel library if it is different from the default. The default
        name is ``libplumedKernel`` with its extension determined automatically depending on the operating system and is
        located in ``os.environ['CONDA_PREFIX']/lib``. This variable is useful for users who did not install or want to
        use Plumed from the conda repository."""
        return self._plumed_kernel_library

    @plumed_kernel_library.setter
    def plumed_kernel_library(self, value: str):
        self._plumed_kernel_library = value

    @property
    def topology(self) -> str:
        """os.PathLike: The topology file to simulate."""
        return self._topology

    @topology.setter
    def topology(self, value: os.PathLike):
        self._topology = value

    @property
    def coordinates(self) -> str:
        """os.PathLike: The coordinates of the MD system."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value: str):
        self._coordinates = value

    @property
    def prefix(self) -> str:
        """
        The prefix for file names generated from this simulation.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, new_prefix: str):
        self._prefix = new_prefix

    @property
    def phase(self) -> str:
        """str: Which ``phase`` of the calculation to simulate.

        .. note ::
            This could probably be combined with the ``window`` attribute for simplicity.
        """
        return self._phase

    @phase.setter
    def phase(self, value: str):
        self._phase = value

    @property
    def window(self) -> str:
        """str: Which ``window`` of the calculation to simulate."""
        return self._window

    @window.setter
    def window(self, value: str):
        self._window = value

    @property
    def converged(self) -> bool:
        """bool: Whether the simulation is converged.

        .. warning ::
            This is just a placeholder and is not implemented yet.
        """
        return self._converged

    @converged.setter
    def converged(self, value: bool):
        self._converged = value

    def __init__(self):
        self.title = "PBC MD Simulation"
        self._path = "."
        self._executable = "sander"
        self._n_threads = 1
        self._gpu_devices = None
        self._plumed_file = None
        self._plumed_kernel_library = None
        self._topology = "prmtop"
        self._coordinates = "rst7"
        self._prefix = "md"
        self._phase = None
        self._window = None
        self.converged = False

    def _set_plumed_kernel(self):
        # Set the Plumed kernel library path
        if self.plumed_kernel_library is None:
            # using "startswith" as recommended from https://docs.python.org/3/library/sys.html#sys.platform
            if platform.startswith("linux"):
                self.plumed_kernel_library = os.path.join(
                    os.environ["CONDA_PREFIX"], "lib", "libplumedKernel.so"
                )
            elif platform.startswith("darwin"):
                self.plumed_kernel_library = os.path.join(
                    os.environ["CONDA_PREFIX"], "lib", "libplumedKernel.dylib"
                )
            else:
                logger.warning(
                    f"Plumed is not supported with the '{platform}' platform."
                )

        # Add Plumed kernel library to env dict
        if os.path.isfile(self.plumed_kernel_library):
            os.environ = dict(os.environ, PLUMED_KERNEL=self.plumed_kernel_library)
        else:
            logger.warning(
                f"Plumed kernel library {self.plumed_kernel_library} does not exist, "
                "PLUMED_KERNEL environment variable is not set."
            )

    def _config_min(self):
        """Place holder function."""
        raise NotImplementedError()

    def _config_md(self):
        """Place holder function."""
        raise NotImplementedError()

    def config_gb_min(self):
        """Place holder function."""
        raise NotImplementedError()

    def config_gb_md(self):
        """Place holder function."""
        raise NotImplementedError()

    def config_pbc_min(self):
        """Place holder function."""
        raise NotImplementedError()

    def config_pbc_md(self):
        """Place holder function."""
        raise NotImplementedError()

    def _write_input_file(self):
        """Place holder function."""
        raise NotImplementedError()

    def run(self):
        """Place holder function."""
        raise NotImplementedError()

    def check_complete(self):
        """Place holder function."""
        raise NotImplementedError()

    def check_converged(self):
        """Place holder function"""
        raise NotImplementedError()


class Amber(BaseSimulation, abc.ABC):
    """
    A wrapper for setting parameters and running MD with AMBER.
    """

    @property
    def restraint_file(self) -> str:
        """os.PathLike: The file containing NMR-style restraints for AMBER.

        .. note ::
            When running AMBER simulations, you can only use either an AMBER NMR-style
            restraints or a Plumed-style restraints and not both. If both are specified,
            an ``Exception`` will be thrown.
        """
        return self._restraint_file

    @restraint_file.setter
    def restraint_file(self, value: str):
        self._restraint_file = value

    @property
    def cntrl(self) -> dict:
        """
        An ordered dictionary of simulation "control" namelist parameters. The main purpose of this class attribute is
        to make it easy to override certain simulation parameters, such as positional restraints on dummy atoms, and
        the inclusion of exclusion of the NMR-style APR restraints.

        As of AMBER20, these are described, in part, in chapter 20 on ``pmemd``.

        .. note ::
            I can't recall why we wanted an ``OrderedDict`` here.

        .. note ::
            **This is fragile** and this could be hardened by making these ``ENUM``s and doing much more type-checking.
            These must be valid AMBER keywords or else the simulation will crash. The default keywords that are set when
            this class is initialized are partially overridden by the specific "config" functions detailed below.

        The default dictionary keys and values are as follows:

            - ``imin``       : 0
            - ``ntx``        : 1
            - ``irest``      : 0
            - ``maxcyc``     : 0
            - ``ncyc``       : 0
            - ``dt``         : 0.002
            - ``nstlim``     : 5000
            - ``ntpr``       : 500
            - ``ntwe``       : 500
            - ``ntwr``       : 5000
            - ``ntwx``       : 500
            - ``ntxo``       : 1
            - ``ioutfm``     : 1
            - ``ntf``        : 2
            - ``ntc``        : 2
            - ``cut``        : 8
            - ``igb``        : 0
            - ``tempi``      : 298.15
            - ``tempo``      : 298.15
            - ``ntt``        : 3
            - ``gamma_ln``   : 1.0
            - ``ig``         : -1
            - ``ntp``        : 1
            - ``barostat``   : 2
            - ``ntr``        : ``None``
            - ``restraint_wt``  : ``None``
            - ``restraintmask`` : ``None``
            - ``nmropt``     : 1
            - ``pencut``     : -1
        """
        return self._cntrl

    @cntrl.setter
    def cntrl(self, value: dict):
        self._cntrl = value

    @property
    def ewald(self) -> dict:
        """dict: Additional Ewald summation settings."""
        return self._ewald

    @ewald.setter
    def ewald(self, value: dict):
        self._ewald = value

    @property
    def wt(self) -> dict:
        """dict: "Weight" assigned to various simulation parameters. Written to the "&wt" line in the input file."""
        return self._wt

    @wt.setter
    def wt(self, value: dict):
        self._wt = value

    @property
    def group(self) -> dict:
        """dict: Group specification for restraints applied to multiple atoms.

        .. note::
            I don't recall ever having to this option and this seems distinct from applying restraints to a collection
            of atoms.
        """
        return self._group

    @group.setter
    def group(self, value):
        self._group = value

    @property
    def prefix(self) -> str:
        """
        str: The prefix for file names generated from this simulation.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, new_prefix: str):
        self._prefix = new_prefix
        self.input = new_prefix + ".in"
        self.inpcrd = new_prefix + ".inpcrd"
        self.ref = new_prefix + ".inpcrd"
        self.output = new_prefix + ".out"
        self.restart = new_prefix + ".rst7"
        self.mdinfo = new_prefix + ".mdinfo"
        self.mdcrd = new_prefix + ".nc"
        self.mden = new_prefix + ".mden"

    def __init__(self):

        super().__init__()

        # I/O
        self._restraint_file = None

        # File names
        self.input = self._prefix + ".in"
        self.inpcrd = self._prefix + ".inpcrd"
        self.ref = self._prefix + ".inpcrd"
        self.output = self._prefix + ".out"
        self.restart = self._prefix + ".rst7"
        self.mdinfo = self._prefix + ".mdinfo"
        self.mdcrd = self._prefix + ".nc"
        self.mden = self._prefix + ".mden"

        # Input file cntrl settings (Default = NTP)
        self._cntrl = OrderedDict()
        self._cntrl["imin"] = 0
        self._cntrl["ntx"] = 1
        self._cntrl["irest"] = 0
        self._cntrl["maxcyc"] = 0
        self._cntrl["ncyc"] = 0
        self._cntrl["dt"] = 0.002
        self._cntrl["nstlim"] = 5000
        self._cntrl["ntpr"] = 500
        self._cntrl["ntwe"] = 500
        self._cntrl["ntwr"] = 5000
        self._cntrl["ntwx"] = 500
        self._cntrl["ntxo"] = 1
        self._cntrl["ioutfm"] = 1
        self._cntrl["ntf"] = 2
        self._cntrl["ntc"] = 2
        self._cntrl["cut"] = 8.0
        self._cntrl["igb"] = 0
        self._cntrl["tempi"] = 298.15
        self._cntrl["temp0"] = 298.15
        self._cntrl["ntt"] = 3
        self._cntrl["gamma_ln"] = 1.0
        self._cntrl["ig"] = -1
        self._cntrl["ntp"] = 1
        self._cntrl["barostat"] = 2
        self._cntrl["ntr"] = None
        self._cntrl["restraint_wt"] = None
        self._cntrl["restraintmask"] = None
        self._cntrl["nmropt"] = 1
        self._cntrl["pencut"] = -1

        # Other input file sections
        self._ewald = None
        self._wt = None  # or []
        self._group = None  # or []

    def _config_min(self):
        """
        Configure input settings for minimization (without periodic boundary conditions).
        """
        self.cntrl["imin"] = 1
        self.cntrl["ntx"] = 1
        self.cntrl["irest"] = 0
        self.cntrl["maxcyc"] = 5000
        self.cntrl["ncyc"] = 1000
        self.cntrl["dt"] = 0.0
        self.cntrl["nstlim"] = 0
        self.cntrl["ntpr"] = 100
        self.cntrl["ntwr"] = 5000
        self.cntrl["ntwx"] = 0
        self.cntrl["ntwe"] = 0
        self.cntrl["ntxo"] = 1
        self.cntrl["ntf"] = 1
        self.cntrl["ntc"] = 1
        self.cntrl["ntt"] = 0
        self.cntrl["gamma_ln"] = 0.0
        self.cntrl["ig"] = 0
        self.cntrl["ntp"] = 0
        self.cntrl["barostat"] = 0
        self.mdcrd = None
        self.mden = None

    def _config_md(self):
        """
        Configure input settings for MD.
        """
        self.cntrl["imin"] = 0
        self.cntrl["ntx"] = 1
        self.cntrl["irest"] = 0
        self.cntrl["maxcyc"] = 0
        self.cntrl["ncyc"] = 0
        self.cntrl["dt"] = 0.002
        self.cntrl["nstlim"] = 5000
        self.cntrl["ntpr"] = 500
        self.cntrl["ntwe"] = 500
        self.cntrl["ntwr"] = 5000
        self.cntrl["ntwx"] = 500
        self.cntrl["ntxo"] = 1
        self.cntrl["ioutfm"] = 1
        self.cntrl["ntf"] = 2
        self.cntrl["ntc"] = 2
        self.cntrl["ntt"] = 3
        self.cntrl["gamma_ln"] = 1.0
        self.cntrl["ig"] = -1

    def config_gb_min(self):
        """
        Configure input settings for minimization in continuum solvent.
        """

        self._config_min()
        self.title = "GB Minimization"
        self.cntrl["cut"] = 999.0
        self.cntrl["igb"] = 1

    def config_pbc_min(self):
        """
        Configure input settings for minimization in periodic boundary conditions.
        """
        self._config_min()
        self.title = "PBC Minimization"
        self.cntrl["cut"] = 9.0
        self.cntrl["igb"] = 0

    def config_gb_md(self):
        """
        Configure input settings for MD in default GB.
        """

        self._config_md()
        self.title = "GB MD Simulation"
        self.cntrl["cut"] = 999.0
        self.cntrl["igb"] = 1
        self.cntrl["ntp"] = 0
        self.cntrl["barostat"] = 0

    def config_pbc_md(self):
        """
        Configure input settings for default NTP.
        """

        self._config_md()
        self.title = "PBC MD Simulation"
        self.cntrl["cut"] = 9.0
        self.cntrl["igb"] = 0
        self.cntrl["iwrap"] = 1
        self.cntrl["ntp"] = 1
        self.cntrl["barostat"] = 2

    def _write_dict_to_mdin(self, f, dictionary):
        """
        Write dictionary to file, following AMBER format.

        Parameters
        ----------
        f : TextIO
            File where the dictionary should be written
        dictionary : dict
            Dictionary of values

        """

        for key, val in dictionary.items():
            if val is not None:
                f.write("  {:15s} {:s},\n".format(key + " =", str(val)))

        # Write PLUMED option if preferred over Amber NMR restraints
        if self.plumed_file:
            f.write("  {:15s} {:s},\n".format("plumed = ", str(1)))
            f.write("  {:15s} {:s},\n".format("plumedfile =", f"'{self.plumed_file}'"))

        f.write(" /\n")

    def _write_input_file(self):
        """
        Write the input file specification to file.
        """
        logger.debug("Writing {}".format(self.input))
        with open(os.path.join(self.path, self.input), "w") as mdin:
            mdin.write("{}\n".format(self.title))
            mdin.write(" &cntrl\n")
            self._write_dict_to_mdin(mdin, self.cntrl)

            if self.ewald is not None:
                mdin.write(" &ewald\n")
                self._write_dict_to_mdin(mdin, self.ewald)

            if self.cntrl["nmropt"] == 1:
                if self.wt is not None:
                    for line in self.wt:
                        mdin.write(" " + line + "\n")
                mdin.write(" &wt type = 'END', /\n")

                # Specify Amber NMR file
                if self.restraint_file is not None:
                    mdin.write("DISANG = {}\n".format(self.restraint_file))
                    mdin.write("LISTOUT = POUT\n\n")

            if self.group is not None:
                mdin.write("{:s}".format(self.group))

    def run(self, soft_minimize=False, overwrite=False, fail_ok=False):
        """
        Method to run Molecular Dynamics simulation with AMBER.

        Parameters
        ----------
        soft_minimize: bool, optional, default=False
            Whether to slowly turn on non-bonded interactions so that restraints get enforced first.
        overwrite: bool, optional, default=False
            Whether to overwrite simulation files.
        fail_ok: bool, optional, default=False
            Whether a failing simulation should stop execution of ``pAPRika``.

        """

        if overwrite or not self.check_complete():

            # These settings hardcoded at the moment ... possibly expose for
            # editing in the future
            if soft_minimize:
                # Set a burn in value that is 25% of the way between ncyc and
                # maxcyc
                ncyc = self.cntrl["ncyc"]
                maxcyc = self.cntrl["maxcyc"]
                burn_in = int(float(ncyc) + 0.20 * (float(maxcyc) - float(ncyc)))
                # If the burn_in value is nuts, then just set it to zero
                if burn_in < 0 or burn_in >= maxcyc:
                    burn_in = 0
                # Set an end_soft value that is 75% of way between ncyc and
                # maxcyc
                end_soft = int(float(ncyc) + 0.60 * (float(maxcyc) - float(ncyc)))
                self.wt = [
                    "&wt type = 'NB', istep1=0, istep2={:.0f}, value1 = 0.0, value2=0.0, IINC=50, /".format(
                        burn_in
                    ),
                    "&wt type = 'NB', istep1={:.0f}, istep2={:.0f}, value1 = 0.0, value2=1.0, IINC=50, /".format(
                        burn_in, end_soft
                    ),
                ]

            # Check restraints file
            if self.restraint_file and self.plumed_file:
                raise Exception(
                    "Cannot use both NMR-style and Plumed-style restraints at the same time."
                )

            # Add CUDA_VISIBLE_DEVICES variable to env dict
            if self.gpu_devices is not None:
                os.environ = dict(
                    os.environ, CUDA_VISIBLE_DEVICES=str(self.gpu_devices)
                )

            # Set Plumed kernel library to env dict
            if self.plumed_file is not None:
                self._set_plumed_kernel()

            # _amber_write_input_file(self.path+'/'+self.input, self.min, title='GB Minimization.')
            self._write_input_file()

            if self.cntrl["imin"] == 1:
                logger.info("Running Minimization at {}".format(self.path))
            else:
                logger.info("Running MD at {}".format(self.path))

            # Create executable list for subprocess
            exec_list = self.executable.split() + ["-O", "-p", self.topology]
            if self.ref is not None:
                exec_list += ["-ref", self.ref]
            exec_list += [
                "-c",
                self.inpcrd,
                "-i",
                self.input,
                "-o",
                self.output,
                "-r",
                self.restart,
            ]
            if self.mdcrd is not None:
                exec_list += ["-x", self.mdcrd]
            if self.mdinfo is not None:
                exec_list += ["-inf", self.mdinfo]
            if self.mden is not None:
                exec_list += ["-e", self.mden]

            logger.debug("Exec line: " + " ".join(exec_list))

            # Execute
            amber_output = sp.Popen(
                exec_list,
                cwd=self.path,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                env=os.environ,
            )

            amber_stdout = amber_output.stdout.read().splitlines()
            amber_stderr = amber_output.stderr.read().splitlines()

            # Report any stdout/stderr which are output from execution
            if amber_stdout:
                logger.info("STDOUT received from AMBER execution")
                for line in amber_stdout:
                    logger.info(line)

            if amber_stderr:
                logger.info("STDERR received from AMBER execution")
                for line in amber_stderr:
                    logger.info(line)

            # Check completion status
            if self.cntrl["imin"] == 1 and self.check_complete():
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
        Check for the string "TIMINGS" in ``self.output`` file. If "TIMINGS" is found, then the simulation completed
        successfully, as of AMBER20.

        Parameters
        ----------
        alternate_file : os.PathLike, optional, default=None
            If present, check for "TIMINGS" in this file rather than ``self.output``.

        Returns
        -------
        timings : bool
            True if "TIMINGS" is found in file. False, otherwise.
        """

        # Assume not completed
        complete = False

        if alternate_file:
            output_file = alternate_file
        else:
            output_file = os.path.join(self.path, self.output)

        if os.path.isfile(output_file):
            with open(output_file, "r") as f:
                strings = f.read()
                if " TIMINGS" in strings:
                    complete = True
        if complete:
            logger.debug("{} has TIMINGS".format(output_file))
        else:
            logger.debug("{} does not have TIMINGS".format(output_file))

        return complete


class Gromacs(BaseSimulation, abc.ABC):
    """
    A wrapper that can be used to set GROMACS simulation parameters.

    .. todo ::
        possibly modify this module to use the official python wrapper of GROMACS.

    """

    @property
    def index_file(self) -> str:
        """os.PathLike: GROMACS index file that specifies ``groups`` in the system. This is optional for a GROMACS
        simulation."""
        return self._index_file

    @index_file.setter
    def index_file(self, value: str):
        self._index_file = value

    @property
    def checkpoint(self) -> str:
        """os.PathLike: GROMACS checkpoint file (extension is .cpt) that is parsed to ``GROMACS``."""
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

            gmx mdrun -deffnm ``prefix`` -nt ``nthreads`` -gpu_id ``gpu_devices`` -plumed ``plumed.dat``

        This is useful depending on how GROMACS was compiled, e.g. if GROMACS is compiled with the MPI library the
        you will need to use the command below:

        .. code::

            mpirun -np 6 gmx_mpi mdrun -deffnm ``prefix`` -ntomp 1 -gpu_id 0 -plumed ``plumed.dat``
        """
        return self._custom_mdrun_command

    @custom_mdrun_command.setter
    def custom_mdrun_command(self, value: str):
        self._custom_mdrun_command = value

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
        self._nb_method["nstlist"] = 20
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
        self.thermostat["ref_t"] = (
            self.temperature if self.temperature is not None else 298.15
        )
        self.thermostat["tau_t"] = 0.1
        if self.tc_groups:
            self.thermostat["tc-grps"] = self.tc_groups

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
            self.thermostat["gen_vel"] = "yes"
            self.thermostat["gen_temp"] = self.thermostat["ref_t"]
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

    def run(self, overwrite=False, fail_ok=False):
        """
        Method to run Molecular Dynamics simulation with GROMACS.

        Parameters
        ----------
        overwrite: bool, optional, default=False
            Whether to overwrite simulation files.
        fail_ok: bool, optional, default=False
            Whether a failing simulation should stop execution of ``pAPRika``.
        """
        if overwrite or not self.check_complete():
            # Write MDF input file
            self._write_input_file()

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
            # gmx grompp -f npt.mdp -c nvt.gro -p bcd-hep.top -t nvt.cpt -o npt.tpr -n index.ndx
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
                logger.debug("STDOUT received from GROMACS execution")
                for line in grompp_stdout:
                    logger.debug(line)

            if grompp_stderr:
                logger.debug("STDERR received from GROMACS execution")
                for line in grompp_stderr:
                    logger.debug(line)

            # create executable list for MDRUN
            # gmx_mpi mdrun -v -deffnm nvt -nt 6 -gpu_id 0 -plumed plumed.dat
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
                    mdrun_list += ["-nt", f"{self.nthreads}"]

                # Add gpu id if not already specified in custom
                if (
                    self.gpu_devices is not None
                    and "-gpu_id" not in self.custom_mdrun_command
                ):
                    mdrun_list += ["-gpu_id", f"{self.gpu_devices}"]

                # Add plumed file if not already specified in custom
                if self.plumed_file and "-plumed" not in self.custom_mdrun_command:
                    mdrun_list += ["-plumed", self.plumed_file]

            else:
                mdrun_list += [self.executable, "mdrun", "-deffnm", self.prefix]

                # Add number of threads
                mdrun_list += ["-nt", f"{self.nthreads}"]

                # Add gpu id
                if self.gpu_devices is not None:
                    mdrun_list += ["-gpu_id", f"{self.gpu_devices}"]

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
                logger.debug("STDOUT received from MDRUN execution")
                for line in mdrun_out:
                    logger.debug(line)

            if mdrun_err:
                logger.debug("STDERR received from MDRUN execution")
                for line in mdrun_err:
                    logger.debug(line)

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
