import logging
import os
import subprocess as sp
from collections import OrderedDict
from sys import platform

logger = logging.getLogger(__name__)


class Simulation(object):
    """
    A wrapper that can be used to set AMBER simulation parameters.
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
        """The AMBER executable that will be used.

        .. note ::
            This could be made safer by making an ``ENUM`` of ``sander``, ``pmemd`` and ``pmemd.cuda``.
        """
        return self._executable

    @executable.setter
    def executable(self, value):
        self._executable = value

    @property
    def gpu_devices(self):
        """A wrapper around the environmental variable ``CUDA_VISIBLE_DEVICES``."""
        return self._gpu_devices

    @gpu_devices.setter
    def gpu_devices(self, value):
        self._gpu_devices = value

    @property
    def phase(self):
        """str: Which phase of the calculation to simulate.

        .. note ::
            This could probably be combined with the ``window`` attribute for simplicity.
        """
        return self._phase

    @phase.setter
    def phase(self, value):
        self._phase = value

    @property
    def window(self):
        """str: Which window of the calculation to simulate."""
        return self._window

    @window.setter
    def window(self, value):
        self._window = value

    @property
    def topology(self):
        """os.PathLike: The topology file to simulate."""
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    @property
    def restraint_file(self):
        """os.PathLike: The file containing NMR-style restraints for AMBER.

        .. note ::
            When running AMBER simulations, you can only use either an AMBER NMR-style
            restraints or a Plumed-style restraints and not both.
        """
        return self._restraint_file

    @restraint_file.setter
    def restraint_file(self, value):
        self._restraint_file = value

    @property
    def plumed_file(self) -> str:
        """os.PathLike: The name of the Plumed-style restraints file for AMBER.

        .. note ::
            When running AMBER simulations, you can only use either an AMBER NMR-style
            restraints or a Plumed-style restraints and not both.
        """
        return self._plumed_file

    @plumed_file.setter
    def plumed_file(self, value: str):
        self._plumed_file = value

    @property
    def plumed_kernel_library(self) -> str:
        """os.PathLike: The file path of the Plumed kernel library if it is different from the default. The default
        name is `libplumedKernel` with its extension determined automatically depending on the operating system and is
        located in `os.environ['CONDA_PREFIX']/lib`. This variable is useful for users who did not install or want to
        use Plumed from the conda repository."""
        return self._plumed_kernel_library

    @plumed_kernel_library.setter
    def plumed_kernel_library(self, value: str):
        self._plumed_kernel_library = value

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

    @property
    def prefix(self):
        """
        The prefix for file names generated from this simulation.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, new_prefix):
        self._prefix = new_prefix
        self.input = new_prefix + ".in"
        self.inpcrd = new_prefix + ".inpcrd"
        self.ref = new_prefix + ".inpcrd"
        self.output = new_prefix + ".out"
        self.restart = new_prefix + ".rst7"
        self.mdinfo = new_prefix + ".mdinfo"
        self.mdcrd = new_prefix + ".nc"
        self.mden = new_prefix + ".mden"

    @property
    def cntrl(self):
        """
        An ordered dictionary of simulation "control" namelist parameters. The main puprose of this class attribute is
        to make it easy to override certain simulation parameters, such as positional restraints on dummy atoms, and
        the inclusion of exclusion of the NMR-style APR restraints.

        As of AMBER18, these are described, in part, in chapter 18 on ``pmemd``.

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
    def cntrl(self, value):
        self._cntrl = value

    @property
    def ewald(self):
        """Additional Ewald summation settings."""
        return self._ewald

    @ewald.setter
    def ewald(self, value):
        self._ewald = value

    @property
    def wt(self):
        """ "Weight" assigned to various simulation parameters. Written to the "&wt" line in the input file."""
        return self._wt

    @wt.setter
    def wt(self, value):
        self._wt = value

    @property
    def group(self):
        """Group specification for restraints applied to multiple atoms.

        .. note::
            I don't recall ever having to this option and this seems distinct from applying restraints to a collection
            of atoms.
        """
        return self._group

    @group.setter
    def group(self, value):
        self._group = value

    def __init__(self):

        # I/O
        self._path = "."
        self._executable = "sander"
        self._gpu_devices = None
        self._phase = None
        self._window = None
        self._topology = "prmtop"
        self._restraint_file = None
        self._plumed_file = None
        self._plumed_kernel_library = None
        self.title = "PBC MD Simulation"
        self.converged = False

        # File names
        self._prefix = "md"
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

    def config_pbc_min(self):
        """
        Configure input settings for minimization in periodic boundary conditions.
        """
        self._config_min()
        self.title = "PBC Minimization"
        self.cntrl["cut"] = 8.0
        self.cntrl["igb"] = 0

    def config_gb_min(self):
        """
        Configure input settings for minimization in continuum solvent.
        """

        self._config_min()
        self.title = "GB Minimization"
        self.cntrl["cut"] = 999.0
        self.cntrl["igb"] = 1

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
        self.cntrl["cut"] = 8.0
        self.cntrl["igb"] = 0
        self.cntrl["iwrap"] = 1
        self.cntrl["ntp"] = 1
        self.cntrl["barostat"] = 2

    def _write_dict_to_mdin(self, f, dictionary):
        """
        Write dictionary to file, following AMBER format.

        Parameters
        ----------
        f : os.PathLike
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

    def _amber_write_input_file(self):
        """
        Write the input file specification to file.

        """
        logger.debug("Writing {}".format(self.input))
        with open(os.path.join(self.path, self.input), "w") as f:
            f.write("{}\n".format(self.title))
            f.write(" &cntrl\n")
            self._write_dict_to_mdin(f, self.cntrl)

            if self.ewald is not None:
                f.write(" &ewald\n")
                self._write_dict_to_mdin(f, self.ewald)

            if self.cntrl["nmropt"] == 1:
                if self.wt is not None:
                    for line in self.wt:
                        f.write(" " + line + "\n")
                f.write(" &wt type = 'END', /\n")

                # Specify Amber NMR file
                if self.restraint_file is not None:
                    f.write("DISANG = {}\n".format(self.restraint_file))
                    f.write("LISTOUT = POUT\n\n")

            if self.group is not None:
                f.write("{:s}".format(self.group))

    def run(self, soft_minimize=False, overwrite=False, fail_ok=False):
        """
        Run minimization.

        Parameters
        ----------
        soft_minimize: bool
            Whether to slowly turn on non-bonded interactions so that restraints get enforced first.
        overwrite: bool
            Whether to overwrite simulation files.
        fail_ok: bool
            Whether a failing simulation should stop execution of ``pAPRika``.
        """

        if overwrite or not self.has_timings():

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

            # _amber_write_input_file(self.path+'/'+self.input, self.min, title='GB Minimization.')
            self._amber_write_input_file()

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

            # Add CUDA_VISIBLE_DEVICES variable to env dict
            if self.gpu_devices is not None:
                os.environ = dict(
                    os.environ, CUDA_VISIBLE_DEVICES=str(self.gpu_devices)
                )

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
            if self.cntrl["imin"] == 1 and self.has_timings():
                logger.info("Minimization completed...")
            elif self.has_timings():
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

    def has_timings(self, alternate_file=None):
        """
        Check for the string "TIMINGS" in ``self.output`` file. If "TIMINGS" is found, then the simulation completed
        successfully, as of AMBER18.

        Parameters
        ----------
        alternate_file : os.PathLike
            If present, check for "TIMINGS" in this file rather than ``self.output``. Default: None

        Returns
        -------
        timings : bool
            True if "TIMINGS" is found in file. False, otherwise.

        """

        # Assume not completed
        timings = False

        if alternate_file:
            output_file = alternate_file
        else:
            output_file = os.path.join(self.path, self.output)

        if os.path.isfile(output_file):
            with open(output_file, "r") as f:
                strings = f.read()
                if " TIMINGS" in strings:
                    timings = True
        if timings:
            logger.debug("{} has TIMINGS".format(output_file))
        else:
            logger.debug("{} does not have TIMINGS".format(output_file))

        return timings
