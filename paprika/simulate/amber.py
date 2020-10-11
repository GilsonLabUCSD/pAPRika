import abc
import logging
import os
import subprocess as sp
from collections import OrderedDict
from enum import Enum

from .base_class import BaseSimulation

logger = logging.getLogger(__name__)


class AMBER(BaseSimulation, abc.ABC):
    """
    A wrapper for setting parameters and running MD with AMBER.
    """

    class Thermostat(Enum):
        """
        An enumeration of the different thermostat implemented in AMBER (option for ``ntt``).
        """

        Off = 0
        Berendsen = 1
        Andersen = 2
        Langevin = 3
        OptIsoNoseHoover = 9
        StochIsoNoseHoover = 10
        Bussi = 11

    class Barostat(Enum):
        """
        An enumeration of the different barostat implemented in AMBER (option for ``barostat``).
        """

        Off = 0
        Berendsen = 1
        MonteCarlo = 2

    class GBModel(Enum):
        """
        An enumeration of the different Generalized Born Implicit Solvent model implemented in AMBER (option for
        ``igb``).
        """

        Off = 0
        HCT = 1
        OBC1 = 2
        OBC2 = 5
        GBn = 7
        GBn2 = 8
        Vacuum = 6

    class BoxScaling(Enum):
        """
        An enumeration of the different PBC scaling options when running constant pressure simulations in AMBER (
        option for ``ntp``).
        """

        Off = 0
        Isotropic = 1
        Anisotropic = 2
        Semiisotropic = 3

    class Periodicity(Enum):
        """
        An enumeration of the different periodicity options in AMBER (option for ``ntb``).
        """

        Off = 0
        ConstVolume = 1
        ConstPressure = 1

    class SHAKE(Enum):
        """
        An enumeration of the different bond constraint options in AMBER (option for ``ntc``).
        """

        Off = 1
        HBonds = 2
        AllBonds = 3

    @property
    def restraint_file(self) -> str:
        """os.PathLike: The file containing NMR-style restraints for `AMBER`.

        .. note ::
            When running `AMBER` simulations, you can only use either an `AMBER` NMR-style
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

        As of `AMBER20`, these are described, in part, in chapter 20 on ``pmemd``.

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
        self._cntrl["ntc"] = self.SHAKE.HBonds.value
        self._cntrl["cut"] = 8.0
        self._cntrl["igb"] = self.GBModel.Off.value
        self._cntrl["tempi"] = self.temperature
        self._cntrl["temp0"] = self.temperature
        self._cntrl["ntt"] = self.Thermostat.Langevin.value
        self._cntrl["gamma_ln"] = 1.0
        self._cntrl["ig"] = -1
        self._cntrl["ntp"] = self.BoxScaling.Isotropic.value
        self._cntrl["barostat"] = self.Barostat.MonteCarlo.value
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
        self.cntrl["ntc"] = self.SHAKE.Off.value
        self.cntrl["ntt"] = self.Thermostat.Off.value
        self.cntrl["gamma_ln"] = 0.0
        self.cntrl["ig"] = 0
        self.cntrl["ntp"] = self.BoxScaling.Off.value
        self.cntrl["barostat"] = self.Barostat.Off.value
        self.mdcrd = None
        self.mden = None

    def _config_md(self, thermostat):
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
        self.cntrl["ntc"] = self.SHAKE.HBonds.value
        self.cntrl["ntt"] = thermostat.value
        self.cntrl["gamma_ln"] = 1.0
        self.cntrl["ig"] = -1

    def config_gb_min(self, gb_model=GBModel.HCT):
        """
        Configure a reasonable input setting for an energy minimization run with implicit solvent. `Users
        can override the parameters set by this method.`

        Parameters
        ----------
        gb_model: :class:`AMBER.GBModel`, default=HCT
            Option to choose different implicit solvent model.
        """

        self._config_min()
        self.title = "GB Minimization"
        self.cntrl["cut"] = 999.0
        self.cntrl["igb"] = gb_model.value

    def config_pbc_min(self):
        """
        Configure a reasonable input setting for an energy minimization run with periodic boundary conditions. `Users
        can override the parameters set by this method.`
        """
        self._config_min()
        self.title = "PBC Minimization"
        self.cntrl["cut"] = 9.0
        self.cntrl["igb"] = self.GBModel.Off.value

    def config_gb_md(self, gb_model=GBModel.HCT, thermostat=Thermostat.Langevin):
        """
        Configure a reasonable input settings for MD with implicit solvent. `Users can override the parameters
        set by this method.`

        Parameters
        ----------
        gb_model: :class:`AMBER.GBModel`, default=HCT
            Option to choose different implicit solvent model.
        thermostat: :class:`AMBER.Thermostat`, default=Langevin
            Option to choose one of six thermostat implemented in NAMD: (1) `Langevin`, (2) `tCouple`,
            (3) `tRescale`, (4) `vRescale`, (5) `tReassign`, (6) `Lowe-Anderson`.
        """

        self._config_md(thermostat)
        self.title = "GB MD Simulation"
        self.cntrl["cut"] = 999.0
        self.cntrl["igb"] = gb_model.value
        self.cntrl["ntp"] = self.BoxScaling.Off.value
        self.cntrl["barostat"] = self.Barostat.Off.value

    def config_pbc_md(
        self,
        ensemble="npt",
        thermostat=Thermostat.Langevin,
        barostat=Barostat.MonteCarlo,
    ):
        """
        Configure a reasonable input setting for a MD run with periodic boundary conditions. `Users can override the
        parameters set by this method.`

        Parameters
        ----------
        ensemble: str, default="npt"
            Configure MD setting to use ``nvt`` or ``npt``.
        thermostat: :class:`AMBER.Thermostat`, default='Langevin'
            Option to choose one of six thermostats implemented in AMBER: (1) `Berendsen`, (2) `Andersen`,
            (3) `Langevin`, (4) `Optimized Isokinetic Nose-Hoover`, (5) `Stochastic Isokinetic Nose-Hoover`,
            (6) `Bussi`.
        barostat: :class:`AMBER.Barostat`, default='MonteCarlo'
            Option to choose one of two barostats implemented in AMBER: (1) `Berendsen` or (2) `Monte Carlo`.
        """
        if ensemble.lower() not in ["nve", "nvt", "npt"]:
            raise Exception(f"Thermodynamic ensemble {ensemble} is not supported.")

        self._config_md(thermostat)
        self.title = "PBC MD Simulation"
        self.cntrl["cut"] = 9.0
        self.cntrl["igb"] = self.GBModel.Off.value
        self.cntrl["iwrap"] = 1

        if ensemble == "nve":
            self.cntrl["ntt"] = self.Thermostat.Off.value
            self.cntrl["ntb"] = self.Periodicity.ConstVolume.value

        elif ensemble == "nvt":
            self.cntrl["ntb"] = self.Periodicity.ConstVolume.value

        elif ensemble == "npt":
            self.cntrl["ntp"] = self.BoxScaling.Isotropic.value
            self.cntrl["ntb"] = self.Periodicity.ConstPressure.value
            self.cntrl["barostat"] = barostat.value

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
        Method to run Molecular Dynamics simulation with `AMBER`.

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
        successfully, as of `AMBER20`.

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
