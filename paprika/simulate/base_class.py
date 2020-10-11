import logging
import os
from sys import platform

logger = logging.getLogger(__name__)


class BaseSimulation(object):
    """
    Base class for a Molecular Dynamics simulation wrapper.
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
        """str: The Molecular Dynamics executable that will be used.

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
    def thermostat(self) -> dict:
        """dict: Dictionary for the thermostat options."""
        return self._thermostat

    @thermostat.setter
    def thermostat(self, value: dict):
        self._thermostat = value

    @property
    def barostat(self) -> dict:
        """dict: Dictionary for the barostat options."""
        return self._barostat

    @barostat.setter
    def barostat(self, value: dict):
        self._barostat = value

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

        self._temperature = 298.15
        self._pressure = 1.01325
        self._thermostat = {}
        self._barostat = {}

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
