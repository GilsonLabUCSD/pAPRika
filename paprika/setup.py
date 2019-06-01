"""
This class is going to contain the simulation setup...
"""

import pkg_resources

from pathlib import Path
from paprika.tleap import System
from paprika.restraints import static_DAT_restraint, DAT_restraint
from paprika.restraints.restraints import create_window_list
from paprika.restraints.read_yaml import read_yaml
from paprika.io import save_restraints

import logging

logger = logging.getLogger(__name__)


def _get_installed_benchmarks():
    _installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        _installed_benchmarks[entry_point.name] = entry_point.load()
    return _installed_benchmarks


class Setup(object):
    """
    The Setup class provides a wrapper function around the preparation of the host-guest system and the application of restraints.
    """

    def __init__(self, host, guest, backend="openmm"):
        self.host = host
        self.guest = guest
        self.backend = backend

        self.directory = Path("benchmarks").joinpath(self.host).joinpath(self.guest)
        self.directory.mkdir(parents=True, exist_ok=True)

        if self.backend == "amber":
            raise NotImplementedError

        installed_benchmarks = get_benchmarks()
        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks)
        self.benchmark_path = host_yaml.parent
        self.host_yaml = read_yaml(host_yaml)
        self.guest_yaml = read_yaml(guest_yaml)

        self.build_bound()
        self.restraints = self.initialize_restraints()
        self.initialize_calculation()
        self.build_windows()

    def parse_yaml(self, installed_benchmarks):
        """
        Read the YAML recipe for the host and guest.

        Returns
        -------

        """
        try:
            host_yaml = installed_benchmarks["host_guest_pairs"][self.host]["yaml"]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for host: {self.host}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        try:
            guest_yaml = installed_benchmarks["host_guest_pairs"][self.host][self.guest]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        return host_yaml, guest_yaml

    def build_bound(self):
        """
        This is going to change to accommodate SMIRNOFF99Frosst.

        Returns
        -------
        """
        system = System()
        system.output_path = self.directory
        system.output_prefix = f"{self.host}-{self.guest}"
        system.pbc_type = None
        system.neutralize = True
        system.template_lines = [
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
            f"loadamberparams {self.benchmark_path.joinpath(self.host_yaml['frcmod'])}",
            f"loadamberparams {self.benchmark_path.joinpath('dummy.frcmod')}",
            f"{self.host.upper()} = loadmol2 {self.benchmark_path.joinpath(self.host_yaml['structure'])}",
            f"{self.guest.upper()} = loadmol2 {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['structure'])}",
            f"DM1 = loadmol2 {self.benchmark_path.joinpath('dm1.mol2')}",
            f"DM2 = loadmol2 {self.benchmark_path.joinpath('dm2.mol2')}",
            f"DM3 = loadmol2 {self.benchmark_path.joinpath('dm3.mol2')}",
            f"model = loadpdb {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['complex'])}",
        ]
        system.build()

    def initialize_restraints(self):

        windows = [
            self.host_yaml["calculation"]["windows"]["attach"],
            self.host_yaml["calculation"]["windows"]["pull"],
            None,
        ]
        structure = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

        static_restraints = []
        for restraint in self.host_yaml["calculation"]["restraints"]["static"]:
            static = static_DAT_restraint(
                restraint_mask_list=restraint["restraint"]["atoms"].split(),
                num_window_list=windows,
                ref_structure=str(structure),
                force_constant=restraint["restraint"]["force_constant"],
                amber_index=False if self.backend == "openmm" else True,
            )
            static_restraints.append(static)

        conformational_restraints = []
        if self.host_yaml["calculation"]["restraints"]["conformational"]:
            for conformational in self.host_yaml["calculation"]["restraints"][
                "conformational"
            ]:
                raise NotImplementedError
        else:
            logger.debug("Skipping conformational restraints...")

        wall_restraints = []
        if self.host_yaml["calculation"]["restraints"]["wall"]:
            for wall in self.host_yaml["calculation"]["restraints"]["wall"]:
                raise NotImplementedError
        else:
            logger.debug("Skipping wall restraints...")

        guest_restraints = []
        for restraint in self.guest_yaml["restraints"]:
            mask = restraint["restraint"]["atoms"].split()

            guest_restraint = DAT_restraint()
            guest_restraint.auto_apr = True
            guest_restraint.continuous_apr = True
            guest_restraint.amber_index = False if self.backend == "openmm" else True
            guest_restraint.topology = str(structure)
            guest_restraint.mask1 = mask[0]
            guest_restraint.mask2 = mask[1]
            guest_restraint.mask3 = mask[2] if len(mask) > 2 else None
            guest_restraint.mask4 = mask[3] if len(mask) > 3 else None

            guest_restraint.attach["target"] = restraint["restraint"]["attach"][
                "target"
            ]
            guest_restraint.attach["fc_final"] = restraint["restraint"]["attach"][
                "force_constant"
            ]
            guest_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
                "lambda"
            ]["attach"]

            guest_restraint.pull["target_final"] = self.host_yaml["calculation"][
                "target"
            ]["pull"]
            guest_restraint.pull["num_windows"] = windows[1]

            guest_restraint.initialize()
            guest_restraints.append(guest_restraint)

        return static_restraints + conformational_restraints + wall_restraints + guest_restraints

    def build_windows(self):
        window_list = create_window_list(self.restraints)
        for window in window_list:
            self.directory.joinpath("windows").joinpath(window).mkdir(parents=True, exist_ok=True)
        save_restraints(self.directory.joinpath(self.restraints))

        if self.backend == "amber":
            # Write the restraint file in each window.
            raise NotImplementedError

    def initialize_calculation(self):
        print(f"Initialize calculation with host = {self.host}, guest = {self.guest}")
        # Should the end point of this class be a file that contains restraints as JSON and either an AMBER input or serialized OpenMM system?
        pass


def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks
