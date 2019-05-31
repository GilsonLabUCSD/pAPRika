"""
This class is going to contain the simulation setup...
"""

import pkg_resources

from pathlib import Path
from paprika.restraints import static_DAT_restraint, DAT_restraint
from paprika.restraints.read_yaml import read_yaml
from paprika.restraints import DAT_restraint,static_DAT_restraint

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
        self.conformational_restraints = False


        self.directory = "benchmarks"

        if self.backend == "amber":
            raise NotImplementedError

        installed_benchmarks = get_benchmarks()
        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks)
        self.build()
        self.initialize_restraints(host_yaml, guest_yaml)
        self.initialize_calculation()


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

    def build(self):
        """
        This is going to change to accomodate SMIRNOFF99Frosst.

        Returns
        -------
        """
        pass

    def initialize_restraints(self, host_yaml, guest_yaml):


        self.host_yaml = read_yaml(host_yaml)
        self.guest_yaml = read_yaml(guest_yaml)
        windows = [self.host_yaml["calculation"]["windows"]["attach"],
                   self.host_yaml["calculation"]["windows"]["pull"],
                   self.host_yaml["calculation"]["windows"]["release"]
        ]
        structure = self.guest_yaml["complex"]

        # It appears that ParmEd cannot currently read a PosixPath object, so I am going to
        # convert these into strings.

        static_restraints = []
        for restraint in self.host_yaml["calculation"]["restraints"]["static"]:
            r = static_DAT_restraint(restraint_mask_list=restraint["restraint"]["atoms"].split(),
                                                num_window_list=windows,
                                                ref_structure=str(guest_yaml.parent.joinpath(structure)),
                                                force_constant=restraint["restraint"]["force_constant"],
                                                amber_index=False if self.backend == "openmm" else True
                                     )
            static_restraints.append(r)

        # for conformational in self.restraints_dictionary["calculation"]["restraints"]["conformational"]:
        #     raise NotImplementedError
        #
        # for wall in self.restraints_dictionary["calculation"]["restraints"]["wall"]:
        #     raise NotImplementedError

        # for static in self.restraints_dictionary["calculation"]["restraints"]["static"]:
        #     print(dir(static))
        #
        #     # For the static restraints, I think we should only set force constants, not targets!
        #
        # for static in self.restraints_dictionary["calculation"]["restraints"]["static"]:
        #     windows = [self.restraints_dictionary["calculation"]["windows"]["attach"],
        #                self.restraints_dictionary["calculation"]["windows"]["pull"],
        #                self.restraints_dictionary["calculation"]["windows"]["release"]]
        #
        #     print(static)
        #     restraint = static_DAT_restraint(restraint_mask_list=self.restraints_dictionary["calculation"]["restraints"]["static"]["atoms"],
        #                                      num_window_list=windows,
        #                                      ref_structure="data/cb6-but-apr/vac.pdb",
        #                                      force_constant=self.restraints_dictionary["calculation"]["restraints"]["static"]["force_constant"],
        #                                      amber_index=True)


        if self.backend == "amber":

            # save_restraints()
            # create_windows_directory()
            # solvate and/or prepare each window()
            print("AMBER backend")
            print(f"Initialize restraints with host = {self.host}, guest = {self.guest}")
            pass

        if self.backend == "openmm":
            pass

    def initialize_calculation(self):
        print(f"Initialize calculation with host = {self.host}, guest = {self.guest}")
        if self.backend == "amber":
            pass
        if self.backend == "openmm":
            pass

        # Should the end point of this class be a file that contains restraints as JSON and either an AMBER input or serialized OpenMM system?



def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks
