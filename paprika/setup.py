"""
This class is going to contain the simulation setup...
"""


class Setup(object):
    """
    The Setup class provides a wrapper function around the preparation of the host-guest system and the application of restraints.
    """

    def __init__(self):
        self.host = None
        self.guest = None
        self.conformational_restraints = False
        self.solvate = False
        self.backend = "openmm"

        @property
        def value(self):
            return self._value

        @value.setter
        def value(self, v):
            if not (v > 0): raise Exception("value must be greater than zero")
            self._value = v

    def initialize_restraints(self):
        if self.backend == "amber":
            # read_restraints()
            # setup_restraints()
            # save_restraints()
            # create_windows_directory()
            # solvate and/or prepare each window()

            pass

        if self.backend == "openmm":
            pass

    def initialize_calculation(self):
        if self.backend == "amber":
            pass
        if self.backend == "openmm":
            pass

        # Should the end point of this class be a file that contains restraints as JSON and either an AMBER input or serialized OpenMM system?