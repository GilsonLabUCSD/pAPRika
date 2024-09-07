import logging
from enum import Enum

import numpy
import parmed
import pytraj
from openff.units import unit as openff_unit

from paprika import utils
from paprika.utils import check_unit

logger = logging.getLogger(__name__)


class BiasPotentialType(Enum):
    Harmonic = "restraint"
    UpperWall = "upper_wall"
    LowerWall = "lower_wall"


class RestraintType(Enum):
    Distance = "distance"
    Angle = "angle"
    Torsion = "torsion"


class DAT_restraint(object):
    """
    Distance or angle or torsion restraints on atoms in the simulation.
    """

    @property
    def instances(self):
        """A list of ``DAT_restraints`` that have been initialized.

        .. note::
            This should never be called directly and ought to be private.
        """
        return self._instances

    @instances.setter
    def instances(self, value):
        self._instances = value

    @property
    def restraint_type(self) -> RestraintType:
        """str: The type of restraints for this instance: Distance, Angle or Dihedral."""
        return self._restraint_type

    @property
    def topology(self):
        """str or :class:`parmed.Structure`: The topology file used to initialize the restraints. This can be a PDB file
        parsed as string or a :class:`parmed.Structure` object."""
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    @property
    def mask1(self):
        """str: The first atom mask used for the restraint."""
        return self._mask1

    @mask1.setter
    def mask1(self, value):
        self._mask1 = value

    @property
    def mask2(self):
        """str: The second atom mask used for the restraint."""
        return self._mask2

    @mask2.setter
    def mask2(self, value):
        self._mask2 = value

    @property
    def mask3(self):
        """str: The third atom mask used for the restraint. Two atom masks are the minimum requirement to use this
        class."""
        return self._mask3

    @mask3.setter
    def mask3(self, value):
        self._mask3 = value

    @property
    def mask4(self):
        """str: The fourth atom mask used for the restraint. Two atom masks are the minimum requirement to use this
        class."""
        return self._mask4

    @mask4.setter
    def mask4(self, value):
        self._mask4 = value

    @property
    def custom_restraint_values(self):
        """dict: In the case a non-harmonic restraint is desired, the pre-calculated values (``r1``, ``r2``, and so on)
        can be overridden with ones from this dictionary. These values will be directly written to the AMBER
        restraint file, so the keys must be valid AMBER keywords.

        Specifically, for target distance, R:

        • ``R < r1`` **Linear**, with the slope of the "left-hand" parabola at the point ``R = r1``.
        • ``r1 <= R < r2`` **Parabolic**, with restraint energy ``k2(R − r2)²``.
        • ``r2 <= R < r3`` -> ``E = 0``.
        • ``r3 <= R < r4`` **Parabolic**, with restraint energy ``k3(R − r3)²``.
        • ``r4 <= R`` **Linear**, with the slope of the "right-hand" parabola at the point ``R = r4``.

        In the case of AMBER18, this is covered in section 25.1 of the manual.

        For example, a flat bottom restraint can be specified by setting ``r2`` not equal to ``r3``.

        >>> from paprika.restraints import DAT_restraint
        >>> this = DAT_restraint
        >>> this.custom_restraint_values["r2"] = 0.0
        >>> this.custom_restraint_values["r3"] = 12.0

        The values in this dictionary use default AMBER units for distances and force constants.

        This attribute has **no effect** for OpenMM.
        """
        return self._custom_restraint_values

    @custom_restraint_values.setter
    def custom_restraint_values(self, value):
        self._custom_restraint_values = value

    @property
    def auto_apr(self):
        """bool: If ``True``, **(1)** the force constant during the pulling phase will be set to the final force
        constant from the attach phase, and **(2)** the initial target during the pulling phase will be equal to the
        final target from the attach phase. If ``False``, these values must be set manually.
        """
        return self._auto_apr

    @auto_apr.setter
    def auto_apr(self, value):
        if not isinstance(value, bool):
            raise TypeError("The ``auto_apr`` attribute is a boolean.")
        self._auto_apr = value

    @property
    def continuous_apr(self):
        """bool: If ``True``, **(1)** the final window of the attach phase is used as the first window of the pull
        phase, and **(2)** the final window of the pull phase is used as the first window of the release phase.
        """
        return self._continuous_apr

    @continuous_apr.setter
    def continuous_apr(self, value):
        if not isinstance(value, bool):
            raise TypeError("The ``continuous_apr`` attribute is a boolean.")
        self._continuous_apr = value

    @property
    def attach(self):
        """
        dict: Dictionary specifying the APR parameters during the attach phase. The dictionary keys are as follows:

            - ``target``             : The target value for the restraint (mandatory)
            - ``fc_initial``         : The initial force constant (optional)
            - ``fc_final``           : The final force constant (optional)
            - ``num_windows``        : The number of windows (optional)
            - ``fc_increment``       : The force constant increment (optional)
            - ``fraction_increment`` : The percentage of the force constant increment (optional)
            - ``fraction_list``      : The list of force constant percentages (optional)
            - ``fc_list``            : The list of force constants (will be created if not given)

        .. note ::
            This is fragile and this could be hardened by making these ``ENUM`` and doing much more type-checking.
        """

        return self._attach

    @attach.setter
    def attach(self, value):
        self._attach = value

    @property
    def pull(self):
        """
        dict: Dictionary specifying the APR parameters during the pull phase. The dictionary keys are as follows:

            - ``fc``                 : The force constant for the restraint (mandatory)
            - ``target_initial``     : The initial target value (optional)
            - ``target_final``       : The final target value (optional)
            - ``num_windows``        : The number of windows (optional)
            - ``target_increment``   : The target value increment (optional)
            - ``fraction_increment`` : The percentage of the target value increment (optional)
            - ``fraction_list``      : The list of target value percentages (optional)
            - ``target_list``        : The list of target values (will be created if not given)

        .. note ::
            This is fragile and this could be hardened by making these ``ENUM`` and doing much more type-checking.
        """

        return self._pull

    @pull.setter
    def pull(self, value):
        self._pull = value

    @property
    def release(self):
        """
        dict: Dictionary specifying the APR parameters during the release phase. The dictionary keys are as follows:

            - ``target``             : The target value for the restraint (mandatory)
            - ``fc_initial``         : The initial force constant (optional)
            - ``fc_final``           : The final force constant (optional)
            - ``num_windows``        : The number of windows (optional)
            - ``fc_increment``       : The force constant increment (optional)
            - ``fraction_increment`` : The percentage of the force constant increment (optional)
            - ``fraction_list``      : The list of force constant percentages (optional)
            - ``fc_list``            : The list of force constants (will be created if not given)

        .. note ::
            This is fragile and this could be hardened by making these ``ENUM`` and doing much more type-checking.
        """
        return self._release

    @release.setter
    def release(self, value):
        self._release = value

    @property
    def amber_index(self):
        """
        bool: If ``True``, atom indices starts at **1** instead of **0**.
        """
        return self._amber_index

    @amber_index.setter
    def amber_index(self, value):
        self._amber_index = value

    instances = []

    def __init__(self):
        self._restraint_type = None

        self._topology = None
        self._mask1 = None
        self._mask2 = None
        self._mask3 = None
        self._mask4 = None

        # These indices will be automatically populated during :meth:`paprika.restraints.DAT_restraint.initialize`.
        self.index1 = None
        self.index2 = None
        self.index3 = None
        self.index4 = None

        self.group1 = False
        self.group2 = False
        self.group3 = False
        self.group4 = False

        self._custom_restraint_values = {}

        self._auto_apr = False
        self._continuous_apr = True
        self._amber_index = False

        self._attach = {
            "target": None,
            "fc_initial": None,
            "fc_final": None,
            "num_windows": None,
            "fc_increment": None,
            "fraction_increment": None,
            "fraction_list": None,
            "fc_list": None,
        }
        self._pull = {
            "fc": None,
            "target_initial": None,
            "target_final": None,
            "num_windows": None,
            "target_increment": None,
            "fraction_increment": None,
            "fraction_list": None,
            "target_list": None,
        }
        self._release = {
            "target": None,
            "fc_initial": None,
            "fc_final": None,
            "num_windows": None,
            "fc_increment": None,
            "fraction_increment": None,
            "fraction_list": None,
            "fc_list": None,
        }
        self.phase = {
            "attach": {"force_constants": None, "targets": None},
            "pull": {"force_constants": None, "targets": None},
            "release": {"force_constants": None, "targets": None},
        }

        DAT_restraint.instances.append(self)

    def __eq__(self, other):
        """
        Test whether two ``DAT_restraint`` instances are equivalent.
        """
        self_dictionary = self.__dict__
        other_dictionary = other.__dict__
        for dct in [self_dictionary, other_dictionary]:
            # Skip checking topology.
            # It is difficult to check topology for two reasons.
            # One, there are numerical variances in the numeric parameters.
            # We can use `numpy.allcose()` but it needs to be done manually on the
            # nested ParmEd AmberParm or Structure classes.
            # Two, sometimes the dummy atoms seem to be encoded with index `-1`.
            # This can happen if the structure used to setup the restraint
            # does not have the dummy atoms already included.
            # Previously, I checked if the topology was an AmberParm, and tried
            # to replace that with the string representation of the file name,
            # but for some reason, that attribute was not always available.
            # Thus, I simply do not compare the `topology` attribute.
            dct["topology"] = None
        logger.debug(self_dictionary)
        logger.debug(other_dictionary)
        keys = set(self_dictionary.keys()) & set(other_dictionary.keys())
        for key in keys:
            if key != "phase":
                assert self_dictionary[key] == other_dictionary[key]
            else:
                for phs in ["attach", "pull", "release"]:
                    for value in ["force_constants", "targets"]:
                        if (
                            self_dictionary["phase"][phs][value] is None
                            and other_dictionary["phase"][phs][value] is None
                        ):
                            continue
                        else:
                            assert numpy.allclose(
                                self_dictionary["phase"][phs][value],
                                other_dictionary["phase"][phs][value],
                            )
        return True

    @staticmethod
    def _calc_method(phase, restraint_dictionary, method):
        """
        This helper function figures out which values in the restraint dictionary need to be set.

        Parameters
        ----------
        phase: str
            The restraint phase.
        restraint_dictionary: dict
            The restraint dictionary corresponding whose values will be set.
        method: str
            Which "method" will be used to set the remaining values. The "method" refers to which
            combination of restraint inputs are provided -- for example, initial value, final value,
            and number of windows or initial value, increment size, and number of increments -- which
            are defined in :meth:`paprika.restraints.DAT_restraint.initialize`.

            This appears to be a ``str`` but should be an ``int``.

        .. note::
            This could do with some serious sprucing up.

        """

        force_constants = None
        targets = None

        # Attach/Release, Force Constant Method 1
        if phase in ("a", "r") and method == "1":
            force_constants = numpy.linspace(
                restraint_dictionary["fc_initial"],
                restraint_dictionary["fc_final"],
                restraint_dictionary["num_windows"],
            )

        # Attach/Release, Force Constant Method 1a
        elif phase in ("a", "r") and method == "1a":
            force_constants = numpy.linspace(
                0.0,
                restraint_dictionary["fc_final"],
                restraint_dictionary["num_windows"],
            )

        # Attach/Release, Force Constant Method 2
        elif phase in ("a", "r") and method == "2":
            units = restraint_dictionary["fc_initial"].units
            force_constants = (
                numpy.arange(
                    restraint_dictionary["fc_initial"].magnitude,
                    restraint_dictionary["fc_final"].magnitude
                    + restraint_dictionary["fc_increment"].magnitude,
                    restraint_dictionary["fc_increment"].magnitude,
                )
                * units
            )

        # Attach/Release, Force Constant Method 2a
        elif phase in ("a", "r") and method == "2a":
            units = restraint_dictionary["fc_final"].units
            force_constants = (
                numpy.arange(
                    0.0,
                    restraint_dictionary["fc_final"].magnitude
                    + restraint_dictionary["fc_increment"].magnitude,
                    restraint_dictionary["fc_increment"].magnitude,
                )
                * units
            )

        # Attach/Release, Force Constant Method 3
        elif phase in ("a", "r") and method == "3":
            units = restraint_dictionary["fc_final"].units
            force_constants = numpy.asarray(
                [
                    fraction * restraint_dictionary["fc_final"].magnitude
                    for fraction in restraint_dictionary["fraction_list"]
                ]
            )
            force_constants *= units

        # Attach/Release, Force Constant Method 4
        elif phase in ("a", "r") and method == "4":
            fractions = numpy.arange(
                0,
                1.0 + restraint_dictionary["fraction_increment"],
                restraint_dictionary["fraction_increment"],
            )
            # Delete if last element is greater than 1.0
            if fractions[-1] > 1.0:
                fractions = numpy.delete(fractions, -1)

            # Round up last element if it does not equal to 1.0
            if not numpy.isclose(fractions[-1], 1.0):
                fractions[-1] = 1.0

            units = restraint_dictionary["fc_final"].units
            force_constants = numpy.asarray(
                [
                    fraction * restraint_dictionary["fc_final"].magnitude
                    for fraction in fractions
                ]
            )
            force_constants *= units

        # Attach/Release, Force Constant Method 5
        elif phase in ("a", "r") and method == "5":
            units = restraint_dictionary["fc_list"][0].units
            force_constants = numpy.asarray(
                [x.magnitude for x in restraint_dictionary["fc_list"]]
            )
            force_constants *= units

        # Attach/Release, Target Method
        if phase in ("a", "r"):
            units = restraint_dictionary["target"].units
            targets = (
                numpy.asarray(
                    [restraint_dictionary["target"].magnitude] * len(force_constants)
                )
                * units
            )

        # Pull, Target Method 1
        if phase == "p" and method == "1":
            targets = numpy.linspace(
                restraint_dictionary["target_initial"],
                restraint_dictionary["target_final"],
                restraint_dictionary["num_windows"],
            )

        # Pull, Target Method 1a
        elif phase == "p" and method == "1a":
            targets = numpy.linspace(
                0.0,
                restraint_dictionary["target_final"],
                restraint_dictionary["num_windows"],
            )

        # Pull, Target Method 2
        elif phase == "p" and method == "2":
            units = restraint_dictionary["target_initial"].units
            targets = (
                numpy.arange(
                    restraint_dictionary["target_initial"].magnitude,
                    (
                        restraint_dictionary["target_final"]
                        + restraint_dictionary["target_increment"]
                    ).magnitude,
                    restraint_dictionary["target_increment"].magnitude,
                )
                * units
            )

        # Pull, Target Method 2a
        elif phase == "p" and method == "2a":
            units = restraint_dictionary["target_final"].units
            targets = (
                numpy.arange(
                    0.0,
                    restraint_dictionary["target_final"].magnitude
                    + restraint_dictionary["target_increment"].magnitude,
                    restraint_dictionary["target_increment"].magnitude,
                )
                * units
            )

        # Pull, Target Method 3
        elif phase == "p" and method == "3":
            units = restraint_dictionary["target_final"].units
            targets = (
                numpy.asarray(
                    [
                        fraction * restraint_dictionary["target_final"].magnitude
                        for fraction in restraint_dictionary["fraction_list"]
                    ]
                )
                * units
            )

        # Pull, Target Method 4
        elif phase == "p" and method == "4":
            fractions = numpy.arange(
                0,
                1.0 + restraint_dictionary["fraction_increment"],
                restraint_dictionary["fraction_increment"],
            )
            units = restraint_dictionary["target_final"].units
            targets = (
                numpy.asarray(
                    [
                        fraction * restraint_dictionary["target_final"].magnitude
                        for fraction in fractions
                    ]
                )
                * units
            )

        # Pull, Target Method 5
        elif phase == "p" and method == "5":
            units = restraint_dictionary["target_list"][0].units
            targets = (
                numpy.asarray(
                    [x.magnitude for x in restraint_dictionary["target_list"]]
                )
                * units
            )

        # Pull, Force Constant Method
        if phase == "p":
            units = restraint_dictionary["fc"].units
            force_constants = (
                numpy.asarray([restraint_dictionary["fc"].magnitude] * len(targets))
                * units
            )

        if force_constants is None and targets is None:
            logger.error("Unsupported Phase/Method: {} / {}".format(phase, method))
            raise Exception("Unexpected phase/method combination passed to _calc_meth")

        return force_constants, targets

    def initialize(self):
        """
        Automatically set remaining force constants and targets.

        Depending on which values are provided for each phase, a different method will
        be used to determine the list of force constants and targets (below).

        For attach and release, a ``target`` value is required and the method is determined if the
        following values are not ``None``:

            - **Method 1** :  num_windows, fc_initial, fc_final
            - **Method 1a**:  num_windows, fc_final
            - **Method 2** :  fc_increment, fc_initial, fc_final
            - **Method 2a**:  fc_increment, fc_final
            - **Method 3** :  fraction_list, fc_final
            - **Method 4** :  fraction_increment, fc_final
            - **Method 5** :  fc_list

        For pull, a ``fc`` value is required and the method is determined if the
        following values are not ``None``:

            - **Method 1**:   num_windows, target_initial, target_final
            - **Method 1a**:  num_windows, target_final
            - **Method 2**:   target_increment, target_initial, target_final
            - **Method 2a**:  target_increment, target_final
            - **Method 3**:   fraction_list, target_final
            - **Method 4**:   fraction_increment, target_final
            - **Method 5**:   target_list

        .. note ::
            This is unnecessary overengineering.
        """

        # Set default units (Based on Amber)
        energy_unit = openff_unit.kcal / openff_unit.mole
        target_unit = openff_unit.angstrom
        force_constant_unit = energy_unit / openff_unit.angstrom**2
        if self.mask3 or self.mask4:
            target_unit = openff_unit.degrees
            force_constant_unit = energy_unit / openff_unit.radians**2

        # Check attach/release units
        for phase in [self._attach, self._release]:
            for key in ["target", "fc_initial", "fc_final", "fc_increment", "fc_list"]:
                if phase[key] is not None:
                    phase[key] = check_unit(
                        phase[key],
                        base_unit=(
                            target_unit if key == "target" else force_constant_unit
                        ),
                    )

        # Check pull units
        for key in [
            "target_initial",
            "target_final",
            "target_increment",
            "target_list",
            "fc",
        ]:
            if self._pull[key] is not None:
                self._pull[key] = check_unit(
                    self._pull[key],
                    base_unit=force_constant_unit if key == "fc" else target_unit,
                )

        # Check custom restraint units
        if self._custom_restraint_values:
            for key in ["r1", "r2", "r3", "r4", "rk2", "rk3"]:
                if key in self._custom_restraint_values:
                    self._custom_restraint_values[key] = check_unit(
                        self._custom_restraint_values[key],
                        base_unit=(
                            force_constant_unit
                            if key in ["rk2", "rk3"]
                            else target_unit
                        ),
                    )
                else:
                    self._custom_restraint_values[key] = None

        # ------------------------------------ ATTACH ------------------------------------ #
        logger.debug("Calculating attach targets and force constants...")

        # Temporary variables to improve readability
        force_constants = None
        targets = None

        if (
            self.attach["num_windows"] is not None
            and self.attach["fc_final"] is not None
        ):
            if self.attach["fc_initial"] is not None:
                logger.debug("Attach, Method #1")
                force_constants, targets = self._calc_method("a", self.attach, "1")
            else:
                logger.debug("Attach, Method #1a")
                force_constants, targets = self._calc_method("a", self.attach, "1a")

        elif (
            self.attach["fc_increment"] is not None
            and self.attach["fc_final"] is not None
        ):
            if self.attach["fc_initial"] is not None:
                logger.debug("Attach, Method #2")
                force_constants, targets = self._calc_method("a", self.attach, "2")
            else:
                logger.debug("Attach, Method #2a")
                force_constants, targets = self._calc_method("a", self.attach, "2a")

        elif (
            self.attach["fraction_list"] is not None
            and self.attach["fc_final"] is not None
        ):
            logger.debug("Attach, Method #3")
            force_constants, targets = self._calc_method("a", self.attach, "3")

        elif (
            self.attach["fraction_increment"] is not None
            and self.attach["fc_final"] is not None
        ):
            logger.debug("Attach, Method #4")
            force_constants, targets = self._calc_method("a", self.attach, "4")

        elif self.attach["fc_list"] is not None:
            logger.debug("Attach, Method #5")
            force_constants, targets = self._calc_method("a", self.attach, "5")

        elif all(v is None for k, v in self.attach.items()):
            logger.debug("No restraint info set for the attach phase! Skipping...")

        else:
            logger.error(
                "Attach restraint input did not match one of the supported methods..."
            )
            for k, v in self.attach.items():
                logger.debug("{} = {}".format(k, v))
            raise Exception(
                "Attach restraint input did not match one of the supported methods..."
            )

        if force_constants is not None and targets is not None:
            self.phase["attach"]["force_constants"] = force_constants
            self.phase["attach"]["targets"] = targets

        # ------------------------------------ PULL ------------------------------------ #
        logger.debug("Calculating pull targets and force constants...")

        force_constants = None
        targets = None

        if self.auto_apr and self.pull["target_final"] is not None:
            self.pull["fc"] = self.phase["attach"]["force_constants"][-1]
            self.pull["target_initial"] = self.phase["attach"]["targets"][-1]

        if (
            self.pull["num_windows"] is not None
            and self.pull["target_final"] is not None
        ):
            if self.pull["target_initial"] is not None:
                logger.debug("Pull, Method #1")
                force_constants, targets = self._calc_method("p", self.pull, "1")
            else:
                logger.debug("Pull, Method #1a")
                force_constants, targets = self._calc_method("p", self.pull, "1a")

        elif (
            self.pull["target_increment"] is not None
            and self.pull["target_final"] is not None
        ):
            if self.pull["target_initial"] is not None:
                logger.debug("Pull, Method #2")
                force_constants, targets = self._calc_method("p", self.pull, "2")
            else:
                logger.debug("Pull, Method #2a")
                force_constants, targets = self._calc_method("p", self.pull, "2a")

        elif (
            self.pull["fraction_list"] is not None
            and self.pull["target_final"] is not None
        ):
            logger.debug("Pull, Method #3")
            force_constants, targets = self._calc_method("p", self.pull, "3")

        elif (
            self.pull["fraction_increment"] is not None
            and self.pull["target_final"] is not None
        ):
            logger.debug("Pull, Method #4")
            force_constants, targets = self._calc_method("p", self.pull, "4")

        elif self.pull["target_list"] is not None:
            logger.debug("Pull, Method #5")
            force_constants, targets = self._calc_method("p", self.pull, "5")

        elif all(v is None for k, v in self.pull.items()):
            logger.debug("No restraint info set for the pull phase! Skipping...")

        else:
            logger.error(
                "Pull restraint input did not match one of the supported methods..."
            )
            for k, v in self.pull.items():
                logger.debug("{} = {}".format(k, v))
            raise Exception(
                "Pull restraint input did not match one of the supported methods..."
            )

        if force_constants is not None and targets is not None:
            self.phase["pull"]["force_constants"] = force_constants
            self.phase["pull"]["targets"] = targets

        # ------------------------------------ RELEASE ------------------------------------ #
        logger.debug("Calculating release targets and force constants...")

        force_constants = None
        targets = None

        # I don't want auto_apr to make release restraints, unless I'm sure the user wants them.
        # I'm gonna assume that specifying self.attach['fc_final'] indicates you want it,
        # although this weakens the whole purpose of auto_apr.

        if self.auto_apr and self.release["fc_final"] is not None:
            self.release["target"] = self.phase["pull"]["targets"][-1]
            for key in [
                "fc_final",
                "fc_initial",
                "num_windows",
                "fc_increment",
                "fraction_increment",
                "fraction_list",
                "fc_list",
            ]:
                if self.attach[key] is not None and self.release[key] is None:
                    self.release[key] = self.attach[key]

        if (
            self.release["num_windows"] is not None
            and self.release["fc_final"] is not None
        ):
            if self.release["fc_initial"] is not None:
                logger.debug("Release, Method #1")
                force_constants, targets = self._calc_method("r", self.release, "1")
            else:
                logger.debug("Release, Method #1a")
                force_constants, targets = self._calc_method("r", self.release, "1a")

        elif (
            self.release["fc_increment"] is not None
            and self.release["fc_final"] is not None
        ):
            if self.release["fc_initial"] is not None:
                logger.debug("Release, Method #2")
                force_constants, targets = self._calc_method("r", self.release, "2")
            else:
                logger.debug("Release, Method #2a")
                force_constants, targets = self._calc_method("r", self.release, "2a")

        elif (
            self.release["fraction_list"] is not None
            and self.release["fc_final"] is not None
        ):
            logger.debug("Release, Method #3")
            force_constants, targets = self._calc_method("r", self.release, "3")

        elif (
            self.release["fraction_increment"] is not None
            and self.release["fc_final"] is not None
        ):
            logger.debug("Release, Method #4")
            force_constants, targets = self._calc_method("r", self.release, "4")

        elif self.release["fc_list"] is not None:
            logger.debug("Release, Method #5")
            force_constants, targets = self._calc_method("r", self.release, "5")

        elif all(v is None for k, v in self.release.items()):
            logger.debug("No restraint info set for the release phase! Skipping...")

        else:
            logger.error(
                "Release restraint input did not match one of the supported methods..."
            )
            for k, v in self.release.items():
                logger.debug("{} = {}".format(k, v))
            raise Exception(
                "Release restraint input did not match one of the supported methods..."
            )

        if force_constants is not None and targets is not None:
            self.phase["release"]["force_constants"] = force_constants
            self.phase["release"]["targets"] = targets

        # ----------------------------------- WINDOWS ------------------------------------ #

        for phase in ["attach", "pull", "release"]:
            if self.phase[phase]["targets"] is not None:
                window_count = len(self.phase[phase]["targets"])
                logger.debug("Number of {} windows = {}".format(phase, window_count))
            else:
                logger.debug(
                    "This restraint will be skipped in the {} phase".format(phase)
                )

        # ---------------------------------- ATOM MASKS ---------------------------------- #
        logger.debug("Assigning atom indices...")
        self.index1 = utils.index_from_mask(self.topology, self.mask1, self.amber_index)
        self.index2 = utils.index_from_mask(self.topology, self.mask2, self.amber_index)
        self._restraint_type = RestraintType.Distance
        if self.mask3:
            self.index3 = utils.index_from_mask(
                self.topology, self.mask3, self.amber_index
            )
            self._restraint_type = RestraintType.Angle
        else:
            self.index3 = None
        if self.mask4:
            self.index4 = utils.index_from_mask(
                self.topology, self.mask4, self.amber_index
            )
            self._restraint_type = RestraintType.Torsion
        else:
            self.index4 = None
        # If any `index` has more than one atom, mark it as a group restraint.
        if self.mask1 and len(self.index1) > 1:
            self.group1 = True

        if self.mask2 and len(self.index2) > 1:
            self.group2 = True

        if self.mask3 and len(self.index3) > 1:
            self.group3 = True

        if self.mask4 and len(self.index4) > 1:
            self.group4 = True


def static_DAT_restraint(
    restraint_mask_list,
    num_window_list,
    ref_structure,
    force_constant,
    continuous_apr=True,
    amber_index=False,
):
    """
    Create a restraint whose value does not change during a calculation.

    Parameters
    ----------
    restraint_mask_list: list
        A list of masks for which this restraint applies.
    num_window_list: list
        A list of windows during which this restraint will be applied, which should be in the form: [attach windows,
        pull windows, release windows].
    ref_structure: os.PathLike or :class:`parmed.Structure`
        The reference structure that is used to determine the initial, **static** value for this restraint.
    force_constant: float or openff.units.unit.Quantity
        The force constant for this restraint. If float, the number will be transformed to kcal/mol/[Angstrom,radians].
    continuous_apr: bool, optional
        Whether this restraint uses ``continuous_apr``. This must be consistent with existing restraints.
    amber_index: bool, optional
        Whether the atom indices for the restraint should be AMBER-style (starts at 1) or not.

    Returns
    -------
    rest: :class:`paprika.restraints.DAT_restraint`
        A static restraint.

    """

    # Check num_window_list
    if len(num_window_list) != 3:
        raise ValueError(
            "The num_window_list needs to contain three integers corresponding to the number of windows in the "
            "attach, pull, and release phase, respectively "
        )

    rest = DAT_restraint()
    rest.continuous_apr = continuous_apr
    rest.amber_index = amber_index

    # Amber format
    if isinstance(ref_structure, parmed.amber._amberparm.AmberParm):
        reference_trajectory = pytraj.load_parmed(ref_structure, traj=True)
        rest.topology = ref_structure
    # String
    elif isinstance(ref_structure, str):
        reference_trajectory = pytraj.iterload(ref_structure, traj=True)
        rest.topology = parmed.load_file(ref_structure, structure=True)
    else:
        raise TypeError(
            "static_DAT_restraint does not support the type associated with ref_structure:"
            + type(ref_structure)
        )

    # Get atom masks
    rest.mask1 = restraint_mask_list[0]
    rest.mask2 = restraint_mask_list[1]
    if len(restraint_mask_list) >= 3:
        rest.mask3 = restraint_mask_list[2]
    if len(restraint_mask_list) == 4:
        rest.mask4 = restraint_mask_list[3]

    # Force constant - convert to openff.units.unit.Quantity
    force_constant = check_unit(
        force_constant,
        base_unit=(
            openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
            if len(restraint_mask_list) == 2
            else openff_unit.kcal / openff_unit.mole / openff_unit.radians**2
        ),
    )

    # Target value
    rest_type = "distance"
    mask_string = " ".join(restraint_mask_list)
    if len(restraint_mask_list) == 2:
        # Distance restraint
        if reference_trajectory.top.has_box():
            target = pytraj.distance(reference_trajectory, mask_string, image=True)[0]
            logger.debug("Calculating distance with 'image = True' ...")
        else:
            target = pytraj.distance(reference_trajectory, mask_string, image=False)[0]
            logger.debug("Calculating distance with 'image = False' ...")
        target = openff_unit.Quantity(target, units=openff_unit.angstrom)

    elif len(restraint_mask_list) == 3:
        # Angle restraint
        target = pytraj.angle(reference_trajectory, mask_string)[0]
        target = openff_unit.Quantity(target, units=openff_unit.degrees)

        rest_type = "angle"

    elif len(restraint_mask_list) == 4:
        # Dihedral restraint
        target = pytraj.dihedral(reference_trajectory, mask_string)[0]
        target = openff_unit.Quantity(target, units=openff_unit.degrees)

        rest_type = "dihedral"

    else:
        raise IndexError(
            f"The number of masks -- {len(restraint_mask_list)} -- is not 2, 3, or 4 and thus is not one of the "
            f"supported types: distance, angle, or dihedral."
        )

    if numpy.isnan(target):
        raise ValueError(
            f"Target for {rest_type} restraint is a NaN. Check topology file."
        )

    # Attach phase
    if num_window_list[0] is not None and num_window_list[0] != 0:
        rest.attach["target"] = target
        rest.attach["fc_initial"] = force_constant
        rest.attach["fc_final"] = force_constant
        rest.attach["num_windows"] = num_window_list[0]

    # Pull phase
    if num_window_list[1] is not None and num_window_list[1] != 0:
        rest.pull["fc"] = force_constant
        rest.pull["target_initial"] = target
        rest.pull["target_final"] = target
        rest.pull["num_windows"] = num_window_list[1]

    # Release phase
    if num_window_list[2] is not None and num_window_list[2] != 0:
        rest.release["target"] = target
        rest.release["fc_initial"] = force_constant
        rest.release["fc_final"] = force_constant
        rest.release["num_windows"] = num_window_list[2]

    rest.initialize()

    return rest


def check_restraints(restraint_list, create_window_list=False):
    """
    Do basic tests to ensure a list of ``DAT_restraint``s are consistent.
    This function is also overloaded to create window lists as well, which is not ideal.

    Parameters
    ----------
    restraint_list: list
        A list of restraints.
    create_window_list: bool, optional, default=False
        Whether to use the restraints to create windows for the calculation.

    """

    if all(restraint.continuous_apr is True for restraint in restraint_list):
        logger.debug('All restraints are "continuous_apr" style.')
        all_continuous_apr = True
    elif all(restraint.continuous_apr is False for restraint in restraint_list):
        logger.debug('All restraints are not "continuous_apr" style.')
        all_continuous_apr = False
    else:
        raise ValueError(
            "All restraints must have the same setting for ``.continuous_apr``."
        )

    window_list = []
    phases = ["attach", "pull", "release"]
    for phase in phases:
        win_counts = []
        for restraint in restraint_list:
            if restraint.phase[phase]["targets"] is not None:
                win_counts.append(len(restraint.phase[phase]["targets"]))
            else:
                win_counts.append(0)
        max_count = numpy.max(win_counts)

        if max_count > 999:
            logger.info("Window name zero padding only applied up to 999.")

        # For each restraint, make sure the number of windows is either 0 (the restraint
        # is not active) or equal to the maximum number of windows for any
        # restraint.
        if all(count == 0 or count == max_count for count in win_counts):
            if max_count > 0:
                # `continuous_apr` during attach means that the final attach window
                # should be skipped and replaced with `p000`. `continuous_apr` during
                # release means that `r000` should be skipped and replaced with the
                # final pull window.

                if phase == "attach" and all_continuous_apr:
                    window_list += [
                        phase[0] + str("{:03.0f}".format(val))
                        for val in numpy.arange(0, max_count - 1, 1)
                    ]
                elif phase == "attach" and not all_continuous_apr:
                    window_list += [
                        phase[0] + str("{:03.0f}".format(val))
                        for val in numpy.arange(0, max_count, 1)
                    ]
                elif phase == "pull":
                    window_list += [
                        phase[0] + str("{:03.0f}".format(val))
                        for val in numpy.arange(0, max_count, 1)
                    ]
                elif phase == "release" and all_continuous_apr:
                    window_list += [
                        phase[0] + str("{:03.0f}".format(val))
                        for val in numpy.arange(1, max_count, 1)
                    ]
                elif phase == "release" and not all_continuous_apr:
                    window_list += [
                        phase[0] + str("{:03.0f}".format(val))
                        for val in numpy.arange(0, max_count, 1)
                    ]
        else:
            logger.error(
                "Restraints have unequal number of windows during the {} phase.".format(
                    phase
                )
            )
            logger.debug("Window counts for each restraint are as follows:")
            logger.debug(win_counts)
            raise Exception(
                "Restraints have unequal number of windows during the {} "
                "phase.".format(phase)
            )

    # Check that each restraint have at least two atoms
    for restraint in restraint_list:
        if not restraint.index1 or not restraint.index2:
            raise Exception("There must be at least two atoms in a restraint.")

    logger.info("Restraints appear to be consistent")

    if create_window_list:
        return window_list
