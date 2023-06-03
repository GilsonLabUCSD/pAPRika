import base64
import importlib
import json
import logging
import os
import traceback
from enum import Enum

import numpy as np
import pytraj as pt
from openff.units import unit as openff_unit
from parmed import Structure
from parmed.amber import AmberParm

from paprika.restraints import BiasPotentialType, DAT_restraint, RestraintType

# https://stackoverflow.com/questions/27909658/json-encoder-and-decoder-for-complex-numpy-arrays
# https://stackoverflow.com/a/24375113/901925
# https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
# https://github.com/openforcefield/openff-evaluator/blob/main/openff/evaluator/utils/serialization.py
logger = logging.getLogger(__name__)


def _type_string_to_object(type_string):
    last_period_index = type_string.rfind(".")

    if last_period_index < 0 or last_period_index == len(type_string) - 1:
        raise ValueError(
            "The type string is invalid - it should be of the form "
            "module_path.class_name: {}".format(type_string)
        )

    type_string_split = type_string.split(".")

    class_object = None
    module_path = None

    while len(type_string_split) > 0:
        class_name = type_string_split.pop(0)

        try:
            if module_path is None:
                module_path = class_name
            else:
                module_path = module_path + "." + class_name

            # First try and treat the current string as a module
            module = importlib.import_module(module_path)
            class_object = module

        except ImportError:
            # If we get an import error, try then to treat the string
            # as the name of a nested class.
            class_object = getattr(class_object, class_name)

    return class_object


class PaprikaEncoder(json.JSONEncoder):
    """Save :class:`DAT_restraint` as JSON by re-encoding :class:`numpy` arrays."""

    def default(self, obj):
        """If input object is a ndarray it will be converted into a dict
        holding dtype, shape and the data, base64 encoded.

        Parameters
        ----------
        obj: str or :class:`parmed.amber.AmberParm`
            Input object to encode.

        Returns
        -------
        json.JSONEncoder
            Encoded object.
        """
        if isinstance(obj, AmberParm):
            logging.info("Encountered AmberParm, returning name.")
            return obj.name
        if isinstance(obj, Structure):
            logging.warning("Encountered Structure, which does not store filename.")
            return ""

        if isinstance(obj, np.ndarray):
            if obj.flags["C_CONTIGUOUS"]:
                obj_data = obj.data
            else:
                cont_obj = np.ascontiguousarray(obj)
                assert cont_obj.flags["C_CONTIGUOUS"]
                obj_data = cont_obj.data
            data_b64 = base64.b64encode(obj_data)
            # obj_data = obj.tolist()
            return dict(
                __ndarray__=data_b64.decode("utf-8"),
                dtype=str(obj.dtype),
                shape=obj.shape,
            )
        elif isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        elif isinstance(obj, openff_unit.Quantity):
            return serialize_quantity(obj)
        elif isinstance(obj, Enum):
            return serialize_enum(obj)

        # Let the base class default method raise the TypeError
        # return json.JSONEncoder(self, obj)
        return super(PaprikaEncoder, self).default(obj)


class PaprikaDecoder(json.JSONDecoder):
    """
    Decodes a previously encoded :class:`numpy.ndarray` with proper shape and `dtype`,
    and for `openff.units.unit.Quantity` variables.
    """

    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(
            self, object_hook=self.custom_object_hook, *args, **kwargs
        )

    def custom_object_hook(self, obj):
        if "__ndarray__" in obj:
            data = base64.b64decode(obj["__ndarray__"])
            return np.frombuffer(data, obj["dtype"]).reshape(obj["shape"])

        if "@type" in obj:
            if obj["@type"] == "openff.units.unit.Quantity":
                return deserialize_quantity(obj)
            elif "paprika.restraints" in obj["@type"]:
                return deserialize_enum(obj)

        return obj


def serialize_quantity(quantity):
    """Serializes a `openff.units.unit.Quantity` into a dictionary of the form
    `{'value': quantity.value_in_unit(quantity.unit), 'unit': quantity.unit}`
    This copied from openff-evaluator.

    Parameters
    ----------
    quantity : openff.units.unit.Quantity
        The quantity to serialize

    Returns
    -------
    dict of str and str
        A dictionary representation of a `openff.units.unit.Quantity`
        with keys of {"value", "unit"}
    """

    value = quantity.magnitude

    return {
        "@type": "openff.units.unit.Quantity",
        "value": value,
        "unit": str(quantity.units),
    }


def deserialize_quantity(serialized):
    """Deserialize a `openff.units.unit.Quantity` from a dictionary.
    This copied from openff-evaluator.

    Parameters
    ----------
    serialized : dict of str and str
        A dictionary representation of a `openff.units.unit.Quantity`
        which must have keys {"value", "unit"}

    Returns
    -------
    openff.units.unit.Quantity
        The deserialized quantity.
    """

    if "@type" in serialized:
        serialized.pop("@type")

    value_unit = openff_unit.dimensionless

    if serialized["unit"] is not None:
        value_unit = openff_unit(serialized["unit"])

    return serialized["value"] * value_unit


def serialize_enum(enum):
    if not isinstance(enum, Enum):
        raise ValueError("{} is not an Enum".format(type(enum)))

    object_type = type(enum)
    qualified_name = object_type.__qualname__
    enum_type = "{}.{}".format(object_type.__module__, qualified_name)

    return {"@type": enum_type, "value": enum.value}


def deserialize_enum(enum_dictionary):
    if "@type" not in enum_dictionary:
        raise ValueError(
            "The serialized enum dictionary must include which type the enum is."
        )

    if "value" not in enum_dictionary:
        raise ValueError("The serialized enum dictionary must include the enum value.")

    enum_type_string = enum_dictionary["@type"]
    enum_value = enum_dictionary["value"]

    enum_class = _type_string_to_object(enum_type_string)

    if not issubclass(enum_class, Enum):
        raise ValueError("<{}> is not an Enum".format(enum_class))

    return enum_class(enum_value)


def save_restraints(restraint_list, filepath="restraints.json"):
    """Save a list of :class:`paprika.restraints.DAT_restraint` to a JSON file.

    Parameters
    ----------
    restraint_list: list
        List of :class:`paprika.restraints.DAT_restraint`.
    filepath: os.PathLike
        The name of the JSON file to write to.
    """
    logging.debug("Saving restraint information as JSON.")
    with open(os.path.join(filepath), "w") as f:
        for restraint in restraint_list:
            dumped = json.dumps(restraint.__dict__, cls=PaprikaEncoder)
            f.write(dumped)
            f.write("\n")


def load_restraints(filepath="restraints.json"):
    """Load pAPRika (`DAT_restraint`) restraints from a JSON file.

    Parameters
    ----------
    filepath: os.PathLike
        The name of the JSON file to load.

    Returns
    -------
    restraints: list
        List of :class:`paprika.restraints.DAT_restraint`.
    """
    logging.debug("Loading restraint information from JSON.")
    with open(os.path.join(filepath), "r") as f:
        json_data = f.read()

    restraint_json = json_data.split("\n")
    restraints = []

    for restraint in restraint_json:
        if restraint == "":
            continue

        loaded = json.loads(restraint, cls=PaprikaDecoder)
        tmp = DAT_restraint()
        tmp.__dict__ = loaded

        properties = [
            "mask1",
            "mask2",
            "mask3",
            "mask4",
            "topology",
            "instances",
            "custom_restraint_values",
            "auto_apr",
            "continuous_apr",
            "attach",
            "pull",
            "release",
            "amber_index",
        ]

        for class_property in properties:
            if f"_{class_property}" in tmp.__dict__.keys():
                tmp.__dict__[class_property] = tmp.__dict__[f"_{class_property}"]

        restraints.append(tmp)

    return restraints


def load_trajectory(window, trajectory, topology, single_topology=False):
    """Load a trajectory (or trajectories) and return a pytraj ``trajectory`` object.

    Parameters
    ----------
    window: str
        The simulation window to analyze
    trajectory: str or list
        The name or names of the trajectory
    topology: str or :class:`parmed.Structure`
        The topology the simulation
    single_topology: bool
        Whether a single topology is read for all windows

    Returns
    -------
    traj: pytraj.Trajectory
        The trajectory of stored as a pytraj object.
    """

    logger.debug("Load trajectories from {}/{}...".format(window, trajectory))
    if isinstance(trajectory, str):
        trajectory_path = os.path.join(window, trajectory)
    elif isinstance(trajectory, list):
        trajectory_path = [os.path.join(window, i) for i in trajectory]
        logger.debug("Received list of trajectories: {}".format(trajectory_path))
    else:
        raise RuntimeError("Trajectory path should be a `str` or `list`.")

    traj = None

    if isinstance(topology, str) and not single_topology:
        if not os.path.isfile(os.path.join(window, topology)):
            raise FileNotFoundError(
                f"Cannot find `topology` file: {os.path.join(window, topology)}"
            )
        logger.debug(f"Loading {os.path.join(window, topology)} and {trajectory_path}")
        try:
            traj = pt.iterload(trajectory_path, os.path.join(window, topology))
        except ValueError as e:
            formatted_exception = traceback.format_exception(None, e, e.__traceback__)
            logger.info(
                f"Failed trying to load {os.path.join(window, topology)} and {trajectory_path}: "
                f"{formatted_exception}"
            )
    elif isinstance(topology, str) and single_topology:
        traj = pt.iterload(trajectory_path, os.path.join(topology))
    else:
        try:
            traj = pt.iterload(trajectory_path, topology)
        except BaseException:
            raise Exception("Tried to load `topology` object directly and failed.")

    logger.debug("Loaded {} frames...".format(traj.n_frames))

    return traj


def read_restraint_data(
    trajectory,
    restraint,
    distance_unit=openff_unit.angstrom,
    angle_unit=openff_unit.degree,
):
    """Given a trajectory and restraint, read the restraint and return the DAT values.

    Parameters
    ----------
    trajectory: :class:`pytraj.trajectory`
        A trajectory, probably loaded by load_trajectory
    restraint: :class:`DAT_restraint`
        The restraint to analyze
    distance_unit: openff.unit.Quantity
        The unit for the returned distance values
    angle_unit: openff.unit.Quantity
        The unit for the returned angle values

    Returns
    -------
    data: :class:`np.array`
        The values for this restraint in this window
    """

    data = None

    if (
        restraint.mask1
        and restraint.mask2
        and not restraint.mask3
        and not restraint.mask4
    ):
        data = openff_unit.Quantity(
            pt.distance(
                trajectory, " ".join([restraint.mask1, restraint.mask2]), image=True
            ),
            units=openff_unit.angstrom,
        ).to(distance_unit)

    elif (
        restraint.mask1 and restraint.mask2 and restraint.mask3 and not restraint.mask4
    ):
        data = openff_unit.Quantity(
            pt.angle(
                trajectory,
                " ".join([restraint.mask1, restraint.mask2, restraint.mask3]),
            ),
            units=openff_unit.degrees,
        ).to(angle_unit)

    elif restraint.mask1 and restraint.mask2 and restraint.mask3 and restraint.mask4:
        data = openff_unit.Quantity(
            pt.dihedral(
                trajectory,
                " ".join(
                    [restraint.mask1, restraint.mask2, restraint.mask3, restraint.mask4]
                ),
            ),
            units=openff_unit.degrees,
        ).to(angle_unit)

    return data
