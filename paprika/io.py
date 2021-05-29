import base64
import json
import logging as log
import os

import numpy as np
from openff.units import unit
from parmed import Structure
from parmed.amber import AmberParm

from paprika.restraints import DAT_restraint

# https://stackoverflow.com/questions/27909658/json-encoder-and-decoder-for-complex-numpy-arrays
# https://stackoverflow.com/a/24375113/901925
# https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113


def serialize_quantity(quantity):
    """Serializes a openff.units.unit.Quantity into a dictionary of the form
    `{'value': quantity.value_in_unit(quantity.unit), 'unit': quantity.unit}`
    This copied from openff-evaluator.

    Parameters
    ----------
    quantity : openff.evaluator.unit.Quantity
        The quantity to serialize

    Returns
    -------
    dict of str and str
        A dictionary representation of a openff.units.unit.Quantity
        with keys of {"value", "unit"}
    """

    value = quantity.magnitude
    return {
        "@type": "openff.units.unit.Quantity",
        "value": value,
        "unit": str(quantity.units),
    }


def deserialize_quantity(serialized):
    """Deserialize a openff.units.unit.Quantity from a dictionary.
    This copied from openff-evaluator.

    Parameters
    ----------
    serialized : dict of str and str
        A dictionary representation of a openff.units.unit.Quantity
        which must have keys {"value", "unit"}

    Returns
    -------
    openff.units.unit.Quantity
        The deserialized quantity.
    """

    if "@type" in serialized:
        serialized.pop("@type")

    value_unit = unit.dimensionless

    if serialized["unit"] is not None:
        value_unit = unit(serialized["unit"])

    return serialized["value"] * value_unit


class PaprikaEncoder(json.JSONEncoder):
    """Save :class:`DAT_restraint` as JSON by re-encoding :class:`numpy` arrays."""

    def default(self, obj):
        """If input object is an ndarray it will be converted into a dict
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
            log.info("Encountered AmberParm, returning name.")
            return obj.name
        if isinstance(obj, Structure):
            log.warning("Encountered Structure, which does not store filename.")
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
        elif isinstance(obj, unit.Quantity):
            return serialize_quantity(obj)

        # Let the base class default method raise the TypeError
        # return json.JSONEncoder(self, obj)
        return super(PaprikaEncoder, self).default(obj)


class PaprikaDecoder(json.JSONDecoder):
    """
    Decodes a previously encoded :class:`numpy.ndarray` with proper shape and `dtype`,
    and for unit.Quantity variables.
    """

    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        if "__ndarray__" in obj:
            data = base64.b64decode(obj["__ndarray__"])
            return np.frombuffer(data, obj["dtype"]).reshape(obj["shape"])

        if "@type" in obj:
            return deserialize_quantity(obj)

        return obj


def save_restraints(restraint_list, filepath="restraints.json"):
    """Save a list of :class:`paprika.restraints.DAT_restraint` to a JSON file.

    Parameters
    ----------
    restraint_list: list
        List of :class:`paprika.restraints.DAT_restraint`.
    filepath: os.PathLike
        The name of the JSON file to write to.
    """
    log.debug("Saving restraint information as JSON.")
    with open(os.path.join(filepath), "w") as f:
        for restraint in restraint_list:
            dumped = json.dumps(restraint.__dict__, cls=PaprikaEncoder)
            f.write(dumped)
            f.write("\n")


def load_restraints(filepath="restraints.json"):
    """Load restraints from a JSON file.

    Parameters
    ----------
    filepath: os.PathLike
        The name of the JSON file to load.

    Returns
    -------
    restraints: list
        List of :class:`paprika.restraints.DAT_restraint`.
    """
    log.debug("Loading restraint information from JSON.")
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
