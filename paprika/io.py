import logging as log
import os
import base64
import json
import numpy as np
import parmed as pmd
from paprika.restraints import DAT_restraint
from parmed.amber import AmberParm


# https://stackoverflow.com/questions/27909658/json-encoder-and-decoder-for-complex-numpy-arrays
# https://stackoverflow.com/a/24375113/901925
# https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113


class NumpyEncoder(json.JSONEncoder):
    """Save DAT_restraints as JSON by re-encoding `numpy` arrays."""

    def default(self, obj):
        """If input object is an ndarray it will be converted into a dict
        holding dtype, shape and the data, base64 encoded.
        """
        if isinstance(obj, AmberParm):
            print("Encountered AmberParm, returning prmtop name.")
            return obj.name
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


        # Let the base class default method raise the TypeError
        # return json.JSONEncoder(self, obj)
        return super(NumpyEncoder, self).default(obj)


def json_numpy_obj_hook(dct):
    """Decodes a previously encoded numpy ndarray with proper shape and dtype.

    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray
    """
    if isinstance(dct, dict) and "__ndarray__" in dct:
        data = base64.b64decode(dct["__ndarray__"])
        return np.frombuffer(data, dct["dtype"]).reshape(dct["shape"])
        # return dct['__ndarray__']
    return dct


def save_restraints(restraint_list, filepath="restraints.json"):
    log.debug("Saving restraint information as JSON.")
    with open(os.path.join(filepath), "w") as f:
        for restraint in restraint_list:
            dumped = json.dumps(restraint.__dict__, cls=NumpyEncoder)
            f.write(dumped)
            f.write("\n")


def load_restraints(filepath="restraints.json"):
    log.debug("Loading restraint information from JSON.")
    with open(os.path.join(filepath), "r") as f:
        json_data = f.read()
    restraint_json = json_data.split("\n")
    restraints = []
    for restraint in restraint_json:
        if restraint == "":
            continue
        loaded = json.loads(restraint, object_hook=json_numpy_obj_hook)
        tmp = DAT_restraint()
        tmp.__dict__ = loaded
        try:
            log.debug("Setting topology from file name.")
            tmp.topology = pmd.load_file(loaded["topology"], structure=True)
        except IOError:
            log.debug(
                "Unable to set topology information after loading from JSON.")
            log.debug("Topology is set to the file name of the topology file.")
            tmp.topology = loaded["topology"]
        restraints.append(tmp)
    return restraints
