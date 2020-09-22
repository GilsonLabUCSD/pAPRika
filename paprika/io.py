import logging as log
import os
import base64
import json
import numpy as np
from paprika.restraints import DAT_restraint
from paprika.utils import index_from_mask
from parmed.amber import AmberParm
from parmed import Structure


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

        properties = ["mask1", "mask2", "mask3", "mask4", "topology", "instances", "custom_restraint_values",
                      "auto_apr", "continuous_apr", "attach", "pull", "release", "amber_index"]
        for class_property in properties:
            if f"_{class_property}" in tmp.__dict__.keys():
                tmp.__dict__[class_property] = tmp.__dict__[f"_{class_property}"]
        restraints.append(tmp)
    return restraints


def extract_guest_restraints(structure, guest_resname, restraints):
    """
    Utility function to extract the guest restraints from a list of restraints
    and return individual restraints in the form:
        [r, theta, phi, alpha, beta, gamma]

    If there is no restraint applied to a particular reaction coordinate
    a `None` will be inserted.

    This function is useful for parsing guest restraints in analysis when
    computing `ref_state_work`.

    Parameters
    ----------
    structure : parmed.Structure
        parmed structure of the system.
    guest_resname : str
        Residue name of the guest molecule.
    restraints : list
        list of restraints.

    Returns
    -------
    list
        list of guest-specific DAT_restraint().

    Examples
    --------

        >>> free_energy = analysis.fe_calc()
        >>> free_energy.restraint_list = restraints
            ...
        >>> free_energy.compute_ref_state_work(extract_guest_restraints(structure, "BEN", restraints))
        >>> print(free_energy["ref_state_work"])

    """
    guest_resname = guest_resname.upper()

    r = None
    theta = None
    phi = None
    alpha = None
    beta = None
    gamma = None

    for restraint in restraints:

        mask2_residue_name = structure[restraint.mask2].residues[0].name

        # Distance
        if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and not restraint.mask3 and not \
                restraint.mask4:
            r = restraint

        # Angle
        if restraint.mask3 and not restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name

            if "DM2" in restraint.mask1 and "DM1" in restraint.mask2 and guest_resname in mask3_residue_name:
                theta = restraint

            if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and guest_resname in \
                    mask3_residue_name:
                beta = restraint

        # Dihedral
        if restraint.mask4:
            mask3_residue_name = structure[restraint.mask3].residues[0].name
            mask4_residue_name = structure[restraint.mask4].residues[0].name

            if "DM3" in restraint.mask1 and "DM2" in restraint.mask2 and "DM1" in restraint.mask3 and guest_resname \
                    in mask4_residue_name:
                phi = restraint

            if "DM2" in restraint.mask1 and "DM1" in restraint.mask2 and guest_resname in mask3_residue_name \
                    and guest_resname in mask4_residue_name:
                alpha = restraint

            if "DM1" in restraint.mask1 and guest_resname in mask2_residue_name and guest_resname in \
                    mask3_residue_name and guest_resname in mask4_residue_name:
                gamma = restraint

    return [r, theta, phi, alpha, beta, gamma]


def extract_dummy_atoms(structure, serial=True):
    """
    Extract information about dummy atoms from a parmed structure and
    returns the information as a dictionary.

    Parameters
    ----------
    structure : :class:`parmed.structure.Structure`
        The parmed structure object we want to extract from
    serial : bool
        Get indices in serial (starts from 1) or index (starts from 0).
        (NOTE: This wording "serial/index" is a convention from VMD)

    Returns
    -------
    dummy_atoms : dict
        Information about dummy atoms in dictionary form

    Examples
    --------
        >>> dummy_atoms = extract_dummy_atoms(structure)

        main keys: {'DM1', 'DM2', 'DM3'}
        sub keys: 'pos'      - cartesian coordinates (x,y,z)
                  'idx'      - atom indices
                  'idx_type' - type of atom index (serial or index)
                  'mass'     - mass of dummy atom

    """

    dummy_atoms = {"DM1": {}, "DM2": {}, "DM3": {}}

    # coordinates
    dummy_atoms["DM1"]["pos"] = structure[":DM1"].coordinates[0]
    dummy_atoms["DM2"]["pos"] = structure[":DM2"].coordinates[0]
    dummy_atoms["DM3"]["pos"] = structure[":DM3"].coordinates[0]

    # mass
    dummy_atoms["DM1"]["mass"] = [atom.mass for atom in structure[":DM1"].atoms][0]
    dummy_atoms["DM2"]["mass"] = [atom.mass for atom in structure[":DM2"].atoms][0]
    dummy_atoms["DM3"]["mass"] = [atom.mass for atom in structure[":DM3"].atoms][0]

    # atom index
    dummy_atoms["DM1"]["idx"] = index_from_mask(structure, ":DM1", amber_index=serial)[0]
    dummy_atoms["DM2"]["idx"] = index_from_mask(structure, ":DM2", amber_index=serial)[0]
    dummy_atoms["DM3"]["idx"] = index_from_mask(structure, ":DM3", amber_index=serial)[0]

    # index type
    dummy_atoms["DM1"]["idx_type"] = "serial" if serial else "index"
    dummy_atoms["DM2"]["idx_type"] = "serial" if serial else "index"
    dummy_atoms["DM3"]["idx_type"] = "serial" if serial else "index"

    return dummy_atoms
