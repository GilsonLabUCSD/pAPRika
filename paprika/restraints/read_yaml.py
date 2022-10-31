import logging
import re

import yaml
from openff.units import unit as openff_unit

logger = logging.getLogger(__name__)


def read_taproom_yaml(file):
    """
    Read `Taproom <https://github.com/slochower/host-guest-benchmarks>`_ -style YAML-formatted instructions for
    preparing host-guest systems.

    Parameters
    ----------
    file: os.PathLike
        A YAML-formatted file.

    Returns
    -------
    yaml_data: dict
        Dictionary containing simulation setup parameters.

    """

    # Read YAML file
    with open(file, "r") as f:
        yaml_data = yaml.safe_load(f)
    logger.debug(yaml_data)

    # Convert aliases to atom masks
    if "aliases" in yaml_data.keys():
        logger.debug("Dealiasing atom masks...")
        yaml_data = de_alias(yaml_data)

    # Convert all string to OpenFF Quantity
    logger.debug("Converting string to unit.Quantity...")
    convert_string_to_quantity(yaml_data)

    return yaml_data


def multiple_replace(dct, text):
    """
    Create a regular expression to do multiple find and replace.
    """

    # Create a regular expression from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dct.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dct[mo.string[mo.start() : mo.end()]], text)


def de_alias(yaml_data):
    """
    Replace aliased atoms in a ``taproom`` recipe.
    """

    mapping_list = yaml_data["aliases"]
    mapping_dictionary = {}

    for atom_pair in mapping_list:
        mapping_dictionary.update(atom_pair)
    logger.debug(f"Found mapping: {mapping_dictionary}")

    for restraint_type, restraint_type_list in yaml_data["restraints"].items():
        for restraint in restraint_type_list:
            atoms = restraint["restraint"]["atoms"]
            mapped_atoms = multiple_replace(mapping_dictionary, atoms)
            logger.info(f"{atoms} → {mapped_atoms}")
            restraint["restraint"]["atoms"] = mapped_atoms

    if "symmetry_correction" in yaml_data.keys():
        for restraint in yaml_data["symmetry_correction"]["restraints"]:
            atoms = restraint["restraint"]["atoms"]
            mapped_atoms = multiple_replace(mapping_dictionary, atoms)
            logger.info(f"{atoms} → {mapped_atoms}")
            restraint["restraint"]["atoms"] = mapped_atoms

    return yaml_data


def convert_string_to_quantity(yaml_data):
    """
    Convert strings for 'force_constant' and 'targets' to `unit.Quantity`
    """

    def _to_quantity(string, key):
        if string in key:
            value = key[string]
            key[string] = openff_unit.Quantity(value)

    for restraint_type, restraint_type_list in yaml_data["restraints"].items():
        for restraint in restraint_type_list:
            _to_quantity("force_constant", restraint["restraint"])
            _to_quantity("target", restraint["restraint"])
            if "attach" in restraint["restraint"]:
                _to_quantity("force_constant", restraint["restraint"]["attach"])
                _to_quantity("target", restraint["restraint"]["attach"])
            if "pull" in restraint["restraint"]:
                _to_quantity("force_constant", restraint["restraint"]["pull"])
                _to_quantity("target", restraint["restraint"]["pull"])

    if "symmetry_correction" in yaml_data.keys():
        for restraint in yaml_data["symmetry_correction"]["restraints"]:
            _to_quantity("force_constant", restraint["restraint"])
            _to_quantity("target", restraint["restraint"])
