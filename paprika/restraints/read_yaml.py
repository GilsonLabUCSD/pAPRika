import logging
import re

import yaml

logger = logging.getLogger(__name__)

def read_yaml(file):

    with open(file, "r") as f:
        yaml_data = yaml.safe_load(f)
    logger.debug(yaml_data)

    if "aliases" in yaml_data.keys():
        logger.debug("Dealiasing atom masks...")
        yaml_data = de_alias(yaml_data)

    return yaml_data


def multiple_replace(dict, text):
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)


def de_alias(yaml_data):
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
            atoms = restraint["atoms"]
            mapped_atoms = multiple_replace(mapping_dictionary, atoms)
            logger.info(f"{atoms} → {mapped_atoms}")
            restraint["atoms"] = mapped_atoms

    return yaml_data
