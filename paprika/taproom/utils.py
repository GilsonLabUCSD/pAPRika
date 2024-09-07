import logging

from openff.units import unit as openff_unit

from paprika.utils import multiple_replace

logger = logging.getLogger(__name__)


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
