import importlib.metadata
import logging

import yaml

from paprika.taproom.utils import convert_string_to_quantity, de_alias

logger = logging.getLogger(__name__)


def get_benchmarks():
    """
    Determine the installed ``taproom`` benchmarks.
    """
    installed_benchmarks = {}

    for entry_point in importlib.metadata.entry_points(group="taproom.benchmarks"):
        installed_benchmarks[entry_point.name] = entry_point.load()

    return installed_benchmarks


def read_yaml_schema(file):
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
