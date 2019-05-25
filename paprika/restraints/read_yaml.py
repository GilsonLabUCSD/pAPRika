import yaml
from pathlib import Path
import os
import logging

logger = logging.getLogger(__name__)

def read_yaml(file):

    with open(file, "r") as f:
        yaml_data = yaml.safe_load(f)
    logger.debug(yaml_data)

    return yaml_data