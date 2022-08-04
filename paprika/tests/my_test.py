"""
Tests the restraints utilities.
"""
import logging
from importlib import reload

reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.DEBUG,
)
import os

import numpy as np

try:
    import openmm
    import openmm.app as app
    import openmm.unit as openmm_unit
except ImportError:
    import simtk.openmm as openmm
    import simtk.openmm.app as app
    import simtk.unit as openmm_unit

import parmed as pmd
import pytest
from openff.units import unit as pint_unit

from paprika.restraints.openmm import apply_dat_restraint, apply_positional_restraints
from paprika.restraints.restraints import DAT_restraint, FECalcType, create_window_list
from paprika.restraints.utils import (
    extract_guest_restraints,
    get_bias_potential_type,
    get_restraint_values,
)
from paprika.tests.test_tleap import clean_files

rest = DAT_restraint()
rest.fe_method = FECalcType.DDM
rest.amber_index = True
rest.continuous_apr = False
rest.auto_apr = False
rest.topology = os.path.join(
    os.path.dirname(__file__), "../data/cb6-but/cb6-but-notcentered.pdb"
)
rest.mask1 = ":CB6@O,O2,O4,O6,O8,O10"
rest.mask2 = ":BUT@C3"
rest.attach["target"] = 3.0
rest.attach["num_windows"] = 4
rest.attach["fc_initial"] = 0.0
rest.attach["fc_final"] = 3.0
rest.decouple["fc"] = rest.attach["fc_final"]
rest.decouple["target"] = rest.attach["target"]
rest.decouple["electrostatics"] = {
    "lambda_initial": 1.0,
    "lambda_final": 0.0,
    "num_windows": 6,
}
# rest.decouple["sterics"] = {"lambda_final": 1.0, "num_windows": 6}
rest.release["target"] = 6.0
rest.release["num_windows"] = rest.attach["num_windows"]
rest.release["fc_initial"] = rest.attach["fc_initial"]
rest.release["fc_final"] = rest.attach["fc_final"]
rest.initialize()

print(rest.phase["decouple"]["force_constants"])
print(rest.phase["decouple"]["targets"])
print(rest.phase["decouple"]["electrostatics"])
print(rest.phase["decouple"]["sterics"])

# assert np.allclose(rest.phase["decouple"]["electrostatics"], np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
# assert np.allclose(rest.phase["decouple"]["sterics"], np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
# print(rest.phase["decouple"]["force_constants"].to(force_constant_units).magnitude)
