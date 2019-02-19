import os
import re
import shutil

import numpy as np
import parmed as pmd


import logging
logging.basicConfig(level=logging.DEBUG)

from paprika import amber
from paprika import restraints
from paprika import tleap
from paprika import analysis
from paprika import align
from paprika.io import load_restraints, NumpyEncoder
import json
from paprika.tests import addons



restraints = load_restraints("k-cl-restraints.json")
restraint = restraints[0]

free_energy = analysis.fe_calc()
free_energy.prmtop = "k-cl-sol.prmtop"
free_energy.trajectory = 'production*.nc'
free_energy.path = "tmp/windows"
free_energy.restraint_list = restraints
free_energy.collect_data()
free_energy.methods = ['ti-block']
free_energy.ti_matrix = "full"
free_energy.bootcycles = 1000
free_energy.compute_free_energy()

free_energy.compute_ref_state_work([restraint, None, None, None, None, None])

with open("k-cl-results.json", "w") as f:
    dumped = json.dumps(free_energy.results, cls=NumpyEncoder)
    f.write(dumped)
