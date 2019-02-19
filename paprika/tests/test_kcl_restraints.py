import os
import re
import shutil
import json


import numpy as np
import parmed as pmd

from paprika import amber
from paprika import restraints
from paprika import tleap
from paprika import analysis
from paprika import io
from paprika import align
from paprika.tests import addons


# For jupyter
from pathlib import Path

cwd = Path().resolve()
k_cl_pdb = os.path.abspath(os.path.join(cwd, "../data/k-cl/k-cl.pdb"))

# Build the model in vacuum

sys = tleap.System()
sys.template_lines = [
    "source leaprc.water.tip3p",
    "loadamberparams frcmod.ionsjc_tip3p",
    f"model = loadpdb {k_cl_pdb}",
]

sys.output_path = "tmp"
sys.output_prefix = "k-cl"
sys.pbc_type = None
sys.target_waters = None
sys.neutralize = False
sys.build()

# Setup the windows

attach_fractions = np.linspace(0, 1.0, 25)
initial_distance = 2.65
pull_distances = np.linspace(0 + initial_distance, 16.0 + initial_distance, 40)

# Setup the single distance restraint

restraint = restraints.DAT_restraint()
restraint.continuous_apr = True
restraint.amber_index = True
restraint.topology = k_cl_pdb
restraint.mask1 = "@K+"
restraint.mask2 = "@Cl-"

restraint.attach["target"] = initial_distance
restraint.attach["fraction_list"] = attach_fractions
restraint.attach["fc_final"] = 10.0
restraint.pull["fc"] = restraint.attach["fc_final"]
restraint.pull["target_list"] = pull_distances
restraint.initialize()

# Add wall restraint during attachment

wall = restraints.DAT_restraint()
wall.auto_apr = False
wall.amber_index = True
wall.topology = k_cl_pdb
wall.mask1 = "@K+"
wall.mask2 = "@Cl-"

wall.attach["fc_initial"] = 1.0
wall.attach["fc_final"] = 1.0

wall.custom_restraint_values["rk2"] = 1.0
wall.custom_restraint_values["rk3"] = 1.0
wall.custom_restraint_values["r1"] = 0.0
wall.custom_restraint_values["r2"] = 3.5
wall.custom_restraint_values["r3"] = 3.5
wall.custom_restraint_values["r4"] = 999

wall.attach["target"] = 3.5
wall.attach["num_windows"] = len(attach_fractions)

wall.initialize()

from paprika.io import save_restraints

save_restraints([restraint], "k-cl-restraints.json")