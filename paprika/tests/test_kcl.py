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

# Create the windows
window_list = restraints.create_window_list([restraint])
for window in window_list:
    if os.path.exists(f"tmp/windows/{window}"):
        # shutil.rmtree(f"tmp/windows/{window}")
        continue
    else:
        os.makedirs(f"tmp/windows/{window}")

for window in window_list:
    with open(f"tmp/windows/{window}/disang.rest", "a") as file:
        if window[0] == "a":
            for r in [restraint, wall]:
                string = restraints.amber_restraint_line(r, window)
                if string is not None:
                    file.write(string)
        else:
            string = restraints.amber_restraint_line(restraint, window)
            file.write(string)

for window in window_list:
    if window[0] == "a":
        structure = pmd.load_file("tmp/k-cl.prmtop", "tmp/k-cl.rst7", structure=True)
        for atom in structure.atoms:
            if atom.name == "Cl-":
                atom.xz = 2.65
        structure.save(f"tmp/windows/{window}/k-cl.prmtop", overwrite=True)
        structure.save(f"tmp/windows/{window}/k-cl.rst7", overwrite=True)

    elif window[0] == "p":
        structure = pmd.load_file("tmp/k-cl.prmtop", "tmp/k-cl.rst7", structure=True)
        target_difference = (
            restraint.phase["pull"]["targets"][int(window[1:])]
            - restraint.phase["pull"]["targets"][0]
        )
        print(
            f"In window {window} we will translate the guest {target_difference:0.1f} Angstroms."
        )
        for atom in structure.atoms:
            if atom.name == "Cl-":
                atom.xz += target_difference
        structure.save(f"tmp/windows/{window}/k-cl.prmtop", overwrite=True)
        structure.save(f"tmp/windows/{window}/k-cl.rst7", overwrite=True)

# Adjust K/Cl charge from +/- 1.0 to +/- 1.3

for window in window_list:
    structure = pmd.load_file(
        f"tmp/windows/{window}/k-cl.prmtop",
        f"tmp/windows/{window}/k-cl.rst7",
        structure=True,
    )
    for atom in structure.atoms:
        if atom.name == "Cl-":
            atom.charge = -1.3
        elif atom.name == "K+":
            atom.charge = 1.3
    structure.save(f"tmp/windows/{window}/k-cl.prmtop", overwrite=True)
    structure.save(f"tmp/windows/{window}/k-cl.rst7", overwrite=True)

# Solvate in each window...

for window in window_list:
    print(f"Solvating window {window}...")

    if os.path.exists(f"tmp/windows/{window}/k-cl-sol.prmtop"):
        print("Skipping...")
        continue


    structure = pmd.load_file(
        f"tmp/windows/{window}/k-cl.prmtop", f"tmp/windows/{window}/k-cl.rst7"
    )

    if not os.path.exists(f"tmp/windows/{window}/k-cl.pdb"):
        structure.save(f"tmp/windows/{window}/k-cl.pdb")

    system = tleap.System()
    system.output_path = os.path.join("tmp", "windows", window)
    system.output_prefix = "k-cl-sol"

    system.target_waters = 2000
    system.neutralize = False
    system.template_lines = ["source leaprc.water.tip3p", "model = loadpdb k-cl.pdb"]
    system.build()


# Minimize

for window in window_list:
    simulation = amber.Simulation()
    simulation.executable = "pmemd.cuda"

    simulation.path = f"tmp/windows/{window}/"
    simulation.prefix = "minimize"

    simulation.inpcrd = "k-cl-sol.rst7"
    simulation.ref = "k-cl-sol.rst7"
    simulation.topology = "k-cl-sol.prmtop"
    simulation.restraint_file = "disang.rest"

    simulation.config_pbc_min()
    simulation.cntrl["ntr"] = 1
    simulation.cntrl["restraint_wt"] = 50.0
    simulation.cntrl["restraintmask"] = "':1-2'"
    print(f"Running minimization in window {window}...")
    simulation.run()

# Simulate

for window in window_list:
    simulation = amber.Simulation()
    simulation.executable = "pmemd.cuda"

    simulation.path = f"tmp/windows/{window}/"
    simulation.prefix = "production"

    simulation.inpcrd = "minimize.rst7"
    simulation.ref = "k-cl-sol.rst7"
    simulation.topology = "k-cl-sol.prmtop"
    simulation.restraint_file = "disang.rest"

    simulation.config_pbc_md()
    # 20 ns per window
    simulation.cntrl["nstlim"] = 10000000

    print(f"Running production in window {window}...")
    # simulation.run()
    simulation._amber_write_input_file()
    
# free_energy = analysis.fe_calc()
# free_energy.prmtop = "k-cl-sol.prmtop"
# free_energy.trajectory = "production*.nc"
# free_energy.path = "tmp/windows"
# free_energy.restraint_list = [restraint]
# free_energy.collect_data()
# free_energy.methods = ["ti-block", "mbar-block"]
# free_energy.ti_matrix = "full"
# free_energy.bootcycles = 100
# free_energy.compute_free_energy()
#
# # free_energy.results
# free_energy.compute_ref_state_work([restraint, None, None, None, None, None])
#
# with open("./tmp/results.json", "w") as f:
#     dumped = json.dumps(free_energy.results, cls=io.NumpyEncoder)
#     f.write(dumped)
#
# binding_affinity = -1 * (
#     free_energy.results["attach"]["ti-block"]["fe"]
#     + free_energy.results["pull"]["ti-block"]["fe"]
#     + free_energy.results["ref_state_work"]
# )
#
# sem = np.sqrt(
#     free_energy.results["attach"]["ti-block"]["sem"] ** 2
#     + free_energy.results["pull"]["ti-block"]["sem"] ** 2
# )
#
# print(free_energy.results["attach"]["ti-block"]["fe"])
# print(free_energy.results["pull"]["ti-block"]["fe"])
# print(free_energy.results["ref_state_work"])
#
# print(
#     f"The binding affinity for K+ (+1.3) and Cl- (-1.3) = {binding_affinity:0.2f} +/- {sem:0.2f} kcal/mol"
# )
