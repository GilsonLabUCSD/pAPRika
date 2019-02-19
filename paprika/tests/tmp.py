import os
import re
import shutil

import numpy as np
import parmed as pmd

from paprika import amber
from paprika import restraints
from paprika import tleap
from paprika import analysis
from paprika import align
from paprika.tests import addons


k_cl_pdb = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/k-cl/k-cl.pdb"))

# Build the model in vacuum

sys = tleap.System()
sys.template_lines = [
    "source leaprc.water.tip3p",
    "loadamberparams frcmod.ionsjc_tip3p",
    f"model = loadpdb {k_cl_pdb}"
]

sys.output_path = "tmp"
sys.output_prefix = "k-cl"
sys.pbc_type = "octahedral"
sys.target_waters = 2000
sys.neutralize = False
sys.build()

# Setup the windows

attach_string = "0.00 0.01 0.02 0.04 0.07 0.11 0.15 0.25 0.35 0.55 0.75 1.15 1.60 2.50 5.00 7.50 10.00 12.50 15.00 17.50 20.00 22.50 25.00 27.50 30.00 32.50 35.00 37.50 40.00 42.50 45.00 47.50 50.00 52.50 55.00 57.50 60.00 62.50 65.00 67.50 70.00 72.50 75.00 77.50 80.00 82.50 85.00 87.50 90.00 92.50 95.00 97.50"
attach_fractions = [float(i) / 100 for i in attach_string.split()]
initial_distance = 2.65
pull_distances = np.arange(0.0 + initial_distance, 16.0 + initial_distance, 0.5)

# Setup the single distance restraint

restraint = restraints.DAT_restraint()
restraint.continuous_apr = True
restraint.amber_index = True
restraint.topology = k_cl_pdb
restraint.mask1 = "@K+"
restraint.mask2 = "@Cl-"

restraint.attach["target"] = initial_distance
restraint.attach["fraction_list"] = attach_fractions
restraint.attach["fc_final"] = 5.0
restraint.pull["fc"] = restraint.attach["fc_final"]
restraint.pull["target_list"] = pull_distances
restraint.initialize()

# Create the windows
window_list = restraints.create_window_list([restraint])
for window in window_list:
    os.makedirs(f"tmp/windows/{window}")

for window in window_list:
    with open(f"tmp/windows/{window}/disang.rest", "a") as file:
        string = restraints.amber_restraint_line(restraint, window)
        if string is not None:
            file.write(string)

# Create the coordinates in each window

for window in window_list:
    if window[0] == "a":
        shutil.copy("tmp/k-cl.prmtop", f"tmp/windows/{window}/k-cl.prmtop")
        shutil.copy("tmp/k-cl.rst7", f"tmp/windows/{window}/k-cl.rst7")
    elif window[0] == "p":
        structure = pmd.load_file("tmp/k-cl.prmtop", "tmp/k-cl.rst7",
                          structure = True)
        target_difference = restraint.phase['pull']['targets'][int(window[1:])] - restraint.phase['pull']['targets'][0]
        print(f"In window {window} we will translate the guest {target_difference:0.1f} Angstroms.")
        for atom in structure.atoms:
            if atom.name == "Cl-":
                atom.xz += target_difference
        structure.save(f"tmp/windows/{window}/k-cl.prmtop")
        structure.save(f"tmp/windows/{window}/k-cl.rst7")

# Adjust K/Cl charge from +/- 1.0 to +/- 1.3

for window in window_list:
    structure = pmd.load_file(f"tmp/windows/{window}/k-cl.prmtop",
                              f"tmp/windows/{window}/k-cl.rst7",
                              structure=True)
    for atom in structure.atoms:
        if atom.name == "Cl-":
            atom.charge = -1.3
        elif atom.name == "K+":
            atom.charge = 1.3
    structure.save(f"tmp/windows/{window}/k-cl.prmtop", overwrite=True)
    structure.save(f"tmp/windows/{window}/k-cl.rst7", overwrite=True)

# Minimize

for window in window_list:
    simulation = amber.Simulation()
    simulation.executable = "pmemd.cuda"

    simulation.path = f"tmp/windows/{window}/"
    simulation.prefix = "minimize"

    simulation.inpcrd = "k-cl.rst7"
    simulation.ref = "k-cl.rst7"
    simulation.topology = "k-cl.prmtop"
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
    simulation.ref = "k-cl.rst7"
    simulation.topology = "k-cl.prmtop"
    simulation.restraint_file = "disang.rest"

    simulation.config_pbc_md()
    print(f"Running production in window {window}...")
    simulation.run()


free_energy = analysis.fe_calc()
free_energy.prmtop = "k-cl.prmtop"
free_energy.trajectory = 'production*.nc'
free_energy.path = "tmp/windows"
free_energy.restraint_list = [restraint]
free_energy.collect_data()
free_energy.methods = ['ti-block']
free_energy.ti_matrix = "full"
free_energy.bootcycles = 100
free_energy.compute_free_energy()

free_energy.compute_ref_state_work([
        restraint, None, None, None,
        None, None
])

binding_affinity = -1 * (
free_energy.results["attach"]["ti-block"]["fe"] + \
free_energy.results["pull"]["ti-block"]["fe"] + \
free_energy.results["ref_state_work"]

)

sem = np.sqrt(
free_energy.results["attach"]["ti-block"]["sem"]**2 + \
free_energy.results["pull"]["ti-block"]["sem"]**2
)

print(free_energy.results["attach"]["ti-block"]["fe"])
print(free_energy.results["pull"]["ti-block"]["fe"])
print(free_energy.results["ref_state_work"])

print(f"The binding affinity for K+ (+1.3) and Cl- (-1.3) = {binding_affinity:0.2f} +/- {sem:0.2f} kcal/mol")


def plot_pmf(attach, attach_sem,
             pull, pull_sem,
             release_to_std,
             pull_initial, pull_final):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import colorConverter
    import seaborn as sns

    fig = plt.figure(figsize=(6 * 1.2, 6))
    gs = GridSpec(1, 1, wspace=0.2, hspace=0.5)
    ax1 = plt.subplot(gs[0, 0])

    attach = np.asarray(attach)
    pull = np.asarray(pull)

    attach_range = np.arange(len(attach))
    pull_range = np.arange(attach_range[-1], attach_range[-1] + len(pull))
    analytic_range = [pull_range[-1], pull_range[-1]]

    final_fe = attach[-1] + pull[-1] + release_to_std
    final_sem = np.sqrt(attach_sem[-1] ** 2 + pull_sem[-1] ** 2)

    ax1.errorbar(attach_range, attach, yerr=attach_sem, marker="o",
                 ms=8, markeredgecolor='k', markeredgewidth=1, lw=3,
                 label="Attach")
    ax1.errorbar(pull_range, attach[-1] + pull, yerr=pull_sem, marker="o", ms=8,
                 markeredgecolor='k', markeredgewidth=1, lw=3,
                 label="Pull")
    ax1.errorbar(analytic_range, [attach[-1] + pull[-1], final_fe],
                 yerr=[pull_sem[-1], pull_sem[-1]],
                 label="Analytic")

    ax1.scatter(pull_range[-1], final_fe, c='w', edgecolor='k', lw=2, s=80, zorder=10)
    ax1.annotate(r'${0:2.2f} \pm {1:2.2f}$'.format(final_fe, final_sem),
                 xy=(pull_range[-1] + 2, final_fe), xycoords='data')

    ax1.set_xticks([0, len(attach) - 1,
                    len(attach) - 1, len(attach) - 1 + len(pull) - 1,
                    ])

    ax1.set_xticklabels([0, 1, pull_initial, pull_final])
    va = [0, 0, -0.05, -0.05]
    for t, y in zip(ax1.get_xticklabels(), va):
        t.set_y(y)

    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r'$\lambda$ or distance')
    ax1.set_ylabel('Work (kcal/mol)')
    fig.savefig("tmp.png", bbox_inches="tight")


plot_pmf(free_energy.results["attach"]["ti-block"]["fe_matrix"], free_energy.results["attach"]["ti-block"]["sem_matrix"],
         free_energy.results["pull"]["ti-block"]["fe_matrix"], free_energy.results["pull"]["ti-block"]["sem_matrix"],
         free_energy.results["ref_state_work"],
             5.0, 5.0 + 18.0)
