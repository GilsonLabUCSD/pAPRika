import os as os
import numpy as np

import parmed as pmd
import pytraj as pt

import logging
import datetime as dt

d_date = dt.datetime.now()
logging.basicConfig(
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.DEBUG,
)
logging.info("Started logging...")

import paprika

print(paprika.__version__)

from paprika.restraints import static_DAT_restraint
from paprika.restraints import DAT_restraint
from paprika.restraints import create_window_list
from paprika.analysis import fe_calc

prefix = os.path.join(
    "/home",
    "dslochower",
    "kirkwood",
    "projects",
    "smirnoff-host-guest-simulations",
    "systems",
    "a-bam-p",
    "confirm-original",
)
hg = pmd.load_file(
    os.path.join(prefix, "a000", "full.hmr.topo"),
    os.path.join(prefix, "a000", "full.crds"),
    structure=True,
)

dummy_anchors = [":1", ":2", ":3"]
host_anchors = [":4@O3", ":6@C1", ":8@C6"]
guest_anchors = [":10@C4", ":10@N1"]

attach_string = (
    "0.00 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"
)
attach_fractions = [float(i) / 100 for i in attach_string.split()]

pull_string = (
    "0.00 0.40 0.80 1.20 1.60 2.00 2.40 2.80 3.20 3.60 4.00 4.40 4.80 5.20 5.60 6.00 6.40 "
    "6.80 7.20 7.60 8.00 8.40 8.80 9.20 9.60 10.00 10.40 10.80 11.20 11.60 12.00 12.40 12.80 "
    "13.20 13.60 14.00 14.40 14.80 15.20 15.60 16.00 16.40 16.80 17.20 17.60 18.00"
)
pull_distances = [float(i) + 6.00 for i in pull_string.split()]

release_fractions = attach_fractions[::-1]

windows = [len(attach_fractions), len(pull_distances), len(release_fractions)]
print(f"There are {windows} windows in this attach-pull-release calculation.")

static_restraint_atoms = [
    [dummy_anchors[0], host_anchors[0]],
    [dummy_anchors[1], dummy_anchors[0], host_anchors[0]],
    [dummy_anchors[2], dummy_anchors[1], dummy_anchors[0], host_anchors[0]],
    [dummy_anchors[0], host_anchors[0], host_anchors[1]],
    [dummy_anchors[1], dummy_anchors[0], host_anchors[0], host_anchors[1]],
    [dummy_anchors[0], host_anchors[0], host_anchors[1], host_anchors[2]],
]

static_restraint_distance_fc = 5.0
static_restraint_angle_fc = 100.0

guest_restraint_atoms = [
    [dummy_anchors[0], guest_anchors[0]],
    [dummy_anchors[1], dummy_anchors[0], guest_anchors[0]],
    [dummy_anchors[0], guest_anchors[0], guest_anchors[1]],
]

guest_restraint_targets = [6.0, 180.0, 180.0]
guest_restraint_target_final = [24.0, 180.0, 180.0]
guest_restraint_distance_fc = 5.0
guest_restraint_angle_fc = 100.0

host_conformational_template = [["O5", "C1", "O1", "C4"], ["C1", "O1", "C4", "C5"]]

host_residues = len(hg[":MGO"].residues)
first_host_residue = hg[":MGO"].residues[0].number + 1
conformational_restraint_atoms = []
conformational_restraint_targets = []
conformational_restraint_fc = 6.0

for n in range(first_host_residue, host_residues + first_host_residue):
    if n + 1 < host_residues + first_host_residue:
        next_residue = n + 1
    else:
        next_residue = first_host_residue
    conformational_restraint_atoms.append(
        [
            f":{n}@{host_conformational_template[0][0]}",
            f":{n}@{host_conformational_template[0][1]}",
            f":{n}@{host_conformational_template[0][2]}",
            f":{next_residue}@{host_conformational_template[0][3]}",
        ]
    )
    conformational_restraint_targets.append(104.30)
    conformational_restraint_atoms.append(
        [
            f":{n}@{host_conformational_template[1][0]}",
            f":{n}@{host_conformational_template[1][1]}",
            f":{next_residue}@{host_conformational_template[1][2]}",
            f":{next_residue}@{host_conformational_template[1][3]}",
        ]
    )
    conformational_restraint_targets.append(-108.8)

guest_wall_template = [["O2", guest_anchors[0]], ["O6", guest_anchors[0]]]

guest_wall_restraint_atoms = []
guest_wall_restraint_targets = []
guest_wall_restraint_angle_fc = 500.0
guest_wall_restraint_distance_fc = 50.0

for n in range(first_host_residue, host_residues + first_host_residue):
    guest_wall_restraint_atoms.append(
        [f":{n}@{guest_wall_template[0][0]}", f"{guest_wall_template[0][1]}"]
    )
    guest_wall_restraint_targets.append(11.3)
    guest_wall_restraint_atoms.append(
        [f":{n}@{guest_wall_template[1][0]}", f"{guest_wall_template[1][1]}"]
    )
    guest_wall_restraint_targets.append(13.3)

guest_wall_restraint_atoms.append(
    [dummy_anchors[1], guest_anchors[0], guest_anchors[1]]
)
guest_wall_restraint_targets.append(80.0)

static_restraints = []
for index, atoms in enumerate(static_restraint_atoms):
    this = static_DAT_restraint(
        restraint_mask_list=atoms,
        num_window_list=windows,
        ref_structure=hg,
        force_constant=static_restraint_angle_fc
        if len(atoms) > 2
        else static_restraint_distance_fc,
        amber_index=True,
    )
    static_restraints.append(this)

guest_restraints = []
for index, atoms in enumerate(guest_restraint_atoms):
    if len(atoms) > 2:
        angle = True
    else:
        angle = False
    this = DAT_restraint()
    this.auto_apr = True
    this.amber_index = True
    this.topology = hg
    this.mask1 = atoms[0]
    this.mask2 = atoms[1]
    if angle:
        this.mask3 = atoms[2]
        this.attach["fc_final"] = guest_restraint_angle_fc
        this.release["fc_final"] = guest_restraint_angle_fc
    else:
        this.attach["fc_final"] = guest_restraint_distance_fc
        this.release["fc_final"] = guest_restraint_distance_fc
    this.attach["target"] = guest_restraint_targets[index]
    this.attach["fraction_list"] = attach_fractions

    this.pull["target_final"] = guest_restraint_target_final[index]
    this.pull["num_windows"] = windows[1]

    this.release["target"] = guest_restraint_targets[index]
    # Keep the guest restraints on during release.
    this.release["fraction_list"] = [1.0] * windows[2]

    this.initialize()
    guest_restraints.append(this)

conformational_restraints = []
for index, atoms in enumerate(conformational_restraint_atoms):
    this = DAT_restraint()
    this.auto_apr = True
    this.amber_index = True
    this.topology = hg
    this.mask1 = atoms[0]
    this.mask2 = atoms[1]
    this.mask3 = atoms[2]
    this.mask4 = atoms[3]

    this.attach["fraction_list"] = attach_fractions
    this.attach["target"] = conformational_restraint_targets[index]
    this.attach["fc_final"] = conformational_restraint_fc
    this.pull["target_final"] = conformational_restraint_targets[index]
    this.pull["num_windows"] = windows[1]

    this.release["fraction_list"] = release_fractions
    this.release["target"] = conformational_restraint_targets[index]
    this.release["fc_final"] = conformational_restraint_fc

    this.initialize()
    conformational_restraints.append(this)

wall_restraints = []
for index, atoms in enumerate(guest_wall_restraint_atoms):
    if len(atoms) > 2:
        angle = True
    else:
        angle = False

    this = DAT_restraint()
    this.auto_apr = True
    this.amber_index = True
    this.topology = hg
    this.mask1 = atoms[0]
    this.mask2 = atoms[1]
    if angle:
        this.mask3 = atoms[2]
        this.attach["fc_initial"] = guest_wall_restraint_angle_fc
        this.attach["fc_final"] = guest_wall_restraint_angle_fc
        this.custom_restraint_values["rk2"] = 500.0
        this.custom_restraint_values["rk3"] = 0.0
    else:
        this.attach["fc_initial"] = guest_wall_restraint_distance_fc
        this.attach["fc_final"] = guest_wall_restraint_distance_fc
        this.custom_restraint_values["rk2"] = 50.0
        this.custom_restraint_values["rk3"] = 50.0
        this.custom_restraint_values["r1"] = 0.0
        this.custom_restraint_values["r2"] = 0.0

    this.attach["target"] = guest_wall_restraint_targets[index]
    this.attach["num_windows"] = len(attach_fractions)

    this.initialize()
    wall_restraints.append(this)

window_list = create_window_list(guest_restraints)

# Analysis

structure = pt.load(
    os.path.join(prefix, "a000", "full.crds"),
    os.path.join(prefix, "a000", "full.hmr.topo"),
)

stripped = structure.strip(":WAT,:Na+,:Cl-")

analyze = fe_calc()
analyze.prmtop = stripped.topology
analyze.trajectory = "prod.*.nc"
analyze.path = prefix
analyze.restraint_list = guest_restraints + conformational_restraints
analyze.collect_data()
analyze.methods = ["ti-block"]
analyze.quicker_ti_matrix = True
analyze.bootcycles = 1000
analyze.compute_free_energy()
analyze.compute_ref_state_work(
    [guest_restraints[0], guest_restraints[1], None, None, guest_restraints[2], None]
)

attach_fe = analyze.results["attach"]["ti-block"]["fe"]
pull_fe = analyze.results["pull"]["ti-block"]["fe"]
release_fe = analyze.results["release"]["ti-block"]["fe"]
std_state_fe = analyze.results["ref_state_work"]

print(attach_fe)
print(pull_fe)
print(release_fe)
print(std_state_fe)
print(attach_fe + pull_fe + -1 * release_fe + std_state_fe)

print(analyze.results)
