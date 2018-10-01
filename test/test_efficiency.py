import parmed as pmd
import pytraj as pt
import numpy as np
import logging as log

from paprika import restraints
from paprika import analysis

logger = log.getLogger()
logger.setLevel(log.DEBUG)
log.basicConfig(format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %I:%M:%S %p")

inputpdb = pmd.load_file("cb6_but_gb_apr_ref_data/vac.pdb")

# Distance restraint
rest1 = restraints.DAT_restraint()
rest1.continuous_apr = True
rest1.amber_index = True
rest1.topology = inputpdb
rest1.mask1 = ":CB6@O"
rest1.mask2 = ":BUT@C1"
rest1.attach["target"] = 4.5
rest1.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
rest1.attach["fc_final"] = 5.0
rest1.pull["fc"] = rest1.attach["fc_final"]
rest1.pull["target_initial"] = rest1.attach["target"]
rest1.pull["target_final"] = 18.5
rest1.pull["num_windows"] = 19
rest1.initialize()

# Angle restraint
rest2 = restraints.DAT_restraint()
rest2.continuous_apr = True
rest2.amber_index = True
rest2.topology = inputpdb
rest2.mask1 = ":CB6@O1"
rest2.mask2 = ":CB6@O"
rest2.mask3 = ":BUT@C1"
rest2.attach["target"] = 8.0
rest2.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
rest2.attach["fc_final"] = 50.0
rest2.pull["fc"] = rest2.attach["fc_final"]
rest2.pull["target_initial"] = rest2.attach["target"]
rest2.pull["target_final"] = rest2.attach["target"]
rest2.pull["num_windows"] = 19
rest2.initialize()

# Dihedral restraint
rest3 = restraints.DAT_restraint()
rest3.continuous_apr = True
rest3.amber_index = True
rest3.topology = inputpdb
rest3.mask1 = ":CB6@O11"
rest3.mask2 = ":CB6@O1"
rest3.mask3 = ":CB6@O"
rest3.mask4 = ":BUT@C1"
rest3.attach["target"] = -60.0
rest3.attach["fraction_list"] = [0.00, 0.04, 0.181, 0.496, 1.000]
rest3.attach["fc_final"] = 50.0
rest3.pull["fc"] = rest3.attach["fc_final"]
rest3.pull["target_initial"] = rest3.attach["target"]
rest3.pull["target_final"] = rest3.attach["target"]
rest3.pull["num_windows"] = 19
rest3.initialize()

# Create window directories
window_list = restraints.create_window_list([rest1, rest2, rest3])

# Phase abbreviations
phase_dict = {"a": "attach", "p": "pull", "r": "release"}

np.random.seed(12345)

# fecalc = analysis.fe_calc()
# fecalc.prmtop = 'cb6_but_gb_apr_ref_data/vac.prmtop'
# fecalc.trajectory = '*.nc'
# fecalc.path = 'cb6_but_gb_apr_ref_data/'
# fecalc.restraint_list = [rest1, rest2, rest3]
# fecalc.methods = ['mbar-block']
# fecalc.bootcycles = 100
# fecalc.quick_ti_matrix = True
# fecalc.collect_data(single_prmtop=True)
# fecalc.compute_free_energy()

# print(fecalc.results['attach']['mbar-block']['fe'], fecalc.results['attach']['mbar-block']['sem'])
# print(fecalc.results['pull']['mbar-block']['fe'], fecalc.results['pull']['mbar-block']['sem'])

# fecalc = analysis.fe_calc()
# fecalc.prmtop = 'cb6_but_gb_apr_ref_data/vac.prmtop'
# fecalc.trajectory = '*.nc'
# fecalc.path = 'cb6_but_gb_apr_ref_data/'
# fecalc.restraint_list = [rest1, rest2, rest3]
# fecalc.methods = ['mbar-block']
# fecalc.bootcycles = 100
# fecalc.quick_ti_matrix = True
# fecalc.collect_data(single_prmtop=True)
# fecalc.fractions = [1.0]
# fecalc.compute_free_energy()

# print(fecalc.results['attach']['mbar-block']['fe'], fecalc.results['attach']['mbar-block']['sem'])
# print(fecalc.results['pull']['mbar-block']['fe'], fecalc.results['pull']['mbar-block']['sem'])

fecalc = analysis.fe_calc()
fecalc.prmtop = "cb6_but_gb_apr_ref_data/vac.prmtop"
fecalc.trajectory = "*.nc"
fecalc.path = "cb6_but_gb_apr_ref_data/"
fecalc.restraint_list = [rest1, rest2, rest3]
fecalc.methods = ["mbar-block"]
fecalc.bootcycles = 100
fecalc.quick_ti_matrix = True
fecalc.collect_data(single_prmtop=True)
fecalc.fractions = [0.1, 1.0]
fecalc.compute_free_energy()

for key, value in fecalc.results.items():
    print(key, value)
print(
    fecalc.results["attach"]["mbar-block"]["fe"],
    fecalc.results["attach"]["mbar-block"]["sem"],
)
print(
    fecalc.results["pull"]["mbar-block"]["fe"],
    fecalc.results["pull"]["mbar-block"]["sem"],
)
