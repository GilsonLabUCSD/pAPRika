"""
This class contains a simulation analysis wrapper for use with the Property Estimator.
"""

import pkg_resources
import shutil
import parmed as pmd
import os as os
import simtk.openmm as openmm
import simtk.unit as unit

from pathlib import Path
from paprika.io import load_restraints
from paprika.analysis import fe_calc

import logging
logger = logging.getLogger(__name__)

class Analyze(object):
    """
    The Analyze class provides a wrapper function around the analysis of simulations.
    """
    def __init__(self, host, guest, restraint_file="restraints.json", directory_path="benchmarks",
                 phases=["attach", "pull"]):
        self.host = host
        self.guest = guest
        self.directory = Path(directory_path).joinpath(self.host).joinpath(self.guest)

        restraints = load_restraints(self.directory.joinpath(restraint_file))

        analysis = fe_calc()
        analysis.prmtop = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")
        analysis.trajectory = "*.dcd"
        analysis.path = self.directory
        analysis.restraint_list = restraints
        analysis.methods = ["ti-block"]
        analysis.bootcycles = 1
        analysis.collect_data(single_prmtop=True)
        analysis.compute_free_energy(phases=phases)
