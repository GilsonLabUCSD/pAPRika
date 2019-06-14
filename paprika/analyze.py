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

    def __init__(self, host, guest=None, restraint_file="restraints.json",
                 topology_file='coordinates.pdb', trajectory_mask='*.dcd', directory_path="benchmarks"):

        self.host = host
        self.guest = guest if guest is not None else "release"
        self.directory = Path(directory_path).joinpath(self.host).joinpath(self.guest)

        self.restraints = load_restraints(self.directory.joinpath(restraint_file))
        self.results = self.analyze(topology_file, trajectory_mask).results

    def analyze(self, topology_file, trajectory_mask):

        analysis = fe_calc()
        analysis.prmtop = topology_file   # str(self.directory.joinpath(f"{self.host}-{self.guest}.pdb"))
        analysis.trajectory = trajectory_mask
        analysis.path = self.directory.joinpath('windows')
        analysis.restraint_list = self.restraints
        analysis.methods = ["ti-block"]
        analysis.bootcycles = 100
        analysis.collect_data(single_prmtop=False)
        if self.guest != "release":
            analysis.compute_free_energy(phases=["attach", "pull"])

            for restraint in self.restraints:
                if "DM" in restraint.mask1 and self.guest.upper() in restraint.mask2 and not restraint.mask3 and \
                        not restraint.mask4:
                    r_restraint = restraint
                if "DM" in restraint.mask1 and "DM" in restraint.mask2 and self.guest.upper() in restraint.mask3 \
                        and not restraint.mask4:
                    theta_restraint = restraint
                if "DM" in restraint.mask1 and self.guest.upper() in restraint.mask2 and self.guest.upper() in \
                        restraint.mask3 and not restraint.mask4:
                    beta_restraint = restraint

            analysis.compute_ref_state_work(
                [r_restraint, theta_restraint, None, None, beta_restraint, None]
            )

        else:
            analysis.compute_free_energy(phases=["release"])

        return analysis
