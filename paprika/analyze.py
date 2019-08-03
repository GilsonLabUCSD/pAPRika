"""
This class contains a simulation analysis wrapper for use with the Property Estimator.
"""

import logging
from pathlib import Path
import parmed as pmd

import numpy as np
from paprika.analysis import fe_calc
from paprika.io import load_restraints

logger = logging.getLogger(__name__)


class Analyze(object):
    """
    The Analyze class provides a wrapper function around the analysis of simulations.
    """

    def __init__(self, host, guest=None, guest_orientation=None, restraint_file="restraints.json",
                 topology_file='coordinates.pdb', trajectory_mask='*.dcd', directory_path="benchmarks",
                 symmetry_correction=None, guest_residue_name=None):

        self.host = host
        self.guest = guest if guest is not None else "release"
        self.guest_residue_name = self.guest if not guest_residue_name else guest_residue_name

        self.directory = Path(directory_path).joinpath(self.host).joinpath(f"{self.guest}-{guest_orientation}" if
                                                                           guest_orientation is not None else
                                                                           f"{self.guest}")

        self.restraints = load_restraints(self.directory.joinpath(restraint_file))
        self.results = self.analyze(topology_file, trajectory_mask).results

        if symmetry_correction:
            if symmetry_correction["microstates"] == 1:
                logger.info("No microstates specified for the symmetry correction. Assuming separate calculations.")
            elif symmetry_correction["microstates"] > 1:
                self.results["symmetry_correction"] = -0.593 * np.log(symmetry_correction["microstates"])
            else:
                logger.warning("Symmetry correction defined in YAML but unable to determine microstates.")

    def analyze(self, topology_file, trajectory_mask):

        analysis = fe_calc()
        analysis.prmtop = topology_file
        analysis.trajectory = trajectory_mask
        analysis.path = self.directory.joinpath('windows')
        analysis.restraint_list = self.restraints
        analysis.methods = ["ti-block", "mbar-block"]
        analysis.bootcycles = 100
        analysis.collect_data(single_prmtop=False)
        if self.guest != "release":
            analysis.compute_free_energy(phases=["attach", "pull"])

            structure = pmd.load_file(str(analysis.path.joinpath("a000").joinpath(analysis.prmtop)))
            r_restraint = None
            beta_restraint = None
            theta_restraint = None

            for restraint in self.restraints:

                mask2_residue_name = structure[restraint.mask2].residues[0].name

                if "DM" in restraint.mask1 and self.guest_residue_name.upper() in mask2_residue_name and not restraint.mask3 and \
                        not restraint.mask4:
                    r_restraint = restraint
                if restraint.mask3 and not restraint.mask4:
                    mask3_residue_name = structure[restraint.mask3].residues[0].name
                    if "DM" in restraint.mask1 and "DM" in restraint.mask2 and self.guest_residue_name.upper() in mask3_residue_name:
                        theta_restraint = restraint
                    if "DM" in restraint.mask1 and self.guest_residue_name.upper() in mask2_residue_name and self.guest_residue_name.upper() in \
                            mask3_residue_name:
                        beta_restraint = restraint
            if r_restraint and theta_restraint and beta_restraint:
                analysis.compute_ref_state_work(
                    [r_restraint, theta_restraint, None, None, beta_restraint, None]
                )
            else:
                logging.warning("Could not determine the r, θ, or β restraint for computing the analytic release step.")

        else:
            analysis.compute_free_energy(phases=["release"])

        return analysis
