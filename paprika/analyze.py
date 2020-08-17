"""
This class contains a simulation analysis wrapper for use with the OpenFF Evaluator.
"""

import logging
from typing import List

import numpy as np
from paprika.analysis import fe_calc
from paprika.restraints import DAT_restraint

logger = logging.getLogger(__name__)


class Analyze(object):
    """
    The Analyze class provides a wrapper function around the analysis of simulations.
    """

    @classmethod
    def compute_phase_free_energy(
        cls,
        phase: str,
        restraints: List[DAT_restraint],
        windows_directory: str,
        topology_path: str,
        trajectory_mask: str = "*.dcd",
        bootstrap_cycles: int = 1000,
        analysis_method: str = "ti-block"
    ):

        analysis = fe_calc()

        analysis.prmtop = topology_path
        analysis.trajectory = trajectory_mask
        analysis.path = windows_directory

        analysis.restraint_list = restraints

        analysis.methods = [analysis_method]
        analysis.bootcycles = bootstrap_cycles

        analysis.collect_data(single_prmtop=False)
        analysis.compute_free_energy(phases=[phase])

        return analysis.results

    @classmethod
    def compute_ref_state_work(
        cls,
        temperature: float,
        guest_restraints: List[DAT_restraint],
    ) -> float:
        """

        Parameters
        ----------
        temperature
            The temperature at which the calculation was performed in units of kelvin.
        guest_restraints

        Returns
        -------

        """

        # reference_path = str(analysis.path.joinpath("a000").joinpath(analysis.prmtop))

        analysis = fe_calc()
        analysis.temperature = temperature

        distance_restraint = next(
            (
                restraint
                for restraint in guest_restraints
                if "DM" in restraint.mask1
                and restraint.mask2 is not None
                and restraint.mask3 is None
                and restraint.mask4 is None
            ),
            None
        )

        theta_restraint = next(
            (
                restraint
                for restraint in guest_restraints
                if "DM" in restraint.mask1
                and "DM" in restraint.mask2
                and restraint.mask3 is not None
                and restraint.mask4 is None
            ),
            None
        )
        beta_restraint = next(
            (
                restraint
                for restraint in guest_restraints
                if "DM" in restraint.mask1
                and "DM" not in restraint.mask2
                and restraint.mask3 is not None
                and restraint.mask4 is None
            ),
            None
        )

        if not distance_restraint or not theta_restraint or not beta_restraint:

            raise RuntimeError(
                "Could not determine the r, θ, or β restraint for computing the "
                "analytic release step."
            )

        analysis.compute_ref_state_work(
            [distance_restraint, theta_restraint, None, None, beta_restraint, None]
        )

        return analysis.results["ref_state_work"]

    @classmethod
    def symmetry_correction(cls, microstates: int) -> float:

        assert microstates > 0

        if microstates == 1:
            return 0.0

        return -0.593 * np.log(microstates)
