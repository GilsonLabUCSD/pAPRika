"""
This class contains a simulation analysis wrapper for use with the OpenFF Evaluator.
"""

import logging
from typing import Any, Dict, List

import numpy as np

from paprika.analysis import fe_calc
from paprika.restraints import DAT_restraint

logger = logging.getLogger(__name__)


class Analyze:
    """
    The ``Analyze`` class provides a wrapper function around the analysis of simulations.
    """

    @classmethod
    def compute_phase_free_energy(
        cls,
        phase: str,
        restraints: List[DAT_restraint],
        windows_directory: str,
        topology_name: str,
        trajectory_mask: str = "*.dcd",
        bootstrap_cycles: int = 1000,
        analysis_method: str = "ti-block",
    ) -> Dict[str, Any]:
        """Computes the free energy of a specified phase of an APR calculation.

        Parameters
        ----------
        phase
            The phase to analyze. This should be one of ``'attach'``, ``'pull'`` or
            ``'release'``.
        restraints
            A list of the restraints which were used as part of the phase being
            analyzed.
        windows_directory
            The directory which contains the final trajectories and topologies
            from the simulated phase.
        topology_name
            The expected file name (not path) of the coordinate file which contains
            topological information about the system.
        trajectory_mask
            The pattern to use when searching for the simulated trajectories.
        bootstrap_cycles
            The number of bootstrap iterations to perform.
        analysis_method
            The analysis method to use.

        Returns
        -------
            The computed free energies (and their uncertainties) in units of kcal / mol
        """
        analysis = fe_calc()

        analysis.topology = topology_name
        analysis.trajectory = trajectory_mask
        analysis.path = windows_directory

        analysis.restraint_list = restraints

        analysis.methods = [analysis_method]
        analysis.boot_cycles = bootstrap_cycles

        analysis.collect_data(single_topology=False)
        analysis.compute_free_energy(phases=[phase])

        return analysis.results

    @classmethod
    def compute_ref_state_work(
        cls,
        temperature: float,
        guest_restraints: List[DAT_restraint],
    ) -> float:
        """Computes the reference state work of the attach phase.

        Parameters
        ----------
        temperature
            The temperature at which the calculation was performed in units of kelvin.
        guest_restraints
            The guest restraints which were applied during the pull phase.

        Returns
        -------
            The reference state work in units of kcal / mol.
        """

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
            None,
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
            None,
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
            None,
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
    def symmetry_correction(
        cls,
        n_microstates: int,
        temperature: float,
    ) -> float:
        """Computes the free energy corrections to apply to symmetrical guest
        when the guest is restrained to just one of several possible symmetrical
        configurations (e.g. butane restrained in a cyclic host).

        Parameters
        ----------
        temperature
            The temperature at which the calculation was performed in units of kelvin.
        n_microstates
            The number of different symmetrical microstates that the guest can exist
            in (e.g. for butane this is two).

        Returns
        -------
            The symmetry correction to apply in units of kcal / mol.
        """

        assert n_microstates > 0

        if n_microstates == 1:
            return 0.0

        return -temperature * 0.001987204258640832 * np.log(n_microstates)
