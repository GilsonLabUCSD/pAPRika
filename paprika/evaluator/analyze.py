"""
This class contains a simulation analysis wrapper for use with the OpenFF Evaluator.
"""

import logging
from typing import Any, Dict, List, Union

import numpy as np
from openff.units import unit as openff_unit

from paprika.analysis import fe_calc
from paprika.restraints import DAT_restraint
from paprika.utils import check_unit

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
        temperature: Union[float, openff_unit.Quantity],
        guest_restraints: List[DAT_restraint],
    ) -> openff_unit.Quantity:
        """Computes the reference state work of the 'attach' phase.

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

        boresch_restraints = {
            "r": None,
            "theta": None,
            "phi": None,
            "alpha": None,
            "beta": None,
            "gamma": None,
        }

        for restraint in guest_restraints:
            # Distance
            if not restraint.mask3 and not restraint.mask4:
                boresch_restraints["r"] = restraint

            # Angle
            elif restraint.mask3 and not restraint.mask4:
                if "DM2" in restraint.mask1 and "DM1" in restraint.mask2:
                    boresch_restraints["theta"] = restraint
                elif "DM1" in restraint.mask1 and "DM" not in restraint.mask2:
                    boresch_restraints["beta"] = restraint

            # Dihedral restraints
            if restraint.mask4:
                if (
                    "DM3" in restraint.mask1
                    and "DM2" in restraint.mask2
                    and "DM1" in restraint.mask3
                ):
                    boresch_restraints["phi"] = restraint

                elif "DM2" in restraint.mask1 and "DM1" in restraint.mask2:
                    boresch_restraints["beta"] = restraint

                elif "DM1" in restraint.mask1:
                    boresch_restraints["gamma"] = restraint

        if boresch_restraints["r"] is None:
            raise RuntimeError(
                "Could not determine the `r` restraint! Need at least `r` restraint "
                "for computing the analytic release step."
            )

        analysis.compute_ref_state_work(boresch_restraints)

        return analysis.results["ref_state_work"]

    @classmethod
    def symmetry_correction(
        cls,
        n_microstates: int,
        temperature: Union[float, openff_unit.Quantity],
    ) -> openff_unit.Quantity:
        """Computes the free energy corrections to apply to symmetrical guest
        when the guest is restrained to just one of several possible symmetrical
        configurations (e.g. butane restrained in a cyclic host).

        Parameters
        ----------
        n_microstates
            The number of different symmetrical microstates that the guest can exist
            in (e.g. for butane this is two).
        temperature
            The temperature at which the calculation was performed in units of kelvin.

        Returns
        -------
            The symmetry correction to apply in units of kcal / mol.
        """

        assert n_microstates > 0

        if n_microstates == 1:
            return 0.0 * openff_unit.kcal / openff_unit.mol

        k_B = 1.987204118e-3 * openff_unit.kcal / openff_unit.mol / openff_unit.kelvin
        temperature = check_unit(temperature, base_unit=openff_unit.kelvin)

        return -temperature * k_B * np.log(n_microstates)
