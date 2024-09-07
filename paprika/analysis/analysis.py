import json
import logging
import os
import warnings
from itertools import compress

import numpy as np
import pymbar
from openff.units import unit as openff_unit

try:  # pymbar >= 4
    from pymbar.timeseries import detect_equilibration
except ImportError:
    from pymbar.timeseries import (
        detectEquilibration as detect_equilibration,
    )

from paprika.analysis.bootstrap import integrate_bootstraps
from paprika.analysis.utils import (
    get_block_sem,
    get_nearest_max,
    get_subsampled_indices,
)
from paprika.io import (
    PaprikaDecoder,
    PaprikaEncoder,
    load_trajectory,
    read_restraint_data,
)
from paprika.utils import check_unit, is_file_and_not_empty

logger = logging.getLogger(__name__)


class fe_calc(object):
    """
    Computes the free energy for an APR transformation. After calling ``compute_free_energy()``, the
    results are stored in a dictionary called ``results`` in units of kcal/mol.

    To get the binding free energy at standard state (ΔG°), we have used the following sign convention:

    .. math ::
        ΔG° = ΔG_{attach} + ΔG_{pull} - ΔG_{release} + ΔG_{reference}

    .. note ::
        It would be great to add unit support here from ``pint``.

    .. warning ::
        This class really ought to be split into a few smaller classes that sublcass a ``BaseAnalysis`` class.
        There could be separate ``MBARAnalysis`` and ``TIAnalysis`` classes, for example. This would make it much
        more modular and easy to test where TI and MBAR disagree. This could also be used to benchmark
        autocorrelation-based and blocking-analysis-based methods for evaluating the statistical inefficiency.

    """

    @property
    def temperature(self):
        """
        float or openff.units.unit.Quantity: The temperature used during the simulation. This will update β (1/kT) as well.
        """
        return self._temperature

    @temperature.setter
    def temperature(self, new_temperature):
        """Update temperature and β with a new temperature."""
        self._temperature = check_unit(new_temperature, base_unit=openff_unit.kelvin)
        self.beta = 1 / (self.k_B * self._temperature)

    @property
    def topology(self):
        """
        os.PathLike: The topology (prmtop or pdb) file used for the analysis.
        This should match the topology used for the simulation.
        """
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    @property
    def trajectory(self):
        """
        os.PathLike: File name of the trajectories to be analyzed in each window (can include a wildcard).
        """
        return self._trajectory

    @trajectory.setter
    def trajectory(self, value):
        self._trajectory = value

    @property
    def path(self):
        """
        os.PathLike: The parent directory that contains the simulation windows.
        """
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

    @property
    def restraint_list(self):
        """
        list: The list of restraints to be used for analysis.
        """
        return self._restraint_list

    @restraint_list.setter
    def restraint_list(self, value):
        self._restraint_list = value

    @property
    def changing_restraints(self):
        """
        dict: Dictionary containing which restraints change during which phase of the calculation.

        .. note ::
            This should probably be a private attribute, as this property is determined automatically.

        """
        return self._changing_restraints

    @changing_restraints.setter
    def changing_restraints(self, value):
        self._changing_restraints = value

    @property
    def orders(self):
        """
        The sorted order of windows for analysis. In principle, windows could be out-of-order if subsequent
        additional sampling was requested.

        .. note ::
            As far as I know, we have not tested analysis on windows out of order. We imagined that one could add
            additional λ windows through multiple runs by modifying existing restraints (e.g., restraint values of
            [0, 0.5, 1.0, 0.8]) and this module would be able to correctly sort the simulation directories and
            restraint targets.

        .. note ::
            This should probably be a private attribute, as this property is determined automatically.

        """
        return self._orders

    @orders.setter
    def orders(self, value):
        self._orders = value

    @property
    def simulation_data(self):
        """
        dict: Dictionary containing collected trajectory values for the relevant restraints and windows

        .. note ::
            This should probably be a private attribute, as this property is determined automatically.

        """
        return self._simulation_data

    @simulation_data.setter
    def simulation_data(self, value):
        self._simulation_data = value

    @property
    def methods(self):
        """
        list: List of analysis methods to be performed. This is a combination of free energy method (e.g., MBAR, TI,
        ...) and de-correlation method (blocking, autocorrelation, ...).

        Implemented methods are:

            - ``mbar-autoc``
            - ``mbar-block``
            - ``mbar-boot``
            - ``ti-block``

        .. note ::
            This is really fragile. We should definitely use something like an ``ENUM`` here to check for the few
            combinations that we support.

        .. todo ::
            - Add ``ti-autoc``.
            - Add ``wham-block``.
            - Add ``wham-autoc``.
            - Clarify naming.
            - Look into using ``pymbar`` for decorrelation.

        """
        return self._methods

    @methods.setter
    def methods(self, value):
        self._methods = value

    @property
    def conservative_subsample(self):
        """
        bool: Whether the statistical inefficiency is rounded up to the nearest integer. If ``False``, a non-integer
        value
        is used.
        """
        return self._conservative_subsample

    @conservative_subsample.setter
    def conservative_subsample(self, value):
        self._conservative_subsample = value

    @property
    def boot_cycles(self) -> int:
        """
        int: The number of bootstrap iterations for the TI methods.

        **Default**: ``1000``
        """
        return self._boot_cycles

    @boot_cycles.setter
    def boot_cycles(self, value):
        self._boot_cycles = value

    @property
    def compute_roi(self) -> bool:
        """
        bool: Whether to compute the return on investment (ROI) value for each window. The more negative the ROI, the
        better the return on investment of computing more frames for this particular window.
        """
        return self._compute_roi

    @compute_roi.setter
    def compute_roi(self, value):
        self._compute_roi = value

    @property
    def compute_largest_neighbor(self):
        """
        bool: Whether to find and store the maximum SEM to the neighbor windows. This is useful if using the "scale_w"
        approach which is described in Equation (1) in: https://pubs.acs.org/doi/10.1021/acs.jctc.9b00748.
        """
        return self._compute_largest_neighbor

    @compute_largest_neighbor.setter
    def compute_largest_neighbor(self, value):
        self._compute_largest_neighbor = value

    @property
    def ti_matrix(self):
        """
        str: If ``full``, the TI mean and SEM free energy is computed between all windows (i.e., the energy differences
        between all windows and all other windows). If ``diagonal``, the mean and SEM of the free energy is computed
        between the first window and all other windows, as well as between all neighboring windows. If ``endpoints``,
        the mean and SEM free energy is computed between only the first and the last window.

        For most cases, ``diagonal`` should be sufficient and ``full`` is overkill.

        .. note ::
            This is fragile, and we should be checking for the valid strings supported.

        """
        return self._ti_matrix

    @ti_matrix.setter
    def ti_matrix(self, value):
        if value not in ["full", "diagonal", "endpoints"]:
            raise ValueError(f"{value} is not a supported integration scheme.")

        self._ti_matrix = value

    @property
    def exact_sem_each_ti_fraction(self):
        """
        bool: Whether the SEM is computed once for the full data set and then the SEM for each fraction is estimated
        based on the total SEM and the fractional number of uncorrelated data points. If ``True``, the SEM will be
        recomputed each fraction using just that fraction of the raw data.
        """
        return self._exact_sem_each_ti_fraction

    @exact_sem_each_ti_fraction.setter
    def exact_sem_each_ti_fraction(self, value):
        self._exact_sem_each_ti_fraction = value

    @property
    def fractions(self) -> list:
        """
        list: A list of fractions of the total data for which a free energy will be computed. Default is [1.0].
        """
        return self._fractions

    @fractions.setter
    def fractions(self, value):
        for fraction in value:
            if fraction < 0.0 or fraction > 1.0:
                raise ValueError(
                    f"Unable to calculation fraction of the data: {fraction}."
                )

        self._fractions = value

    @property
    def results(self):
        """
        dict: A dictionary containing the results. The ``results`` dictionary is indexed first by phase, and then by
        method. That is, ``results["attach"]["mbar-block"]`` will contain the results from analyzing the attach
        phase with the MBAR free energy estimator and blocking analysis used to estimate the SEM.

        The final free energy and uncertainty estimate is specified in the ``fe`` and ``sem`` keys.

        If multiple fractions of the data are specified, the free energy and SEM from each fraction are stored in
        a nested dictionary under the keys ``fraction_fe`` and ``fraction_sem``. The number of frames analyzed for
        each fraction is stored under ``fraction_n_frames``.

        The full free energy matrix (window-to-window free energy differences) is stored in the ``fe_matrix`` entry.
        Likewise, the full SEM matrix (window-to-window free energy SEMs) is stored in the ``sem_matrix`` entry.

        The work to release the guest to standard concentration is stored under ``ref_state_work``, and is
        negative by convention.

        The format of a typical dictionary will resemble this (where I have a omitted the values in the matrices
        for clarity):

        ::

            {'attach': {'mbar-block': {'fe': 13.019994367423227,
                                       'fe_matrix': array(),
                                       'fraction_fe': {'1.0': 13.019994367423227},
                                       'fraction_fe_matrix': {'1.0': array()},
                                       'fraction_n_frames': {'1.0': 2280000},
                                       'fraction_sem': {'1.0': 0.12783177525365627},
                                       'fraction_sem_matrix': {'1.0': array()},
                                       'n_frames': 2280000,
                                       'sem': 0.12783177525365627,
                                       'sem_matrix': array()},
                        'window_order': array([ 0,  1,  2, ..., 12, 13, 14])},
             'pull': {'mbar-block': {'fe': 5.367885925413715,
                                     'fe_matrix': array(),
                                     'fraction_fe': {'1.0': 5.367885925413715},
                                     'fraction_fe_matrix': {'1.0': array()},
                                     'fraction_n_frames': {'1.0': 2900000},
                                     'fraction_sem': {'1.0': 0.21231171453084993},
                                     'fraction_sem_matrix': {'1.0': array()},
                                     'n_frames': 2900000,
                                     'sem': 0.21231171453084993,
                                     'sem_matrix': array()},
                      'window_order': array([ 0,  1,  2, ..., 43, 44, 45])},
             'ref_state_work': -7.141514582862005}

        .. note ::
            This will be automatically populated with the outcome of the analysis. Do not directly modify these values.


        """
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    @property
    def energy_unit(self):
        """openff.unit: The based unit for energy."""
        return self._energy_unit

    @energy_unit.setter
    def energy_unit(self, value: openff_unit.Quantity):
        self._energy_unit = value

    @property
    def distance_unit(self):
        """openff.unit: The base unit for distance."""
        return self._distance_unit

    @distance_unit.setter
    def distance_unit(self, value: openff_unit.Quantity):
        self._distance_unit = value

    @property
    def angle_unit(self):
        """openff.unit: The base unit for angles."""
        return self._angle_unit

    @angle_unit.setter
    def angle_unit(self, value: openff_unit.Quantity):
        self._angle_unit = value

    @property
    def temperature_unit(self):
        """openff.unit: The base unit for temperature."""
        return self._temperature_unit

    @temperature_unit.setter
    def temperature_unit(self, value: openff_unit.Quantity):
        self._temperature_unit = value

    def __init__(self):
        self._energy_unit = openff_unit.kcal / openff_unit.mole
        self._distance_unit = openff_unit.angstrom
        self._angle_unit = openff_unit.radian
        self._temperature_unit = openff_unit.kelvin

        self._temperature = 298.15 * self._temperature_unit
        self.k_B = 1.987204118e-3 * self._energy_unit / self._temperature_unit
        self.beta = 1 / (self.k_B * self._temperature)
        self._topology = None
        self._trajectory = None
        self._path = None
        self._restraint_list = []
        self._changing_restraints = None
        self._orders = None
        self._simulation_data = None
        self._methods = ["mbar-block"]
        self._conservative_subsample = False
        self._boot_cycles = 1000
        self._compute_roi = False
        self._compute_largest_neighbor = False
        self._ti_matrix = "full"
        self._exact_sem_each_ti_fraction = False
        self._fractions = [1.0]
        self._results = {}

    def collect_data(self, single_topology=False):
        """
        Gather simulation data on the distance, angle, and torsion restraints that change during the simulation.

        Parameters
        ----------
        single_topology: bool
            Whether a single `topology` file is read for all windows
        """

        self.changing_restraints = self.identify_changing_restraints()
        self.orders = self.determine_window_order()
        self.simulation_data = self.read_trajectories(single_topology=single_topology)

    def collect_data_from_json(self, filepath):
        """
        Read in simulation data from a JSON file.

        Parameters
        ----------
        filepath: os.PathLike
            The name of the JSON file.
        """
        with open(filepath, "r") as f:
            json_data = f.read()
            data = json.loads(json_data, cls=PaprikaDecoder)

        self.changing_restraints = data["changing_restraints"]
        self.orders = data["orders"]
        self.simulation_data = data["simulation_data"]

    def identify_changing_restraints(self):
        """Figure out which restraints change during each phase of the calculation.

        Returns
        -------
        changing_restraints : dict
            A dictionary containing which restraints change during which phase of the calculation
        """

        changing_restraints = {"attach": [], "pull": [], "release": []}

        for phase in ["attach", "pull", "release"]:
            if phase == "attach" or phase == "release":
                changing_parameter = "force_constants"
            else:
                changing_parameter = "targets"
            for restraint in self.restraint_list:
                if restraint.phase[phase][changing_parameter] is not None:
                    static = all(
                        np.isclose(x, restraint.phase[phase][changing_parameter][0])
                        for x in restraint.phase[phase][changing_parameter]
                    )
                else:
                    static = True

                changing_restraints[phase].append(not static)

        return changing_restraints

    def determine_window_order(self):
        """Order the trajectories (i.e., simulation windows) in terms of increasing force constants and
        targets for each restraint.

        Returns
        -------
        orders : dict
            The sorted order of windows for analysis
        """

        orders = {"attach": [], "pull": [], "release": []}
        active_attach_restraints = list(
            compress(self.restraint_list, self.changing_restraints["attach"])
        )
        active_pull_restraints = list(
            compress(self.restraint_list, self.changing_restraints["pull"])
        )
        active_release_restraints = list(
            compress(self.restraint_list, self.changing_restraints["release"])
        )

        attach_orders = []
        pull_orders = []
        release_orders = []

        for restraint in active_attach_restraints:
            attach_orders.append(
                np.argsort(restraint.phase["attach"]["force_constants"])
            )
        if not all([np.array_equal(attach_orders[0], i) for i in attach_orders]):
            raise Exception(
                "The order of increasing force constants is not the same in all restraints."
            )
        elif attach_orders:
            orders["attach"] = attach_orders[0]
        else:
            orders["attach"] = np.empty(0)

        for restraint in active_pull_restraints:
            pull_orders.append(np.argsort(restraint.phase["pull"]["targets"]))
        if not all([np.array_equal(pull_orders[0], i) for i in pull_orders]):
            raise Exception(
                "The order of increasing target distances is not the same in all restraints."
            )
        elif pull_orders:
            orders["pull"] = pull_orders[0]
        else:
            orders["pull"] = np.empty(0)

        for restraint in active_release_restraints:
            release_orders.append(
                np.argsort(restraint.phase["release"]["force_constants"])
            )
        if not all([np.array_equal(release_orders[0], i) for i in release_orders]):
            raise Exception(
                "The order of increasing force constants is not the same in all restraints."
            )
        elif release_orders:
            orders["release"] = release_orders[0]
        else:
            orders["release"] = np.empty(0)

        return orders

    def read_trajectories(self, single_topology=False):
        """For each phase and window, and for each non-static restraint, parse the trajectories to
        get the restraint values.

        Parameters
        ----------
        single_topology : bool
            Whether a single `topology` file is read for all windows

        Returns
        -------
        data : dict
            Dictionary containing restraint values for analysis
        """

        data = {"attach": [], "pull": [], "release": []}

        ordered_attach_windows = [
            os.path.join(self.path, "a{:03d}".format(i))
            for i in self.orders["attach"]
            if i is not None
        ]
        ordered_pull_windows = [
            os.path.join(self.path, "p{:03d}".format(i))
            for i in self.orders["pull"]
            if i is not None
        ]
        ordered_release_windows = [
            os.path.join(self.path, "r{:03d}".format(i))
            for i in self.orders["release"]
            if i is not None
        ]

        active_attach_restraints = np.asarray(self.restraint_list)[
            self.changing_restraints["attach"]
        ]
        active_pull_restraints = np.asarray(self.restraint_list)[
            self.changing_restraints["pull"]
        ]
        active_release_restraints = np.asarray(self.restraint_list)[
            self.changing_restraints["release"]
        ]

        # Niel: I'm just checking if *one* restraint is `continuous_apr`,
        # which should be the same value for all restraints.
        if len(active_attach_restraints) > 0:
            if (
                active_attach_restraints[0].continuous_apr
                and self.orders["attach"].size > 0
                and self.orders["pull"].size > 0
            ):
                logger.debug(
                    "Replacing {} with {} in {} for `continuous_apr`...".format(
                        ordered_attach_windows[-1],
                        ordered_pull_windows[0],
                        ordered_attach_windows,
                    )
                )
                ordered_attach_windows[-1] = ordered_pull_windows[0]

        if len(active_release_restraints) > 0:
            if (
                active_release_restraints[0].continuous_apr
                and self.orders["release"].size > 0
                and self.orders["pull"].size > 0
            ):
                logger.debug(
                    "Replacing {} with {} in {} for `continuous_apr`...".format(
                        ordered_release_windows[-1],
                        ordered_pull_windows[-1],
                        ordered_release_windows,
                    )
                )
                ordered_release_windows[-1] = ordered_pull_windows[-1]

        for window_index, window in enumerate(ordered_attach_windows):
            phase = "attach"
            data[phase].append([])
            traj = load_trajectory(
                window, self.trajectory, self.topology, single_topology
            )
            for restraint_index, restraint in enumerate(active_attach_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(
                    traj, restraint, self.distance_unit, self.angle_unit
                )

        for window_index, window in enumerate(ordered_pull_windows):
            phase = "pull"
            data[phase].append([])
            traj = load_trajectory(
                window, self.trajectory, self.topology, single_topology
            )
            for restraint_index, restraint in enumerate(active_pull_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(
                    traj, restraint, self.distance_unit, self.angle_unit
                )

        for window_index, window in enumerate(ordered_release_windows):
            phase = "release"
            data[phase].append([])
            traj = load_trajectory(
                window, self.trajectory, self.topology, single_topology
            )
            for restraint_index, restraint in enumerate(active_release_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(
                    traj, restraint, self.distance_unit, self.angle_unit
                )

        return data

    def prepare_data(self, phase):
        number_of_windows = len(self.simulation_data[phase])
        data_points = [len(np.asarray(x).T) for x in self.simulation_data[phase]]
        max_data_points = max(data_points)
        active_restraints = list(
            compress(self.restraint_list, self.changing_restraints[phase])
        )
        force_constants = [
            np.copy(i.phase[phase]["force_constants"]) for i in active_restraints
        ]
        targets = [i.phase[phase]["targets"] for i in active_restraints]

        ordered_force_constants = [i[self.orders[phase]] for i in force_constants]
        ordered_targets = [i[self.orders[phase]] for i in targets]

        # Convert targets and force constants to proper base units
        for i, rest in enumerate(active_restraints):
            # Distance
            if rest.mask3 is None and rest.mask4 is None:
                ordered_targets[i] = ordered_targets[i].to(self.distance_unit)
                ordered_force_constants[i] = ordered_force_constants[i].to(
                    self.energy_unit / self.distance_unit**2
                )
            # Angles
            elif rest.mask3 is not None or rest.mask4 is not None:
                ordered_targets[i] = ordered_targets[i].to(self.angle_unit)
                ordered_force_constants[i] = ordered_force_constants[i].to(
                    self.energy_unit / self.angle_unit**2
                )

        return (
            number_of_windows,
            data_points,
            max_data_points,
            active_restraints,
            ordered_force_constants,
            ordered_targets,
            self.simulation_data[phase],
        )

    def run_mbar(self, phase, prepared_data, method, verbose=False):
        """
        Compute the free energy matrix for a series of windows. We'll follow the pymbar nomenclature
        for data structures.

        Parameters
        ----------
        phase: str
            The phase of the calculation to analyze.
        prepared_data: :class:`np.array`
            The list of "prepared data" including the number of windows, data points, which restraints are changing,
            their force constants and targets, and well as the order of the windows. This probably ought to be
            redesigned.
        method: str
            The method used to calculate the SEM.
        verbose: bool, optional
            Whether to set the `verbose` option on pyMBAR.
        """

        # Unpack the prepared data
        (
            num_win,
            data_points,
            max_data_points,
            active_rest,
            force_constants,
            targets,
            ordered_values,
        ) = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points)

        # Set up the reduced potential energy array. ie, the potential of each window's
        # coordinates in each window's potential function
        u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)

        # Transpose force_constants and targets into "per window" format, instead of
        # the "per restraint" format.
        target_units = np.array([targets[r][0].units for r in range(len(active_rest))])
        force_units = np.array(
            [force_constants[r][0].units for r in range(len(active_rest))]
        )

        force_constants_T = np.asarray(force_constants).T * force_units
        targets_T = np.asarray(targets).T * target_units

        # Note, the organization of k = coordinate windows, l = potential windows
        # seems to be opposite of the documentation. But I got wrong numbers
        # the other way around.
        for k in range(num_win):  # Coordinate windows
            for l in range(num_win):  # Potential Windows
                for r, rest in enumerate(active_rest):  # Restraints
                    # If this is a dihedral, we need to shift around restraint value
                    # on the periodic axis to make sure the lowest potential is
                    # used.
                    if rest.mask3 is not None and rest.mask4 is not None:
                        target = targets_T[l, r]
                        # Coords from coord window, k
                        bool_list = ordered_values[k][r] < target - (
                            180.0 * openff_unit.degrees
                        ).to(self.angle_unit)
                        ordered_values[k][r][bool_list] += (
                            360.0 * openff_unit.degrees
                        ).to(self.angle_unit)
                        bool_list = ordered_values[k][r] > target + (
                            180.0 * openff_unit.degrees
                        ).to(self.angle_unit)
                        ordered_values[k][r][bool_list] -= (
                            360.0 * openff_unit.degrees
                        ).to(self.angle_unit)

                # Compute the potential ... for each frame, sum the contributions for each restraint
                # Note, we multiply by beta, and do some extra [l,:,None] to
                # get the math operation correct.
                u_kln[k, l, 0 : N_k[k]] = sum(
                    [
                        self.beta * k * (val - eq) ** 2
                        for k, val, eq in zip(
                            force_constants_T[l], ordered_values[k], targets_T[l]
                        )
                    ]
                ).magnitude

        g_k = np.ones([num_win], np.float64)
        # Should I subsample based on the restraint coordinate values? Here I'm
        # doing it on the potential.  Should be pretty close ....
        if method == "mbar-block":
            # We want to use all possible data to get the free energy estimates Deltaf_ij,
            # but for uncertainty estimates we'll subsample to create
            # uncorrelated data.
            for k in range(num_win):
                l = k
                # If the potential is zero everywhere, we can't estimate the uncertainty, so
                # check the next *potential* window which probably had non-zero
                # force constants
                while not u_kln[k, l, 0 : N_k[k]].any():
                    l += 1
                # Now compute statistical inefficiency: g = N*(SEM**2)/variance
                nearest_max = get_nearest_max(N_k[k])
                sem = get_block_sem(u_kln[k, l, 0:nearest_max])
                variance = np.var(u_kln[k, l, 0 : N_k[k]])
                g_k[k] = N_k[k] * (sem**2) / variance

        if method == "mbar-autoc" or method == "mbar-boot":
            for k in range(num_win):
                [t0, g_k[k], Neff_max] = detect_equilibration(
                    N_k[k]
                )  # compute indices of uncorrelated using timeseries

        # Create subsampled indices and count their lengths. If g=1, ie no correlation,
        # then subsampling will return identical indices to original
        # (hopefully)
        ss_indices = []
        N_ss = np.zeros([num_win], np.int32)  # N_subsample
        for k in range(num_win):
            ss_indices.append(
                get_subsampled_indices(
                    N_k[k], g_k[k], conservative=self.conservative_subsample
                )
            )
            N_ss[k] = len(ss_indices[k])

        self.results[phase][method]["fraction_fe_matrix"] = {}
        self.results[phase][method]["fraction_sem_matrix"] = {}
        self.results[phase][method]["fraction_fe_Neffective"] = {}
        self.results[phase][method]["fraction_sem_Neffective"] = {}

        for fraction in self.fractions:
            # Setup mbar calc, and get matrix of free energies, uncertainties
            # To estimate the free energy, we won't do subsampling.  We'll do
            # another MBAR calculation later with subsampling to estimate the
            # uncertainty.
            frac_N_k = np.array([int(fraction * n) for n in N_k], dtype=np.int32)

            mbar = pymbar.MBAR(u_kln, frac_N_k, verbose=verbose)
            try:  # pymbar >= 4
                mbar_results = mbar.compute_free_energy_differences(
                    compute_uncertainty=True
                )
            except AttributeError:
                mbar_results = mbar.getFreeEnergyDifferences(
                    compute_uncertainty=True, return_dict=True
                )

            Deltaf_ij = mbar_results["Delta_f"]
            dDeltaf_ij = mbar_results["dDelta_f"]

            try:  # pymbar >= 4
                Deltaf_ij_N_eff = mbar.compute_effective_sample_number()
            except AttributeError:
                Deltaf_ij_N_eff = mbar.computeEffectiveSampleNumber()

            # Estimate uncertainty from decorrelated samples
            # Create subsampled indices and count their lengths
            frac_N_ss = np.array([int(fraction * n) for n in N_ss], dtype=np.int32)

            # Create a new potential array for the uncertainty calculation
            # (are we using too much memory?)
            u_kln_err = np.zeros([num_win, num_win, np.max(frac_N_ss)], np.float64)

            # Populate the subsampled array, drawing the appropriate
            # fraction of subsamples from the original
            for k in range(num_win):
                for l in range(num_win):
                    u_kln_err[k, l, 0 : frac_N_ss[k]] = u_kln[
                        k, l, ss_indices[k][0 : frac_N_ss[k]]
                    ]

            # We toss junk_Deltaf_ij, because we got a better estimate for it from above using all data.
            # But dDeltaf_ij will replace the previous, because it correctly accounts for the
            # correlation in the data.
            if method == "mbar-boot":
                try:  # pymbar >= 4
                    mbar = pymbar.MBAR(
                        u_kln_err,
                        frac_N_ss,
                        verbose=verbose,
                        n_bootstraps=self.boot_cycles,
                    )
                    mbar_results = mbar.compute_free_energy_differences(
                        compute_uncertainty=True, uncertainty_method="bootstrap"
                    )
                except AttributeError:
                    logger.warning(
                        "MBAR with bootstrapping is not available, this feature is only available with PyMBAR >= 4. "
                        "Reverting to `mbar-autoc` for PyMBAR >=3,<4."
                    )
                    mbar = pymbar.MBAR(u_kln_err, frac_N_ss, verbose=verbose)
                    mbar_results = mbar.getFreeEnergyDifferences(
                        compute_uncertainty=True, return_dict=True
                    )
            else:
                mbar = pymbar.MBAR(u_kln_err, frac_N_ss, verbose=verbose)
                try:  # pymbar >= 4
                    mbar_results = mbar.compute_free_energy_differences(
                        compute_uncertainty=True
                    )
                except AttributeError:
                    mbar_results = mbar.getFreeEnergyDifferences(
                        compute_uncertainty=True, return_dict=True
                    )

            dDeltaf_ij = mbar_results["dDelta_f"]
            try:  # pymbar >= 4
                dDeltaf_ij_N_eff = mbar.compute_effective_sample_number()
            except AttributeError:
                dDeltaf_ij_N_eff = mbar.computeEffectiveSampleNumber()

            # Put back into kcal/mol
            Deltaf_ij /= self.beta
            dDeltaf_ij /= self.beta

            self.results[phase][method]["fraction_fe_matrix"][fraction] = Deltaf_ij
            self.results[phase][method]["fraction_fe_Neffective"][
                fraction
            ] = Deltaf_ij_N_eff
            self.results[phase][method]["fraction_sem_matrix"][fraction] = dDeltaf_ij
            self.results[phase][method]["fraction_sem_Neffective"][
                fraction
            ] = dDeltaf_ij_N_eff

    def run_ti(self, phase, prepared_data, method):
        """
        Compute the free energy using the TI method.

        We compute the partial derivative of the potential (i.e., forces), for each frame, with respect to the
        changing parameter, either a lambda or target value. The force constants are scaled by the λ parameter which
        controls their strength: ``0`` to ``fc_max``.

        Potential:
          U = λ × fc_max × (values - target)²
        Forces during attach:
          dU/dλ = fc_max × (values - target)²
        Forces during pull:
          dU/d(target) = 2 × λ × fc_max × (values - target)
        Forces during release:
          (same as attach)

        Then we integrate over the interval covered by λ or target.

        Parameters
        ----------
        phase: str
            The phase of the calculation to analyze.
        prepared_data: :class:`np.array`
            The list of "prepared data" including the number of windows, data points, which restraints are changing,
            their force constants and targets, and well as the order of the windows. This probably ought to be
            redesigned.
        method: str
            The method used to calculate the SEM.

        .. note ::
            ``phase`` and ``method`` should be `enum` types and tied to class attributes.

        .. warning ::
            We have only considered whether this will work for the case where the pull phase is a single distance
            restraint with a changing target value. This has not been tested for a changing angle or dihedral.

        """

        # Unpack the prepared data
        (
            num_win,
            data_points,
            max_data_points,
            active_rest,
            force_constants,
            targets,
            ordered_values,
        ) = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points)

        # The dU array to store the partial derivative of the potential with respect lambda or target,
        # depending on the whether attach/release or pull. Data stored for each frame.  This just a
        # temporary storage space.
        dU = np.zeros([num_win, max_data_points], np.float64)

        # The mean, SEM, standard deviation, and number of uncorrelated dU values for each window.
        dU_avgs = np.zeros([num_win], np.float64)
        dU_sems = np.zeros([num_win], np.float64)
        dU_stdv = np.zeros([num_win], np.float64)
        dU_Nunc = np.zeros([num_win], np.float64)
        # The statistical inefficiency
        g = np.zeros([num_win], np.float64)

        # Array for values of the changing coordinate (x-axis), either lambda or target.
        # I'll name them dl_vals for dlambda values.
        dl_vals = np.zeros([num_win], np.float64)

        # Setup interpolation array for the dLambda (dl) coordinate. We're gonna create
        # this progressively by appending ...
        dl_intp = np.zeros([0], np.float64)

        # Get units
        target_units = np.array([targets[r][0].units for r in range(len(active_rest))])
        force_units = np.array(
            [force_constants[r][0].units for r in range(len(active_rest))]
        )

        # Store the max force constant value for each restraint.
        max_force_constants = (
            np.array(
                [np.max(force_constants[r]).magnitude for r in range(len(active_rest))]
            )
            * force_units
        )

        # Transpose force_constants and targets into "per window" format, instead of
        # the "per restraint" format.
        # print(targets)
        force_constants_T = np.asarray(force_constants).T * force_units
        targets_T = np.asarray(targets).T * target_units

        # For each window: do dihedral wrapping, compute forces, append dl_intp
        for k in range(num_win):  # Coordinate windows
            # Wrap dihedrals so we get the right potential
            for r, rest in enumerate(active_rest):  # Restraints
                # If this is a dihedral, we need to shift around restraint value
                # on the periodic axis to make sure the lowest potential is used.

                if rest.mask3 is not None and rest.mask4 is not None:
                    target = targets_T[k, r]
                    bool_list = ordered_values[k][r] < target - (
                        180.0 * openff_unit.degrees
                    ).to(self.angle_unit)
                    ordered_values[k][r][bool_list] += (360.0 * openff_unit.degrees).to(
                        self.angle_unit
                    )
                    bool_list = ordered_values[k][r] > target + (
                        180.0 * openff_unit.degrees
                    ).to(self.angle_unit)
                    ordered_values[k][r][bool_list] -= (360.0 * openff_unit.degrees).to(
                        self.angle_unit
                    )

            # Compute forces and store the values of the changing coordinate,
            # either lambda or target
            if phase == "attach" or phase == "release":
                dU[k, 0 : N_k[k]] = (
                    sum(
                        [
                            (k * (val - eq) ** 2)
                            for k, val, eq in zip(
                                max_force_constants, ordered_values[k], targets_T[k]
                            )
                        ]
                    )
                    .to(self.energy_unit)
                    .magnitude
                )

                # this is lambda. assume the same scaling for all restraints
                dl_vals[k] = force_constants_T[k, 0] / max_force_constants[0]
            else:
                dU[k, 0 : N_k[k]] = (
                    sum(
                        [
                            2.0 * k * (val - eq)
                            for k, val, eq in zip(
                                max_force_constants, ordered_values[k], targets_T[k]
                            )
                        ]
                    )
                    .to(self.energy_unit / self.distance_unit)
                    .magnitude
                )

                # Currently assuming a single distance restraint
                dl_vals[k] = targets_T[k, 0].to(self.distance_unit).magnitude

            # Compute standard deviations and SEMs, unless we're going to do
            # exact_sem_each_ti_fraction
            dU_avgs[k] = np.mean(dU[k, 0 : N_k[k]])
            dU_stdv[k] = np.std(dU[k, 0 : N_k[k]])
            if method == "ti-block":
                nearest_max = get_nearest_max(N_k[k])
                dU_sems[k] = get_block_sem(dU[k, 0:nearest_max])
                # Rearrange SEM = StdDev/sqrt(N) to get N_uncorrelated
                dU_Nunc[k] = (dU_stdv[k] / dU_sems[k]) ** 2
            elif method == "ti-nocor":
                dU_sems[k] = dU_stdv[k] / np.sqrt(N_k[k])
                dU_Nunc[k] = N_k[k]
            g[k] = N_k[k] / dU_Nunc[k]

            # Create the interpolation by appending 100 points between each window.
            # Start with k=1 so we don't double count.
            if k > 0:
                dl_intp = np.append(
                    dl_intp,
                    np.linspace(dl_vals[k - 1], dl_vals[k], num=100, endpoint=False),
                )

        # Tack on the final value to the dl interpolation
        dl_intp = np.append(dl_intp, dl_vals[-1])

        logger.debug("Running bootstrap calculations...")

        # Setup fractions. For simplicity, we'll always do this, even
        # if we're doing the total data, ie self.fractions=[1.0].
        self.results[phase][method]["fraction_fe_matrix"] = {}
        self.results[phase][method]["fraction_sem_matrix"] = {}

        for fraction in self.fractions:
            logger.debug("Working on fraction ... {}".format(fraction))

            # Compute means for this fraction.
            frac_dU_avgs = np.array(
                [np.mean(dU[k, 0 : int(fraction * n)]) for k, n in enumerate(N_k)]
            )

            # If self.exact_sem_each_ti_fraction, we're gonna recompute the SEM for each fraction
            # rather than estimating it from the standard deviation (dU_stdv) and number of
            # uncorrelated data points (dU_Nunc) from the total data set.
            if method == "ti-block" and self.exact_sem_each_ti_fraction:
                frac_dU_sems = np.zero([k], np.float64)
                for k in range(num_win):
                    nearest_max = get_nearest_max(int(fraction * N_k[k]))
                    frac_dU_sems[k] = get_block_sem(dU[k, 0:nearest_max])
            elif method == "ti-nocor" and self.exact_sem_each_ti_fraction:
                frac_dU_sems = np.zero([k], np.float64)
                for k in range(num_win):
                    frac_dU_sems[k] = np.std(
                        dU[k, 0 : int(fraction * N_k[k])]
                    ) / np.sqrt(int(fraction * N_k[k]))
            else:
                frac_dU_sems = dU_stdv / np.sqrt(fraction * dU_Nunc)

            dU_samples = np.random.normal(
                frac_dU_avgs, frac_dU_sems, size=(self.boot_cycles, frac_dU_avgs.size)
            )

            # Run bootstraps
            (
                self.results[phase][method]["fraction_fe_matrix"][fraction],
                self.results[phase][method]["fraction_sem_matrix"][fraction],
            ) = integrate_bootstraps(
                dl_vals, dU_samples, x_intp=dl_intp, matrix=self.ti_matrix
            )

            # Put units back on
            self.results[phase][method]["fraction_fe_matrix"][
                fraction
            ] *= self.energy_unit
            self.results[phase][method]["fraction_sem_matrix"][
                fraction
            ] *= self.energy_unit

            # The attach/release work (integration) yields appropriately positive work, but
            # the pull work needs a negative multiplier. Think W = −Force × distance type thing.
            if phase == "pull":
                self.results[phase][method]["fraction_fe_matrix"][fraction] *= -1.0

        if self.compute_roi:
            logger.info(phase + ": computing ROI for " + method)
            # Do ROI calc
            max_fraction = np.max(self.fractions)
            # If we didn't compute fe/sem for fraction 1.0 already, do it now
            dU_samples = np.random.normal(
                dU_avgs, dU_sems, size=(self.boot_cycles, dU_avgs.size)
            )
            if not np.isclose(max_fraction, 1.0):
                junk_fe, total_sem_matrix = integrate_bootstraps(
                    dl_vals, dU_samples, x_intp=dl_intp, matrix=self.ti_matrix
                )
            else:
                total_sem_matrix = self.results[phase][method]["fraction_sem_matrix"][
                    max_fraction
                ].magnitude
            self.results[phase][method]["roi"] = np.zeros([num_win], np.float64)

            for k in range(num_win):
                # Compute overall integrated SEM with 10% smaller SEM for dU[k]
                cnvg_dU_samples = np.array(dU_samples)
                cnvg_dU_samples[:, k] = np.random.normal(
                    dU_avgs[k], 0.9 * dU_sems[k], self.boot_cycles
                )
                junk_fe, cnvg_sem_matrix = integrate_bootstraps(
                    dl_vals, cnvg_dU_samples, x_intp=dl_intp, matrix=self.ti_matrix
                )

                #         d( dG_sem )      d( dUdl_sem )
                # ROI = --------------- * ---------------
                #        d( dUdl_sem )     d( n_frames )
                #
                # Deriv1----^---^---^        ^---^---^----Deriv2

                # Deriv1:
                deriv1 = (cnvg_sem_matrix[0, -1] - total_sem_matrix[0, -1]) / (
                    -0.1 * dU_sems[k]
                )

                # Deriv2:
                #
                # dUdl_sem = dUdl_stddev / sqrt(n_frames/g)
                #
                # d( dUdl_sem )              dUdl_stddev
                # -------------- =  - --------------------------
                # d( n_frames )          2g * (n_frames/g)**3/2

                deriv2 = (
                    -1.0 * dU_stdv[k] / (2.0 * g[k] * (N_k[k] / g[k]) ** (3.0 / 2.0))
                )

                # ROI
                self.results[phase][method]["roi"][k] = deriv1 * deriv2

    def compute_free_energy(self, phases=["attach", "pull", "release"], seed=None):
        """
        Compute the free energy of binding from a simulation. This function populates the ``results`` dictionary
        of the :class:`fe_calc` object.

        Parameters
        ----------
        phases: list
            Which phases of the calculation to analyze.
        seed: int
            Random number seed.
        """

        for fraction in self.fractions:
            if fraction <= 0.0 or fraction > 1.0:
                raise Exception(
                    "The fraction of data to analyze must be 0 < fraction ≤ 1.0."
                )

        for phase in phases:
            self.results[phase] = {}
            self.results[phase]["window_order"] = self.orders[phase]

            for method in self.methods:
                if seed is not None:
                    np.random.seed(seed)
                    logger.debug(f"Setting random number seed = {seed}")

                self.results[phase][method] = {}

                # Prepare data
                if sum(self.changing_restraints[phase]) == 0:
                    logger.debug("Skipping free energy calculation for %s" % phase)
                    continue
                prepared_data = self.prepare_data(phase)
                self.results[phase][method]["n_frames"] = np.sum(prepared_data[1])

                logger.debug(
                    "Running {} analysis on {} phase ...".format(method, phase)
                )

                if (
                    method == "mbar-block"
                    or method == "mbar-autoc"
                    or method == "mbar-boot"
                ):
                    self.run_mbar(phase, prepared_data, method)
                elif method == "ti-block":
                    self.run_ti(phase, prepared_data, method)
                else:
                    raise NotImplementedError(
                        f"Method ({method}) is not implemented yet."
                    )

                # Store endpoint free energy and SEM for each fraction
                self.results[phase][method]["fraction_n_frames"] = {}
                self.results[phase][method]["fraction_fe"] = {}
                self.results[phase][method]["fraction_sem"] = {}

                for fraction in self.fractions:
                    self.results[phase][method]["fraction_n_frames"][fraction] = int(
                        fraction * self.results[phase][method]["n_frames"]
                    )

                    self.results[phase][method]["fraction_fe"][fraction] = self.results[
                        phase
                    ][method]["fraction_fe_matrix"][fraction][0, -1]

                    self.results[phase][method]["fraction_sem"][fraction] = (
                        self.results[phase][method]["fraction_sem_matrix"][fraction][
                            0, -1
                        ]
                    )

                # Set these higher level (total) values, which will be slightly
                # easier to access
                max_fraction = np.max(self.fractions)
                self.results[phase][method]["fe_matrix"] = self.results[phase][method][
                    "fraction_fe_matrix"
                ][max_fraction]
                self.results[phase][method]["sem_matrix"] = self.results[phase][method][
                    "fraction_sem_matrix"
                ][max_fraction]
                self.results[phase][method]["fe"] = self.results[phase][method][
                    "fe_matrix"
                ][0, -1]
                self.results[phase][method]["sem"] = self.results[phase][method][
                    "sem_matrix"
                ][0, -1]

                if self.compute_largest_neighbor:
                    # Store convergence values, which are helpful for running
                    # simulations
                    windows = len(self.results[phase][method]["sem_matrix"])
                    self.results[phase][method]["largest_neighbor"] = (
                        openff_unit.Quantity(
                            np.ones([windows], np.float64) * -1.0,
                            units=self.energy_unit,
                        )
                    )
                    logger.info(f"{phase}: computing largest_neighbor for {method}...")

                    for i in range(windows):
                        if i == 0:
                            self.results[phase][method]["largest_neighbor"][i] = (
                                self.results[phase][method]["sem_matrix"][i][i + 1]
                            )
                        elif i == windows - 1:
                            self.results[phase][method]["largest_neighbor"][i] = (
                                self.results[phase][method]["sem_matrix"][i][i - 1]
                            )
                        else:
                            left = self.results[phase][method]["sem_matrix"][i][i - 1]
                            right = self.results[phase][method]["sem_matrix"][i][i + 1]
                            if left > right:
                                max_val = left
                            elif right > left:
                                max_val = right
                            else:
                                max_val = right
                            self.results[phase][method]["largest_neighbor"][i] = max_val

    def compute_ref_state_work(self, restraints, state="final"):
        """
        Compute the work to place a molecule at standard reference state conditions starting from a state defined by
        up to six restraints. These are Boresch-style restraints.

        Parameters
        ----------
        restraints : list or dict
            A list or dict of :class:`paprika.restraints.DAT_restraint` objects. If a list is specified, the restraints
            must be in order of the six translational and orientational restraints (Boresch-style). The six restraints
            are: [r, theta, phi, alpha, beta, gamma] and they should be passed to this function in that order. If any
            of these coordinates is not being restrained, use a `None` in place of a
            :class:`paprika.restraints.DAT_restraint` object. For a dictionary, the restraints should have `r`, `theta`,
            `phi`, `alpha`, `beta`, `gamma`, as the dictionary keys.

            See :meth:`paprika.analysis.ref_state_work` for details on the calculation.
        state : str, optional, default="final"
            Option to estimate the reference work using the `initial` state (DDM) or `final` pulling state (APR) for
            `r0` when calculating the reference work of releasing guest restraints in bulk water. See Equation 9 in ref
            Henriksen et al. (2015).
        """
        # Check input argument
        state = state.lower()
        if state not in ["initial", "final"]:
            raise ValueError(
                'Function argument `state` can only take values "initial" or "final" as input.'
            )

        # Convert list to dictionary
        if isinstance(restraints, list):
            warnings.warn(
                "Converting restraint list to dictionary, make sure the restraints are listed "
                "in the order of [`r`, `theta`,`phi`, `alpha`, `beta`, `gamma`]."
            )
            restraints = {
                "r": restraints[0],
                "theta": restraints[1],
                "phi": restraints[2],
                "alpha": restraints[3],
                "beta": restraints[4],
                "gamma": restraints[5],
            }

        # Raise error if no distance restraint exist
        if not restraints or restraints["r"] is None:
            raise ValueError(
                "At minimum, a single distance restraint is necessary to compute the work of releasing "
                "the guest to standard state."
            )

        fcs = []
        targs = []

        # Distance restraint - special care when extracting `r0`, depends on whether calculating with APR or DDM.
        for colvar in ["r", "theta", "phi", "alpha", "beta", "gamma"]:
            restraint = restraints[colvar]

            target_index = -1
            if colvar == "r" and state == "initial":
                target_index = 0

            if restraint is None:
                fcs.append(None)
                targs.append(None)

            else:
                target_and_force_exist = False

                phases = ["release", "pull"]
                if state == "initial":
                    phases = ["attach", "release"]

                for phase in phases:
                    if restraint.phase[phase]["force_constants"] is not None:
                        force_index = -1
                        if phase == "release":
                            force_index = 0

                        fcs.append(
                            np.sort(restraint.phase[phase]["force_constants"])[
                                force_index
                            ]
                        )
                        targs.append(
                            np.sort(restraint.phase[phase]["targets"])[target_index]
                        )

                        target_and_force_exist = True

                        break

                if not target_and_force_exist:
                    raise ValueError(
                        f"The `{colvar}` restraints should have at least attach/release values (for DDM) or "
                        "pull/release values (for APR) initialized in order to `compute_ref_state_work`."
                    )

        # Store reference work of release guest restraints
        self.results["ref_state_work"] = ref_state_work(
            self.temperature,
            fcs[0],
            targs[0],
            fcs[1],
            targs[1],
            fcs[2],
            targs[2],
            fcs[3],
            targs[3],
            fcs[4],
            targs[4],
            fcs[5],
            targs[5],
            self.energy_unit,
        )

    def save_results(self, filepath="results.json", overwrite=False):
        """
        Save the analysis results to a JSON file.

        Parameters
        ----------
        filepath: str
            The name of the JSON file to write to.
        overwrite: bool
            Option to whether overwrite file if already exist.
        """
        if not overwrite and is_file_and_not_empty(filepath):
            raise FileExistsError(f"File `{filepath}` exists, will not overwrite.")

        with open(filepath, "w") as f:
            dumped = json.dumps(self.results, cls=PaprikaEncoder)
            f.write(dumped)

    def load_results(self, filepath):
        """
        Read in a JSON file for the results.

        Parameters
        ----------
        filepath: str
            The name of the JSON file to read.
        """
        with open(filepath, "r") as f:
            data = f.read()

        self.results = json.loads(data, cls=PaprikaDecoder)

    def save_data(self, filepath="simulation_data.json", overwrite=False):
        """
        Save the simulation data (DAT values) to a JSON file.

        Parameters
        ----------
        filepath: str
            The name of the JSON file to write to.
        overwrite: bool
            Option to whether overwrite file if already exist.
        """
        if not overwrite and is_file_and_not_empty(filepath):
            raise FileExistsError(f"File `{filepath}` exists, will not overwrite.")

        with open(filepath, "w") as f:
            dumped = json.dumps(
                {
                    "simulation_data": self.simulation_data,
                    "changing_restraints": self.changing_restraints,
                    "orders": self.orders,
                },
                cls=PaprikaEncoder,
            )
            f.write(dumped)

    def load_data(self, filepath):
        """
        Load the simulation data (DAT values) from a JSON file.

        Parameters
        ----------
        filepath: str
            The name of the JSON file to read.
        """
        with open(filepath, "r") as f:
            data = f.read()

        simulation_data = json.loads(data, cls=PaprikaDecoder)

        self.simulation_data = simulation_data["simulation_data"]
        self.changing_restraints = simulation_data["changing_restraints"]
        self.orders = simulation_data["orders"]


def ref_state_work(
    temperature,
    r_fc,
    r_tg,
    th_fc,
    th_tg,
    ph_fc,
    ph_tg,
    a_fc,
    a_tg,
    b_fc,
    b_tg,
    g_fc,
    g_tg,
    energy_unit=openff_unit.kcal / openff_unit.mole,
):
    """
    Computes the free energy to release a molecule from some restrained translational
    and orientational configuration (relative to another molecule or lab frame) into
    the reference configuration: standard concentration (1.0/1660.5392 Å³) and
    unrestrained orientational freedom (8π²).

    The Wikipedia entry on Euler angles is useful.

    Assume two molecules (H and G). Three translational and three orientational degrees
    of freedom define their relative configuration. In order to match experimentally
    reported free energies, we often need to compute the work (free energy) of moving
    a molecule (G) from a restrained configuration relative to H into the experimental
    reference state (usually defined as 1 M, or 1 molecule per 1660.5 Å³)

    ::

        H3
          \        [a1]                [a2]
           H2-------H1-------<d1>-------G1-------G2
              {t1}           {t2}           {t3}   \
                                                   G3

        Degrees of Freedom
        -----------------------------------------------------------------
         id   atoms        type, spherical coordinate/Euler angle
        -----------------------------------------------------------------
        <d1>: H1-G1         distance, r
        [a1]: H2-H1-G1      angle, theta
        {t1}: H3-H2-H1-G1   torsion, phi
        {t2}: H2-H1-G1-G2   torsion, alpha
        [a2]: H1-G1-G2      angle, beta
        {t3}: H1-G1-G2-G3   torsion, gamma

    Parameters
    ----------
    temperature : openff.units.unit.Quantity
        The temperature (in Kelvin) at which the reference state calculation will take place.
    r_fc: openff.units.unit.Quantity
        The distance, :math:`r`, restraint force constant (kcal/mol-Å²).
    r_tg : openff.units.unit.Quantity
        The distance, :math:`r`, restraint target values (Å). The target range is 0 to infinity.
    th_fc : openff.units.unit.Quantity
        The angle, :math:`θ`, restraint force constant (kcal/mol-radian²).
    th_tg : openff.units.unit.Quantity
        The angle, :math:`θ`, restraint target values (radian). The target range is 0 to π.
    ph_fc : openff.units.unit.Quantity
        The torsion, :math:`\\phi`, restraint force constant (kcal/mol-radian²).
    ph_tg : openff.units.unit.Quantity
        The torsion, :math:`\\phi`, restraint target values (radian). The target range is 0 to 2π.
    a_fc : openff.units.unit.Quantity
        The torsion, :math:`α`, restraint force constant (kcal/mol-radian²).
    a_tg: openff.units.unit.Quantity
        The torsion, :math:`α`, restraint target values (radian). The target range is 0 to 2π.
    b_fc : openff.units.unit.Quantity
        The angle, :math:`β`, restraint force (kcal/mol-radian²).
    b_tg : openff.units.unit.Quantity
        The angle, :math:`β`, restraint target values (radian). The target range is 0 to π.
    g_fc: openff.units.unit.Quantity
        The angle, :math:`γ`, restraint force (kcal/mol-radian²).
    g_tg: openff.units.unit.Quantity
        The angle, :math:`γ`, restraint target values (radian). The target range is 0 to 2π.
    energy_unit: openff.units.unit.Quantity
        The unit for the output free energy.

    Returns
    -------
    RT * np.log(trans * orient): openff.units.unit.Quantity
        The free energy associated with releasing the restraints (in kcal/mol openff units).
    """

    # Convert all distances and angles to units
    distance_unit = openff_unit.angstrom
    angle_unit = openff_unit.radian

    R = 1.987204118e-3 * openff_unit.kcal / openff_unit.mole / openff_unit.kelvin
    RT = (R * temperature).to(energy_unit)

    # Distance Integration Function
    def dist_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return (arange**2) * np.exp((-1.0 / RT) * fc * (arange - targ) ** 2)

        targ = targ.to(distance_unit)
        fc = fc.to(energy_unit / distance_unit**2)
        arange = (np.arange(0.0, 100.0, 0.0001) * openff_unit.angstrom).to(
            distance_unit
        )

        return np.trapz(potential(arange, RT, fc, targ), arange)

    # Angle Integration Function
    def ang_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return np.sin(arange) * np.exp((-1.0 / RT) * fc * (arange - targ) ** 2)

        targ = targ.to(angle_unit)
        fc = fc.to(energy_unit / angle_unit**2)
        arange = (np.arange(0.0, np.pi, 0.00005) * openff_unit.radians).to(angle_unit)

        return np.trapz(potential(arange, RT, fc, targ), arange)

    # Torsion Integration Function
    def tors_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return np.exp((-1.0 / RT) * fc * (arange - targ) ** 2)

        # Note, because of periodicity, I'm gonna wrap +/- pi around target for integration.
        targ = targ.to(angle_unit)
        fc = fc.to(energy_unit / angle_unit**2)
        arange = (
            np.arange(targ.magnitude - np.pi, targ.magnitude + np.pi, 0.00005)
            * openff_unit.radians
        ).to(angle_unit)

        return np.trapz(potential(arange, RT, fc, targ), arange)

    # Distance restraint, r
    if None in [r_fc, r_tg]:
        raise Exception("Distance restraint info (r_fc, r_tg) must be specified")
    else:
        r_int = dist_int(RT, r_fc, r_tg)

    # Angle restraint, theta
    if None in [th_fc, th_tg]:
        th_int = 2.0
    else:
        th_int = ang_int(RT, th_fc, th_tg)

    # Torsion restraint, phi
    if None in [ph_fc, ph_tg]:
        ph_int = 2.0 * np.pi * openff_unit.radians
    else:
        ph_int = tors_int(RT, ph_fc, ph_tg)

    # Torsion restraint, alpha
    if None in [a_fc, a_tg]:
        a_int = 2.0 * np.pi * openff_unit.radians
    else:
        a_int = tors_int(RT, a_fc, a_tg)

    # Angle restraint, beta
    if None in [b_fc, b_tg]:
        b_int = 2.0
    else:
        b_int = ang_int(RT, b_fc, b_tg)

    # Torsion restraint, gamma
    if None in [g_fc, g_tg]:
        g_int = 2.0 * np.pi * openff_unit.radians
    else:
        g_int = tors_int(RT, g_fc, g_tg)

    # Concentration term
    V0 = 1660.5392 * openff_unit.angstrom**3
    translational = r_int * th_int * ph_int * (1.0 / V0)  # C^o = 1/V^o

    # Orientational term
    rotational_volume = 8.0 * np.pi**2
    orientational = a_int * b_int * g_int / rotational_volume

    # Return the free energy
    return RT * np.log(translational * orientational)
