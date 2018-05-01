import logging as log
import os as os
from itertools import compress
import numpy as np
import pytraj as pt
import pymbar


class fe_calc(object):
    """
    Computes the free energy for an APR transformation. After calling `compute_free_energy()`, the 
    results are stored in a dictionary called `results` in kcal/mol.

    Attributes
    ----------
    temperature : {float}
        The simulation temperature (Kelvin)
    k_B : {float}
        Boltzmann's constant (kcal/mol-K)
    beta : {float}
        1 / (k_B * T)  (kcal/mol)
    prmtop : {str} or ParmEd AmberParm
        Simulation parameters
    trajectory : {str}
        File name of the trajectories (can probably include a wildcard)
    path : {str}
        The parent directory that contains the simulation windows
    restraint_list : {list}
        List of restraints to be analyzed (this can now be a list of all simulation restraints)
    changing_restraints : {dict}
        Dictionary containing which restraints change during which phase of the calculation
    orders : {dict}
        The sorted order of windows for analysis
    simulation_data : {dict}
        Dictionary containing collected trajectory values for the relevant restraints and windows
    methods : {list}
        List of analysis methods to be performed (e.g., MBAR, TI, ...)
    results : {dict}
        TODO: description
    """

    def __init__(self):

        self._temperature = 298.15
        self.k_B = 0.0019872041
        self.beta = 1 / (self.k_B * self._temperature)

        self.prmtop = None
        self.trajectory = None
        self.path = None

        self.restraint_list = []
        self.changing_restraints = None
        self.orders = None
        self.simulation_data = None

        self.methods = ['mbar-block']  # mbar-autoc, mbar-none, ti-block, ti-autoc, ti-none
        # TODO: Add check that fe_methods and subsample_methods have correct keywords

        self.results = {}

    @property
    def temperature(self):
        """Allow updating of beta with temperature."""
        return self._temperature

    @temperature.setter
    def temperature(self, new_temperature):
        """Set updating of beta with new temperature."""
        self.beta = 1 / (self.k_B * new_temperature)

    def collect_data(self, single_prmtop=False, fraction=1.0):
        """Gather simulation data on the distance, angle, and torsion restraints that change during the simulation.

        """

        self.changing_restraints = self.determine_static_restraints()
        self.orders = self.determine_window_order()
        self.simulation_data = self.read_trajectories(single_prmtop=single_prmtop, fraction=fraction)

    def determine_static_restraints(self):
        """Figure out which restraints change during each phase of the calculation.
        
        Returns
        -------
        changing_restraints : {dict}
            A dictionary containing which restraints change during which phase of the calculation
        """

        changing_restraints = {'attach': [], 'pull': [], 'release': []}

        for phase in ['attach', 'pull', 'release']:
            if phase == 'attach' or phase == 'release':
                changing_parameter = 'force_constants'
            else:
                changing_parameter = 'targets'
            for restraint in self.restraint_list:
                if restraint.phase[phase][changing_parameter] is not None:
                    static = all(
                        np.isclose(x, restraint.phase[phase][changing_parameter][0])
                        for x in restraint.phase[phase][changing_parameter])
                else:
                    static = True

                changing_restraints[phase].append(not static)

        return changing_restraints

    def determine_window_order(self):
        """Order the trajectories (i.e., simulation windows) in terms of increasing force constants and
        targets for each restraint.
        
        Returns
        -------
        orders : {dict}
            The sorted order of windows for analysis
        """

        orders = {'attach': [], 'pull': [], 'release': []}
        active_attach_restraints = np.asarray(self.restraint_list)[self.changing_restraints['attach']]
        active_pull_restraints = np.asarray(self.restraint_list)[self.changing_restraints['pull']]
        active_release_restraints = np.asarray(self.restraint_list)[self.changing_restraints['release']]

        attach_orders = []
        pull_orders = []
        release_orders = []

        for restraint in active_attach_restraints:
            attach_orders.append(np.argsort(restraint.phase['attach']['force_constants']))
        if not all([np.array_equal(attach_orders[0], i) for i in attach_orders]):
            raise Exception('The order of increasing force constants is not the same in all restraints.')
        elif attach_orders:
            orders['attach'] = attach_orders[0]
        else:
            orders['attach'] = []

        for restraint in active_pull_restraints:
            pull_orders.append(np.argsort(restraint.phase['pull']['targets']))
        if not all([np.array_equal(pull_orders[0], i) for i in pull_orders]):
            raise Exception('The order of increasing target distances is not the same in all restraints.')
        elif pull_orders:
            orders['pull'] = pull_orders[0]
        else:
            orders['pull'] = []

        for restraint in active_release_restraints:
            release_orders.append(np.argsort(restraint.phase['release']['force_constants']))
        if not all([np.array_equal(release_orders[0], i) for i in release_orders]):
            raise Exception('The order of increasing force constants is not the same in all restraints.')
        elif release_orders:
            orders['release'] = release_orders[0]
        else:
            orders['release'] = []

        return orders

    def read_trajectories(self, single_prmtop=False, fraction=1.0):
        """For each each phase and window, and for each non-static restraint, parse the trajectories to 
        get the restraint values.

        Parameters
        ----------
        single_prmtop : {bool}
            Whether a single `prmtop` is read for all windows
        fraction : {float}
            Fraction of data to read, to check free energy convergence

        Returns
        -------
        data : {dict}
            Dictionary containing restraint values for analysis 
        """

        data = {'attach': [], 'pull': [], 'release': []}

        ordered_attach_windows = [
            os.path.join(self.path, 'a{:03d}'.format(i)) for i in self.orders['attach'] if i is not None
        ]
        ordered_pull_windows = [
            os.path.join(self.path, 'p{:03d}'.format(i)) for i in self.orders['pull'] if i is not None
        ]
        ordered_release_windows = [
            os.path.join(self.path, 'r{:03d}'.format(i)) for i in self.orders['release'] if i is not None
        ]

        active_attach_restraints = np.asarray(self.restraint_list)[self.changing_restraints['attach']]
        active_pull_restraints = np.asarray(self.restraint_list)[self.changing_restraints['pull']]
        active_release_restraints = np.asarray(self.restraint_list)[self.changing_restraints['release']]

        # Niel: I'm just checking if *one* restraint is `continuous_apr`,
        # which should be the same value for all restraints.
        if active_attach_restraints[0].continuous_apr and self.orders['attach'].any and self.orders['pull'].any:
            log.debug('Replacing {} with {} in {} for `continuous_apr`...'.format(
                ordered_attach_windows[-1], ordered_pull_windows[0], ordered_attach_windows))
            ordered_attach_windows[-1] = ordered_pull_windows[0]

        # This is inefficient and slow.
        # I am going to separately loop through the attach, then pull, then release windows.
        # Niel: I'm sure you can think of a better solution.

        for window_index, window in enumerate(ordered_attach_windows):
            phase = 'attach'
            data[phase].append([])
            for restraint_index, restraint in enumerate(active_attach_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(restraint, window, self.trajectory,
                                                                                 self.prmtop, single_prmtop, fraction)

        for window_index, window in enumerate(ordered_pull_windows):
            phase = 'pull'
            data[phase].append([])
            for restraint_index, restraint in enumerate(active_pull_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(restraint, window, self.trajectory,
                                                                                 self.prmtop, single_prmtop, fraction)

        for window_index, window in enumerate(ordered_release_windows):
            phase = 'release'
            data[phase].append([])
            for restraint_index, restraint in enumerate(active_release_restraints):
                data[phase][window_index].append([])
                data[phase][window_index][restraint_index] = read_restraint_data(restraint, window, self.trajectory,
                                                                                 self.prmtop, single_prmtop, fraction)

        return data

    def _prepare_data(self, phase):
        number_of_windows = len(self.simulation_data[phase])
        data_points = [len(np.asarray(x).T) for x in self.simulation_data[phase]]
        max_data_points = max(data_points)
        active_restraints = list(compress(self.restraint_list, self.changing_restraints[phase]))
        force_constants = [i.phase[phase]['force_constants'] for i in active_restraints]
        targets = [i.phase[phase]['targets'] for i in active_restraints]

        return number_of_windows, data_points, max_data_points, active_restraints, force_constants, targets, self.simulation_data[
            phase]

    def _run_mbar(self, prepared_data, verbose=False):
        """
        Compute the free energy matrix for a series of windows. We'll follow the pymbar nomenclature for data structures.
        """

        # Unpack the prepared data
        num_win, data_points, max_data_points, active_rest, force_constants, targets, ordered_values = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points)

        # Setup the reduced potential energy array. ie, the potential of each window's
        # coordinates in each window's potential function
        u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)

        # Note, the organization of k = coordinate windows, l = potential windows
        # seems to be opposite of the documentation. But I got wrong numbers the other way around.

        for r, rest in enumerate(active_rest):
            if rest.mask3 is not None:
                force_constants[r] *= (np.pi / 180.0)**2

        for k in range(num_win):  # Coordinate windows
            for l in range(num_win):  # Potential Windows
                force_constants_T = np.asarray(force_constants).T[l, :, None]
                targets_T = np.asarray(targets).T[l, :, None]

                for r, rest in enumerate(active_rest):  # Restraints

                    # If this is a dihedral, we need to shift around restraint value
                    # on the periodic axis to make sure the lowest potential is used.
                    if rest.mask3 is not None and rest.mask4 is not None:
                        target = np.asarray(targets).T[l][r]  # Taken from potential window, l
                        bool_list = ordered_values[k][r] < target - 180.0  # Coords from coord window, k
                        ordered_values[k][r][bool_list] += 360.0
                        bool_list = ordered_values[k][r] > target + 180.0
                        ordered_values[k][r][bool_list] -= 360.0

                # Compute the potential ... for each frame, sum the contributions for each restraint
                # Note, we multiply by beta, and do some extra [l,:,None] to get the math operation correct.

                u_kln[k, l, 0:N_k[k]] = np.sum(
                    self.beta * force_constants_T * (ordered_values[k] - targets_T)**2, axis=0)

        # Setup mbar calc, and get matrix of free energies, uncertainties
        mbar = pymbar.MBAR(u_kln, N_k, verbose=verbose)
        Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)

        # Should I subsample based on the restraint coordinate values? Here I'm
        # doing it on the potential.  Should be pretty close ....
        if 'mbar-block' in self.methods:
            # We want to use all possible data to get the free energy estimates Deltaf_ij,
            # but for uncertainty estimates we'll subsample to create uncorrelated data.
            g_k = np.zeros([num_win], np.float64)
            ss_indices = []
            N_ss = np.zeros([num_win], np.int32)  # N_subsample
            for k in range(num_win):
                l = k
                # If the potential is zero everywhere, we can't estimate the uncertainty, so
                # check the next *potential* window which probably had non-zero force constants
                while not u_kln[k, l, 0:N_k[k]].any():
                    l += 1
                # Now compute statistical inefficiency: g = N*(SEM**2)/variance
                nearest_max = get_nearest_max(N_k[k])
                sem = get_block_sem(u_kln[k, l, 0:nearest_max])
                variance = np.var(u_kln[k, l, 0:N_k[k]])
                g_k[k] = (N_k[k] * (sem**2) / variance)
                # Create subsampled indices and count their lengths
                ss_indices.append(get_subsampled_indices(N_k[k], g_k[k]))
                N_ss[k] = len(ss_indices[k])

            # Create a new potential array for the uncertainty calculation (are we using too much memory?)
            u_kln_err = np.zeros([num_win, num_win, np.max(N_ss)], np.float64)

            # Populate the subsampled array, drawing values from the original
            for k in range(num_win):
                for l in range(num_win):
                    u_kln_err[k, l, 0:N_ss[k]] = u_kln[k, l, ss_indices[k]]

            mbar = pymbar.MBAR(u_kln_err, N_ss, verbose=verbose)
            tmp_Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)

        # Put back into kcal/mol
        Deltaf_ij /= self.beta
        dDeltaf_ij /= self.beta

        # Return Matrix of free energies and uncertainties
        return Deltaf_ij, dDeltaf_ij

    def _run_ti(self, phase, prepared_data):
        """
        Compute the free energy using the TI method.
        """

        # Unpack the prepared data
        num_win, data_points, max_data_points, active_rest, force_constants, targets, ordered_values = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points)

        # Setup the dU array. ie, the partial derivative of the potential with respect force constant or target,
        # depending on the whether attach/release or pull.
        dU = np.zeros([max_data_points], np.float64)

        # The mean/sem dU value for each window
        dU_avgs = np.zeros([num_win], np.float64)
        dU_sems = np.zeros([num_win], np.float64)

        # Convert the force_constants to degrees
        for r, rest in enumerate(active_rest):
            if rest.mask3 is not None:
                force_constants[r] *= (np.pi / 180.0)**2

        # Deal with dihedral wrapping
        for k in range(num_win):  # Coordinate windows
            force_constants_T = np.asarray(force_constants).T[k, :, None]
            targets_T = np.asarray(targets).T[k, :, None]

            for r, rest in enumerate(active_rest):  # Restraints

                # If this is a dihedral, we need to shift around restraint value
                # on the periodic axis to make sure the lowest potential is used.
                if rest.mask3 is not None and rest.mask4 is not None:
                    target = np.asarray(targets).T[k][r]
                    bool_list = ordered_values[k][r] < target - 180.0
                    ordered_values[k][r][bool_list] += 360.0
                    bool_list = ordered_values[k][r] > target + 180.0
                    ordered_values[k][r][bool_list] -= 360.0

            # Compute the partial derivative (ie forces) for each frame. The force constants
            # are scaled by a lambda (l) parameter which controls their strength: 0 to max.
            # Potential:         U = l*fc_max*(values - targets)**2
            # Forces Attach: dU/dl = 
            
            if phase == 'attach' or phase == 'release':
                dU[0:N_k[k]] = np.sum(force_constants_T * (ordered_values[k] - targets_T)**2, axis=0)
            else:
                dU[0:N_k[k]] = np.sum(2.0 * force_constants_T * (ordered_values[k] - targets_T), axis=0)

            dU_avgs[k] = np.mean( dU[0:N_k[k]] )
            dU_sems[k] = get_block_sem( dU[0:N_k[k]] )



    def compute_free_energy(self):
        """
        Do free energy calc.
        """

        for phase in ['attach', 'pull', 'release']:
            self.results[phase] = {}
            for method in self.methods:
                # Initialize some values that we will compute
                # The matrix gives all possible fe/sem for any window to any other window
                self.results[phase][method] = {}
                self.results[phase][method]['fe'] = None
                self.results[phase][method]['sem'] = None
                self.results[phase][method]['fe_matrix'] = None
                self.results[phase][method]['sem_matrix'] = None

                # mbar with blocking are currently supported.
                if method == 'mbar-block':
                    if sum(self.changing_restraints[phase]) == 0:
                        log.debug('Skipping free energy calculation for %s' % phase)
                        break
                    prepared_data = self._prepare_data(phase)
                    self.results[phase][method]['n_frames'] = np.sum(prepared_data[1])
                    self.results[phase][method]['fe_matrix'],self.results[phase][method]['sem_matrix']\
                        = self._run_mbar(prepared_data)
                    self.results[phase][method]['fe'] = self.results[phase][method]['fe_matrix'][0, -1]
                    self.results[phase][method]['sem'] = self.results[phase][method]['sem_matrix'][0, -1]

                    windows = len(self.results[phase][method]['sem_matrix'])
                    self.results[phase][method]['convergence'] = np.ones([windows], np.float64) * -1.0
                    self.results[phase][method]['ordered_convergence'] = np.ones([windows], np.float64) * -1.0
                    log.info(phase + ': computing convergence for mbar-blocking')
                    for i in range(windows):
                        if i == 0:
                            self.results[phase][method]['ordered_convergence'][i]\
                                = self.results[phase][method]['sem_matrix'][i][i+1]
                        elif i == windows - 1:
                            self.results[phase][method]['ordered_convergence'][i]\
                                = self.results[phase][method]['sem_matrix'][i][i-1]
                        else:
                            left = self.results[phase][method]['sem_matrix'][i][i - 1]
                            right = self.results[phase][method]['sem_matrix'][i][i + 1]
                            if left > right:
                                max_val = left
                            elif right > left:
                                max_val = right
                            else:
                                max_val = right
                            self.results[phase][method]['ordered_convergence'][i] = max_val

                    self.results[phase][method]['convergence'] = \
                        [self.results[phase][method]['ordered_convergence'][i] for i in self.orders[phase]]


    def compute_ref_state_work(self, restraints):
        """
        Compute the work to place a molecule at standard reference state conditions
        starting from a state defined by up to six restraints. (see ref_state_work for
        details)

        Parameters
        ----------
        restraints : list [r, theta, phi, alpha, beta, gamma]
            A list of paprika DAT_restraint objects in order of the six translational and
            orientational restraints needed to describe the configuration of one molecule
            relative to another. The six restraints are: r, theta, phi, alpha, beta, gamma.
            If any of these coordinates is not being restrained, use a None in place of a
            DAT_restraint object. (see ref_state_work for details on the six restraints)
        """

        if not restraints or restraints[0] is None:
            raise Exception('At minimum, a single distance restraint must be supplied to compute_ref_state_work')

        fcs = []
        targs = []

        for restraint in restraints:
            if restraint is None:
                fcs.append(None)
                targs.append(None)
            elif restraint.phase['release']['force_constants'] is not None:
                fcs.append( np.sort(restraint.phase['release']['force_constants'])[-1] )
                targs.append( np.sort(restraint.phase['release']['targets'])[-1] )
            elif restraint.phase['pull']['force_constants'] is not None:
                fcs.append( np.sort(restraint.phase['pull']['force_constants'])[-1] )
                targs.append( np.sort(restraint.phase['pull']['targets'])[-1] )
            else:
                raise Exception('Restraints should have pull or release values initialized in order to compute_ref_state_work')

        # Convert degrees to radians for theta, phi, alpha, beta, gamma
        for i in range(1,5):
            if targs[i] is not None:
                targs[i] = np.radians(targs[i])


        self.results['ref_state_work'] = ref_state_work(self.temperature,
                                                        fcs[0], targs[0],
                                                        fcs[1], targs[1],
                                                        fcs[2], targs[2],
                                                        fcs[3], targs[3],
                                                        fcs[4], targs[4],
                                                        fcs[5], targs[5])


def get_factors(n):
    """
    Return a list of integer factors for a number.
    """
    factors = []
    sqrt_n = int(round(np.sqrt(n) + 0.5))
    i = 1
    while i <= sqrt_n:
        if n % i == 0:
            factors.append(int(i))
            j = n / i
            if j != i:
                factors.append(int(j))
        i += 1
    return sorted(factors, key=int)


def get_nearest_max(n):
    """
    Return the number with the largest number of factors between n-100 and n.
    """
    num_factors = []
    max_factors = 0
    if n % 2 == 0:
        beg = n - 100
        end = n
    else:
        beg = n - 101
        end = n - 1
    if beg < 0:
        beg = 0
    for i in range(beg, end + 2, 2):
        num_factors = len(get_factors(i))
        if num_factors >= max_factors:
            max_factors = num_factors
            most_factors = i
    return most_factors


def get_block_sem(data_array):
    """
    Compute the standard error of the mean (SEM) for a data_array using the blocking method."
    """
    # Get the integer factors for the number of data points. These
    # are equivalent to the block sizes we will check.
    block_sizes = get_factors(len(data_array))

    # An array to store means for each block ... make it bigger than we need.
    block_means = np.zeros([block_sizes[-1]], np.float64)

    # Store the SEM for each block size, except the last two size for which
    # there will only be two or one blocks total and thus very noisy.
    sems = np.zeros([len(block_sizes) - 2], np.float64)

    # Check each block size except the last two.
    for size_idx in range(len(block_sizes) - 2):
        # Check each block, the number of which is conveniently found as
        # the other number of the factor pair in block_sizes
        num_blocks = block_sizes[-size_idx - 1]
        for blk_idx in range(num_blocks):
            # Find the index for beg and end of data points for each block
            data_beg_idx = blk_idx * block_sizes[size_idx]
            data_end_idx = (blk_idx + 1) * block_sizes[size_idx]
            # Compute the mean of this block and store in array
            block_means[blk_idx] = np.mean(data_array[data_beg_idx:data_end_idx])
        # Compute the standard deviation across all blocks, devide by num_blocks-1 for SEM
        sems[size_idx] = np.std(block_means[0:num_blocks], ddof=0) / np.sqrt(num_blocks - 1)
        # Hmm or should ddof=1? I think 0, see Flyvbjerg -----^

    # Return the max SEM found ... this is a conservative approach.
    return np.max(sems)


def get_subsampled_indices(N, g, conservative=False):
    """ Get subsampling indices. Adapted from pymbar's implementation. """

    # g should not be less than 1.0
    if g < 1.0:
        g = 1.0

    # if conservative, assume integer g and round up
    if conservative:
        g = np.ceil(g)

    # initialize
    indices = [0]
    g_idx = 1.0
    int_step = int(np.round(g_idx * g))

    while int_step < N:
        indices.append(int_step)
        g_idx += 1
        int_step = int(np.round(g_idx * g))

    return indices


def read_restraint_data(restraint, window, trajectory, prmtop, single_prmtop=False, fraction=1.0):
    """Given a trajectory (or trajectories) and restraint, read the restraint values.

    Note this is *slow* because it will load the trajectory for *each* restraint. This is done on purpose,
    so it is easier to debug and follow the code. It also makes it easier to package the data for MBAR by 
    logically separating simulation windows and simulation restraints.
    
    Parameters:
    ----------
    restraint : {DAT_restraint}
        The restraint to analyze
    window : {str}
        The simulation window to analyze
    trajectory : {str} or {list}
        The name or names of the trajectory
    prmtop : {str} or ParmEd AmberParm
        The parameters for the simulation
    single_prmtop : {bool}
        Whether a single `prmtop` is read for all windows
    fraction : {float}
        Fraction of data to read, to check free energy convergence
    Returns
    -------
    data : {np.array}
        The values for this restraint in this window
    """

    log.debug('Reading restraint data for {}...'.format(window, trajectory))
    if isinstance(trajectory, str):
        trajectory_path = os.path.join(window, trajectory)
    elif isinstance(trajectory, list):
        trajectory_path = [os.path.join(window, i) for i in trajectory]
        log.debug('Received list of trajectories: {}'.format(trajectory_path))

    if isinstance(prmtop, str) and not single_prmtop:
        traj = pt.iterload(trajectory_path, os.path.join(window, prmtop))
    elif isinstance(prmtop, str) and single_prmtop:
        traj = pt.iterload(trajectory_path, os.path.join(prmtop))
    else:
        try:
            traj = pt.iterload(trajectory_path, prmtop)
        except:
            raise Exception('Tried to load `prmtop` object directly and failed.')

    if fraction > 1:
        raise Exception('The fraction of data to analyze cannot be greater than 1.')
    elif np.isclose(fraction, 1):
        pass
    else:
        log.debug('Loaded {} frames...'.format(traj.n_frames))
        traj = traj[0:int(fraction * traj.n_frames)]
        log.debug('Analyzing {} frames...'.format(traj.n_frames))

    if restraint.mask1 and restraint.mask2 and \
            not restraint.mask3 and not restraint.mask4:
        data = pt.distance(traj, ' '.join([restraint.mask1, restraint.mask2]))
    elif restraint.mask1 and restraint.mask2 and \
            restraint.mask3 and not restraint.mask4:
        data = pt.angle(traj, ' '.join([restraint.mask1, restraint.mask2, restraint.mask3]))
    elif restraint.mask1 and restraint.mask2 and \
            restraint.mask3 and restraint.mask4:
        data = pt.dihedral(traj, ' '.join([restraint.mask1, restraint.mask2, \
                                           restraint.mask3, restraint.mask4]))
    return data


def ref_state_work(temperature,
                   r_fc,  r_tg,
                   th_fc, th_tg,
                   ph_fc, ph_tg,
                   a_fc,  a_tg,
                   b_fc,  b_tg,
                   g_fc,  g_tg):
    """
    Computes the free energy to release a molecule from some restrained translational
    and orientational configuration (relative to another molecule or lab frame) into
    the reference configuration: standard concentration (1.0/1660.5392 Angstom^3) and
    unrestrained orientational freedom (8*pi^2).  (See Euler angles)

    Assume two molecules (H and G). Three translational and three orientational degrees
    of freedom define their relative configuration. In order to match experimentally
    reported free energies, we often need to compute the work (free energy) of moving
    a molecule (G) from a restrained configuration relative to H into the experimental
    reference state (usually defined as 1 M, or 1 molecule per 1660.5 Angstrom^3)

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
    temperature : float
        temperature (K) at which the reference state calculation will take place
    r_fc, r_tg : float, float
        The distance, r, restraint force constant and target values (Angstrom, kcal/mol-Angstom^2). The target range lies within 0 to infinity.
    th_fc, th_tg : float, float
        The angle, theta, restraint force constant and target values (radian, kcal/mol-radian^2). The target range is 0 to pi.
    ph_fc, ph_tg : float, float
        The torsion, phi, restraint force constant and target values (radian, kcal/mol-radian^2). The target range is 0 to 2*pi.
    a_fc, a_tg : float, float
        The torsion, alpha, restraint force constant and target values (radian, kcal/mol-radian^2). The target range is 0 to 2*pi.
    b_fc, b_tg : float, float
        The angle, beta, restraint force constant and target values (radian, kcal/mol-radian^2). The target range is 0 to pi.
    g_fc, g_tg : float, float
        The angle, gamma, restraint force constant and target values (radian, kcal/mol-radian^2). The target range is 0 to 2*pi.

    """

    R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
    RT = R*temperature

    # Distance Integration Function
    def dist_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return (arange**2) * np.exp( (-1.0/RT)*fc*(arange - targ)**2 )
        arange = np.arange(0.0, 100.0, 0.0001)
        return np.trapz( potential(arange, RT, fc, targ), arange )

    # Angle Integration Function
    def ang_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return np.sin(arange) * np.exp( (-1.0/RT)*fc*(arange - targ)**2 )
        arange = np.arange(0.0, np.pi, 0.00005)
        return np.trapz( potential(arange, RT, fc, targ), arange )

    # Torsion Integration Function
    def tors_int(RT, fc, targ):
        def potential(arange, RT, fc, targ):
            return np.exp( (-1.0/RT)*fc*(arange - targ)**2 )
        # Note, because of periodicity, I'm gonna wrap +/- pi around targ for integration
        arange = np.arange(targ-np.pi, targ+np.pi, 0.00005)
        return np.trapz( potential(arange, RT, fc, targ), arange )

    # Distance restraint, r
    if None in [r_fc, r_tg]:
        raise Exception('Distance restraint info (r_fc, r_tg) must be specified')
    else:
        r_int  = dist_int( RT, r_fc,  r_tg  )

    # Angle restraint, theta
    if None in [th_fc, th_tg]:
        th_int = 2.0
    else:
        th_int = ang_int(  RT, th_fc, th_tg )

    # Torsion restraint, phi
    if None in [ph_fc, ph_tg]:
        ph_int = 2.0*np.pi
    else:
        ph_int = tors_int( RT, ph_fc, ph_tg )

    # Torsion restraint, alpha
    if None in [a_fc, a_tg]:
        a_int = 2.0*np.pi
    else:
        a_int  = tors_int( RT, a_fc,  a_tg  )

    # Angle restraint, beta
    if None in [b_fc, b_tg]:
        b_int = 2.0
    else:
        b_int  = ang_int(  RT, b_fc,  b_tg  )

    # Torsion restraint, gamma
    if None in [g_fc, g_tg]:
        g_int = 2.0*np.pi
    else:
        g_int  = tors_int( RT, g_fc,  g_tg  )

    # Concentration term
    trans = r_int * th_int * ph_int * (1.0/1660.5392)   # C^o = 1/V^o
    # Orientational term
    orient = a_int * b_int * g_int / (8.0*(np.pi**2))

    # Return the free energy
    return RT*np.log( trans * orient )

