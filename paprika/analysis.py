import logging as log
import os as os
import numpy as np
import pytraj as pt
import pymbar

class fe_calc(object):
    """
    Computes the free energy for an APR transformation.
    """

    def __init__(self):

        self.temperature = 298.15
        self.k_b = 0.0019872041   # kcal/mol-K
        self.beta = 1/(self.k_b*self.temperature)  # Add auto updating for beta

        # List of DAT_restraints
        self.restraint_list = None

        # All raw restraint values. Organized as:
        # simulation_values[phase][0:num_win][0:num_rest][np.array]
        self.simulation_values = None  # Add a function to create this automatically

        # FE/Uncertainty methods.
        self.methods = ['mbar-block'] # mbar-autoc, mbar-none, ti-block, ti-autoc, ti-none
        # TODO: Add check that fe_methods and subsample_methods have correct keywords

        # Keep track of the order of increasing force_constants/targets for each
        # phase. This helps in cases where reordering is required due to post-hoc
        # window additions.
        self.ordered_index = {}

        # FE calculation results will be stored here
        self.results = {}

    def _identify_active_rest(self, phase, change_param, restraint_list):
        """ Identify the restraints which are changing in the specified phase."""

        # This is a verbose approach, but it got messy when trying to be more pythonic.

        # The list of restraints which are changing
        active_rest = []
        # A boolean of which restriants in restraint_list are changing
        active_rest_bool = []

        # Iterate through the restraints and, for attach/release, check if force_constants are changing,
        # or, for pull, check if targets are changing
        for rest in restraint_list:
            if rest.phase[phase][change_param] is not None:
                # These is probably a better way to do this ... but
                # create an array in which the first element is repeated
                test_val = rest.phase[phase][change_param][0]
                test_len = rest.phase[phase][change_param].size
                test_arr = test_val*np.ones([test_len])
                # Then check if it is identical to the real array ...
                if np.allclose(rest.phase[phase][change_param], test_arr):
                    # If yes, then it's not changing, and thus not active
                    active_rest_bool.append(False)
                else:
                    # Add to active restraint list
                    active_rest.append(rest)
                    active_rest_bool.append(True)
            else:
                # If nothing is set for this phase, it's definitely not active
                active_rest_bool.append(False)

        return active_rest, active_rest_bool


    def _prepare_active_rest_data(self, phase, restraint_list, simulation_values):
        """ 
        Identifies the restraints which are changing in the specified phase
        and reorders windows numerically, in case additional
        windows are added out of order after an initial run.
        """

        # Identify the restraint parameter which is expected to change for this phase
        if phase in ['attach', 'release']:
            change_param = 'force_constants'
        elif phase == 'pull':
            change_param = 'targets'
        else:
            raise Exception('Invalid phase specified: '+phase+'. Should be either attach, pull, or release')
    
        # First we identify the restraints which are changing in the specified phase.
        # active_rest is the list of active restraints, active_rest_bool is bool of
        # active restraints in the full restraint_list.
        active_rest, active_rest_bool = self._identify_active_rest(phase, change_param, restraint_list)

        # If we didn't find any active restraints, we're done
        if len(active_rest) == 0:
            return None

        # Create a window index list which matches numerical ascending order for the change_param
        reorder = np.argsort(active_rest[0].phase[phase][change_param])
        self.ordered_index[phase] = reorder

        # If continuous_apr and attach/release, we need to know the order of pull too
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            pull_active_rest, pull_active_rest_bool = self._identify_active_rest('pull', 'targets', restraint_list)
            if len(pull_active_rest) == 0:
                raise Exception('There should be active pull restraints, but none were found!')
            pull_reorder =  np.argsort(pull_active_rest[0].phase['pull']['targets'])
            pull_first_active_idx = next((i for i, rbool in enumerate(pull_active_rest_bool) if rbool is True), None)

        # Count number of windows
        num_win = len(reorder)

        # Count the number of active restraints
        num_rest = len(active_rest)

        # Identify the index of the first active restraint in the full
        # restraint_list by checking active_rest_bool. Probably don't need
        # to do this but it ensures we're checking the size of an active restraint
        first_active_idx = next((i for i, rbool in enumerate(active_rest_bool) if rbool is True), None)

        # Now count the size of data_points for each np.array in simulation_values.
        # We only need to check size for first_active_idx, since they should all be same.
        # If continuous_apr and attach/release, we'll need to take the data
        # for the final attach/release window from the first/last pull window.
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            data_points = np.zeros([num_win], np.int32)
            for win_idx in reorder[:-1]:
                data_points[win_idx] = simulation_values[phase][win_idx][first_active_idx].size
            # If attach, then take data from first pull window
            if phase == 'attach':
                pull_idx = 0
            # Else, must be release, so take from last pull window
            else:
                pull_idx = -1
            data_points[-1] = simulation_values['pull'][pull_reorder[pull_idx]][pull_first_active_idx].size
        # If not continuous_apr or pull, then it's more straightforward
        else:
            data_points = np.array([simulation_values[phase][i][first_active_idx].size for i in reorder])

        # Get max_data_points ... will be needed for MBAR
        max_data_points = int(np.max(data_points))

        # Now create properly ordered force_constants and targets.
        force_constants = np.zeros([num_win, num_rest], np.float64)
        targets = np.zeros([num_win, num_rest], np.float64)
        for win_idx in reorder:
            for rest_idx in range(num_rest):
                force_constants[win_idx, rest_idx] = active_rest[rest_idx].phase[phase]['force_constants'][win_idx]
                if active_rest[rest_idx].mask3 is not None:
                    # Convert force constants with radians into degrees
                    force_constants[win_idx, rest_idx] *= ( (np.pi/180.0)**2 )
                targets[win_idx, rest_idx] = active_rest[rest_idx].phase[phase]['targets'][win_idx]

        # Also reorder the data to match the force_constants, targets.
        # Again, if attach/release, must take last window from pull
        ordered_values = [[] for i in range(num_win)]
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            for win_idx in reorder[:-1]:
                for rest_idx,rest_bool in enumerate(active_rest_bool):
                    if rest_bool:
                        ordered_values[win_idx].append(simulation_values[phase][win_idx][rest_idx])
            for rest_idx,rest_bool in enumerate(active_rest_bool):
                if rest_bool:
                    # Same logic as above ... must figure out which pull window to use as
                    # final attach/release window.
                    if phase == 'attach':
                        pull_idx = 0
                    else:
                        pull_idx = -1
                    ordered_values[-1].append(simulation_values['pull'][pull_reorder[pull_idx]][rest_idx])
        else:
            for win_idx in reorder:
                for rest_idx,rest_bool in enumerate(active_rest_bool):
                    if rest_bool:
                        ordered_values[win_idx].append(simulation_values[phase][win_idx][rest_idx])

        # Return all this in a list. Will pass to _run_mbar, which unpacks.
        # That seems inefficient, but we'll keep it for now.
        return [num_win, num_rest, data_points, max_data_points, active_rest, force_constants, targets, ordered_values]


    def _run_mbar(self, prepared_data, verbose=False):
        """
        Compute the free energy matrix for a series of windows. We'll follow the pymbar nomenclature for data structures.
        """

        # Unpack the prepared data
        num_win, num_rest, data_points, max_data_points, active_rest, force_constants, targets, ordered_values = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points, np.int32)

        # Setup the reduced potential energy array. ie, the potential of each window's
        # coordinates in each window's potential function
        u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)

        # Note, the organization of k = coordinate windows, l = potential windows
        # seems to be opposite of the documentation. But I got wrong numbers the other way around.
        for k in range(num_win):  # Coordinate windows
            for l in range(num_win): # Potential Windows
                for r,rest in enumerate(active_rest): # Restraints
                    # If this is a dihedral, we need to shift around restraint value
                    # on the periodic axis to make sure the lowest potential is used.
                    if rest.mask3 is not None and rest.mask4 is not None:
                        target = targets[l][r] # Taken from potential window, l
                        bool_list = ordered_values[k][r] < target - 180.0 # Coords from coord window, k
                        ordered_values[k][r][bool_list] += 360.0
                        bool_list = ordered_values[k][r] > target + 180.0
                        ordered_values[k][r][bool_list] -= 360.0

                # Compute the potential ... for each frame, sum the contributions for each restraint
                # Note, we multiply by beta, and do some extra [l,:,None] to get the math operation correct.
                u_kln[k,l,0:N_k[k]] = np.sum(self.beta*force_constants[l,:,None]*(ordered_values[k] - targets[l,:,None])**2, axis=0)

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
            N_ss = np.zeros([num_win], np.int32) # N_subsample
            for k in range(num_win):
                l = k
                # If the potential is zero everywhere, we can't estimate the uncertainty, so
                # check the next *potential* window which probably had non-zero force constants
                while not u_kln[k,l,0:N_k[k]].any():
                    l += 1
                # Now compute statistical inefficiency: g = N*(SEM**2)/variance
                nearest_max = get_nearest_max(N_k[k]) 
                sem = get_block_sem(u_kln[k,l,0:nearest_max])
                variance = np.var( u_kln[k,l,0:N_k[k]] )
                g_k[k] = ( N_k[k]*(sem**2)/variance )
                # Create subsampled indices and count their lengths
                ss_indices.append( get_subsampled_indices(N_k[k], g_k[k]) )
                N_ss[k] = len(ss_indices[k])

            # Create a new potential array for the uncertainty calculation (are we using too much memory?)
            u_kln_err = np.zeros([num_win, num_win, np.max(N_ss)], np.float64)

            # Populate the subsampled array, drawing values from the original
            for k in range(num_win):
                for l in range(num_win):
                    u_kln_err[k,l,0:N_ss[k]] = u_kln[k,l,ss_indices[k]]

            mbar = pymbar.MBAR(u_kln_err, N_ss, verbose=verbose)
            tmp_Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)
                    
        # Put back into kcal/mol
        Deltaf_ij /= self.beta
        dDeltaf_ij /= self.beta

        # Return Matrix of free energies and uncertainties
        return Deltaf_ij, dDeltaf_ij

    def compute_free_energy(self):
        """
        Do free energy calc.
        """

        for phase in ['attach', 'pull', 'release']:
            self.results[phase] = {}
            for method in self.methods:
                self.results[phase][method] = {}
                self.results[phase][method]['fe'] = None
                self.results[phase][method]['sem'] = None
                self.results[phase][method]['fe_matrix'] = None
                self.results[phase][method]['sem_matrix'] = None
                    
                # mbar with blocking are currently supported.
                if method == 'mbar-block':
                    prepared_data = self._prepare_active_rest_data(phase, self.restraint_list, self.simulation_values)
                    if prepared_data:
                        self.results[phase][method]['fe_matrix'],self.results[phase][method]['sem_matrix']\
                            = self._run_mbar(prepared_data)
                        self.results[phase][method]['fe'] = self.results[phase][method]['fe_matrix'][0,-1]
                        self.results[phase][method]['sem'] = self.results[phase][method]['sem_matrix'][0,-1]

                        windows = len(self.results[phase][method]['sem_matrix'])
                        self.results[phase][method]['convergence'] = np.ones([windows], np.float64)*-1.0
                        self.results[phase][method]['ordered_convergence'] = np.ones([windows], np.float64)*-1.0
                        log.info(phase+': computing convergence for mbar-blocking')
                        for i in range(windows):
                            if i == 0:
                                self.results[phase][method]['ordered_convergence'][i]\
                                    = self.results[phase][method]['sem_matrix'][i][i+1]
                            elif i == windows-1:
                                self.results[phase][method]['ordered_convergence'][i]\
                                    = self.results[phase][method]['sem_matrix'][i][i-1]
                            else:
                                left = self.results[phase][method]['sem_matrix'][i][i-1]
                                right = self.results[phase][method]['sem_matrix'][i][i+1]
                                if left > right:
                                    max_val = left
                                elif right > left:
                                    max_val = right
                                else:
                                    max_val = right
                                self.results[phase][method]['ordered_convergence'][i]\
                                    = max_val

                        # Un-reorder so that convergence easily matches up with original window order
                        unreorder = np.argsort(self.ordered_index[phase])
                        self.results[phase][method]['convergence'] =\
                            self.results[phase][method]['ordered_convergence'][unreorder]

                            

### Additional Functions

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
        j = n/i
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
        beg = n-101
        end = n-1
    if beg < 0:
        beg = 0
    for i in range(beg, end + 2, 2):
        num_factors = len( get_factors(i) )
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
    block_sizes = get_factors( len(data_array) )

    # An array to store means for each block ... make it bigger than we need.
    block_means = np.zeros( [block_sizes[-1]], np.float64 )

    # Store the SEM for each block size, except the last two size for which
    # there will only be two or one blocks total and thus very noisy.
    sems = np.zeros( [len(block_sizes) - 2], np.float64 )

    # Check each block size except the last two.
    for size_idx in range(len(block_sizes) - 2):
        # Check each block, the number of which is conveniently found as
        # the other number of the factor pair in block_sizes
        num_blocks = block_sizes[-size_idx - 1]
        for blk_idx in range(num_blocks):
            # Find the index for beg and end of data points for each block
            data_beg_idx = blk_idx*block_sizes[size_idx]
            data_end_idx = (blk_idx+1)*block_sizes[size_idx]
            # Compute the mean of this block and store in array
            block_means[blk_idx] = np.mean( data_array[ data_beg_idx : data_end_idx ] )
        # Compute the standard deviation across all blocks, devide by num_blocks-1 for SEM
        sems[size_idx] = np.std( block_means[0:num_blocks], ddof=0 ) / np.sqrt( num_blocks - 1 )
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
    int_step = int( np.round( g_idx*g ) )

    while int_step < N:
        indices.append(int_step)
        g_idx += 1
        int_step = int( np.round( g_idx*g ) )

    return indices



def collect_data(restraint_list):
    # determine_static_restraints(restraint_list)
    # determine_window_order(restraint_list)
    # read_trajectories(restraint_list)
    pass

def determine_static_restraints(restraint_list):

    # Actually, first identify the changing restraints...
    changing_restraints = {
        'attach' : [],
        'pull' : [],
        'release' : []
    }

    for phase in ['attach', 'pull', 'release']:
        if phase == 'attach' or phase == 'release':
            changing_parameter = 'force_constants'
        else:
            changing_parameter = 'targets'
        for restraint in restraint_list:
            if restraint.phase[phase][changing_parameter] is not None:
                static = all(np.isclose(x, restraint.phase[phase][changing_parameter][0]) for x in
                    restraint.phase[phase][changing_parameter]
                    )
            else:
                static = True

            changing_restraints[phase].append(not static)

    return changing_restraints

def determine_window_order(restraint_list, changing_restraints):
    orders = {
        'attach' : [],
        'pull' : [],
        'release' : []
    }
    active_attach_restraints = np.asarray(restraint_list)[changing_restraints['attach']]
    active_pull_restraints = np.asarray(restraint_list)[changing_restraints['pull']]
    active_release_restraints = np.asarray(restraint_list)[changing_restraints['release']]

    attach_orders = []
    pull_orders = []
    release_orders = []

    for restraint in active_attach_restraints:
        attach_orders.append(np.argsort(restraint.phase['attach']['force_constants']))
    if not all([np.array_equal(attach_orders[0], i) for i in attach_orders]):
        raise Exception
    elif attach_orders:
        orders['attach'] = attach_orders[0]
    else:
        orders['attach'] = []
    for restraint in active_pull_restraints:
        pull_orders.append(np.argsort(restraint.phase['pull']['targets']))
    if not all([np.array_equal(pull_orders[0], i) for i in pull_orders]):
        raise Exception
    elif pull_orders:
        orders['pull'] = pull_orders[0]
    else:
        orders['pull'] = []

    # Still a problem here, with no release windows.
    for restraint in active_release_restraints:
        release_orders.append(np.argsort(restraint.phase['release']['force_constants']))
    if not all([np.array_equal(release_orders[0], i) for i in release_orders]):
        raise Exception
    elif release_orders:
        orders['release'] = release_orders[0]
    else:
        orders['release'] = []
    # Make sure each of these is a list with the same order -- i.e., that each restraint has the order the same.
    # all([np.array_equal(attach_orders[0], i) for i in attach_orders])
    # Then, return the ordering, perhaps concatenated together.
    return orders


def read_trajectories(restraint_list, changing_restraints, order, prmtop, traj, path='./',
                      strip=None, inpcrd=None):

    # Map between the ordering and the windows.
    simulation_data = {
        'attach' : [],
        'pull'   : [],
        'release': []
    }

    ordered_windows = \
        [os.path.join('./', 'a{:03d}'.format(i)) for i in order['attach'] if i] + \
        [os.path.join('./', 'p{:03d}'.format(i)) for i in order['pull'] if i] + \
        [os.path.join('./', 'r{:03d}'.format(i)) for i in order['release'] if i]

    active_attach_restraints = np.asarray(restraint_list)[changing_restraints['attach']]
    active_pull_restraints = np.asarray(restraint_list)[changing_restraints['pull']]
    active_release_restraints = np.asarray(restraint_list)[changing_restraints['release']]

    active_restraints = [active_attach_restraints, active_pull_restraints, active_release_restraints]

    for restraint_index, restraint in enumerate(active_restraints):
        log.debug('Restraint = {}'.format(restraint_index))
        for window_index, window in enumerate(ordered_windows):
            # Each relevant window gets its own list inside `simulation_data[phase]`
            if restraint in active_attach_restraints:
                phase = 'attach'
            elif restraint in active_pull_restraints:
                phase = 'pull'
            elif restraint in active_release_restraints:
                phase = 'release'
            else:
                raise Exception
            simulation_data[phase].append([])
            log.debug('Window = {}'.format(window_index))
            simulation_data[phase][window_index].append([])

            if strip:
                structure = pt.load(os.path.join(window, inpcrd),
                                    os.path.join(window, prmtop))

                stripped = structure.strip(':WAT,:Na+,:Cl-')

                traj = pt.iterload(os.path.join(window, traj),
                                   top=stripped.topology)
            else:
                traj = pt.iterload(os.path.join(window, traj),
                                   os.path.join(window, prmtop))

            if restraint.mask1 and restraint.mask2 and \
                    not restraint.mask3 and not restraint.mask4:
                data = pt.distance(traj, ' '.join([restraint.mask1, restraint.mask2]))
            elif restraint.mask1 and restraint.mask2 and \
                    restraint.mask3 and not restraint.mask4:
                data = pt.angle(traj, ' '.join([restraint.mask1, restraint.mask2, restraint.mask3]))
            elif restraint.mask1 and restraint.mask2 and \
                    restraint.mask3 and restraint.mask4:
                data = pt.dihedral(traj, ' '.join([restraint.mask1, restraint.mask2, restraint.mask3, restraint.mask4]))

            simulation_data[phase][window_index][restraint_index] = data

    return simulation_data