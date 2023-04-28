import logging
import os
import traceback

import numpy as np
import pytraj as pt
from openff.units import unit as openff_unit
from scipy import stats as statistics

logger = logging.getLogger(__name__)


def get_factors(n):
    """
    Return a list of integer factors for a number.

    Parameters
    ----------
    n: int or float
        Number to factor

    Returns
    -------
    sorted: list
        A list of sorted factors.

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
    Return the number with the largest number of factors between n − 100 and n.

    Parameters
    ----------
    n: int
        Desired number to factor.

    Returns
    -------
    most_factors: int
        The number with the most factors.

    """
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
    Compute the standard error of the mean (SEM) using the blocking method.

    Note
    ----
        This is a conservative approach. Here, we report the maximum SEM determined from blocking analysis (cf. the
        "plateau" on the blocking curve).

    Parameters
    ----------
    data_array: :class:`np.array`
        Array containing data values.

    Returns
    -------
    np.max(sems): float
        The maximum SEM obtained from te blocking curve.

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
        # Compute the standard deviation across all blocks, devide by
        # num_blocks-1 for SEM
        sems[size_idx] = np.std(block_means[0:num_blocks], ddof=0) / np.sqrt(
            num_blocks - 1
        )
        # Hmm or should ddof=1? I think 0, see Flyvbjerg -----^

    return np.max(sems)


def get_subsampled_indices(N, g, conservative=False):
    """Get the indices of independent (subsampled) frames. This is adapted from the implementation in `pymbar`.

    Parameters
    ----------
    N: int
        The length of the array to be indexed.
    g: int
        The statistical inefficiency of the data.
    conservative: bool, optional, default=False
        Whether `g` should be rounded up to the nearest integer.

    Returns
    -------
    indices: list
        A list of indices that can be used to pull out de-correlated frames from a time series.

    """

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
        g_idx += 1.0
        int_step = int(np.round(g_idx * g))

    return indices


def load_trajectory(window, trajectory, topology, single_topology=False):
    """Load a trajectory (or trajectories) and return a pytraj ``trajectory`` object.

    Parameters
    ----------
    window: str
        The simulation window to analyze
    trajectory: str or list
        The name or names of the trajectory
    topology: str or :class:`parmed.Structure`
        The topology the simulation
    single_topology: bool
        Whether a single topology is read for all windows

    Returns
    -------
    traj: pytraj.trajectory
        The trajectory of stored as a pytraj object.
    """

    logger.debug("Load trajectories from {}/{}...".format(window, trajectory))
    if isinstance(trajectory, str):
        trajectory_path = os.path.join(window, trajectory)
    elif isinstance(trajectory, list):
        trajectory_path = [os.path.join(window, i) for i in trajectory]
        logger.debug("Received list of trajectories: {}".format(trajectory_path))
    else:
        raise RuntimeError("Trajectory path should be a `str` or `list`.")

    if isinstance(topology, str) and not single_topology:
        if not os.path.isfile(os.path.join(window, topology)):
            raise FileNotFoundError(
                f"Cannot find `topology` file: {os.path.join(window, topology)}"
            )
        logger.debug(f"Loading {os.path.join(window, topology)} and {trajectory_path}")
        try:
            traj = pt.iterload(trajectory_path, os.path.join(window, topology))
        except ValueError as e:
            formatted_exception = traceback.format_exception(None, e, e.__traceback__)
            logger.info(
                f"Failed trying to load {os.path.join(window, topology)} and {trajectory_path}: "
                f"{formatted_exception}"
            )
    elif isinstance(topology, str) and single_topology:
        traj = pt.iterload(trajectory_path, os.path.join(topology))
    else:
        try:
            traj = pt.iterload(trajectory_path, topology)
        except BaseException:
            raise Exception("Tried to load `topology` object directly and failed.")

    logger.debug("Loaded {} frames...".format(traj.n_frames))

    return traj


def read_restraint_data(
    traj, restraint, distance_unit=openff_unit.angstrom, angle_unit=openff_unit.degree
):
    """Given a trajectory and restraint, read the restraint and return the DAT values.

    Parameters
    ----------
    traj: :class:`pytraj.trajectory`
        A trajectory, probably loaded by load_trajectory
    restraint: :class:`DAT_restraint`
        The restraint to analyze
    distance_unit: openff.unit.Quantity
        The unit for the returned distance values
    angle_unit: openff.unit.Quantity
        The unit for the returned angle values

    Returns
    -------
    data: :class:`np.array`
        The values for this restraint in this window
    """

    if (
        restraint.mask1
        and restraint.mask2
        and not restraint.mask3
        and not restraint.mask4
    ):
        data = openff_unit.Quantity(
            pt.distance(traj, " ".join([restraint.mask1, restraint.mask2]), image=True),
            units=openff_unit.angstrom,
        ).to(distance_unit)

    elif (
        restraint.mask1 and restraint.mask2 and restraint.mask3 and not restraint.mask4
    ):
        data = openff_unit.Quantity(
            pt.angle(
                traj, " ".join([restraint.mask1, restraint.mask2, restraint.mask3])
            ),
            units=openff_unit.degrees,
        ).to(angle_unit)

    elif restraint.mask1 and restraint.mask2 and restraint.mask3 and restraint.mask4:
        data = openff_unit.Quantity(
            pt.dihedral(
                traj,
                " ".join(
                    [restraint.mask1, restraint.mask2, restraint.mask3, restraint.mask4]
                ),
            ),
            units=openff_unit.degrees,
        ).to(angle_unit)

    return data


def summarize_statistics(x, y):
    """
    Summarize statistics for linear regression of "y vs x".

    Parameters
    ----------
    x: np.array
        X data
    y: np.array
        Y data

    Returns
    -------
    summary_statistics: List
        * slope
        * intercept
        * R - Pearson's correlation coefficient
        * R^2 - coefficient of determination
        * RMSE - Root-Mean-Squared-Error
        * MSE - Mean Signed Error
        * MUE - Mean Unsigned Error or Mean Absolute Error
        * Tau - Kendall's Tau
    """
    summary_statistics = np.empty(8)
    # Slope, intercept, R - Pearson correlation coefficient
    (
        summary_statistics[0],
        summary_statistics[1],
        summary_statistics[2],
        pval,
        stderr,
    ) = statistics.linregress(x, y)

    # R^2 - coefficient of determination
    summary_statistics[3] = summary_statistics[2] ** 2

    # RMSE - Root-Mean-Squared-Error
    summary_statistics[4] = np.sqrt(np.mean((y - x) ** 2))

    # MSE - Mean Signed Error
    summary_statistics[5] = np.mean(y - x)

    # MUE - Mean Unsigned Error
    summary_statistics[6] = np.mean(np.absolute(y - x))

    # Tau - Kendall's Tau
    summary_statistics[7], prob = statistics.kendalltau(x, y)

    return summary_statistics
