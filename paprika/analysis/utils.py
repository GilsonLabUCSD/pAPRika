import numpy


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
    sqrt_n = int(round(numpy.sqrt(n) + 0.5))
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

    most_factors = 0
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
    data_array: :class:`numpy.array`
        Array containing data values.

    Returns
    -------
    numpy.max(sems): float
        The maximum SEM obtained from te blocking curve.

    """
    # Get the integer factors for the number of data points. These
    # are equivalent to the block sizes we will check.
    block_sizes = get_factors(len(data_array))

    # An array to store means for each block ... make it bigger than we need.
    block_means = numpy.zeros([block_sizes[-1]], numpy.float64)

    # Store the SEM for each block size, except the last two size for which
    # there will only be two or one blocks total and thus very noisy.
    sems = numpy.zeros([len(block_sizes) - 2], numpy.float64)

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
            block_means[blk_idx] = numpy.mean(data_array[data_beg_idx:data_end_idx])
        # Compute the standard deviation across all blocks, devide by
        # num_blocks-1 for SEM
        sems[size_idx] = numpy.std(block_means[0:num_blocks], ddof=0) / numpy.sqrt(
            num_blocks - 1
        )
        # Hmm or should ddof=1? I think 0, see Flyvbjerg -----^

    return numpy.max(sems)


def get_subsampled_indices(N, g, conservative=False):
    """Get the indices of independent (subsampled) frames. This is adapted from the implementation in `pymbar`.

    Parameters
    ----------
    N: int
        The length of the array to be indexed.
    g: float
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
        g = numpy.ceil(g)

    # initialize
    indices = [0]
    g_idx = 1.0
    int_step = int(numpy.round(g_idx * g))

    while int_step < N:
        indices.append(int_step)
        g_idx += 1.0
        int_step = int(numpy.round(g_idx * g))

    return indices
