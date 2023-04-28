from typing import Dict, List, Union

import numpy as np
import openmm.unit as openmm_unit
from openff.units import unit
from openff.units.openmm import from_openmm
from scipy import stats as statistics
from scipy.interpolate import Akima1DInterpolator

from paprika.utils import check_unit

R_gas = from_openmm(
    openmm_unit.AVOGADRO_CONSTANT_NA * openmm_unit.BOLTZMANN_CONSTANT_kB
)


def integrate_bootstraps(x, ys, x_intp=None, matrix="full"):
    """
    Integrate splines created via bootstrapping.

    Parameters
    ----------
    x: :class:`np.array`
        The x coordinate of the curve to be integrated.
    ys: :class:`np.array`
        Two dimensional array in which the first dimension is boot_cycles and the second
        dimension contains the arrays of y values which correspond to the x values and will
        be used for integration. The shape of this is :code:`(boot_cycles, len(x))`.
    x_intp: :class:`np.array`, optional, default=None
        An array which finely interpolates the x values. If not provided, it will be generated
        by adding 100 evenly spaced points between each x value. Default: None.
    matrix: str, optional, default='full`
        If ``full``, the mean and SEM integration is computed between x values. If ``diagonal``,
        the mean and SEM integration is computed between the first value and all other values,
        as well as the neighboring values to each value. If ``endpoints``, the integration
        is computed between only the first and last x value.

    Returns
    -------
    avg_matrix: :class:`np.array`
        Matrix of the integration mean between each x value (as specified by 'matrix')
    sem_matrix: :class:`np.array`
        Matrix of the uncertainty (SEM) between each x value (as specified by 'matrix')

    """

    num_x = len(x)

    # Prepare to store the index location of the x values in the x_intp array
    x_idxs = np.zeros([num_x], np.int32)

    # If not provided, generate x interpolation with 100 inpolated points between
    # each x value. Store the index locations of the x values in the x_intp
    # array.
    if x_intp is None:
        x_intp = np.zeros([0], np.float64)
        for i in range(1, num_x):
            x_intp = np.append(
                x_intp, np.linspace(x[i - 1], x[i], num=100, endpoint=False)
            )
            x_idxs = len(x_intp)
        # Tack on the final value onto the interpolation
        x_intp = np.append(x_intp, x[-1])
    # If x_intp is provided, find the locations of x values in x_intp
    else:
        i = 0
        for j in range(len(x_intp)):
            if np.isclose(x[i], x_intp[j]):
                x_idxs[i] = j
                i += 1
        if i != num_x:
            raise Exception(
                "One or more x values seem to be missing in the x_intp array,"
                + " or one of the lists is not monotonically increasing!"
            )

    cycles = len(ys)

    # Setup array to store integration bootstraps
    int_matrix = np.zeros([num_x, num_x, cycles], np.float64)

    # Do the integration bootstraps. Originally, I had matrix=endpoints in the loop
    # below with everything else, but I'll split it out here in case that's faster
    # due to avoiding the if statements.
    if matrix == "endpoints":
        for cycle in range(cycles):
            intp_func = Akima1DInterpolator(x, ys[cycle])
            y_intp = intp_func(x_intp)
            #            for i in range(0, num_x):
            #                for j in range(i+1, num_x):
            #                    int_matrix[i, j, cycle] = np.trapz( y_intp, x_intp )
            int_matrix[0, num_x - 1, cycle] = np.trapz(y_intp, x_intp)
    else:
        for cycle in range(cycles):
            intp_func = Akima1DInterpolator(x, ys[cycle])
            y_intp = intp_func(x_intp)
            for i in range(0, num_x):
                for j in range(i + 1, num_x):
                    if matrix == "diagonal" and i != 0 and j - i > 1:
                        continue
                    beg = x_idxs[i]
                    end = x_idxs[j]
                    int_matrix[i, j, cycle] = np.trapz(y_intp[beg:end], x_intp[beg:end])

    # Setup matrices to store the average/sem values.
    # Is it bad that the default is 0.0 rather than None?
    avg_matrix = np.zeros([num_x, num_x], np.float64)
    sem_matrix = np.zeros([num_x, num_x], np.float64)

    # Second pass to compute the mean and standard deviation.
    for i in range(0, num_x):
        for j in range(i + 1, num_x):
            # If quick_ti_matrix, only populate first row and neighbors in
            # matrix
            if matrix == "diagonal" and i != 0 and j - i > 1:
                continue
            if matrix == "endpoints" and i != 0 and j != num_x - 1:
                continue
            avg_matrix[i, j] = np.mean(int_matrix[i, j])
            avg_matrix[j, i] = -1.0 * avg_matrix[i, j]
            sem_matrix[i, j] = np.std(int_matrix[i, j])
            sem_matrix[j, i] = sem_matrix[i, j]

    return avg_matrix, sem_matrix


def regression_bootstrap(
    x_data,
    x_sem,
    y_data,
    y_sem,
    cycles=1000,
    with_replacement=True,
    with_uncertainty=True,
):
    """
    Estimate the statistics for linear regression given x and y data. This is used to
    compare calculated and experimental values.

    Parameters
    ----------
    x_data
    x_sem
    y_data
    y_sem
    cycles: int
        Number of bootstrap cycles
    with_replacement: bool
        If True, choose bootstrap samples randomly
    with_uncertainty: bool
        If true, generate samples from normal distribution based on mean=data, mu=sem

    Returns
    -------
    Returns the mean, sem, ci_low, ci_high for [slope, intercept, R, R^2, RMSE, MSE, MUE, Kendall's Tau]
    """

    summary_statistics = np.empty((cycles, 8))

    # Bootstrapping
    for cycle in range(cycles):
        new_x = np.empty_like(x_data)
        new_y = np.empty_like(y_data)
        for index in range(len(x_data)):
            if with_replacement:
                j = np.random.randint(len(x_data))
            else:
                j = index
            if with_uncertainty and x_sem is not None:
                new_x[index] = np.random.normal(x_data[j], x_sem[j])
            elif with_uncertainty and x_sem is None:
                new_x[index] = x_data[j]
            if with_uncertainty and y_sem is not None:
                new_y[index] = np.random.normal(y_data[j], y_sem[j])
            elif with_uncertainty and y_sem is None:
                new_y[index] = y_data[j]

        summary_statistics[cycle] = summarize_statistics(new_x, new_y)

    # Confidence interval
    ci = np.empty((8, 2))
    for statistic in range(8):
        sorted_statistic = np.sort(summary_statistics[:, statistic])
        ci[statistic][0] = sorted_statistic[int(0.025 * cycles)]
        ci[statistic][1] = sorted_statistic[int(0.975 * cycles)]

    # Summarize results
    results = {
        "mean": {
            "slope": np.mean(summary_statistics[:, 0]),
            "intercept": np.mean(summary_statistics[:, 1]),
            "R": np.mean(summary_statistics[:, 2]),
            "R**2": np.mean(summary_statistics[:, 3]),
            "RMSE": np.mean(summary_statistics[:, 4]),
            "MSE": np.mean(summary_statistics[:, 5]),
            "MUE": np.mean(summary_statistics[:, 6]),
            "Tau": np.mean(summary_statistics[:, 7]),
        },
        "sem": {
            "slope": np.std(summary_statistics[:, 0]),
            "intercept": np.std(summary_statistics[:, 1]),
            "R": np.std(summary_statistics[:, 2]),
            "R**2": np.std(summary_statistics[:, 3]),
            "RMSE": np.std(summary_statistics[:, 4]),
            "MSE": np.std(summary_statistics[:, 5]),
            "MUE": np.std(summary_statistics[:, 6]),
            "Tau": np.std(summary_statistics[:, 7]),
        },
        "ci_low": {
            "slope": ci[0][0],
            "intercept": ci[1][0],
            "R": ci[2][0],
            "R**2": ci[3][0],
            "RMSE": ci[4][0],
            "MSE": ci[5][0],
            "MUE": ci[6][0],
            "Tau": ci[7][0],
        },
        "ci_high": {
            "slope": ci[0][1],
            "intercept": ci[1][1],
            "R": ci[2][1],
            "R**2": ci[3][1],
            "RMSE": ci[4][1],
            "MSE": ci[5][1],
            "MUE": ci[6][1],
            "Tau": ci[7][1],
        },
    }

    return results


def dG_bootstrap(
    x_data: Union[List, np.array],
    x_sem: Union[List, np.array],
    y_data: Union[List, np.array],
    y_sem: Union[List, np.array],
    cycles: int = 1000,
    temperature: Union[float, unit.Quantity] = 298.15 * unit.kelvin,
    with_uncertainty: bool = True,
    energy_units: unit.Quantity = None,
):
    """
    Combine dG when for multiple binding poses. This is from Eq(A12) and Eq(A13) from
    Henriksen, Fenley and Gilson - 10.1021/acs.jctc.5b00405

    Parameters
    ----------
    x_data: numpy.array or list
        Array of delta G
    x_sem: numpy.array or list
        Array of error for x data
    y_data: numpy.array or list
        Array of delta G
    y_sem: numpy.array or list
        Array of error for y data
    cycles: int
        Number of bootstrap cycles
    temperature: float or openff.units.unit.Quantity
        Temperature of the simulation
    with_uncertainty: bool
        If true, generate samples from normal distribution based on mean=data, mu=sem
    energy_units: unit.Quantity
        If true, return values as openff.units.unit.Quantity. Default is kcal/mol.

    Return
    ------
    results: Dict
        Returns the mean, standard deviation and confidence interval from bootstrapping.
    """

    x_data = check_unit(x_data, base_unit=unit.kilocalorie_per_mole).magnitude
    x_sem = check_unit(x_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    y_data = check_unit(y_data, base_unit=unit.kilocalorie_per_mole).magnitude
    y_sem = check_unit(y_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    temperature = check_unit(temperature, base_unit=unit.kelvin)

    summary_statistics = np.empty((cycles))
    RT = (R_gas * temperature).to(unit.kilocalorie_per_mole).magnitude
    beta = 1.0 / RT

    ci = np.empty((2))
    for cycle in range(cycles):
        new_x = np.empty_like(x_data)
        new_y = np.empty_like(y_data)

        if with_uncertainty and x_sem is not None:
            new_x = np.random.normal(x_data, x_sem)
        elif with_uncertainty and x_sem is None:
            new_x = x_data
        if with_uncertainty and y_sem is not None:
            new_y = np.random.normal(y_data, y_sem)
        elif with_uncertainty and y_sem is None:
            new_y = y_data
        summary_statistics[cycle] = -RT * np.log(
            np.exp(-beta * new_x) + np.exp(-beta * new_y)
        )

        # Get confidence interval
        sorted_statistic = np.sort(summary_statistics)
        ci[0] = sorted_statistic[int(0.025 * cycles)]
        ci[1] = sorted_statistic[int(0.975 * cycles)]

    results = {
        "mean": np.mean(summary_statistics),
        "sem": np.std(summary_statistics),
        "ci": ci,
    }

    if energy_units is not None:
        for stat in results:
            results[stat] *= energy_units

    return results


def dH_bootstrap(
    dH_x_data,
    dH_x_sem,
    dH_y_data,
    dH_y_sem,
    dG_x_data,
    dG_x_sem,
    dG_y_data,
    dG_y_sem,
    cycles=1000,
    temperature=298.15 * unit.kelvin,
    with_uncertainty=True,
    energy_units: unit.Quantity = None,
):
    """
    Combine dH when for multiple binding poses. This is from Eq(A16) and Eq(A17) from
    Henriksen, Fenley and Gilson - 10.1021/acs.jctc.5b00405

    Parameters
    ----------
    dH_x_data: numpy.array or list
        Array of delta H
    dH_x_sem: numpy.array or list
        Array of error for x data
    dH_y_data: numpy.array or list
        Array of delta H
    dH_y_sem: numpy.array or list
        Array of error for y data
    dG_x_data: numpy.array or list
        Array of delta G
    dG_x_sem: numpy.array or list
        Array of error for x data
    dG_y_data: numpy.array or list
        Array of delta G
    dG_y_sem: numpy.array or list
        Array of error for y data
    cycles: int
        Number of bootstrap cycles
    temperature: float or openff.units.unit.Quantity
        Temperature of the simulation
    with_uncertainty: bool
        If true, generate samples from normal distribution based on mean=data, mu=sem
    energy_units: unit.Quantity
        If true, return values as openff.units.unit.Quantity. Default is kcal/mol.

    Return
    ------
    results: Dict
        Returns the mean, standard deviation and confidence interval from bootstrapping.
    """

    dH_x_data = check_unit(dH_x_data, base_unit=unit.kilocalorie_per_mole).magnitude
    dH_x_sem = check_unit(dH_x_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    dH_y_data = check_unit(dH_y_data, base_unit=unit.kilocalorie_per_mole).magnitude
    dH_y_sem = check_unit(dH_y_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    dG_x_data = check_unit(dG_x_data, base_unit=unit.kilocalorie_per_mole).magnitude
    dG_x_sem = check_unit(dG_x_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    dG_y_data = check_unit(dG_y_data, base_unit=unit.kilocalorie_per_mole).magnitude
    dG_y_sem = check_unit(dG_y_sem, base_unit=unit.kilocalorie_per_mole).magnitude
    temperature = check_unit(temperature, base_unit=unit.kelvin)

    summary_statistics = np.empty((cycles))
    RT = (R_gas * temperature).to(unit.kcal / unit.mole).magnitude
    beta = 1.0 / RT

    ci = np.empty((2))
    for cycle in range(cycles):
        new_dH_x = np.empty_like(dH_x_data)
        new_dH_y = np.empty_like(dH_y_data)

        new_dG_x = np.empty_like(dG_x_data)
        new_dG_y = np.empty_like(dG_y_data)

        # Resample dH
        if with_uncertainty and dH_x_sem is not None:
            new_dH_x = np.random.normal(dH_x_data, dH_x_sem)
        elif with_uncertainty and dH_x_sem is None:
            new_dH_x = dH_x_data
        if with_uncertainty and dH_y_sem is not None:
            new_dH_y = np.random.normal(dH_y_data, dH_y_sem)
        elif with_uncertainty and dH_y_sem is None:
            new_dH_y = dH_y_data

        # Resample dG
        if with_uncertainty and dG_x_sem is not None:
            new_dG_x = np.random.normal(dG_x_data, dG_x_sem)
        elif with_uncertainty and dG_x_sem is None:
            new_dG_x = dG_x_data
        if with_uncertainty and dG_y_sem is not None:
            new_dG_y = np.random.normal(dG_y_data, dG_y_sem)
        elif with_uncertainty and dG_y_sem is None:
            new_dG_y = dG_y_data

        summary_statistics[cycle] = (
            new_dH_x * np.exp(-beta * new_dG_x) + new_dH_y * np.exp(-beta * new_dG_y)
        ) / (np.exp(-beta * new_dG_x) + np.exp(-beta * new_dG_y))

    # Confidence interval
    sorted_statistic = np.sort(summary_statistics)
    ci[0] = sorted_statistic[int(0.025 * cycles)]
    ci[1] = sorted_statistic[int(0.975 * cycles)]

    results = {
        "mean": np.mean(summary_statistics),
        "sem": np.std(summary_statistics),
        "ci": ci,
    }

    if energy_units is not None:
        for stat in results:
            results[stat] *= energy_units

    return results


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
