"""
Python module to compute the Mann-Kendall test for trend in time series data.

This module contains a single function 'test' which implements the Mann-Kendall
test for a linear trend in a given time series. 

Introduction to the Mann-Kendall test
-------------------------------------

The Mann-Kendall test is used to determine whether or not there is a linear
monotonic trend in a given time series data. It is a non-parametric trend
closely related to the concept of Kendall's correlation coefficient [1]_. The
null hypothesis, :math:`H_0`, states that there is no monotonic trend, and this
is tested against one of three possible alternative hypotheses, :math:`H_a`:
(i) there is an upward monotonic trend, (ii) there is a downward monotonic
trend, or (iii) there is either an upward monotonic trend or a downward
monotonic trend. It is a robust test for trend detection used widely in
financial, climatological, hydrological, and environmental time series
analysis.

Assumptions underlying the Mann-Kendall test
--------------------------------------------
The Mann-Kendall test involves the following assumptions regarding the given
time series data:
    1. In the absence of a trend, the data are independently and identically
    distributed (iid).
    2. The measurements represent the true states of the observables at the
    times of measurements.
    3. The methods used for sample collection, instrumental measurements and
    data handling are unbiased.

Advantages of the Mann-Kendall test
-----------------------------------
The Mann-Kendall test provides the following advantages:
    1. It does not assume the data to be distributed according to any
    particular rule, i.e., e.g., it does not require that the data be normally
    distributed.
    2. It is not effected by missing data other than the fact the number of
    sample points are reduced and hence might effect the statistical
    significance adversely.
    3. It is not effected by irregular spacing of the time points of
    measurement.
    4. It is not effected by the length of the time series.

Limitations of the Mann-Kendall test
------------------------------------
The following limitations have to be kept in mind:
    1. The Mann-Kendall test is not suited for data with periodicities (i.e.,
    seasonal effects). In order for the test to be effective, it is recommended
    that all known periodic effects be removed from the data in a preprocessing
    step before computing the Mann-Kendall test.
    2. The Mann-Kendall test tends to give more negative results for shorter
    datasets, i.e., the longer the time series the more effective is the trend
    detection computation. 

Formulae
--------
The first step in the Mann-Kendall test for a time series :math:`x_1, x_2,
\dots, x_n` of length :math:`n` is to compute the indicator function
:math:`sgn(x_i - x_j)` such that:

    .. math::
        sgn(x_i - x_j) = \begin{cases}
                            1, if x_i - x_j > 0,
                            0, if x_i - x_j = 0,
                            -1, if x_i - x_j < 0,
                         \end{cases}

"""
# Created: Mon Apr 17, 2017  01:18PM
# Last modified: Mon Apr 17, 2017  05:30PM
# Copyright: Bedartha Goswami <goswami@uni-potsdam.de>


import numpy as np
from scipy.special import ndtri, ndtr
import sys


def test(t, x, precx=None, alpha=None, Ha=None):
    """
    Runs the Mann-Kendall test for trend in time series data.

    Parameters
    ----------
    t : 1D numpy.ndarray
        array of the tmie points of measurements
    x : 1D numpy.ndarray
        array containing the measurements corresponding to entries of 't'
    precx : scalar, float, greater than zero
        least count error of measurements which help determine ties in the data
    alpha : scalar, float, greater than zero
        significance level of the statistical test (Type I error)
    Ha : string, options include 'up', 'down', 'upordown'
        type of test: one-sided ('up' or 'down') or two-sided ('updown')

    Returns
    -------
    MK : string
        result of the statistical test indicating whether or not to accept hte
        alternative hypothesis 'Ha'
    m : scalar, float
        slope of the linear fit to the data
    c : scalar, float
        intercept of the linear fit to the data
    p : scalar, float, greater than zero
        p-value of the obtained Z-score statistic for the Mann-Kendall test

    Raises
    ------
    AssertionError:
        when the least count error of measurements 'precx' is not given
        when significance level of test 'alpha' is not given
        when type of alternative hypothesis 'Ha' is not given

    """
    # assert a least count for the measurements x
    assert precx, "Please provide least count error for measurements 'x'"
    assert alpha, "Please provide significance level 'alpha' for the test"
    assert Ha, "Please provide the alternative hypothesis 'Ha'"

    # estimate sign of all possible (n(n-1)) / 2 differences
    n = len(t)
    sgn = np.zeros((n, n), dtype="int")
    for i in range(n):
        tmp = x - x[i]
        tmp[np.where(np.fabs(tmp) < precx)] = 0.
        sgn[i] = np.sign(tmp)

    # estimate mean of the sign of all possible differences
    S = sgn[np.triu_indices(n, k=1)].sum()

    # estimate variance of the sign of all possible differences
    # 1. Determine no. of tie groups 'p' and no. of ties in each group 'q'
    np.fill_diagonal(sgn, precx * 1E6)
    i, j = np.where(sgn == 0.)
    ties = np.unique(x[i])
    p = len(ties)
    q = np.zeros(len(ties), dtype="int")
    for k in range(p):
        idx =  np.where(np.fabs(x - ties[k]) < precx)[0]
        q[k] = len(idx)
    # 2. Determine the two terms in the variance calculation
    term1 = n * (n - 1) * (2 * n + 5)
    term2 = (q * (q - 1) * (2 * q + 5)).sum()
    # 3. estimate variance
    varS = float(term1 - term2) / 18.

    # Compute the Z-score based on above estimated mean and variance
    if (S - precx) > 0.:
        Zmk = (S - 1) / np.sqrt(varS)
    elif (S - precx) == 0.:
        Zmk = 0.
    elif (S + precx) < 0.:
        Zmk = (S + 1) / np.sqrt(varS)

    # compute test based on given 'alpha' and alternative hypothesis
    # note: for all the following cases, the null hypothesis Ho is:
    # Ho := there is no monotonic trend
    # 
    # Ha := There is an upward monotonic trend
    if Ha == "up":
        Z_ = ndtri(1. - alpha)
        if Zmk >= Z_:
            MK = "accept Ha := upward trend"
        else:
            MK = "reject Ha := upward trend"
    # Ha := There is a downward monotonic trend
    elif Ha == "down":
        Z_ = ndtri(1. - alpha)
        if Zmk <= -Z_:
            MK = "accept Ha := downward trend"
        else:
            MK = "reject Ha := downward trend"
    # Ha := There is an upward OR downward monotonic trend
    elif Ha == "upordown":
        Z_ = ndtri(1. - alpha / 2.)
        if np.fabs(Zmk) >= Z_:
            MK = "accept Ha := upward OR downward trend"
        else:
            MK = "reject Ha := upward OR downward trend"

    # ----------
    # AS A BONUS
    # ----------
    # estimate the slope and intercept of the line
    m = np.corrcoef(t, x)[0, 1] * (np.std(x) / np.std(t))
    c = np.mean(x) - m * np.mean(t)

    # ----------
    # AS A BONUS
    # ----------
    # estimate the p-value fo rthe obtained Z-score Zmk
    if (S - precx) > 0.:
        if Ha == "up":
            p = 1. - ndtr(Zmk)
        elif Ha == "down":
            p = ndtr(Zmk)
        elif Ha == "upordown":
            p = 0.5 * (1. - ndtr(Zmk))
    elif (S - precx) == 0.:
        p = 0.5
    elif (S + precx) < 0.:
        if Ha == "up":
            p = 1. - ndtr(Zmk)
        elif Ha == "down":
            p = ndtr(Zmk)
        elif Ha == "upordown":
            p = 0.5 * (ndtr(Zmk))

    return MK, m, c, p
