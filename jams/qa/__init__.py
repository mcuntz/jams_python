#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Defines common error measures and discharge signatures


    Error measures
    --------------
    def bias(obs,mod):      bias
                            mean(obs) - mean(mod)
    def mae(obs,mod):       mean absolute error
                            mean(abs(obs - mod))
    def mse(obs,mod):       mean square error
                            mean((obs - mod)**2)
    def rmse(obs,mod):      root mean square error
                            sqrt(mean((obs - mod)**2))
    def nse(obs,mod):       Nash-Sutcliffe efficiency
                            1 - sum((obs - mod)**2) / sum((obs - mean(obs))**2)
    def kge(obs,mod):       Kling-Gupta efficiency
                            1 - sqrt((1-r)**2 + (1-a)**2 + (1-b)**2),
                            where r is the Pearson correlation of obs and mod,
                                  a is mean(mod) / mean(obs), and
                                  b is std(mod) / std(obs)
    def pearson(obs,mod):   Pearson's correlation coefficient
                            mean((obs-mean(obs))/stddev(obs) * (mod-mean(mod))/stddev(mod))

    Input
    -----
    obs         ND-array
    mod         ND-array
    quantiles   Scalar or 1D array_like percentages of exceedance


    Output
    ------
    Measure calculated along the first axis.



    Discharge signatures
    --------------------
    def autocorrelation(dat, lags):
                            Autocorrelation of a data series at given lags.
    def flowdurationcurve(dat, quantiles=None, concavity_index=False,
                          mid_segment_slope=False, mhigh_segment_volume=False,
                          high_segment_volume=False, low_segment_volume=False):
                            Flow duration curves for a given data vector. The Flow duration curve at a
                            certain quantile x is the data point p where x% of the data points are above the value p.
                            Optionally, can be calculated:
                                the concavity index CI can be calculated [Zhang2014].
                                the FDC mid-segment slope [Shafii et. al 2014]
                                the FDC medium high-segment volume [Shafii et. al 2014]
                                the FDC high-segment volume [Shafii et. al 2014]
                                the FDC low-segment volume [Shafii et. al 2014]
    def limbdensities(dat):
                            Rising and declinging limb densities,
                            which are the duration of the data increase (decrease) divided by the number of peaks.
    def maximummonthlyflow(date, dat):
                            Maximum of average flows per month
    def moments(dat, mean_data=False, stddev_data=False, median_data=False,
                            Moments of data and log-transformed data
                            Returns several moments of data series given, i.e.
                                * mean               of data
                                * standard deviation of data
                                * median             of data
                                * maximum/ peak      of data
                                * mean               of log-transformed data
                                * standard deviation of log-transformed data
                                * median             of log-transformed data
                                * maximum/ peak      of log-transformed data
     def peakdistribution(dat, quantiles=None, slope_peak_distribution=False):
                            Calculates the peak distribution.
                            Optionally, the slope of the peak distribution can be returned.


    Input
    -----
    dat         1D array_like
    quantiles   Scalar or 1D array_like percentages of exceedance


    Output
    ------
    Discharge signature.


    Examples
    --------
    >>> # Create some data
    >>> obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])

    >>> # bias
    >>> print(np.round(bias(obs, mod),2))
    -0.43

    >>> print(np.round(mae(obs, mod),2))
    0.55

    >>> print(np.round(mse(obs, mod),2))
    0.54

    >>> print(np.round(rmse(obs, mod),2))
    0.73

    >>> print(np.round(nse(obs, mod),2))
    0.71

    >>> print(np.round(kge(obs, mod),2))
    0.58

    >>> print(np.round(pear2(obs, mod),2))
    0.99

    >>> print(np.round(confint(obs, p=0.95)))
    [ 13.  16.]


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2016 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.


    History
    -------
    Written, MC, May 2016
"""

from .quality_assess import bias, mae, mse, rmse, nse, kge, pearson
from .signatures     import autocorrelation, flowdurationcurve, limbdensities
from .signatures     import maximummonthlyflow, moments, peakdistribution

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 2419"
__date__     = 'Date: 16.05.2016'
