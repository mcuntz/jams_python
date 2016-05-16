#!/usr/bin/env python
"""
    Defines common error measures.

    
    Definition
    ----------
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
    obs        ND-array
    mod        ND-array


    Output
    ------
    Measure calculated along the first axis.
                 

    Restrictions
    ------------
    None
    
    
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
    This file is part of the UFZ Python package.

    The UFZ Python package is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The UFZ Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2014 Arndt Piayda


    History
    -------
    Written, AP, Jul 2014
    Modified MC, Dec 2014 - use simple formulas that work with normal and masked arrays but do not deal with NaN
    Modified AP, Sep 2015 - add confidence interval
    Modified ST, Nov 2015 - added KGE
    Modified MC, May 2016 - calc along 1st axis so that measures work for ND-arrays
                          - all 1-liner
"""

from .quality_assess import bias, mae, mse, rmse, nse, kge, pearson

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 2419"
__date__     = 'Date: 16.05.2016'
