#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

__all__ = ['constant_values']

# --------------------------------------------------------------------

def constant_values(juldate,series, length, eps):
    """
        Checks if a given series of data contains consecutive values which are constant over a certain time period.

        Definition
        ----------
        def constant_values(juldate, series, length, eps):


        Input
        -----
        juldate 1D-array of julian dates for data points in series
        series  1D-array of data
        length  length of time period in which values are checked for constancy
                (unit of length is days)
        eps     difference of minimal and maximal data point in window of length <length> has to
                be above 2*eps otherwise all datapoints in window are regarded to be constant

        
        Output
        ------
        1D-array of boolean values indicating if value is in a constant period (.true.) or not (.false.)


        Examples
        --------
        --> see __init__.py for full example of workflow

        >>> dates  = np.array([100.1, 100.6, 101.1, 101.6, 102.1, 102.6, 103.1, 103.6, 104.1, 104.6, 105.1, 105.6])
        >>> series = np.array([1.0, 2.0, 5.0, 3.000, 3.001, 3.010, 3.012, 2.990, 3.100, 4.0, 3.0, 5.0])
        >>> length = 2
        >>> eps    = 0.1
        >>> print(constant_values(dates, series, length, eps))
        [False False False  True  True  True  True  True  True False False False]

        >>> dates  = np.array([100.1, 100.6, 101.1, 101.6, 102.1, 102.6, 103.1, 103.6, 104.1, 104.6, 105.1, 105.6])
        >>> series = np.array([1.0, 2.0, 5.0, 3.000, 3.001, 3.010, 3.012, 2.990, 3.100, 4.0, 3.0, 5.0])
        >>> length = 3.5
        >>> eps    = 0.1
        >>> print(constant_values(dates, series, length, eps))
        [False False False False False False False False False False False False]

        >>> dates  = np.array([100.1, 100.6, 101.1, 101.6, 102.1, 102.6, 103.1, 103.6, 104.1, 104.6, 105.1, 105.6])
        >>> series = np.array([1.0, 2.0, 5.0, 3.000, 3.001, 4.0, 3.012, 2.990, 3.100, 4.0, 3.0, 5.0])
        >>> length = 2
        >>> eps    = 0.1
        >>> print(constant_values(dates, series, length, eps))
        [False False False False False False False False False False False False]

        >>> dates  = np.array([100.1, 100.6, 101.1, 101.6, 102.1, 102.6, 103.1, 103.6, 104.1, 104.6, 105.1, 105.6])
        >>> series = np.array([1.0, 2.0, 5.0, 4.0, 3.0, 5.0, 4.0, 3.000, 3.001, 3.012, 2.990, 3.100])
        >>> length = 2
        >>> eps    = 0.1
        >>> print(constant_values(dates, series, length, eps))
        [False False False False False False False  True  True  True  True  True]


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  JM, Oct 2015
    """
    if (len(series.shape) > 1):
        raise ValueError('constant_values: only allowed for 1d data arrays')

    if (len(juldate.shape) > 1):
        raise ValueError('constant_values: only allowed for 1d julian date arrays')

    if (np.size(juldate,0) != np.size(series,0) ):
        raise ValueError('constant_values: size of dates and data points are not matching')

    nseries = np.size(series,0)
    
    out = np.zeros(nseries, dtype=np.bool)

    # indexes of entries where difference is too small
    diff_idx = np.where(np.abs(np.diff(series))<2.0*eps)[0]

    for idx in diff_idx:

        wind_idx = np.where((juldate>=juldate[idx]) & (juldate<juldate[idx]+length))[0]
        
        min_val = np.amin(series[wind_idx])
        max_val = np.amax(series[wind_idx])

        while (max_val - min_val < 2.0*eps):

            out[wind_idx] = True

            # try now also with next data point
            # but only if last index is not the last one possible
            if (wind_idx[-1] != nseries - 1):
                wind_idx = np.append(wind_idx,np.array(wind_idx[-1]+1))
                min_val = np.amin(series[wind_idx])
                max_val = np.amax(series[wind_idx])
            else:
                # make sure that it exits the while here
                max_val = min_val + 3.0*eps

        # shorten the idexes which still have to be checked
        # e.g. [2,3,4,5,6,10,11,12] --> [10,11,12]
        diff_idx = diff_idx[np.where(diff_idx>wind_idx[-1]-1)[0]]
    
    return out

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
