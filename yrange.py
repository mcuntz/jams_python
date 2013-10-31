#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from around import * # from ufz
import const         # from ufz

def yrange(arr, symmetric=False):
    """
        Calculates plot range from input array


        Definition
        ----------
        def yrange(arr, symmetric=False):


        Input
        -----
        arr        number array


        Optional Input
        --------------
        symmetric  if True, range will be symmetric around 0. if min(arr)<0 and max(arr)>0.


        Output
        ------
        Range to be used as [xyz]range


        Restrictions
        ------------
        Uses around.
        Does not work well for 0<range<1. Use yrange(arr*10.)/10.


        Examples
        --------
        >>> import numpy as np
        >>> print(yrange(np.arange(102)))
        [0.0, 101.0]

        >>> print(yrange(np.arange(102)-10.))
        [-10.0, 91.0]

        >>> print(yrange(np.arange(102)-10., symmetric=True))
        [-91.0, 91.0]


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, Feb 2013 - ported to Python 3
    """
    #
    #eps = tiny
    eps = const.tiny
    # Check input
    if len(arr) == 1:
        return [arr[0],arr[0]]
    maxarr = np.amax(arr)
    minarr = np.amin(arr)
    if maxarr == minarr:
        return [maxarr,maxarr]
    #
    # Round to max difference between adjacent values
    sarr    = np.sort(arr)
    #maxdiff = np.amax(np.abs((sarr-np.roll(sarr,1))[1:]))
    maxdiff = np.amax(np.diff(sarr))
    expom   = np.log10(maxdiff)
    if expom > 0:
        expom = np.int(np.floor(expom+10.*eps*10.))
    else:
        expom = np.int(np.floor(expom-10.*eps))
    mini = around(minarr, expom, floor=True)
    #mini = np.around(minarr, expom)
    maxi = around(maxarr, expom, ceil=True)
    #maxi = np.around(maxarr, expom)
    if symmetric:
        if (mini*maxi < 0.):
            maxmax =  np.amax(np.abs([mini,maxi]))
            maxi   =  maxmax
            mini   = -maxmax
    #
    # Return range
    return [mini,maxi]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

