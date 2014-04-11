#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from around import around # from ufz
import const              # from ufz

def yrange(*args, **kwargs):
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

        >>> print(yrange(np.arange(102), np.arange(1002), np.arange(10002)))
        [0.0, 10001.0]

        >>> print(yrange(-np.arange(102), np.arange(1002), np.arange(10002)))
        [-101.0, 10001.0]

        >>> print(yrange(-np.arange(102), np.arange(1002), np.arange(10002), symmetric=True))
        [-10001.0, 10001.0]

        >>> print(yrange(-np.arange(102), np.arange(1002), np.arange(10002), symmetric=False))
        [-101.0, 10001.0]

        >>> a = np.ma.arange(102)
        >>> a[-1] = np.ma.masked
        >>> print(yrange(a))
        [0.0, 100.0]


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
                  MC, Nov 2013 - accept more than 1 array as input
                  MC, Nov 2013 - masked arrays
                  MC, Apr 2014 - assert
    """
    #
    #eps = tiny
    eps = const.tiny
    # Check input
    assert len(args) > 0, 'no input argument given.'
    counter = 0
    for i in args:
        arr    = np.ma.array(i)
        minarr = np.ma.amin(arr)
        maxarr = np.ma.amax(arr)
        # Round to max difference between adjacent values
        if arr.size == 1:
            expom = 0
        else:
            sarr    = np.ma.sort(arr)
            maxdiff = np.ma.amax(np.diff(sarr))
            expom   = np.ma.log10(maxdiff)
            if expom > 0:
                expom = np.int(np.floor(expom+10.*eps*10.))
            else:
                expom = np.int(np.floor(expom-10.*eps))
        if counter == 0:
            minall   = minarr
            maxall   = maxarr
            expomall = expom
        else:
            minall   = np.minimum(minall, minarr)
            maxall   = np.maximum(maxall, maxarr)
            if expomall > 0:
                expomall = np.maximum(expomall, expom)
            else:
                expomall = np.minimum(expomall, expom)
        counter += 1
        
    # Round range
    mini = around(minall, expomall, floor=True)
    maxi = around(maxall, expomall, ceil=True)

    if 'symmetric' in kwargs:
        if kwargs['symmetric']:
            if (mini*maxi < 0.):
                maxmax =  np.maximum(np.abs(mini),np.abs(maxi))
                maxi   =  maxmax
                mini   = -maxmax

    # Return range
    return [mini,maxi]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    
    # print(yrange(np.arange(102)))
    # #[0.0, 101.0]

    # print(yrange(np.arange(102)-10.))
    # #[-10.0, 91.0]

    # print(yrange(np.arange(102)-10., symmetric=True))
    # #[-91.0, 91.0]

    # print(yrange(np.arange(102), np.arange(1002), np.arange(10002)))
    # #[0.0, 10001.0]

    # print(yrange(-np.arange(102), np.arange(1002), np.arange(10002)))
    # #[-101.0, 10001.0]

    # print(yrange(-np.arange(102), np.arange(1002), np.arange(10002), symmetric=True))
    # #[-10001.0, 10001.0]

    # print(yrange(-np.arange(102), np.arange(1002), np.arange(10002), symmetric=False))
    # #[-101.0, 10001.0]
    
    # a = np.ma.arange(102)
    # a[-1] = np.ma.masked
    # print(yrange(a))
    # #[0.0, 100.0]
