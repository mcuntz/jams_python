#!/usr/bin/env python
import numpy as np
from ufz import around, tiny

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
        symmetric  if True, range will be symmetric around 0. if min(arr) <0 and max(arr)>0.


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
        >>> print yrange(np.arange(102))
        [0.0, 101.0]
        
        >>> print yrange(np.arange(102)-10.)
        [-10.0, 91.0]
        
        >>> print yrange(np.arange(102)-10., symmetric=True)
        [-91.0, 91.0]
        

        History
        -------
        Written, MC, Jan 2012
    """
    #
    eps = tiny
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
    maxdiff = np.amax(np.abs((sarr-np.roll(sarr,1))[1:]))
    expom   = np.log10(maxdiff)
    if expom > 0:
        expom = np.floor(expom+10.*eps*10.)
    else:
        expom = np.floor(expom-10.*eps)
    mini = around(minarr, expom, floor=True)
    maxi = around(maxarr, expom, ceil=True)
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
    doctest.testmod()
