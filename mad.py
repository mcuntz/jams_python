#!/usr/bin/env python
import numpy as np

def mad(datin, z=7, deriv=0):
    """
        Median absolute deviation test, either on raw values, 1st or 2nd derivatives.
        Returns mask with false everywhere except where <(median-MAD*z/0.6745) or >(md+MAD*z/0.6745).

        Definition
        ----------
        def mad(datin, z=7, deriv=0):


        Input
        -----
        datin      1d array


        Optional Input
        --------------
        z          Input is allowed to deviate maximum z standard deviations from the median (default: 7)
        deriv      0: Act on raw input; 1: Use first derivatives; 2: Use 2nd derivatives


        Output
        ------
        mask with false everywhere except where input deviates more than z standard deviations from median


        Restrictions
        ------------
        1st derivative is
            d = datin[1:n]-datin[0:n-1]
        because mean of left and right would give 0 for spikes


        Examples
        --------
        >>> import numpy as np
        >>> y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,\
                          2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,\
                          4.64,5.34,5.42,8.01],dtype=np.float)
        >>> print mad(y)
        [False False False False False False False False False False False False
         False False False False False False False False False False False False
         False False]

        >>> print mad(y,z=4)
        [False False False False False False False False False False False False
         False False False False False False False False False False False False
         False  True]

        >>> print mad(y,z=3)
        [ True False False False False False False False False False False False
         False False False False False False False False False False False False
          True  True]

        >>> print mad(y,z=4,deriv=2)
        [False False False False False False False False False False False False
         False False False False False False False False False False False  True]

        >>> my = np.ma.array(y, mask=mad(y,z=4))
        >>> print my
        [-0.25 0.68 0.94 1.15 2.26 2.35 2.37 2.4 2.47 2.54 2.62 2.64 2.9 2.92 2.92
         2.93 3.21 3.26 3.3 3.59 3.68 4.3 4.64 5.34 5.42 --]


        History
        -------
        Written, MC, Nov 2011
    """
    n      = np.size(datin)
    if deriv == 0:
        m      = n
        d      = datin
    elif deriv == 1:
        m      = n-1
        d      = datin[1:n]-datin[0:n-1]
    elif deriv == 2:
        m      = n-2
        d      = (datin[1:n-1]-datin[0:n-2]) - (datin[2:n]-datin[1:n-1])
    else:
        print 'Unimplemented option in mad'
        import sys
        sys.exit(1)
    # Median
    md     = np.ma.median(d)
    # Median absolute deviation
    MAD    = np.ma.median(np.ma.abs(md-d))
    # Range around median
    thresh = MAD/0.6745 * z
    # Output is false everywhere
    res    = np.zeros(m, dtype=np.bool)
    # Except where outside z-range
    res[:] = np.logical_or(d<(md-thresh), d>(md+thresh))
    
    return res


if __name__ == '__main__':
    import doctest
    doctest.testmod()
