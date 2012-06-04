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
        datin      array; mad acts on axis=0


        Optional Input
        --------------
        z          Input is allowed to deviate maximum z standard deviations from the median (default: 7)
        deriv      0: Act on raw input; 1: Use first derivatives; 2: Use 2nd derivatives


        Output
        ------
        mask with false everywhere except where input deviates more than z standard deviations from median


        Restrictions
        ------------
        If input is an array then it mad checks along the zeroth axis for outlier.

        1st derivative is
            d = datin[1:n]-datin[0:n-1]
        because mean of left and right would give 0 for spikes.


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

        >>> yy = np.transpose(np.array([y,y]))
        >>> print np.transpose(mad(yy,z=4))
        [[False False False False False False False False False False False False
          False False False False False False False False False False False False
          False  True]
         [False False False False False False False False False False False False
          False False False False False False False False False False False False
          False  True]]

        >>> yyy = np.transpose(np.array([y,y,y]))
        >>> print np.transpose(mad(yyy,z=3))
        [[ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]
         [ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]
         [ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]]


        History
        -------
        Written,  MC, Nov 2011
        Modified, MC, May 2012 - act on axis=0 of array
                  MC, Jun 2012 - axis=0 did not always work: spread md and MAD to input dimensions
    """
    sn = list(np.shape(datin))
    n  = sn[0]
    if deriv == 0:
        m      = n
        d      = datin
    elif deriv == 1:
        m      = n-1
        sm     = sn
        sm[0]  = m
        dd     = np.ravel(datin)
        d      = np.reshape(dd[1:]-dd[:-1], sm)
    elif deriv == 2:
        m      = n-2
        sm     = sn
        sm[0]  = m
        dd     = np.ravel(datin)
        d      = np.reshape((dd[1:-1]-dd[:-2]) - (dd[2:]-dd[1:-1]), sm)
    else:
        print 'Unimplemented option in mad'
        import sys
        sys.exit(1)
    # Median
    md     = np.ma.median(d, axis=0)
    md     = np.ma.array([md for i in xrange(m)])
    # Median absolute deviation
    MAD    = np.ma.median(np.ma.abs(d-md), axis=0)
    MAD    = np.ma.array([MAD for i in xrange(m)])
    # Range around median
    thresh = MAD * (z/0.6745)
    # True where outside z-range
    res = (d<(md-thresh)) | (d>(md+thresh))

    return res


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,\
    #               2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,\
    #               4.64,5.34,5.42,8.01],dtype=np.float)
    # yy = np.transpose(np.array([y,y]))
    # print np.transpose(mad(yy,z=4))
