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

        If all(d.mask==True) then return d.mask, which is all True
        

        Examples
        --------
        >>> import numpy as np
        >>> y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,
        ...               2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,
        ...               4.64,5.34,5.42,8.01],dtype=np.float)
        >>> print(mad(y))
        [False False False False False False False False False False False False
         False False False False False False False False False False False False
         False False]

        >>> print(mad(y,z=4))
        [False False False False False False False False False False False False
         False False False False False False False False False False False False
         False  True]

        >>> print(mad(y,z=3))
        [ True False False False False False False False False False False False
         False False False False False False False False False False False False
          True  True]

        >>> print(mad(y,z=4,deriv=2))
        [False False False False False False False False False False False False
         False False False False False False False False False False False  True]

        >>> my = np.ma.array(y, mask=mad(y,z=4))
        >>> print(my)
        [-0.25 0.68 0.94 1.15 2.26 2.35 2.37 2.4 2.47 2.54 2.62 2.64 2.9 2.92 2.92
         2.93 3.21 3.26 3.3 3.59 3.68 4.3 4.64 5.34 5.42 --]

        >>> yy = np.transpose(np.array([y,y]))
        >>> print(np.transpose(mad(yy,z=4)))
        [[False False False False False False False False False False False False
          False False False False False False False False False False False False
          False  True]
         [False False False False False False False False False False False False
          False False False False False False False False False False False False
          False  True]]

        >>> yyy = np.transpose(np.array([y,y,y]))
        >>> print(np.transpose(mad(yyy,z=3)))
        [[ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]
         [ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]
         [ True False False False False False False False False False False False
          False False False False False False False False False False False False
           True  True]]


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

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Nov 2011
        Modified, MC, May 2012 - act on axis=0 of array
                  MC, Jun 2012 - axis=0 did not always work: spread md and MAD to input dimensions
                  MC, Jun 2012 - use np.diff, remove spreads
                  MC, Feb 2013 - ported to Python 3
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
        d      = np.diff(datin, axis=0)
    elif deriv == 2:
        m      = n-2
        sm     = sn
        sm[0]  = m
        d      = np.diff(datin, n=2, axis=0)
    else:
        raise ValueError('Unimplemented option.')
    # Shortcut if all masked
    if type(d) == type(np.ma.empty(1)):
        if np.all(d.mask == True):
            return d.mask
    # Median
    md     = np.ma.median(d, axis=0)
    # Median absolute deviation
    MAD    = np.ma.median(np.ma.abs(d-md), axis=0)
    # Range around median
    thresh = MAD * (z/0.6745)
    # True where outside z-range
    res = (d<(md-thresh)) | (d>(md+thresh))

    return res


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,
    #               2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,
    #               4.64,5.34,5.42,8.01],dtype=np.float)
    # yy = np.transpose(np.array([y,y]))
    # print np.transpose(mad(yy,z=4))

