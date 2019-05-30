#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.division import division
import jams.const as const

def interpol(xout, xin, yin):
    """
        One-dimensional linear interpolation on first dimension.

        If xin and yin are 1D arrays, this function wraps to numpy.interp.
        If yin is an ND array and xin and xout are 1D, then yin is interpolated in
        its first dimension.


        Definition
        ----------
        def interpol(xout, xin, yin):


        Input
        -----
        xout    1D array of the x-coordinates of the interpolated values if yin is ND array.
                array-like of x-coordinates of the interpolated values if yin is 1D array.
        xin     1D array of the x-coordinates of the data points, must be increasing.
        yin     ND array of the y-coordinates of the data points, 1st dim must have same length as xin.


        Output
        ------
        The interpolated values with shape (np.size(xin),)+yin.shape[1:].


        Examples
        --------
        >>> import numpy as np
        >>> xin  = np.arange(360, dtype=np.float)
        >>> yin  = np.sin(xin)
        >>> xout = np.arange(10)*10. + 0.5
        >>> soll = np.interp(xout, xin, yin)
        >>> yout = interpol(xout, xin, yin)
        >>> print(np.any(yout != soll))
        False

        >>> sout = (3,1)
        >>> yin2 = np.transpose(np.tile(yin,sout))
        >>> yout = interpol(xout, xin, yin2)
        >>> for i in range(3):
        ...    if np.any(yout[:,i] != soll):
        ...        print(True)

        >>> sout = (3,2,1)
        >>> yin3 = np.transpose(np.tile(yin,sout))
        >>> yout = interpol(xout, xin, yin3)
        >>> for i in range(3):
        ...    for j in range(2):
        ...        if np.any(yout[:,j,i] != soll):
        ...            print(True)


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jun 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
                  MC, Nov 2016 - const.tiny -> const.eps
    """
    #
    # If yin 1D-array then call immediately np.interp without check
    # If no np.interp wanted uncomment next line
    # isone = False
    if np.ndim(yin) == 1:
        # If no np.interp wanted comment next line and uncomment the following two lines
        return np.interp(xout, xin, yin)
        # isone = True
        # yin = yin[:,np.newaxis]
    #
    # Check input
    assert np.ndim(xin)  == 1, "x input values not 1D array"
    assert np.ndim(xout) <= 1, "x output values not scalar or 1D array"
    #
    # Subscripts
    eps  = const.eps
    s    = np.minimum(np.maximum(np.searchsorted(xin, xout)-1, 0), xin.size-2) # Subscript intervals
    # Distances
    ums1 = xout-xin[s]                                                         # distance from point before
    ums2 = xin[s+1]-xin[s]
    ums  = division(ums1, ums2, 0.)
    ums  = np.where((np.abs(ums1) < eps) | (np.abs(ums2) < eps), 0., ums)      # for numerical stability
    # Blow to output shape
    sout = yin.shape[1:][::-1] + (1,)
    ums  = np.transpose(np.tile(ums,sout))

    # If no np.interp wanted comment next line and uncomment the following five lines
    return yin[s,...] + ums*(yin[s+1,...]-yin[s,...])
    # yout = yin[s,...] + ums*(yin[s+1,...]-yin[s,...])
    # if isone:
    #     return np.squeeze(yout)
    # else:
    #     return yout


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    # xin  = np.arange(360, dtype=np.float)
    # yin  = np.sin(xin)
    # xout = np.arange(10)*10. + 0.5
    # print np.interp(xout, xin, yin)
    # print interpol(xout, xin, yin)

    # sout = (10,1)
    # yin2 = np.transpose(np.tile(yin,sout))
    # yout = interpol(xout, xin, yin2)
    # print yout[:,0]
    # print yout[:,5]

    # sout = (3,2,1)
    # yin3 = np.transpose(np.tile(yin,sout))
    # yout = interpol(xout, xin, yin3)
    # print yout[:,0,0]
    # print yout[:,1,2]

