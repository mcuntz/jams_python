#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from division import division # from ufz
import const                  # from ufz

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
        Written,  MC, Jun 2012
        Modified, MC, Feb 2013 - ported to Python 3
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
    if np.ndim(xin) != 1: raise ValueError("x input values not 1D array")
    if np.ndim(xout) > 1: raise ValueError("x output values not scalar or 1D array")
    #
    # Subscripts
    #tiny = np.finfo(np.float).eps
    tiny = const.tiny
    s    = np.minimum(np.maximum(np.searchsorted(xin, xout)-1, 0), xin.size-2) # Subscript intervals
    # Distances
    ums1 = xout-xin[s]                                                         # distance from point before
    ums2 = xin[s+1]-xin[s]
    ums  = division(ums1, ums2, 0.)
    ums  = np.where((np.abs(ums1) < tiny) | (np.abs(ums2) < tiny), 0., ums)    # for numerical stability
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

