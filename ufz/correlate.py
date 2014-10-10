#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def correlate(x, y, axis=None):
    """
        Computes the cross-correlation function of two series x and y.

        Note that the computations are performed on anomalies (deviations from average).

        Returns the values of the cross-correlation at different lags.
        Lags are given as [0,1,2,...,n,n-1,n-2,...,-2,-1].


        Definition
        ----------
        def correlate(x, y, axis=None):


        Input
        -----
        x          1D or 2D MaskedArray, time series
        y          1D or 2D MaskedArray, time series


        Optional Input
        --------------
        axis       integer (default: None) Axis along which to compute (0 for rows, 1 for cols)
                   If None, the array is flattened first.


        Output
        ------
        Returns the values of the cross-correlation at different lags.
        Lags are given as [0,1,2,...,n,n-1,n-2,...,-2,-1],
        so that the output has a kind of u-shape.


        References
        ----------
        This routine is from Pierre GM on
            http://sourceforge.net/p/matplotlib/mailman/matplotlib-users/?viewmonth=200611&viewday=14
        or
            http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg01487.html



        Examples
        --------
        >>> x = np.sin(np.random.random(10000)) # some data
        >>> y = np.roll(x, 100)                 # shift by 100
        >>> print(np.argmax(correlate(x,y)))          # should be 100
        100

        >>> x = np.sin(np.random.random(2000)).reshape(10,200)
        >>> y = np.roll(x, 100, axis=1)
        >>> ii = np.argmax(correlate(x,y,axis=1),axis=1) # should be 100 and 300
        >>> ii = np.where(ii > 200, ii-200, ii)
        >>> print(ii)
        [100 100 100 100 100 100 100 100 100 100]


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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Apr 2014 - from Pierre GM at matplotlib-users@lists.sourceforge.net
    """
    assert x.shape == y.shape, "Inconsistent shape !"
    if axis is None:
        if x.ndim > 1:
            x = x.ravel()
            y = y.ravel()
        npad = x.size + y.size
        xanom = (x - x.mean(axis=None))
        yanom = (y - y.mean(axis=None))
        Fx = np.fft.fft(xanom, npad, )
        Fy = np.fft.fft(yanom, npad, )
        iFxy = np.fft.ifft(Fx.conj()*Fy).real
        varxy = np.sqrt(np.inner(xanom,xanom) * np.inner(yanom,yanom))
    else:
        npad = x.shape[axis] + y.shape[axis]
        if axis == 1:
            if x.shape[0] != y.shape[0]:
                raise ValueError("Arrays should have the same length!")
            xanom = (x - x.mean(axis=1)[:,None])
            yanom = (y - y.mean(axis=1)[:,None])
            varxy = np.sqrt((xanom*xanom).sum(1) * (yanom*yanom).sum(1))[:,None]
        else:
            if x.shape[1] != y.shape[1]:
                raise ValueError("Arrays should have the same width!")
            xanom = (x - x.mean(axis=0))
            yanom = (y - y.mean(axis=0))
            varxy = np.sqrt((xanom*xanom).sum(0) * (yanom*yanom).sum(0))
        Fx = np.fft.fft(xanom, npad, axis=axis)
        Fy = np.fft.fft(yanom, npad, axis=axis)
        iFxy = np.fft.ifft(Fx.conj()*Fy,n=npad,axis=axis).real
    #
    return iFxy/varxy


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # x = np.sin(np.random.random(10000)) # some data
    # y = np.roll(x, 100)                 # shift by 100
    # print(np.argmax(correlate(x,y)))          # should be 100

    # x = np.sin(np.random.random(2000)).reshape(10,200)
    # y = np.roll(x, 100, axis=1)
    # ii = np.argmax(correlate(x,y,axis=1),axis=1)   # should be 100 and 300
    # ii = np.where(ii > 200, ii-200, ii)
    # print(ii)
