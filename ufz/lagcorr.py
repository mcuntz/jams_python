#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.correlate import correlate

def lagcorr(y1, y2, max=True, plot=False):

    """
        Calculates the time lag between two 1D arrays where the correlation
        of both arrays is maximum or minimum.

        
        Definition
        ----------
        def lagcorr(y1, y2, max=True, plot=False):


        Input
        -----
        y1           np.ma.array(N) or unmasked array of size N
        y2           np.ma.array(N) or unmasked array of size N


        Optional Input
        --------------
        max          bool: if True, returns lag of maximum correlation, else
                           returns lag of minimum correlation
        plot         bool: if True, correlation function is plotted


        Output
        ------
        lag          int: lag of maximum or minimum correlation between y1 and y2


        Examples
        --------
        >>> # Create some data
        >>> y1 = np.sin(np.arange(np.pi,np.pi*3,0.1))[2:]
        >>> y2 = np.sin(np.arange(np.pi,np.pi*3,0.1))[:-2]
        >>> # calculate correlation lag
        >>> lagcorr(y1, y2)
        2

        >>> x = np.sin(np.random.random(1000)) # some data
        >>> y = np.roll(x, -100)                 # shift by 100 to the right
        >>> lagcorr(x, y)
        -100
        >>> y = np.roll(x, 100)                 # shift by 100 to the left
        >>> lagcorr(x, y)
        100


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

        Copyright 2014 Arndt Piayda, Matthias Cuntz


        History
        -------
        Written,  AP, Mar 2014
        Modified, AP & MC, Apr 2014 - use correlate
    """

    # cross correlation with Fast Fourier
    cc = correlate(y1,y2)
    # minimum or maximum
    if max:
        out = np.argmax(cc)
    else:
        out = np.argmin(cc)
    # if negative lag
    if out > y1.size: out = out - 2*y1.size

    # plot
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure('maxcorr')
        sub = fig.add_subplot(111)
        xx = np.arange(cc.size)
        sub.plot(xx,cc,'-')
        sub.plot(xx[out],cc[out],'o')
        plt.show()

    return out

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # # Create some data
    # y1 = np.sin(np.arange(np.pi,np.pi*3,0.1))[2:]
    # y2 = np.sin(np.arange(np.pi,np.pi*3,0.1))[:-2]
    # # calculate correlation lag
    # lagcorr(y1, y2)
