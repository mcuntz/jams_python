#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.correlate import correlate

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
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2014 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
