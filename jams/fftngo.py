#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt
import warnings          as wa
from   jams.const import eps

def fftngo(t, y, nbins=False, plot=False):
    """
        Fast fourier transformation for dummies (like me). Takes time and
        amplitude arrays and calculates frequency spectrum.


        Definition
        ----------
        def fftngo(t, y, nbins=False, plot=False):


        Input
        -----
        t            numpy.ndarray(N), with independent values (e.g. time or space)
        y            numpy.ndarray(N), with dependent values (e.g. amplitude)


        Optional Input
        --------------
        nbins        int if not False, number of bins to equally bin the frequency
                     spectrum on a logatithmic scale (default: False, no binning).
                     if not False, output caontains binned frequency spectrum additionaly
                     to original spectrum.
        plot         bool, creates a plot of input and output if True (default: False)


        Output
        ------
        X            numpy.ndarray(N/2), frequency spectrum
        Y            numpy.ndarray(N/2), amplitude spectrum
        binm         if nbins not False, numpy.ndarray(nbins-1), bin means
        mv           if nbins not False, numpy.ndarray(nbins-1), bin values


        Examples
        --------
        >>> from jams import around
        >>> # generate example data
        >>> Ts = 0.5 # [h]               # sampling interval of half hour
        >>> Fs = 1./Ts # [1/h]           # sampling frequency
        >>> t = np.arange(0, 365*24, Ts) # one year time data

        >>> tt = 24. # [h]           # period time of observed phenomena
        >>> ff = 1./tt # [1/h]       # frequency of observed phenomena
        >>> s = np.sin(2*np.pi*t*ff) # one year of phenomena data

        >>> # calculate frequency spectrum
        >>> X, Y = fftngo(t, s)
        
        >>> # get dominant frequency
        >>> fd = X[np.argmax(Y)]
        
        >>> # period equals 24 h
        >>> print(int(1./fd))
        24
        
        >>> # get binned frequency spectrum additionally to standard output
        >>> X, Y, binm, mv = fftngo(t, s, nbins=50)


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  AP, Oct 2016
                  MC, Nov 2016 - const.tiny -> const.eps
    """
    # input check
    assert isinstance(t, np.ndarray) & isinstance(y, np.ndarray), 'x and y must be numpy arrays'
    assert (t.ndim == 1) & (y.ndim == 1),                         'x and y must be 1d arrays'
    assert (~np.isnan(t).any()) & (~np.isnan(y).any()),           'no nan allowed in x and y'
    assert (~(np.abs(np.diff(t,2)) > eps)).any(),                'x must contain equal time steps'
    assert isinstance(nbins, int) | (~nbins),                     'nbins must be integer'
    
    # get half size of samples
    n = t.size//2
    
    # get sample frequency 
    dt = t[1] - t[0]
    
    # generate x axis: frequencies
    X = np.fft.fftfreq(n*2, d=dt)[:n]
    
    # generate y axis: amplitudes
    Y = np.fft.fft(y)
    
    # use only abslute, real part of transformation and normalize to half frequencies
    Y = 2.*np.abs(Y[:n])/n

    # calculate logarithmic frequency binning
    if nbins:
        # transform to linear
        Xl, Yl = np.log10(X[1:]), np.log10(Y[1:])
        
        # generate bins
        bins = np.linspace(Xl[0], Xl[-1], nbins)
        
        # binning
        digi = np.digitize(Xl, bins)
        
        # get mean values for bins
        with wa.catch_warnings():
            wa.simplefilter("ignore")
            binm = np.array([Xl[digi == i].mean() for i in range(1, len(bins))])
            mv   = np.array([Yl[digi == i].mean() for i in range(1, len(bins))])
        
        # back-transform
        binm, mv = 10**binm, 10**mv
    
    if plot:
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(t,y)
        ax[0].set_xlabel('t')
        ax[0].set_ylabel('y')
        ax[1].plot(X, Y)
        if nbins:
            ax[1].plot(binm, mv, 'ro')
        ax[1].set_xlabel('Frequency [1/Unit of t]')
        ax[1].set_ylabel('Amplitude [Unit of y]')
        ax[1].set_xscale("log", nonposx='clip')
        ax[1].set_yscale("log", nonposx='clip')
        plt.show()
    
    if nbins:
        return X, Y, binm, mv
    else:
        return X, Y


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
