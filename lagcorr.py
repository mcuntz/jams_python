#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def lagcorr(y1, y2, lagmin, lagmax, max=True, plot=False):
    
    """
        Calculates the time lag between two 1D arrays where the correlation
        of both arrays is maximum or minimum.

        Definition
        ----------
        def lagcorr(x, y, lagmin, lagmax, max=True, plot=False):


        Input
        -----
        y1           np.ma.array(N) or unmasked array of size N
        y2           np.ma.array(N) or unmasked array of size N
        lagmin       int: minimum time lag between y1 and y2 to start with
        lagmax       int: maximum time lag between y1 and y2 to end with
        
        
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
        >>> lagcorr(y1, y2, -10, 10)
        -2
        
        
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

        Copyright 2014 Arndt Piayda, Matthias Cuntz


        History
        -------
        Written,  AP, Mar 2014
    """
    
    # check input
    if not isinstance(lagmin, int) or not isinstance(lagmax, int):
        raise ValueError('maxcorr: lagmin and lagmax have to be integer!')
    if lagmax<lagmin:
        raise ValueError('maxcorr: max has to be bigger than min!')
    y1_size, y2_size = y1.size, y2.size
    if y1_size!=y2_size:
        raise ValueError('maxcorr: y1 and y2 have to be of equal size!')
    
    # check if masked or not
    try:
        temp = y1.mask
        temp = y2.mask
    except AttributeError:
        y1=np.ma.array(y1, mask=np.isnan(y1))
        y2=np.ma.array(y2, mask=np.isnan(y2))
    
    # generate span
    span = np.arange(lagmin,lagmax+1)
    if span.size>=y2_size:
        raise ValueError('maxcorr: lagmin-lagmax span is larger than array length!')
    correlations = []
    
    # rolling over span
    for i in span:
        if i<0:
            y1_roll = y1[:i]
            y2_roll = y2[-i:]
        elif i==0:
            y1_roll = y1
            y2_roll = y2
        elif i>0:
            y1_roll = y1[i:]
            y2_roll = y2[:-i]
        
        # mask matching
        total_mask  = y1_roll.mask | y2_roll.mask
        y1_roll.mask = total_mask
        y2_roll.mask = total_mask
        
        # calculate correlations for each lag
        correlations += [np.correlate(y1_roll.compressed(),y2_roll.compressed())[0]]        
        
    # get index of maximum or minimum correlation
    index=np.argmax(correlations) if max else np.argmin(correlations) 
    
    # plot
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure('maxcorr')
        sub = fig.add_subplot(111)
        sub.plot(span,correlations,'-')
        sub.plot(span[index],correlations[index],'o')
        plt.show()
    
    # return lag of maximum or minimum correlation
    return span[index]

# http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg01487.html
def ccf(x, y, axis=None):    
    """Computes the cross-correlation function of two series `x` and `y`.
Note that the computations are performed on anomalies (deviations from 
average).
Returns the values of the cross-correlation at different lags.
Lags are given as [0,1,2,...,n,n-1,n-2,...,-2,-1].
 
:Parameters:
    `x` : 1D MaskedArray
        Time series.
    `y` : 1D MaskedArray
        Time series.
    `axis` : integer *[None]*
        Axis along which to compute (0 for rows, 1 for cols).
        If `None`, the array is flattened first.
    """
    assert x.ndim == y.ndim, "Inconsistent shape !"
#    assert(x.shape == y.shape, "Inconsistent shape !")   
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
                raise ValueError, "Arrays should have the same length!"
            xanom = (x - x.mean(axis=1)[:,None])
            yanom = (y - y.mean(axis=1)[:,None])
            varxy = np.sqrt((xanom*xanom).sum(1) * (yanom*yanom).sum(1))[:,None]
        else:
            if x.shape[1] != y.shape[1]:
                raise ValueError, "Arrays should have the same width!"
            xanom = (x - x.mean(axis=0))
            yanom = (y - y.mean(axis=0))
            varxy = np.sqrt((xanom*xanom).sum(0) * (yanom*yanom).sum(0))
        Fx = np.fft.fft(xanom, npad, axis=axis)
        Fy = np.fft.fft(yanom, npad, axis=axis)
        iFxy = np.fft.ifft(Fx.conj()*Fy,n=npad,axis=axis).real
    #    
    return iFxy/varxy


def ccf1d(x,y):
    """Computes the crosscorrelation of two flat arrays `x` and `y`, with the
numpy.correlate function.
Note that the computations are performed on anomalies (deviations from 
average).
    """
    if x.ndim > 1:
        x = x.ravel()
    if y.ndim > 1:
        y = y.ravel()
    (xanom, yanom) = (x-x.mean(), y-y.mean())
    corxy = np.correlate(xanom, yanom, 'full')
    n = min(x.size, y.size)
    #    return np.r_[ xc[len(yf)-1:], 0, xc[:len(yf)-1] ]
    corxy = np.r_[ corxy[:n][::-1], 0, corxy[n:][::-1] ]
    varxy = np.sqrt(np.inner(xanom,xanom) * np.inner(yanom,yanom))
    return corxy/varxy

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    '''
    # %paste in iPython
    from lagcorr import lagcorr, ccf, ccf1d
    x = np.sin(np.random.random(10000)) # some data
    y = np.roll(x, 100)                 # shift by 100
    %timeit ccf(x,y)                    # with FFT
    %timeit ccf1d(x,y)                  # with correlate 'full'
    %timeit lagcorr(x, y, 50, 200)      # with correlate 'valid'
    print(lagcorr(x, y, 50, 200))          # should be 100
    print(np.argmax(ccf(x,y)))                 # should be 100
    '''

