#!/usr/bin/env python
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()