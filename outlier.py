#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.stats as s
#import pdb

def outlier(y,alpha=0.01,k=-1,quiet=True):
    """
        Estimates the number of outliers in an approximately normal distributed sample
        according to Rossner''s (1983) generalized Extreme Studentized Deviate (ESD).
        It returns the indeces of the outliers.

        Definition
        ----------
        def outlier(y, alpha=alpha, k=k, quiet=quiet)


        Input
        -----
        y            Approx. normal distributed data


        Options
        -------
        alpha        Significance level (default: 0.01)
        k            Numbers of outliers to be checked (default: size(data)/2)
        quiet        True: Perform quietly (default: True)


        Output
        ------
        Indeces of outliers in input array
        

        Restrictions
        ------------
        Assumptions: The input sample is assumed to be approximately normal distributed


        References
        ----------
        Rosner B, Percentage Points for a Generalized ESD Many-Outlier Procedure,
            Technometrics, 25(2), 165-172, 1983
        http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm (July 2010)
        Yu RC, Teh HW, Jaques PA, Sioutas C, Froines JR,
            Quality control of semi-continuous mobility size-fractionated
            particle number concentration data; atmospheric Environment 38, 3341- 3348, 2004
        

        Examples
        --------
        >>> import numpy as np
        >>> y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,
        ...               2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,
        ...               4.64,5.34,5.42,8.01],dtype=np.float)
        >>> print(outlier(y,0.05,2,quiet=False))
        # of outliers
        1
        outliers
        [ 8.01]
        indeces
        [25]
        [25]
        >>> #
        >>> # Example from Rosner paper
        >>> y = np.array([-0.25,0.68,0.94,1.15,1.20,1.26,1.26,1.34,1.38,1.43,1.49,
        ...               1.49,1.55,1.56,1.58,1.65,1.69,1.70,1.76,1.77,1.81,1.91,
        ...               1.94,1.96,1.99,2.06,2.09,2.10,2.14,2.15,2.23,2.24,2.26,
        ...               2.35,2.37,2.40,2.47,2.54,2.62,2.64,2.90,2.92,2.92,2.93,
        ...               3.21,3.26,3.30,3.59,3.68,4.30,4.64,5.34,5.42,6.01],dtype=np.float)
        >>> print(outlier(y,0.05,10,quiet=False))
        # of outliers
        3
        outliers
        [ 5.34  5.42  6.01]
        indeces
        [51 52 53]
        [51 52 53]
        >>> print(outlier(y,0.05,quiet=False))
        # of outliers
        3
        outliers
        [ 5.34  5.42  6.01]
        indeces
        [51 52 53]
        [51 52 53]
        >>> print(outlier(y,quiet=False))
        outlier: Found no outliers.
        [-1]


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

        Copyright 2010-2013 Maren Goehler, Matthias Cuntz


        History
        -------
        Written,  MG,   Jul 2010
        Modified, MC,   Mar 2011 - formulation was wrong: after removal of a point, the 
                                   absolute deviation was not recalculated with the new mean
                                   it gives exactly the same as the example in the above NIST-webpage
                                 - calc mean and standard deviation with cummulative formula
                                   so that there is only three sums of whole array
                                 - include case with no outliers: returns -1.
               MC & MG, Oct 2011 - break loop if new variance < 0 in iteration
               MC,      Feb 2013 - ported to Python 3
    """
    # check user input
    if ((alpha<=0) or (alpha>1)):
        raise ValueError('alpha must between 0 and 1.')

    # Student''s t-distribution
    nn    = np.ma.count(y)
    if k == -1:
        k = nn//2
    if k > nn//2:
        print('WARNING: k > size(data)/2 -> set k=size(data)/2')
        k = nn-1
    out      = 0
    outlrs   = np.zeros(k,dtype=np.float)
    iioutlrs = np.zeros(k,dtype=np.int)
    i     = np.arange(k,dtype=np.int) + 1
    nnip1 = np.array(nn-i+1,dtype=np.float)
    nnim1 = np.array(nn-i-1,dtype=np.float)
    pcrit = 1.-(alpha/(2.*nnip1))
    t     = s.t.ppf(pcrit,nnim1)
    lamb  = ((nnim1+1.)*t)/np.sqrt((nnim1+t*t)*nnip1)

    # Calc mean and stddev ourselves (2-pass formula) so that we understand
    # the cummulative way to calculate mean and stddev below
    #mm = np.ma.mean(y[0:nn])
    #ss  = np.std(y[0:nn],ddof=1)
    mm       = np.ma.sum(y)/nn
    dev      = y-mm
    ep       = np.ma.sum(dev)
    tt       = np.ma.sum(dev*dev)
    ss       = np.ma.sqrt((tt-ep/nn)/(nn-1))
    
    adev        = np.ma.abs(dev)
    thisii      = np.ma.argmax(adev)
    iioutlrs[0] = thisii
    outlrs[0]   = y[thisii]
    # normalised deviation R
    R   = (np.ma.abs(y[thisii]-mm)) / ss
    # Test R against t-distribution
    if R > lamb[0]:
        out  = 1
    
    if k > 1:
        for i in range(k-1):
            # remove the highest value from above and redo the procedure
            mm1 = mm
            # Cummulative way to calculate mean and stddev
            mm  = (mm*(nn-i) - y[thisii])/(nn-1-i)
            ep  = ep + mm1*(nn-i) - y[thisii] - mm*(nn-1-i)
            tt  = tt - y[thisii]*y[thisii] + (nn-i)*mm1*mm1 - (nn-1-i)*mm*mm
            if ((tt-ep/(nn-1-i)) < 0.):
                break
            ss     = np.sqrt((tt-ep/(nn-1-i))/(nn-1-i-1))
            adev   = np.ma.masked_array(np.ma.abs(y-mm), mask=np.ma.mask_or(adev.mask,y==y[thisii]))
            thisii = np.ma.argmax(adev)
            iioutlrs[i+1] = thisii
            outlrs[i+1]   = y[thisii]
            R  = (np.abs(y[thisii]-mm)) / ss
            if R > lamb[i+1]:
                out = i+2

    
    if out > 0:
        # report results
        iiout = iioutlrs[0:out]
        iiout = iiout[::-1] # reverse
        if (quiet == False):
            print('# of outliers')
            print(out)
            print('outliers')
            print(y[iiout])
            print('indeces')
            print(iiout)
        return iiout
    else:
        if (quiet == False):
            print('outlier: Found no outliers.')
        return np.array([-1])

if __name__ == '__main__':
    import doctest
    doctest.testmod()
