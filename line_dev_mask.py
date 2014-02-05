#!/usr/bin/env python
import numpy as np
from scipy.optimize import fmin
from plotngo import plotngo

def line_dev_mask(x, y, z, p_guess=[1.,0.], plot=False):
    
    """
        Checks an np.ma.array y for outliers deviating from a straight line
        fitted to y(x). Outliers are detected when they deviate z*std(y) from
        the fitted line. Detection loops the fitting and std calculation until
        no outlier is detected anymore. 

        Definition
        ----------
        def line_dev_mask(x, y, z, p_guess=[1.,0.], plot=False):


        Input
        -----
        x            np.ma.array(N)
        y            np.ma.array(N) or unmasked array of size N
        z            int: outlier threshold (line fit +- z*std(y))
        
        
        Optional Input
        --------------
        p_guess      list(2): initial parameters for line fit [slope, offset]
        plot         bool: if True, iterative plotting detection, if False, no
                           plotting


        Output
        ------
        new_mask     mask of y where additionally to the originally masked 
                     elements, all elements are masked that are detected as
                     outliers deviating from the line fit by z*std(y)


        Examples
        --------
        >>> # Create some data
        >>> x = np.ma.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
        >>> y = np.ma.array([1.,1.5,2.,2.5,3.,100.,8.,4.5,5.,5.5], mask=[0,1,0,0,0,0,0,0,0,0])
        
        >>> # detect outliers
        >>> new_mask = line_dev_mask(x,y,0.5)
        >>> # apply new mask to x
        >>> np.ma.array(y.data, mask=new_mask)
        masked_array(data = [1.0 -- 2.0 2.5 3.0 -- -- 4.5 5.0 5.5],
                     mask = [False  True False False False  True  True False False False],
               fill_value = 1e+20)
        <BLANKLINE>

        
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

        Copyright 2009-2012 Matthias Cuntz


        History
        -------
        Written, AP, Feb 2014
    """
    
    if plot:
        import matplotlib.pyplot as plt
        
    # create working array
    y_new = np.ma.copy(y)
    
    # loop until no outliers are detected anymore
    go_on = True
    while go_on:
        # fit line to data
        p_opt = fmin(absdif, p_guess, args=(x, y_new), disp=0)
        
        # calculate maximum and minimum deviation limit depending on std
        max_dev = line(x,[p_opt[0],p_opt[1]+np.ma.std(y_new)*z])
        min_dev = line(x,[p_opt[0],p_opt[1]-np.ma.std(y_new)*z])
        
        if plot:
            fig = plt.figure('line_dev_mask')
            sub = fig.add_subplot(111)
            sub.plot(x,y_new,'ko')
            sub.plot(x,line(x,p_opt),'k-')
            sub.plot(x,min_dev,'r-')
            sub.plot(x,max_dev,'r-')
            plt.show()
        
        # check for limit trespass
        new_mask = (y_new>max_dev) | (y_new<min_dev)
        y_new.mask += new_mask
        
        # go on when outliers are detected
        go_on = True if np.sum(new_mask)>0 else False
    
    # return mask where outliers are masked
    return y_new.mask
    
def line(x, p):
    
    return p[0] * x + p[1]

def absdif(p,x,y):
    return np.ma.sum(np.ma.abs(y-line(x, p)))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
