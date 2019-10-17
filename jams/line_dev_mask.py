#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.optimize as opt
import jams.functions as functions

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
        >>> new_mask = line_dev_mask(x,y,0.5,plot=False)
        >>> # apply new mask to x
        >>> print(np.ma.array(y.data, mask=new_mask))
        [1.0 -- 2.0 2.5 3.0 -- -- 4.5 5.0 5.5]

        
        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

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
        Written,  AP, Feb 2014
        Modified, MC, Feb 2014 - use jams.functions
    """
    
    if plot:
        import matplotlib.pyplot as plt
        
    # create working array
    y_new = np.ma.copy(y)
    
    # loop until no outliers are detected anymore
    go_on = True
    while go_on:
        # fit line to data
        p_opt = opt.fmin(functions.cost_abs, p_guess, args=(functions.line_p,x,y_new), disp=False)
        
        # calculate maximum and minimum deviation limit depending on std
        max_dev = functions.line_p(x,[p_opt[0]+np.ma.std(y_new)*z,p_opt[1]])
        min_dev = functions.line_p(x,[p_opt[0]-np.ma.std(y_new)*z,p_opt[1]])
        
        if plot:
            fig = plt.figure('line_dev_mask')
            sub = fig.add_subplot(111)
            sub.plot(x,y_new,'ko')
            sub.plot(x,functions.line_p(x,p_opt),'k-')
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()
