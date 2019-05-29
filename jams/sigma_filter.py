#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.optimize as opt
import jams.functions as functions

def sigma_filter(x, y, z=3, func=functions.line_p, p=[0.,1.], plot=False, popt=False):
    
    """
        Checks a 1D-array for points deviating more than z standard deviations
        from a given function.
        The default function is a line and a 3-sigma filter is applied.
        If the data is very noisy then the filter profits from good starting parameters.


        Definition
        ----------
        def sigma_filter(x, y, z=3, func=functions.line_p, p=[0.,1.], plot=False):


        Input
        -----
        x            1D-array, or ND-array, which is filtered in the last dimension, i.e. x.shape[-1] == y.size
        y            1D-(masked) array
        
        
        Optional Input
        --------------
        z            multiplier to standard deviation used so that all deviations
                     abs(y-func(x,p)) > z*std(abs(y-func(x,p))) are masked
                     Default: 3
        func         func object. Default: functions.line_p
        p            initial function parameters, Default: [0,1]
        plot         True: plot iterative fits; Default: False
        popt         True: output also optimal parameters after filtering; False: output only mask
                     Default: False


        Output
        ------
        if popt:
            new_mask, popt    mask of y where y deviates more than z*std from fitted function, optimal parameters
        else:
            new_mask          mask of y where y deviates more than z*std from fitted function


        Examples
        --------
        >>> # Create some data
        >>> np.random.seed(1)
        >>> x = np.arange(40)
        >>> y = np.ma.array(2.+1.5*x+np.random.normal(0,1,40)) # line
        >>> y[2] = np.ma.masked
        >>> y[10] *= 10.
        >>> y[20] *= 10.
        >>> y[30] /= 10.
        >>> # detect outliers
        >>> new_mask = sigma_filter(x, y, z=5, plot=False)
        >>> # apply new mask to x
        >>> from autostring import astr
        >>> print(astr(np.ma.array(y, mask=new_mask), 1)[0:12])
        [' 3.6' ' 2.9' '--  ' ' 5.4' ' 8.9' ' 7.2' '12.7' '11.7' '14.3' '15.3' '--  ' '16.4']

        >>> def fun(x, p):
        ...     return(p[0]*np.sin(p[1]*x))
        >>> y = np.ma.array(fun(x, [2.,1.5])+np.random.normal(0,0.1,40)) # sine wave
        >>> y[2] = np.ma.masked
        >>> y[10] *= 10.
        >>> y[20] *= 10.
        >>> # detect outliers - line
        >>> print(sigma_filter(x, y, z=3, plot=False)[0:12])
        [False  True  True  True False  True False  True False False  True  True]
        >>> # detect outliers - sine
        >>> print(sigma_filter(x, y, z=5, p=[2.1,1.4], func=fun, plot=False)[0:12])
        [False False  True False False False False False False False  True False]
        
        >>> def fun2(x, p):
        ...     return(p[0]+p[1]*np.sin(p[2]*x))
        >>> y = np.ma.array(fun2(x, [2.,1.5,3.])+np.random.normal(0,0.1,40)) # sin wave
        >>> y[2] = np.ma.masked
        >>> y[10] *= 10.
        >>> y[20] *= 10.
        >>> # detect outliers
        >>> s, p = sigma_filter(x, y, z=5, p=[2.1,1.4,3.1], func=fun2, plot=False, popt=True)
        >>> print(s[0:12])
        [False False  True False False False False False False False  True False]
        >>> print(astr(p, 1))
        ['2.0' '1.5' '3.0']

        
        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2018 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2014-2018 Matthias Cuntz


        History
        -------
        Written,  MC, Feb 2014 - changed line_dev_mask
        Modified, MC, Dec 2014 - x.compressed error because x is non-masked
                  MC, May 2018 - pout, x ND-array
    """
    if plot:
        import matplotlib.pyplot as plt
        
    # create masked work array
    y_new  = np.ma.array(y.copy())
    if np.ma.count_masked(y_new) == 0: y_new.mask = np.zeros(y_new.shape)
    # Initialise fitted  parameters
    p_opt = np.copy(p)
    
    # loop until no outliers are detected anymore
    go_on = True
    while go_on:
        # fit function to data
        mm = np.where(~y_new.mask)[0]
        if x.ndim == 1:
            xx = x[mm]
        else:
            xx = x[...,mm]
        yy = y_new.compressed()
        p_opt = opt.fmin(functions.cost_abs, p_opt, args=(func,xx,yy), disp=False)

        # look for anomalies > z*standard dev
        dev  = np.abs(yy - func(xx,p_opt))
        sdev = np.std(dev, ddof=1)
        ii   = np.where(dev > (z*sdev))[0]

        if plot:
            fig = plt.figure('sigma_filter')
            sub = fig.add_subplot(111)
            sub.plot(x, y_new, 'ko-')
            iix = np.argsort(x)
            sub.plot(x[iix], func(x[iix],p_opt), 'b--')
            sub.plot(x[iix], func(x[iix],p_opt)-z*sdev, 'r--')
            sub.plot(x[iix], func(x[iix],p_opt)+z*sdev, 'r--')
            plt.show()

        if ii.size>0:
            y_new[mm[ii]] = np.ma.masked
        else:
            go_on = False

    # return mask where outliers are masked
    if popt:
        return y_new.mask, p_opt
    else:
        return y_new.mask


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # # Create some data
    # np.random.seed(1)
    # x = np.arange(40)
    # y = np.ma.array(2.+1.5*x+np.random.normal(0,1,40)) # line
    # y[2] = np.ma.masked
    # y[10] *= 10.
    # y[20] *= 10.
    # y[30] /= 10.
    # # detect outliers
    # new_mask = sigma_filter(x, y, z=5, plot=True)
    # # apply new mask to x
    # from autostring import astr
    # print(astr(np.ma.array(y, mask=new_mask), 1))
    # # [' 3.6' ' 2.9' '--  ' ' 5.4' ' 8.9' ' 7.2'
    # #  '12.7' '11.7' '14.3' '15.3' '--  ' '16.4'
    # #  '19.7' '21.1' '24.1' '23.4' '25.8' '26.6'
    # #  '29.0' '31.1' '--  ' '34.6' '35.9' '37.0'
    # #  '38.9' '38.8' '40.9' '41.6' '43.7' '46.0'
    # #  '--  ' '48.1' '49.3' '50.7' '52.3' '54.5'
    # #  '54.9' '57.7' '60.7' '61.2']

    # def fun(x, p):
    #     return(p[0]*np.sin(p[1]*x))
    # y = np.ma.array(fun(x, [2.,1.5])+np.random.normal(0,0.1,40)) # sine wave
    # y[2] = np.ma.masked
    # y[10] *= 10.
    # y[20] *= 10.
    # # detect outliers - line
    # print(sigma_filter(x, y, z=3, plot=True))
    # # [False  True  True  True False  True False  True False False  True  True
    # #   True False False  True  True False False False  True False False False
    # #   True False False False  True  True False False  True  True False False
    # #   True  True False False]
    # # detect outliers - sine
    # print(sigma_filter(x, y, z=5, p=[2.1,1.4], func=fun, plot=True))
    # # [False False  True False False False False False False False  True
    # #  False False False False False False False False False  True False
    # #  False False False False False False False False False False False
    # #  False False False False False False False]
    
    # def fun2(x, p):
    #     return(p[0]+p[1]*np.sin(p[2]*x))
    # y = np.ma.array(fun2(x, [2.,1.5,3.])+np.random.normal(0,0.1,40)) # sin wave
    # y[2] = np.ma.masked
    # y[10] *= 10.
    # y[20] *= 10.
    # # detect outliers
    # print(sigma_filter(x, y, z=5, p=[2.1,1.4,3.1], func=fun2, plot=True))
    # # [False False  True False False False False False False False  True False
    # #  False False False False False False False False  True False False False
    # #  False False False False False False False False False False False False
    # #  False False False False]
