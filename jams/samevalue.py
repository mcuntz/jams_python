#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
Created on Sun Aug 30 18:35:32 2015

@author: wiedeman
"""
import numpy as np

"""
        Checks if abs. differences of array values within a certain window are smaller than threshold.
        Return mask. 


        Definition
        ----------
        def samevalue(x,tol,window)


        Input
        -----
        x            1D-array
        tol          tolerance for difference of two subsequent values.
                     If absolute difference between two values is smaller than tol array entries will be masked
        window       size of window where all abs. differences of subsequent values has to be smaller than tol to get masked
                     if at least one difference is larger than tol all array entries within window stay unmasked
        

        Output
        ------
        new_mask     mask of x where abs. difference of two subsequent values within a certain window are smaller that tol


        Examples
        --------
        >>> # Create some data
        >>> x                = np.arange(0.,10.,0.01)
        >>> x_curve          = (np.sin(x))
        >>> x_curve[700:800] = x_curve[699]
        >>> x_mask           = samevalue(x_curve,0.001,50)
        >>> xx               = ma.array(x_curve, mask=x_mask)

        >>> # apply new mask to x
        >>> print(x_mask[650:850])
        [False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True,  True,  True,  True,
        True,  True,  True,  True,  True,  True, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015 Andreas Wiedemann

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
        Written,  AW, Aug 2015

"""
def samevalue(x,tol,window):
    # define mask where input values smaller than tolerance
    x         = np.ma.array(x)        
    diff      = np.abs(np.diff(x))
    tol       = np.float(tol)
    x_mask    = np.ma.masked_where(diff > tol, diff)
    
    # mask x only if values are smaller than tolerance in a range as long or longer than the window-size  
    count = 0    
    for i in range(np.shape(x)[0]):
        if i > window and sum(x_mask.mask[count:][i-window:i]) < 1:
            first_unique_index                      = i-window+1
            last_unique_index                       = np.argmax(diff[first_unique_index:] > tol)
            x[first_unique_index:first_unique_index + last_unique_index] = np.ma.masked
            count += last_unique_index
        else:
            count += 1
        
    return x.mask

if __name__ == '__main__':
    import doctest
    doctest.testmod()
