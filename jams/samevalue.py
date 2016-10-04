# -*- coding: utf-8 -*-
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

        Copyright 2015 Andreas Wiedemann


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
