#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

def rolling(x, win):
    """
        Uses binary strides to reshape an array in a "rolling window" style
        to perform any function on a rolling window without for loop


        Definition
        ----------
        def rolling(x, win):


        Input
        -----
        x         1D np.array(N), array to be rehaped
        win       int, window size


        Output
        ------
        x2D,      2D np.array(N-(win+1),win), x reshaped in windows


        Restrictions
        ------------
        The method is only applicable for 1D numpy arrays.


        References
        ----------
        This routine is taken from Alex Rogozhnikov 
        (http://arogozhnikov.github.io/)


        Example
        --------
        >>> x = np.arange(10)

        # get windowed array with window size 3
        >>> win = 3
        >>> print(rolling(x, win).shape)
        (8, 3)


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

        Copyright 2016 Arndt Piayda


        History
        -------
        Written,  AP, Dec 2016
    """
    # test input
    assert x.ndim == 1,          'x must be 1D'
    assert isinstance(win, int), 'win must be integer'

    s   = x.strides[0]
    x2D = np.lib.stride_tricks.as_strided(x, shape=[x.size - win+ 1, win], strides=[s, s])
    return x2D

    
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
