#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
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
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2016 Arndt Piayda

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
