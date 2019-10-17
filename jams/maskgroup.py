#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def maskgroup(x, n):
    """
        Masks elements in masked 1d-array gathered in groups of n elements.
        Retains any existing mask.

        For example for n=1: if an unmasked value is surrounded by two masked values, it will be masked as well.
        Or for n=2: two consecutive unmasked elements surrounded by masked elements will be masked.


        Definition
        ----------
        def maskgroup(x, n):


        Input
        -----
        x            1d-array
        n            int: any number <= n of consecutive unmasked data points surrounded by masked elements
                          should get masked


        Output
        ------
        new_mask     mask of x where additionally to the originally masked
                     elements, all elements are masked that are in a consecutive
                     group of data points with a length smaller equal to n


        Examples
        --------
        >>> # Create some data
        >>> x=np.ma.array([1,2,3,4,5,6,7,8,9,10,11,12], mask=[1,0,0,0,0,1,0,0,0,1,0,0])

        >>> # mask all elements which are gathered in groups smaller then 3
        >>> new_mask = maskgroup(x,2)
        >>> # apply new mask to x
        >>> print(np.ma.array(x.data, mask=new_mask))
        [-- 2 3 4 5 -- 7 8 9 -- -- --]

        >>> # mask all elements with are gathered in groups smaller equal to 3
        >>> new_mask = maskgroup(x,3)
        >>> # apply new mask to x
        >>> print(np.ma.array(x.data, mask=new_mask))
        [-- 2 3 4 5 -- -- -- -- -- -- --]


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
        Modified, MC, Feb 2014 - call it maskgroup instead of small_kickout1d, mask <=n instead of <n
                  MC, Apr 2014 - removed enumerate and zip for simplification
    """
    # create running lists
    index  = []
    length = []
    count  = 0

    # count indices and lengths of unmasked data groups
    for i in range(x.mask.size):
        if i==0:
            if not x.mask[i]:
                index += [i]
                count  = 1
        if i>0:
            if not x.mask[i] and x.mask[i-1]:
                index += [i]
                count  = 1
            elif not x.mask[i]:
                count += 1
            elif x.mask[i] and not x.mask[i-1]:
                length += [count]
                count = 0
            else:
                pass
    if count>0:
        length += [count]

    # create new mask
    new_mask = np.copy(x.mask)
    for i in range(len(index)):
        if length[i]<=n:
            new_mask[index[i]:index[i]+length[i]] = True

    return new_mask

if __name__ == '__main__':
    import doctest
    doctest.testmod()
