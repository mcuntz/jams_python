#!/usr/bin/env python
from __future__ import print_function
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
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Arndt Piayda


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
