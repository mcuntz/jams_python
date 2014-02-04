#!/usr/bin/env python
import numpy as np

def small_kickout1d(x, n):
    
    """
        Generates a new mask from an already masked array x. In the new mask,
        all elements are masked which are masked in the original array x.
        Additionally, all originally unmasked elements, which gather in groups
        smaller than n are masked.

        Definition
        ----------
        def small_kickout1d(x, n):


        Input
        -----
        x            np.ma.array(N)
        n            int: maximum number of consecutive unmasked data points
                          that should be masked


        Output
        ------
        new_mask     mask of x where additionally to the originally masked 
                     elements, all elements are masked that are in a consecutive
                     group of data points with a length smaller than n


        Examples
        --------
        >>> # Create some data
        >>> x=np.ma.array([1,2,3,4,5,6,7,8,9,10,11,12], mask=[1,0,0,0,0,1,0,0,0,1,0,0])
        
        >>> # mask all elements with are gathered in groups smaller then 3
        >>> new_mask = small_kickout1d(x,3)
        >>> # apply new mask to x
        >>> np.ma.array(x.data, mask=new_mask)
        masked_array(data = [-- 2 3 4 5 -- 7 8 9 -- -- --],
                     mask = [ True False False False False  True False False False  True  True  True],
               fill_value = 999999)
        <BLANKLINE>
        
        >>> # mask all elements with are gathered in groups smaller then 4
        >>> new_mask = small_kickout1d(x,4)
        >>> # apply new mask to x
        >>> np.ma.array(x.data, mask=new_mask)
        masked_array(data = [-- 2 3 4 5 -- -- -- -- -- -- --],
                     mask = [ True False False False False  True  True  True  True  True  True  True],
               fill_value = 999999)
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
    
    # create running lists
    index = []
    length = []
    count = 0
    
    # count indices and lengths of unmasked data groups
    for i, item in enumerate(x.mask):
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
    for item in zip(index,length):
        if item[1]<n:
            new_mask[item[0]:item[0]+item[1]] = True
    
    return new_mask

if __name__ == '__main__':
    import doctest
    doctest.testmod()