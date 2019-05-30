#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def unpack(array, mask, value=0.):
    """
        Mimics Fortran unpack intrinsic (without optional field).
        Unpacks the last dimension into several dimensions under a mask.
        The unpacked elements will be set to user-defined value.
        The mask can have any dimensions up to the array dimensions.


        Definition
        ----------
        def unpack(array, mask, value=0.):


        Input
        -----
        array         input array
        mask          boolean array with dimensions <= array dimensions


        Optional input
        -----
        value         value of the new elements (default: 0)


        Output
        ------
        results: array with mask-1 dimensions more than input array.
                 The new elements have the user-defined value.


        Restrictions
        ------------
        All mask values false is undefined.


        Examples
        --------
        # Create some data
        # for example an island in the middle of an ocean
        >>> import numpy as np
        >>> a = np.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        ...               [0., 0., 0., 1., 1., 1., 0., 0., 0., 0.],
        ...               [0., 0., 0., 1., 1., 1., 0., 0., 0., 0.],
        ...               [0., 0., 0., 1., 1., 1., 0., 0., 0., 0.],
        ...               [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])
        >>> nn = list(np.shape(a))
        >>> nn.insert(0,2)
        >>> a3 = np.empty(nn)
        >>> for i in range(2): a3[i,...] = a

        # Pack array to  keep only the island elements
        # Mask
        >>> from pack import pack
        >>> mask = a == 1.0
        >>> b = pack(a, mask)
        >>> b3 = pack(a3, mask)
        >>> c = unpack(b, mask)
        >>> print(np.any(a-c != 0.))
        False
        >>> print(c)
        [[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
         [0. 0. 0. 1. 1. 1. 0. 0. 0. 0.]
         [0. 0. 0. 1. 1. 1. 0. 0. 0. 0.]
         [0. 0. 0. 1. 1. 1. 0. 0. 0. 0.]
         [0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]
        >>> c3 = unpack(b3, mask)
        >>> print(np.any(a3-c3 != 0.))
        False
        >>>
        >>> d = unpack(b, mask, value=2)
        >>> print(d)
        [[2. 2. 2. 2. 2. 2. 2. 2. 2. 2.]
         [2. 2. 2. 1. 1. 1. 2. 2. 2. 2.]
         [2. 2. 2. 1. 1. 1. 2. 2. 2. 2.]
         [2. 2. 2. 1. 1. 1. 2. 2. 2. 2.]
         [2. 2. 2. 2. 2. 2. 2. 2. 2. 2.]]


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2009-2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jul. 2009
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    dmask   = np.shape(mask)
    ndmask  = len(dmask)
    nmask   = mask.size
    darray  = np.shape(array)
    ndarray = len(darray)
    narray  = array.size
    #
    # Check array and mask
    ntmask = mask.sum()
    assert darray[-1] == ntmask, 'Last dimension of input array %s must have same size as true values in mask %s.' % (darray[-1], ntmask)
    #
    # Make multi mask array array
    icount = nmask
    for i in range(ndarray-1):
        icount *= darray[i]
    mask1d = np.ravel(mask)
    masknd = np.empty(icount, dtype='bool')
    for i in range(icount//nmask):
        masknd[i*nmask:(i+1)*nmask] = mask1d[:]
    #
    # Make indeces
    index = np.arange(icount)
    ii = index[masknd]
    if len(ii) != narray:
        raise ValueError('Index creation failed. Index %s Array %s.' % (len(ii), narray))
    #
    # Flat output array
    array1d = np.ravel(array)
    arraynd = np.ones(icount)*value
    arraynd[ii] = array1d[:]
    #
    # Reshaped output array
    newdim = list(darray[0:-1])
    for i in range(ndmask):
        newdim.append(dmask[i])
    out = np.reshape(arraynd, newdim)
    #
    return out

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

