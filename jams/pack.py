#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['pack']


def pack(array, mask):
    """
    Mimics Fortran pack intrinsic (without optional vector).

    Packs the last dimensions of an arbitrary shaped array
    into a one dimensional array under a mask.

    The mask can have any dimensions up to the array dimensions


    Definition
    ----------
    def pack(array, mask):


    Input
    -----
    array         input ND-array
    mask          boolean ND-array with dimensions <= array dimensions


    Output
    ------
    results: array with mask dimensions-1 less dimensions than input array
             Last dimension has only elements that correspond to true
             elements of mask


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
    >>> mask = a == 1.0
    >>> b = pack(a, mask)
    >>> print(sum(b))
    9.0
    >>> print(a.sum())
    9.0
    >>> print(b)
    [1. 1. 1. 1. 1. 1. 1. 1. 1.]
    >>> b3 = pack(a3, mask)
    >>> print(a3.sum())
    18.0
    >>> print(b3.sum())
    18.0
    >>> print(b3)
    [[1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1.]]


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python
    library, Department of Computational Hydrosystems, Helmholtz Centre for
    Environmental Research - UFZ, Leipzig, Germany.

    Copyright (c) 2009-2021 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Matthias Cuntz, Jul 2009
    Modified, Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Apr 2014 - assert
              Matthias Cuntz, Sep 2021 - code refactoring
    """
    dmask   = mask.shape
    ndmask  = np.ndim(mask)
    nmask   = mask.size
    darray  = array.shape
    ndarray = np.ndim(array)
    narray  = array.size
    #
    # Check array and mask
    assert ndarray >= ndmask, 'Input array has less dimensions ' + str(ndarray) + ' then mask ' + str(ndmask)
    k = 0
    while k > -ndmask:
        k -= 1
        assert dmask[k] == darray[k], 'Input array and mask must have the same last dimensions. Array: ' + str(darray) + ' Mask: ' + str(dmask)
    #
    # Make array and mask 1d
    farray = array.ravel()  # flat in Fortran=column-major mode
    fmask  = mask.ravel()
    afmask = np.empty(narray, dtype=bool)
    nn = narray // nmask
    k  = 0
    while k < nn:
        afmask[k*nmask:(k+1)*nmask] = fmask[:]
        k += 1
    #
    # Mask array and reshape
    afarray = farray[afmask]
    dout = list(darray)
    k = 0
    while k < ndmask:
        del dout[-1]
        k += 1
    nnmask = mask.sum()
    dout.append(nnmask)
    out = np.reshape(afarray, dout)
    #
    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # a = np.array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
    #               [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
    #               [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
    #               [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
    #               [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])
    # nn = list(np.shape(a))
    # nn.insert(0,2)
    # a3 = np.empty(nn)
    # for i in range(2): a3[i,...] = a

    # # Pack array to  keep only the island elements
    # # Mask
    # mask = a == 1.0
    # b = pack(a, mask)
    # print sum(b)
    # #9.0
    # print a.sum()
    # #9.0
    # b
    # #array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    # b3 = pack(a3, mask)
    # print a3.sum()
    # #18.0
    # print b3.sum()
    # #18.0
    # b3
    # #array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
    # #       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]])
