#!/usr/bin/env python
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
    
        Example
        -------
        # Create some data
        # for example an island in the middle of an ocean
        >>> import numpy as np
        >>> a = np.array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
        ...        [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
        ...        [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
        ...        [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
        ...        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])
        >>> nn = list(np.shape(a))
        >>> nn.insert(0,2)
        >>> a3 = np.empty(nn)
        >>> for i in range(2): a3[i,...] = a

        # Pack array to  keep only the island elements
        # Mask
        >>> from pack import *
        >>> mask = a == 1.0
        >>> b = pack(a, mask)
        >>> b3 = pack(a3, mask)
        >>> c = unpack(b, mask)
        >>> print np.any(a-c != 0.)
        False
        >>> print c
        [[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  0.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  0.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]
        >>> c3 = unpack(b3, mask)
        >>> print np.any(a3-c3 != 0.)
        False
        >>> 
        >>> d = unpack(b, mask, value=2)
        >>> print d
        [[ 2.  2.  2.  2.  2.  2.  2.  2.  2.  2.]
         [ 2.  2.  2.  1.  1.  1.  2.  2.  2.  2.]
         [ 2.  2.  2.  1.  1.  1.  2.  2.  2.  2.]
         [ 2.  2.  2.  1.  1.  1.  2.  2.  2.  2.]
         [ 2.  2.  2.  2.  2.  2.  2.  2.  2.  2.]]
       
        History
        -------
        Written, MC, Jul. 2009
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
    if darray[-1] < ntmask:
        print 'UNPACK: Last dimension of input array %s must have same size as true values in mask %s.' % (darray[-1], ntmask)
        return None
    #
    # Make multi mask array array
    icount = nmask
    for i in range(ndarray-1):
        icount *= darray[i]
    mask1d = np.ravel(mask)
    masknd = np.empty(icount, dtype='bool')
    for i in range(icount/nmask):
        masknd[i*nmask:(i+1)*nmask] = mask1d[:]
    #
    # Make indeces
    index = np.arange(icount)
    ii = index[masknd]
    if len(ii) != narray:
        print 'UNPACK: Indes creation failed. Index %s Array %s.' % (len(ii), narray)
        return None    
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
    doctest.testmod()
