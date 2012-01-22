#!/usr/bin/env python
import numpy as np

def pack(array, mask):
    """
        Mimics Fortran pack intrinsic (without optional vector).
        Packs the last dimensions of an arbitrary shaped array
        into a one dimension under a mask.
        The mask can have any dimensions up to teh array dimensions

        Definition
        ----------
        def pack(array, mask):


        Input
        -----
        array         input array
        mask          boolean array with dimensions <= array dimensions


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
        >>> mask = a == 1.0
        >>> b = pack(a, mask)
        >>> print sum(b)
        9.0
        >>> print a.sum()
        9.0
        >>> b
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
        >>> b3 = pack(a3, mask)
        >>> print a3.sum()
        18.0
        >>> print b3.sum()
        18.0
        >>> b3
        array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
               [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]])
        

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
    if ndarray < ndmask:
        print 'PACK: Input array has less dimensions %s then mask %s.' % (ndarray, ndmask)
        return None
    k = 0
    while k > -ndmask:
        k -= 1
        if dmask[k] != darray[k]:
            print 'PACK: Input array and mask have the same last dimensions. Array: %s Mask: %s' % (darray, dmask)
            return None
    #
    # Make array and mask 1d
    farray = np.ravel(array) # flat in Fortran=column-major mode
    fmask  = np.ravel(mask)
    afmask = np.empty(narray, dtype='bool')
    nn = narray / nmask
    k = 0
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
    doctest.testmod()
