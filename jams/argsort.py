#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['argsort']

def argsort(a, *args, **kwargs):
    """
        Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables.

        Passes all keywords directly to the individual routines, i.e.
            numpy.argsort(a, axis=-1, kind='quicksort', order=None)
            numpy.ma.argsort(a, axis=None, kind='quicksort', order=None, fill_value=None)
            sorted(iterable[, cmp[, key[, reverse]]])

        Only key cannot be given for Python iterables because the input array is used
        as key in the sorted function.


        Definition
        ----------
        def argsort(a, *args, **kwargs):


        Input
        -----
        a            numpy.ndarray, numpy.ma.MaskedArray or Python iterable
        all arguments of numpy.argsort, numpy.ma.argsort, or sorted.


        Optional Input
        --------------
        all keyword arguments of numpy.argsort, numpy.ma.argsort, or sorted.


        Output
        ------
        Array of indices that sort a along the specified axis. In other words, a[index_array] yields a sorted a.


        Examples
        --------
        >>> import numpy as np
        >>> a = np.array([0,4,6,2,1,5,3,5])
        >>> ii = argsort(a)
        >>> print(a[ii])
        [0 1 2 3 4 5 5 6]

        >>> ii = argsort(a, kind='quicksort')
        >>> print(a[ii])
        [0 1 2 3 4 5 5 6]

        >>> a = np.ma.array([0,4,6,2,1,5,3,5], mask=[0,0,1,1,0,0,0,0])
        >>> ii = argsort(a)
        >>> print(a[ii])
        [0 1 3 4 5 5 -- --]

        >>> ii = argsort(a, fill_value=1)
        >>> print(a[ii])
        [0 -- -- 1 3 4 5 5]

        >>> a = [0,4,6,2,1,5,3,5]
        >>> ii = argsort(a)
        >>> b = [ a[i] for i in ii ]
        >>> print(b)
        [0, 1, 2, 3, 4, 5, 5, 6]

        >>> a = [0,4,6,2,1,5,3,5]
        >>> ii = argsort(a, reverse=True)
        >>> b = [ a[i] for i in ii ]
        >>> print(b)
        [6, 5, 5, 4, 3, 2, 1, 0]


        Notes
        -----
        argsort for iterables from
        http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python
        http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    if isinstance(a, np.ma.MaskedArray):
        return np.ma.argsort(a, *args, **kwargs)
    elif isinstance(a, np.ndarray):
        return np.argsort(a, *args, **kwargs)
    else:
        return iargsort(a, *args, **kwargs)


# same numpy.argsort but for python iterables
def iargsort(seq, *args, **kwargs):
    if 'key' in kwargs:
        raise KeyError('keyword key cannot be given to argsort.')
    return sorted(range(len(seq)), *args, key=seq.__getitem__, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
