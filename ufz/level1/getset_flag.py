#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['get_flag', 'set_flag']

# --------------------------------------------------------------------

def get_flag(flags, n):
    """
        Get the flags at position n from CHS data flag vector.


        Definition
        ----------
        def get_flag(flags, n):


        Input
        -----
        flags   ND-array of integers
        n       position to extract

        
        Output
        ------
        ND-array of integers with flag at position n of flags
        returns -1 if input integer is less than 10**(n+1)


        Examples
        --------
        >>> flags = np.array([3, 30, 31, 300, 301, 3001, 3201, 312121212])
        >>> print(get_flag(flags, 0))
        [3 3 3 3 3 3 3 3]

        >>> print(get_flag(flags, 1))
        [-1  0  1  0  0  0  2  1]

        >>> print(get_flag(flags, 2))
        [-1 -1 -1  0  1  0  0  2]

        >>> print(get_flag(flags, 3))
        [-1 -1 -1 -1 -1  1  1  1]

        >>> print(get_flag(flags, 4))
        [-1 -1 -1 -1 -1 -1 -1  2]


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
    """
    ilog10 = np.log10(flags).astype(np.int)         # how many digits
    out = np.ones(flags.shape, dtype=np.int)*(-1)   # if not enough digits -> -1
    ii = np.where((ilog10-n) >= 0)[0]               # otherwise extract n-th digit
    if ii.size > 0:
        out[ii] = (flags[ii] // 10**(ilog10[ii]-n)) % 10

    return out

# --------------------------------------------------------------------

def set_flag(flags, n, iflag, ii=None):
    """
        Set the flags at position n to iflag at indeces ii of CHS data flag vector.


        Definition
        ----------
        def set_flag(flags, n, iflag, ii=None):


        Input
        -----
        flags   CHS data flag vector, ND-array of integers
        n       position in flag to set, missing positions will be cretaed and filled with 0
        iflag   set flag to this value


        Optional Input
        --------------
        ii      indices at which to set flag (default: None, i.e. each entry)

        
        Output
        ------
        ND-array of integers with flag at position n of flags set to iflag


        Examples
        --------
        >>> flags = np.array([3, 30, 301, 3101, 312121212])
        >>> print(set_flag(flags, 1, 2, [0,1,2]))
        [       32        32       321      3101 312121212]

        >>> print(set_flag(flags, 3, 2, [0,1,2,3]))
        [     3002      3002      3012      3102 312121212]

        >>> print(set_flag(flags, 1, 2))
        [       32        32       321      3201 322121212]


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
    """
    fflags = flags.copy()
    # extend flag vector if needed
    ilog10 = np.log10(fflags).astype(np.int)  # are there enough flag positions (i)
    jj = np.where(ilog10 < n)[0]
    if jj.size > 0:                          # increase number (filled with 0)
        fflags[jj] *= 10**(n-ilog10[jj])
        ilog10 = np.log10(fflags).astype(np.int)
        
    # get entries in flag vector
    isflags = get_flag(fflags, n)

    # set entries in flag vector
    if ii is None:                                  # set all to iflag
        fflags += 10**(ilog10-n) * (iflag-isflags)
    else:
        if np.size(ii) > 0:                         # set ii to iflag
            fflags[ii] += 10**(ilog10[ii]-n) * (iflag-isflags[ii])

    return fflags

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
