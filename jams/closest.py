#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def closest(vec, num, value=False):
    """
        Index in array at which the entry is closest to a given number.


        Definition
        ----------
        def closest(vec, num, value=False):


        Input
        -----
        vec        ND-array
        num        number to compare


        Optional Input
        --------------
        value      If true, give closest array element instead of index


        Output
        ------
        Flat index of element closest to given number num.
        Use np.unravel_index to get index tuple.


        Restrictions
        ------------
        None.


        Examples
        --------
        >>> vec = np.arange(100)/99.*5.
        >>> print(closest(vec, 3.125))
        62
        >>> from autostring import astr
        >>> print(astr(closest(vec, 3.125, value=True),3,pp=True))
        3.131

        >>> vec = np.arange(100).reshape((10,10))/99.*5.
        >>> print(astr(closest(vec, 3.125, value=True),3,pp=True))
        3.131
        >>> print(closest(vec, 3.125))
        62
        >>> print(np.unravel_index(closest(vec, 3.125), vec.shape))
        (6, 2)
        >>> print(astr(vec[np.unravel_index(closest(vec, 3.125), vec.shape)],3,pp=True))
        3.131


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Jul 2015 - 
    """
    out = np.ma.argmin(np.ma.abs(np.ma.array(vec)-num))
    if value:
      return vec.flat[out]
    else:
      return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
