#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import ufz.const as const

def around(num, powten, ceil=False, floor=False):
    """
        Round to the passed power of ten.


        Definition
        ----------
        def around(num, powten=None, ceil=False, floor=False):


        Input
        -----
        num        number array


        Optional Input
        --------------
        powten     Power of ten array
                   If missing, simple round (ceil, floor) is taken.
        ceil       ceil instead of round to the nearest power of ten
        floor      floor instead of round to the nearest power of ten


        Output
        ------
        Rounded values.


        Restrictions
        ------------
        Powten is exactly opposite of decimal keyword of numpy.around.

        From numpy.around documentation:
        'For values exactly halfway between rounded decimal values,
        Numpy rounds to the nearest even value. Thus 1.5 and 2.5 round to 2.0,
        -0.5 and 0.5 round to 0.0, etc. Results may also be surprising due to the
        inexact representation of decimal fractions in the IEEE floating point
        standard and errors introduced when scaling by powers of ten.'


        Examples
        --------
        >>> from autostring import astr
        >>> print(astr(around(np.array([3.5967,345.5967]), -3),3,pp=True))
        ['  3.597' '345.597']

        >>> print(astr(around(np.array([1994344,345.5967]), [3,-3]),3,pp=True))
        ['1.994e+06' '3.456e+02']

        >>> print(astr(around(np.array([1994344,345.5967]), [3,-3], ceil=True),3,pp=True))
        ['1.995e+06' '3.456e+02']

        >>> print(astr(around(np.array([1994344,345.5967]), [3,-3], floor=True),3,pp=True))
        ['1.994e+06' '3.456e+02']

        >>> print(astr(around(np.array([3.5967,345.5967]), 3),3,pp=True))
        ['0.000' '0.000']

        >>> print(astr(around(np.array([3.5967,345.5967]), 3, ceil=True),3,pp=True))
        ['1000.000' '1000.000']


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

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2011
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    #
    # Check input
    assert (ceil+floor) < 2, 'ceil and floor keywords given.'
    if (powten == None):
        ipowten = 0
    else:
        ipowten = np.array(powten)
    nnum = np.size(num)
    npowten = np.size(ipowten)
    assert not ((npowten != nnum) and (npowten != 1)), 'powten must be scalar or have array size of input numbers.'
    #
    # Shift decimal point
    # Does not work, too imprecise: out = num * np.exp(-ipowten*np.log(10.))
    out = num * 10.**(-ipowten)
    # Round/ceil/floor
    #eps = np.MachAr().eps
    #eps = np.finfo(np.float).eps
    eps = const.tiny
    if (ceil):
        out = np.ceil(out-10.*eps)
    elif (floor):
        out = np.floor(out+10.*eps)
    else:
       out = np.around(out)
    # Shift back decimal point
    # Does not work, too imprecise: out *= np.exp(ipowten*np.log(10.))
    out *= 10.**ipowten

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
