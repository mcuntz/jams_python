#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.const import eps

def around(num, powten, ceil=False, floor=False):
    """
        Round to the passed power of ten.


        Definition
        ----------
        def around(num, powten=None, ceil=False, floor=False):


        Input
        -----
        num        ND array of floats


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
        >>> from jams import astr
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
        This file is part of the JAMS Python package.

        Copyright (c) 2011-2016 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2011
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
                  MC, Nov 2016 - const.tiny -> const.eps
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
