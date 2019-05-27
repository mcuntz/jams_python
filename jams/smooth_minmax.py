#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

def smax(x, y, eta=1. - 1.e-5):
    """
        Calculate smooth maximum of two numbers.


        Definition
        ----------
        def around(x, y, eta=1.):


        Input
        -----
        x        ND array of floats
        y        ND array of floats


        Optional Input
        --------------
        eta        Shape of smoothing function (default: 1._dp)


        Output
        ------
        Smoothed maximum


        Restrictions
        ------------
        If eta is set to 1, then this function becomes the standard maximum function.


        Examples
        --------
        >>> from jams.autostring import astr
        >>> print(astr(smax(3.5967, 3.6),3,pp=True))
            3.610


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

        Copyright 2016-2016 Stephan Thober


        History
        -------
        Written,  ST, Mar 2016
        Modified, 
    """
    #
    if (eta > 1.):
        raise ValueError('eta is {:3.4f}, but cannot be larger than 1.'.format(eta))
    #
    out = (x + y + ((x + y)**2. - 4 * eta * x * y)**0.5) / (2. * eta)
    return out


def smin(x, y, eta=1. - 1.e-5):
    """
        Calculate smooth minimum of two numbers.


        Definition
        ----------
        def around(x, y, eta=1. - 1.e-5):


        Input
        -----
        x        ND array of floats
        y        ND array of floats


        Optional Input
        --------------
        eta        Shape of smoothing function (default: 1._dp)


        Output
        ------
        Smoothed minimum


        Restrictions
        ------------
        If eta is set to 1, then this function becomes the standard minimum function.


        Examples
        --------
        >>> from jams.autostring import astr
        >>> print(astr(smin(3.5967, 3.6),3,pp=True))
            3.587


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

        Copyright 2016-2016 Stephan Thober


        History
        -------
        Written,  ST, Mar 2016
        Modified, 
    """
    #
    if (eta > 1.):
        raise ValueError('eta is {:3.4f}, but cannot be larger than 1.'.format(eta))
    #
    out = (x + y - ((x + y)**2. - 4 * eta * x * y)**0.5) / (2. * eta)
    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
