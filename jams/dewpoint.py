#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def dewpoint(pres, Celsius=False):
    """
        Calculates the dew point [K] from ambient humidity [Pa].


        Definition
        ----------
        def dewpoint(pres):


        Input
        -----
        pres       actual vapour pressure [Pa]


        Optional Input
        --------------
        Celsius    If True, return degree Celsius instead of Kelvin


        Output
        ------
        Dew point in Kelvin [K]


        Restrictions
        ------------
        None.


        References
        ---------
        Uses Bucks original vapour pressure formulation based on Tetens formula
        Buck, A. L., New equations for computing vapour pressure and enhancement factor,
                     J. Appl. Meteorol., 20, 1527-1532, 1981.

        Examples
        --------
        >>> Ta = 20. + 273.15
        >>> from jams import esat
        >>> es = esat(Ta, formula='Buck_original')
        >>> from autostring import astr
        >>> print(astr(dewpoint(es),3,pp=True))
        293.150

        >>> print(astr(dewpoint(es, Celsius=True),3,pp=True))
        20.000


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2013 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jan 2012
        Modified, MC, Feb 2013 - ported to Python 3
    """
    pw  = 611.21
    c1  = 240.97
    c2  = 17.502
    w   = np.log(pres/pw)
    out = w*c1/(c2-w)
    if Celsius:
        return out
    else:
        import jams.const as const
        return out+const.T0


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
