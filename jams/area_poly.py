#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def area_poly(x, y):
    """
        Determines the are of a polygon represented as a list of (x,y) vertex coordinates,
        implicitly wrapping around from the last vertex to the first.


        Definition
        ----------
        def area_poly(x, y):


        Input
        -----
        x   np.array, x coordinates of the polygon
        y   np.array, y coordinates of the polygon


        Output
        ------
        area


        References
        ----------
        Code from http://stackoverflow.com/questions/451426/how-do-i-calculate-the-surface-area-of-a-2d-polygon
        From the website:
          It is an application of Green's theorem (http://en.wikipedia.org/wiki/Green%27s_theorem#Area_Calculation)
          for the functions -y and x; exactly in the way a planimeter works.
          More specifically:
            Formula above = integral_perimeter(-y dx + x dy) = integral_area((-(-dy)/dy+dx/dx)dydyx = 2*Area


        Example
        --------
        >>> x = np.array([1.0, 2.0, 2.0, 1.0])
        >>> y = np.array([1.0, 1.0, 2.0, 2.0])
        >>> from autostring import astr
        >>> print(astr(area_poly(x,y),1,pp=True))
        1.0


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2019 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Nov 2012 - stackoverflow.com
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, May 2019 - np.sum needs iterable instead of generator in Python 3
    """

    # Could include some checks here
    p = list(zip(x,y))
    return 0.5 * np.abs(np.sum([ x0*y1 - x1*y0 for ((x0, y0), (x1, y1)) in segments(p) ]))


def segments(p):
    return list(zip(p, p[1:] + [p[0]]))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # x = np.array([1.0, 2.0, 2.0, 1.0])
    # y = np.array([1.0, 1.0, 2.0, 2.0])
    # print(area_poly(x,y))

