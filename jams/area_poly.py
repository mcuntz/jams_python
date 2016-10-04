#!/usr/bin/env python
from __future__ import print_function
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
        Written,  MC, Nov 2012 - stackoverflow.com
        Modified, MC, Feb 2013 - ported to Python 3
    """

    # Could include some checks here
    p = list(zip(x,y))
    return 0.5 * np.abs(np.sum(x0*y1 - x1*y0 for ((x0, y0), (x1, y1)) in segments(p)))


def segments(p):
    return list(zip(p, p[1:] + [p[0]]))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # x = np.array([1.0, 2.0, 2.0, 1.0])
    # y = np.array([1.0, 1.0, 2.0, 2.0])
    # print(area_poly(x,y))

