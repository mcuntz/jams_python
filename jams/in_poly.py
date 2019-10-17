#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def in_poly(P, coord_x, coord_y):
    """
        Determines whether a 2D point falls within a polygon, on a vertex or
        edge, or is located outside of a polygon. The polygon can be convex or
        not.


        Definition
        ----------
        def in_poly(P, coord_x, coord_y):


        Input
        -----
        P         2D list or np.array,
                  x and y coordinates of the point in question in the form [x,y]
        coord_x   np.array, x coordinates of the polygon
        coord_y   np.array, y coordinates of the polygon


        Output
        ------
        integer, 1 = point inside polygon
                 0 = point on vertex/edge
                -1 = point outside polygon


        Restrictions
        ------------
        The method is only applicable for 2D polygons and points.


        References
        ----------
        This routine is re-coded from the JAMS Fortran package.
        Copyright: Juliane Mai, 2012.
        
        The original version of the source code (pnpoly) was implemented by
        W. Randolph Franklin. It had been, however, assigning insufficiently
        vertex/edge points.


        Example
        --------
        >>> coord_x = np.array([2.,2.,5.,7.,5.])
        >>> coord_y = np.array([1.,4.,6.,3.,1.])

        # point inside polygon
        >>> P = [4.,3.]
        >>> print(in_poly(P, coord_x, coord_y))
        1

        # point outside polygon
        >>> P = [8.,6.]
        >>> print(in_poly(P, coord_x, coord_y))
        -1

        # point on edge/vertex of polygon
        >>> P = [2.,2.]
        >>> print(in_poly(P, coord_x, coord_y))
        0


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  AP, Nov 2012
        Modified, MC, Nov 2012 - documentation change, return 0 immediately
                  MC, Feb 2013 - ported to Python 3
                  MC, Oct 2013 - inpoly
                  MC, Apr 2014 - assert
    """

    # ironing :-)
    coord_x, coord_y = coord_x.flatten(), coord_y.flatten()

    # test input sizes
    assert np.size(coord_x) == np.size(coord_y), 'in_poly: coord_x and coord_y must have same size.'

    n = coord_x.size
    # result is outside as long as no other test works
    erg = -1

    # relative coordinates
    X = coord_x - P[0]
    Y = coord_y - P[1]

    # Edge test
    if np.any((X==0.) & (Y==0.)): return 0

    for i in range(n):
        # vertical Vertex test
        j = (i+1) % n
        if (coord_x[i] == coord_x[j]) and (coord_x[i] == P[0]):
            ly = (P[1]-coord_y[j]) / (coord_y[i]-coord_y[j])
            if (ly >= 0.) and (ly <= 1.): return 0

        # horizontal Vertex test
        if (coord_y[i] == coord_y[j]) and (coord_y[i] == P[1]):
            lx = (P[0]-coord_x[j]) / (coord_x[i]-coord_x[j])
            if (lx >= 0.) and (lx <= 1.): return 0

        # Inside test
        MX = X[i] >= 0.
        NX = X[j] >= 0.
        MY = Y[i] >= 0.
        NY = Y[j] >= 0.

        test1 = not((MY or NY) and (MX or NX)) or (MX and NX)
        test2 = not(MY and NY and (MX or NX) and not(MX and NX))

        if (not test1):
            if test2:
                tt = (Y[i]*X[j] - X[i]*Y[j]) / (X[j] - X[i])
                if tt == 0.:
                    return 0
                elif tt > 0.:
                    erg = -erg
            else:
                erg = -erg
    return erg


def inpoly(*args, **kwargs):
    """
        wrapper for in_poly
        def in_poly(P, coord_x, coord_y):


        Example
        --------
        >>> coord_x = np.array([2.,2.,5.,7.,5.])
        >>> coord_y = np.array([1.,4.,6.,3.,1.])

        # point inside polygon
        >>> P = [4.,3.]
        >>> print(inpoly(P, coord_x, coord_y))
        1

        # point outside polygon
        >>> P = [8.,6.]
        >>> print(inpoly(P, coord_x, coord_y))
        -1

        # point on edge/vertex of polygon
        >>> P = [2.,2.]
        >>> print(inpoly(P, coord_x, coord_y))
        0
    """
    return in_poly(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


