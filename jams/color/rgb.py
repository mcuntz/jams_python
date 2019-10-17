#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Interpolate between colours in RGB space; make continuous colour maps.


    Definition
    ----------
    Calculates colour between two given colors in RGB space
        def rgb_blend(col1, col2, fraction=0.5):
    n interpolated colours between two colours in RGB space
        def rgb_range(col1, col2, n=255, cmap=None, pow=1):
    n interpolated colours in RGB space between several colours changing at certain fractions
        def rgb_gradient(colours, fractions, n=255, cmap=None):


    Input
    -----
    rgb_blend
        col1      1st rgb colour tuple
        col2      2nd rgb colour tuple

    rgb_range
        col1      1st rgb colour tuple
        col2      2nd rgb colour tuple

    rgb_gradient
        colours   Nx3 array like of colour tuples
        fractions N array like of fractions for blending between colours



    Optional Input
    --------------
    rgb_blend
        fraction   fraction between 0=col1 and 1=col2; default: 0.5

    rgb_range
        n          number of interpolated colours with first colour=col1 and last colour=col2; default: 255
        cmap       if given, register colour map under that name
        pow        1 (default) is linear interpolation
                   >1 remains longer near col1, i.e. higher values are detailed
                   <1 remains longer near col2, i.e. lower values are detailed

    rgb_gradient
        n          number of interpolated colours with first colour=first colour in colours
                   and last colour=last colour in colours; default: 255
        cmap       if given, register colour map under that name


    Output
    ------
    rgb_blend      rgb tuple
    rgb_range      list of rgb tuples
    rgb_gradient   list of rgb tuples


    Examples
    --------
    >>> r = (1.0,0.0,0.0)
    >>> b = (0.0,0.0,1.0)
    >>> print(rgb_blend(r,b,0.0), rgb_blend(r,b,0.5), rgb_blend(r,b,1.0))
    (1.0, 0.0, 0.0) (0.5, 0.0, 0.5) (0.0, 0.0, 1.0)
    >>> print(rgb_range(r,b,3))
    [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]
    >>> print(rgb_range(r,b,3,pow=2))
    [(1.0, 0.0, 0.0), (0.75, 0.0, 0.25), (0.0, 0.0, 1.0)]
    >>> print(rgb_gradient([r,b],[0.0,1.0],3))
    [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]
    >>> print(rgb_gradient([r,r,b,b],[0.0,0.25,0.75,1.0],5, cmap='MyGradient'))
    [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)]


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2013 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Apr 2013
"""
import numpy as np
import jams.const as const

__all__ = ['rgb_blend', 'rgb_gradient', 'rgb_range']

# ---------------------------------------------------------------------

# http://stackoverflow.com/questions/25007/conditional-formatting-percentage-to-color-conversion
def rgb_blend(col1, col2, fraction=0.5):
    """
        Calculates colour between two colors.

        Definition
        ----------
        def rgb_blend(col1, col2, fraction=0.5):

        Input
        -----
            col1      1st rgb colour tuple
            col2      2nd rgb colour tuple

        Optional Input
        --------------
            fraction  fraction between 0=col1 and 1=col2; default: 0.5

        Output
        ------
        rgb_blend     rgb tuple

        Examples
        --------
        >>> r = (1.0,0.0,0.0)
        >>> b = (0.0,0.0,1.0)
        >>> print(rgb_blend(r,b,0.0), rgb_blend(r,b,0.5), rgb_blend(r,b,1.0))
        (1.0, 0.0, 0.0) (0.5, 0.0, 0.5) (0.0, 0.0, 1.0)

        History
        -------
        Written,  MC, Apr 2013
    """
    return tuple([v1 + (v2-v1)*fraction for (v1, v2) in zip(col1, col2)])

# ---------------------------------------------------------------------

def rgb_range(col1, col2, n=255, cmap=None, pow=1):
    """
        n interpolated colours between two colours.

        Definition
        ----------
        def rgb_range(col1, col2, n=255, cmap=None, pow=1):

        Input
        -----
            col1      1st rgb colour tuple
            col2      2nd rgb colour tuple

        Optional Input
        --------------
            n         number of interpolated colours with first colour=col1 and last colour=col2; default: 255
            cmap      if given, register colour map under that name
            pow       1 (default) is linear interpolation
                      >1 remains longer near col1, i.e. higher values are detailed
                      <1 remains longer near col2, i.e. lower values are detailed

        Output
        ------
        rgb_range     list of rgb tuples

        Examples
        --------
        >>> r = (1.0,0.0,0.0)
        >>> b = (0.0,0.0,1.0)
        >>> print(rgb_range(r,b,3))
        [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]
        >>> print(rgb_range(r,b,3,pow=2))
        [(1.0, 0.0, 0.0), (0.75, 0.0, 0.25), (0.0, 0.0, 1.0)]

        History
        -------
        Written,  MC, Apr 2013
    """
    colr = [rgb_blend(col1, col2, (np.float(i)/np.float(n-1))**pow) for i in range(n)]
    if cmap is not None:
        import matplotlib.colors as col
        import matplotlib.cm as cm
        iscmap = col.ListedColormap(colr, name=cmap, N=n)
        cm.register_cmap(name=cmap,cmap=iscmap)
    return colr

# ---------------------------------------------------------------------

def rgb_gradient(colours, fractions, n=255, cmap=None):
    """
        n interpolated colours between several colours changing at certain fractions.

        Definition
        ----------
        def rgb_gradient(colours, fractions, n=255, cmap=None):

        Input
        -----
            colours   Nx3 array like of colour tuples
            fractions N array like of fractions for blending between colours

        Optional Input
        --------------
            n         number of interpolated colours with first colour=first colour in colours
                      and last colour=last colour in colours; default: 255
            cmap      if given, register colour map under that name

        Output
        ------
        rgb_gradient  list of rgb tuples

        Examples
        --------
        >>> r = (1.0,0.0,0.0)
        >>> b = (0.0,0.0,1.0)
        >>> print(rgb_gradient([r,b],[0.0,1.0],3))
        [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]
        >>> print(rgb_gradient([r,r,b,b],[0.0,0.25,0.75,1.0],5))
        [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)]
        >>> print(rgb_gradient([r,r,b,b],[0.0,0.25,0.75,1.0],5, cmap='MyGradient'))
        [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)]

        History
        -------
        Written,  MC, Apr 2013
    """
    cols  = np.array(colours)
    fracs = np.array(fractions)
    if cols.shape[0] != fracs.size: raise Error('colours.shape[0] != fractions.size')
    colors = []
    for i in range(n):
        frac = np.float(i)/np.float(n-1)
        if frac <= fracs[0]:                         # before first fraction
            colors += [tuple(cols[0,:])]
        elif frac >= fracs[-1]:                      # after last fraction
            colors += [tuple(cols[-1,:])]
        else:
            ii = np.where(fracs >= frac)[0][0]
            if np.abs(fracs[ii]-frac) > const.eps:   # exactly a fraction
                colors += [rgb_blend(cols[ii-1,:], cols[ii,:], frac)]
            else:
                colors += [tuple(cols[ii,:])]

    if cmap is not None:
        import matplotlib.colors as col
        import matplotlib.cm as cm
        iscmap = col.ListedColormap(colors, name=cmap, N=n)
        cm.register_cmap(name=cmap,cmap=iscmap)

    return colors

# ---------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # r = (1.0,0.0,0.0)
    # b = (0.0,0.0,1.0)
    # print(rgb_blend(r,b,0), rgb_blend(r,b,0.1), rgb_blend(r,b,0.2), rgb_blend(r,b,0.3),
    #       rgb_blend(r,b,0.4), rgb_blend(r,b,0.5), rgb_blend(r,b,0.6), rgb_blend(r,b,0.7),
    #       rgb_blend(r,b,0.8), rgb_blend(r,b,0.9), rgb_blend(r,b,1))
    # print(rgb_range(r,b,11))
    # print(rgb_gradient([r,b],[0,1],11))
    # print(rgb_gradient([r,r,b,b],[0,0.4,0.6,1],11))
    # print(rgb_gradient([r,r,b,b],[0,0.3,0.7,1],11))

    # print(rgb_range(r,b,11))
    # print(rgb_range(r,b,11, pow=3))
    # print(rgb_range(r,b,11, pow=0.3))
