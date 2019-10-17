#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy.spatial.distance import cdist

def get_nearest(xy, xyz):
    """
        Returns the z value for each point in an xy coordinate array 
        which is equal to the z value at the nearest point in a given
        xyz field.

        Definition
        ----------
        get_nearest(xy, xyz):


        Input
        -----
        xy        2D np.array (n,2),
                  x and y coordinates of the points where a z value is desired
        xyz       2D np.array (m,3), x, y and z vaules of the field


        Output
        ------
        z         1D np.array (n,), z values at all points of xy


        Restrictions
        ------------
        Does not work with NAN values or masked arrays. You need to provide only
        valid data.


        Example
        --------
        # make some data
        >>> np.random.seed(seed=1)
        >>> xy = np.random.random(6).reshape(-1,2)
        >>> np.random.seed(seed=2)
        >>> xyz = np.random.random(150).reshape(-1,3)

        # get z values at the points xy
        >>> z = get_nearest(xy, xyz)
        >>> print(np.round(z,2))
        [0.25 0.7  0.28]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014 Arndt Piayda

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
        Written,  AP, May 2014
        Modified, ST, Jun 2014 - minor change to documentation
    """

    dist   = cdist(xy, xyz[:,:2], metric='euclidean')
    dist   = np.ma.array(dist, mask=dist==0.)
    indmin = np.ma.argmin(dist, axis=1)
    
    return xyz[indmin,2]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
