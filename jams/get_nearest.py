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

        Copyright 2014 Arndt Piayda


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
