#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def cellarea(lat, lon, globe=False):
    """
        Calculates the area of Earth grid cells in metre square


        Definition
        ----------
        def cellarea(lon, lat, globe=False):


        Input
        -----
        lat          latitudes in degrees N
        lon          longitudes in degrees E


        Options
        -------
        globe        assumes that latitudes span the globe
                     i.e. they are bounded by 90/-90 degrees

        Output
        ------
        array with area in m^2


        Restrictions
        ------------
        No irregular spacings.


        Examples
        --------
        # Gaussian latitudes
        >>> import numpy as np
        >>> lat = np.array([ 12.98898858, 9.27785325, 5.56671363])
        >>> lon = np.array([ 0., 3.75, 7.5])
        >>> from autostring import astr
        >>> print(astr(cellarea(lat,lon)[0,:],3,pp=True))
        ['1.676e+11' '1.676e+11' '1.676e+11']
        >>> print(astr(cellarea(lat,lon)[1,:],3,pp=True))
        ['1.698e+11' '1.698e+11' '1.698e+11']
        >>> print(astr(cellarea(lat,lon)[2,:],3,pp=True))
        ['1.712e+11' '1.712e+11' '1.712e+11']


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

        Copyright 2009-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    nlat = len(lat)
    assert nlat >= 2, 'at least 2 latitudes must be given'
    nlon = len(lon)
    assert nlon >= 2, 'at least 2 longitudes must be given'
    #
    # cell sizes in degrees # still + or -
    dlat = lat-np.roll(lat,1)
    if globe:
        if lat[0] > 0: # descending lats
            l0 = 90.
        else:
            l0 = -90. # ascending lats
        dlat[0]  = (lat[0]-l0)+0.5*(lat[1]-lat[0])
        dlat[-1] = (-l0-lat[-1])-0.5*(lat[-2]-lat[-1])
    else:
        dlat[0]  = lat[1]-lat[0]
        dlat[-1] = lat[-1]-lat[-2]
    loni = lon
    # check if meridian in lon range -> shift to -180,180
    if np.any(np.abs(np.diff(loni)) > 360./nlon):
        loni = np.where(loni > 180., loni-360., loni)
    # check if -180,180 longitude in lon range -> shift to 0,360
    if np.any(np.abs(np.diff(loni)) > 360./nlon):
        loni = np.where(loni < 0, loni+360., loni)
    dlon = np.abs(loni-np.roll(loni,1))
    dlon[0] = np.abs(loni[1]-loni[0])
    #
    # Northern latitude of grid cell edges
    nlat = np.size(lat)
    n_lat = lat[0] + dlat[0]/2. + np.cumsum(dlat)
    #
    # Area of grid cells in m^2 with lat/lon in degree
    ae  = 6.371e6    # radius of Earth
    d2r = np.pi/180. # degree to radian
    #
    area = np.empty([nlat,nlon])
    dlat = np.abs(dlat[:])
    for i in range(nlon):
        area[:,i] = 2.*(ae**2) * dlon[i]*d2r * np.sin(0.5*dlat*d2r) * np.cos(lat*d2r)
    #
    return area

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
