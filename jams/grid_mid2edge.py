#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

__all__ = ['grid_mid2edge']

def grid_mid2edge(lon, lat):
    """
        2D arrays of longitude and latitude grid edges from grid midpoints.


        Definition
        ----------
        def grid_mid2edge(lon, lat):


        Input
        -----
        lon        1D (nlon) or 2D array (nlat,nlon) of longitude grid midpoints
        lat        1D (nlat) or 2D array (nlat,nlon) of latitude grid midpoints


        Optional Input
        --------------
        None


        Output
        ------
        lon2, lat2 = 2D arrays (nlat+1,nlon+1) of longitude and latitude grid edges.


        Restrictions
        ------------
        Longitudes and latitudes are interpolated as on a plane, not a sphere.
        However, rotated grids should be possible.


        Examples
        --------
        >>> from jams.autostring import astr
        >>> nlon = 36
        >>> nlat = 18

        Lon -180 - +180, Lat -90 - +90, i.e. S-N
        >>> lon = np.arange(nlon)*(360./nlon) - 175. # -180 - +180
        >>> lat = np.arange(nlat)*(180./nlat) - 85.  # S-N
        >>> lonh, lath = grid_mid2edge(lon, lat)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['-180' ' 180']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        ['-90' ' 90']

        Lon -180 - +180, Lat +90 - -90, i.e. N-S
        >>> lat = lat[::-1]                          # N-S
        >>> lonh, lath = grid_mid2edge(lon, lat)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['-180' ' 180']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        [' 90' '-90']

        Lon 0 - 360, Lat +90 - -90, i.e. N-S
        >>> lon += 180.                              # 0-360
        >>> lonh, lath = grid_mid2edge(lon, lat)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['  0' '360']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        [' 90' '-90']

        Rolled axis for different map projections, Lat +90 - -90, i.e. N-S
        >>> lon  = np.roll(lon, nlon//3)             # 245-235
        >>> lonh, lath = grid_mid2edge(lon, lat)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['240' '240']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        [' 90' '-90']
        >>> lon -= 180.                              # 65-55
        >>> lonh, lath = grid_mid2edge(lon, lat)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['60' '60']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        [' 90' '-90']

        2D input
        >>> lon = np.arange(nlon)*(360./nlon) - 175. # -180 - +180
        >>> lat = np.arange(nlat)*(180./nlat) - 85.  # S-N
        >>> lon2, lat2 = np.meshgrid(lon,lat)
        >>> lonh, lath = grid_mid2edge(lon2, lat2)
        >>> print(astr([lonh[0,0],lonh[-1,-1]]))
        ['-180' ' 180']
        >>> print(astr([lath[0,0],lath[-1,-1]]))
        ['-90' ' 90']


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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Oct 2014
    """
    #
    # Check input
    do2 = np.ndim(lon)
    da2 = np.ndim(lat)
    assert (do2 >=1) & (do2 <= 2), 'longitudes have to be either 1D or 2D arrays.'
    assert (da2 >=1) & (da2 <= 2), 'latitudes have to be either 1D or 2D arrays.'
    assert do2 == da2, 'longitudes and latitudes have to be either 1D or 2D arrays; no mixture possible.'
    #
    # Make 2D
    if do2 == 1:
        lon2, lat2 = np.meshgrid(lon,lat) # have sizes 
    else:
        lon2 = lon # no copy
        lat2 = lat
    #
    # N-S or S-N
    isns = True                                 # descending latitudes
    if np.min(np.diff(lat2)) > 0.: isns = False # ascending latitudes
    #
    # 0-360 0r -180-+180
    is360 = True          #    0 -  360
    if np.any(lon2 < 0.): # -180 - +180
        lon2  = np.copy(lon2) + 180. # do all calculations in 0-360
        is360 = False
    #
    # out arrays
    nlat = lon2.shape[0]
    nlon = lon2.shape[1]
    lonh = np.empty((nlat+1, nlon+1))
    lath = np.empty((nlat+1, nlon+1))
    # Edge points in interior = mean of 4 surrounding grid boxes
    lonh[1:-1,1:-1] = 0.25*(lon2[0:-1,0:-1] + # upper left
                            lon2[0:-1,1:] +   # upper right
                            lon2[1:,0:-1] +   # lower left
                            lon2[1:,1:])      # lower right
    lath[1:-1,1:-1] = 0.25*(lat2[0:-1,0:-1] +
                            lat2[0:-1,1:] +   # same as lon
                            lat2[1:,0:-1] +
                            lat2[1:,1:])
    # left column w/o corners = left - half distance
    lonh[1:-1,0]  = 0.75*lon2[0:-1,0] + 0.75*lon2[1:,0] - 0.25*lon2[0:-1,1] - 0.25*lon2[1:,1]
    lath[1:-1,0]  = 0.5*(lat2[0:-1,0] + lat2[1:,0])
    # right column w/o corners = right + half distance
    lonh[1:-1,-1] = 0.75*lon2[0:-1,-1] + 0.75*lon2[1:,-1] - 0.25*lon2[0:-1,-2] - 0.25*lon2[1:,-2]
    lath[1:-1,-1] = 0.5*(lat2[0:-1,-1] + lat2[1:,-1])
    # upper row w/o corners = up + half distance
    lonh[0,1:-1]  = 0.5*(lon2[0,0:-1] + lon2[0,1:])
    lath[0,1:-1]  = 0.75*lat2[0,0:-1] + 0.75*lat2[0,1:] - 0.25*lat2[1,0:-1] - 0.25*lat2[1,1:]
    # lower row w/o corners = low - half distance
    lonh[-1,1:-1]  = 0.5*(lon2[-1,0:-1] + lon2[-1,1:])
    lath[-1,1:-1]  = 0.75*lat2[-1,0:-1] + 0.75*lat2[-1,1:] - 0.25*lat2[-2,0:-1] - 0.25*lat2[-2,1:]
    # corners = midpoint plus or minus dist to opposite corner
    lonh[0,0]   = lon2[0,0]   - (lonh[1,1] - lon2[0,0])     # upper left
    lath[0,0]   = lat2[0,0]   - (lath[1,1] - lat2[0,0])
    lonh[0,-1]  = lon2[0,-1]  - (lonh[1,-2] - lon2[0,-1])   # upper right
    lath[0,-1]  = lat2[0,-1]  - (lath[1,-2] - lat2[0,-1])
    lonh[-1,0]  = lon2[-1,0]  - (lonh[-2,1] - lon2[-1,0])   # lower left
    lath[-1,0]  = lat2[-1,0]  - (lath[-2,1] - lat2[-1,0])
    lonh[-1,-1] = lon2[-1,-1] - (lonh[-2,-2] - lon2[-1,-1]) # lower right
    lath[-1,-1] = lat2[-1,-1] - (lath[-2,-2] - lat2[-1,-1])

    # lon can be > 360
    # this makes 360=0 and 180=-180: stay with 360 and 180
    #    lonh %= 360
    if np.any(lonh > 360.):
        lonh = np.where(lonh > 360., lonh % 360., lonh)
    
    # return to -180-+180
    if not is360:
        lonh -= 180.

    return lonh, lath


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from jams.autostring import astr
    # nlon = 36
    # nlat = 18
    # lon = np.arange(nlon)*(360./nlon) - 175. # -180 - +180
    # lat = np.arange(nlat)*(180./nlat) - 85.  # S-N
    # lonh, lath = grid_mid2edge(lon, lat)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['-180' ' 180']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    ['-90' ' 90']
    # lat = lat[::-1]                          # N-S
    # lonh, lath = grid_mid2edge(lon, lat)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['-180' ' 180']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    [' 90' '-90']
    # lon += 180.                              # 0-360
    # lonh, lath = grid_mid2edge(lon, lat)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['  0' '360']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    [' 90' '-90']
    # lon2, lat2 = np.meshgrid(lon,lat)
    # lonh, lath = grid_mid2edge(lon2, lat2)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['  0' '360']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    [' 90' '-90']
    # lon  = np.roll(lon, nlon//3)             # 245-235
    # lonh, lath = grid_mid2edge(lon, lat)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['240' '240']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    [' 90' '-90']
    # lon -= 180.                              # 65 - 55
    # lonh, lath = grid_mid2edge(lon, lat)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # #    ['60' '60']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # #    [' 90' '-90']
    # # 2D
    # lon = np.arange(nlon)*(360./nlon) - 175. # -180 - +180
    # lat = np.arange(nlat)*(180./nlat) - 85.  # S-N
    # lon2, lat2 = np.meshgrid(lon,lat)
    # lonh, lath = grid_mid2edge(lon2, lat2)
    # print(astr([lonh[0,0],lonh[-1,-1]]))
    # # ['-180' ' 180']
    # print(astr([lath[0,0],lath[-1,-1]]))
    # # ['-90' ' 90']
