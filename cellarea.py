#!/usr/bin/env python
import numpy as np

def cellarea(lat, lon, globe=False):
    """
        Calculates the area of grid cells in metre square

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
    
        Example
        -------
        # Gaussian latitudes
        >>> import numpy as np
        >>> lat = np.array([ 12.98898858, 9.27785325, 5.56671363])
        >>> lon = np.array([ 0., 3.75, 7.5])
        >>> cellarea(lat,lon)
        array([[  1.67639084e+11,   1.67639084e+11,   1.67639084e+11],
               [  1.69790428e+11,   1.69790428e+11,   1.69790428e+11],
               [  1.71229889e+11,   1.71229889e+11,   1.71229889e+11]])

        History
        -------
        Written, MC, Jul. 2009
    """
    nlat = len(lat)
    if nlat < 2:
        print 'CELLAREA: at least 2 latitudes must be given'
        return None
    nlon = len(lon)
    if nlon < 2:
        print 'CELLAREA: at least 2 longitudes must be given'
        return None
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
    if np.any(abs(loni-np.roll(loni,1)) > 360./nlon):
        loni = np.where(loni > 180., loni-360., loni)
    # check if -180,180 longitude in lon range -> shift to 0,360
    if np.any(abs(loni-np.roll(loni,1)) > 360./nlon):
        loni = np.where(loni < 0, loni+360., loni)
    dlon = abs(loni-np.roll(loni,1))
    dlon[0] = abs(loni[1]-loni[0])
    #
    # Northern latitude of grid cell edges
    nlat = len(lat)
    n_lat = lat[0] + dlat[0]/2. + np.cumsum(dlat)
    # 
    # Area fo grid cells in m^2 with lat/lon in degree
    ae  = 6.371e6    # radius of Earth
    d2r = np.pi/180. # degree to radian
    #
    area = np.empty([nlat,nlon])
    dlat = np.abs(dlat[:])
    for i in range(nlon):
        area[:,i] = (2.*(ae**2) * dlon[i]*d2r 
                     * np.sin(0.5*dlat*d2r) * np.cos(lat*d2r))
    #
    return area

if __name__ == '__main__':
    import doctest
    doctest.testmod()
