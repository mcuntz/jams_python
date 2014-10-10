#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.date2dec import date2dec

def t2sap(date, data, swd=None, undef=-9999.):
    """
        Conversion of temperature difference measured with sap flow sensors
        (Granier type) in (mV) to sap flux density (cm^3 cm^-2 h^-1) (Granier, 1987).

        In addition, the correction according to Clearwater (1999) is possible.


        Definition
        ----------
        def t2sap(date, data, swd=None, undef=-9999.):


        Input
        -----
        date    1D array (n,) of ascii times in format DD.MM.YYYY hh:mm:ss
        data    ND array (n,...) of the raw data in mV. First dimension is time


        Options
        --------------
        undef   values are excluded from calculations (default: -9999)
        swd     if None:     sapflux calculation according to original Granier
                             calibration function.
                if not None: sapwood depth in cm for Clearwater correction
                             T = (T - b*Tmax)/a
                             with a = active sapwood   = swd/2
                                  b = inactive sapwood = 1-a


        Output
        ------
        SFD     2D array (n,m) with sap flux density in [cm^3 cm^-2 h^-1] according
                to the original Granier calibration function:
                SFD = 0.0119*K**1.231*3600 [Granier, 1985]


        Restrictions
        ------------
        None


        References
        ----------
        Clearwater, M. J., Meinzer, F. C., Andrade, J. L., Goldstein, G., Holbrook, N. M.,
            Potential errors in measurement of nonuniform sap flow using heat dissipation probes,
            Tree Physiology 19, 681-687, 1999
        Granier, A.,
            Evaluation of transpiration in a Douglas-fir stand by means of sap flow measurements,
            Tree Physiology 3, 309-320, 1987


        Examples
        --------
        # normal sapflux conversion
        >>> data    = np.array([0.434, 0.433, 0.432, 0.431, 0.431, 0.432])
        >>> date  = np.array(['18.05.2013 08:00', '18.05.2013 08:10', '18.05.2013 08:20',
        ...                   '18.05.2013 08:30', '18.05.2013 08:40', '18.05.2013 08:50'])
        >>> SFD     = t2sap(date, data)
        >>> print(np.round(SFD,3))
        [ 0.     0.024  0.057  0.095  0.095  0.057]


        >>> # sapflux conversion including clearwater correction
        >>> data    = np.array([0.434, 0.433, 0.432, 0.431, 0.431, 0.432])
        >>> date  = np.array(['18.05.2013 08:00', '18.05.2013 08:10', '18.05.2013 08:20',
        ...                   '18.05.2013 08:30', '18.05.2013 08:40', '18.05.2013 08:50'])
        >>> SFD     = t2sap(date, data, swd=1.5)
        >>> print(np.round(SFD,3))
        [ 0.     0.035  0.082  0.135  0.135  0.082]


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
        along with The UFZ Python package.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Andreas Wiedemann


        History
        -------
        Written,  AW, Jun 2014
        Modified, MC, Jun 2014 - sap_app -> t2sap, incl. undef, first date then data
                               - squeeze 1D array, fill undef if not masked
    """

    isvec = False
    if data.ndim == 1:
        isvec = True
        data = data.reshape((-1, 1))
    # mask data
    isnotmask = True
    if type(data) == np.ma.core.MaskedArray:
        isnotmask = False
    data_masked = np.ma.array(data, mask=(data==undef))

    # julian day
    jd        = np.array(date2dec(ascii = date))
    # integer julian day 
    jd_int    = jd.astype(np.int)
    # unique days
    jd_uni    = np.unique(jd_int)
    # Tmax per day
    tmax_day  = np.ma.ones((jd_uni.size, data.shape[1]))*undef
    # Time of Tmax per day
    jdmax_day = np.ma.zeros((jd_uni.size, data.shape[1]), dtype=np.int)

    # Determine Tmax per day
    for count, i in enumerate(jd_uni):
        # where is given day
        ii = np.where(jd_int == i)[0]
        # index of Tmax of day
        jj = np.ma.argmax(data_masked[ii,:], axis=0)
        # time at Tmax of day
        jdmax_day[count,:] = jd[ii[jj]]
        # Tmax per day
        tmax_day[count,:]  = data_masked[ii[jj],:]

    # Tmax for every record
    tmax = np.empty(data.shape)
    for i in range(data.shape[1]):
        # time points at Tmax per day
        xx = jdmax_day[:,i]
        # Tmax per day
        yy = tmax_day[:,i]
        ii = ~yy.mask
        # Tmax per time point
        tmax[:,i] = np.interp(jd, xx[ii], yy[ii])

    if swd == None:
        Tsw = data_masked
    else:
        # sapwood depth in cm / 2cm for needle lenght = portion of sensor in sapwood
        a = 0.5*swd
        # portion of sensor in inactive xylem
        b = 1.0-a
        # define weighted mean of T in the sapwood (a) and T in the inactive xylem (b)
        Tsw  = (data_masked - (b*tmax)) / a
        
    # converts raw data [mV] into Sap Flux density [cm3*cm-2*h-1]
    SFD  = (0.0119* ((tmax - Tsw)/Tsw)**1.231) * 3600.

    # if it was not masked then the original undef will be undef
    if isnotmask:
        SFD  = SFD.filled(undef)
    
    # 1D in = 1D out
    if isvec:
        SFD = np.squeeze(SFD)

    return SFD


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
