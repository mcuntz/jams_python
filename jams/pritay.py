#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def pritay(T, Rg, elev, a=1.12):
    '''
        Daily reference evapotranspiration after Priestley & Taylor
        

        Definition
        ----------
        def pritay(T, Rg, elev, a=1.12):


        Input
        -----
        T            np.array, np.ma.array or pd.dataframe with daily mean
                     temperature [degC]
        Rg           np.array, np.ma.array or pd.dataframe with daily sum of 
                     global radiation [kJ m-2]
        elev         float elevation above sea level [m]
            

        Optional Input
        --------------
        a            int coefficient for conditions (wet or dry) [-] 
                     (default: wet: 1.12)


        Output
        ------
        ETref        np.array, np.ma.array or pd.dataframe with daily sum of
                     reference evapotranspiration [mm]


        Examples
        --------
        >>> import numpy as np
        >>> T    = np.array([13.4, 17.7])  # [degC]
        >>> Rg   = np.array([8143., 24296.]) # [kJ m-2]
        >>> elev = 89.5 # [m] above sea level
        >>> print(np.round(pritay(T, Rg, elev), 2))
        [2.24 7.3 ]


        Literature
        -----
        Priestley, C.H.B., Taylor, R.J., 1972. On the assessment of surface
        heat flux and evaporation using large-scale parameters. Monthly
        Weather Review 100, 81-92.


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
        Written,  AP, Okt 2016
    '''
    # slope of vapor pressure curve [kPa degC-1]
    dv = (4098.*(0.6108*np.exp(17.27*T/(T + 237.3))))/(T + 237.3)**2

    # atmospheric pressure [kPa]    
    pair = 101.3*((T+273.16 - 0.0065*elev)/(T+273.16))**5.26                              
    
    # psychrometric constant [kPa C-1]    
    gamma = 0.665*10**-3*pair

    # reference evapotranspiration [???]
    ETref = a*(dv/(dv + gamma))*((Rg/1000)/2.45)
    
    return ETref

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
