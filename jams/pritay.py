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
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2016 Arndt Piayda

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
