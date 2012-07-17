#!/usr/bin/env python
import numpy as np

def calcvpd(rel_hum,Tair):
    '''
    PURPOSE:
    --------
    Calculates Vapour Pressure Deficit (VPD) according to Buck, A. L.
    (1981): New equations for computing vapour pressure and 
    enhancement factor, Journal of Applied Meteorology,20,1507-1520.
    mentioned in Campbell and Norman derived from Tetens formula.
    
    e_s = 6.1121 * e**(17.502 * T / (T + 240.97))
    vpd = e_s - (rel_hum / 100.)
    
    DEFINITION:
    -----------
    def calcvpd(rel_hum,Tair): 
    
    INPUT:
    ------
    Relatvie Air Humidity (0 to 1)
    Air Temperature in Kelvin

    OUTPUT:
    ------ 
    Vapour Pressure Deficit in Pa


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2011 Toni Rau

        
    WRITTEN BY TR MAY 2011
    '''
    print 'CALCVPD: obsolete function. Use esat: calcvpd(rh,T) = (1-rh)*esat(t)'
    print '         This function is precisely: calcvpd(rh,T) = (1-rh)*esat(T, formula="Buck_original")'
    sat_pres = 6.1121*np.exp(17.502 * (Tair - 273.15) / (Tair - 273.15 + 240.97))
    vap_def = sat_pres*(1.- rel_hum) * 100.
    return vap_def
