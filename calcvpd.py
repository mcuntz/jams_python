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

    WRITTEN BY TR MAY 2011
    '''
    print 'CALCVPD: obsolete function. Use esat: calcvpd(rh,T) = (1-rh)*esat(t)'
    print '         This function is precisely: calcvpd(rh,T) = (1-rh)*esat(T, formula="Buck_original")'
    sat_pres = 6.1121*np.exp(17.502 * (Tair - 273.15) / (Tair - 273.15 + 240.97))
    vap_def = sat_pres*(1.- rel_hum) * 100
    return vap_def
