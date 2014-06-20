# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import ufz

def sap_app(data,period,swd=None):
        """
        Conversion of temperature difference measured with sap flow sensors 
        (Granier type) in (mV) to sap flux density (cm^3 cm^-2 h^-1) 
        (Granier, 1985). In addition, the correction according to
        Clearwater (1999) is possible.
        
        
        Definition
        ----------        
        def sap_app(data,period,swd=None):


        Input
        -----        
        data    2D array (n,m) or 1D array (n,) of the raw data in mV, arranged in: records
                in rows & sensors in columns.
        period  1D array of time corresponding to the records (n,)
        
        
        Options
        --------------   
        swd     if None: sapflux calculation according to the original Granier 
                calibration function
                if Clearwater Correction is needed: enter sapwood depth in cm (float)
                -> Tsw = (ΔT - b*Δtmax)/a with a = active sapwood, b = inactive sapwood
        
        
        Output
        ------
        SFD     2D array (n,m) with sap flux density in [cm^3 cm^-2 h^-1] according
                to the original Granier calibration function:
                SFD = 0.0119*K**1.231*3600 [Granier, 1985]
        
        
        Restrictions
        ------------
        Needs masked arrays, nan does not work. Values of -9999 are masked and not
        processed
        
        
        References
        ----------
        Clearwater, M. J., Meinzer, F. C., Andrade, J. L., Goldstein, G., 
        Holbrook, N. M., Potential errors in measurement of nonuniform sap flow
        using heat dissipation probes, Tree physiology (1999).
        Granier, A. Evaluation of transpiration in a Douglas-fir stand by means
        of sap flow measurements. Tree physiology (1987).
        
        
        Examples
        --------
        # normal sapflux conversion
        >>> data    = np.array([0.434,0.433,0.432,0.431,0.431,0.432])
        >>> period  = np.array(['18.05.2013 08:00', '18.05.2013 08:10', '18.05.2013 08:20', '18.05.2013 08:30', '18.05.2013 08:40', '18.05.2013 08:50'])
        >>> SFD     = sap_app(data,period)
        >>> print(np.round(SFD,3))
        [[ 0.   ]
         [ 0.024]
         [ 0.057]
         [ 0.095]
         [ 0.095]
         [ 0.057]]

            
        >>> # sapflux conversion including clearwater correction
        >>> data    = np.array([0.434,0.433,0.432,0.431,0.431,0.432])
        >>> period  = np.array(['18.05.2013 08:00', '18.05.2013 08:10', '18.05.2013 08:20', '18.05.2013 08:30', '18.05.2013 08:40', '18.05.2013 08:50'])
        >>> SFD     = sap_app(data,period,swd=1.5)
        >>> print(np.round(SFD,3))
        [[ 0.   ]
         [ 0.035]
         [ 0.082]
         [ 0.135]
         [ 0.135]
         [ 0.082]]

            
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

        Copyright 2014 Andreas Wiedemann
        
        
        History
        -------
        Written,  AW, Jun 2014
       
        """
    
    #    if error_sign  == None:
    #        es = -9999
    #    else: 
    #    if (days==None) | (days==0):
    #        day = 1
    #    elif (days >= 1):
    #        
    #        day = days
        if data.ndim == 1:
            data = data.reshape((-1, 1))
            
        es             =  -9999                                                                                  ## defines error_sign  
        jd             = np.array(ufz.date2dec(ascii = period))                                                  ## converts date into julian day
        jd_int         = jd.astype(np.int)                                                                       ## converts julian day into integer
        jd_uni         = np.unique(jd_int)                                                                       ## takes only one entry per julian day
        tmax_day       = np.ma.ones((jd_uni.size, data.shape[1]))*es                                             ## Space for Tmax values, as long as given days
        tmax_day_indx  = np.ma.zeros((jd_uni.size, data.shape[1]), dtype=np.int)                                 ## Space for Tmax index, as long as given days
        data_masked    = np.ma.array(data, mask=(data==es))                                                      ## maskes original data where values are not valid
            
        count = 0
        for i in jd_uni:                                                                                         ## loop along the number of days
            ii                     = np.where(jd_int == i)[0]                                                    ## finds index for each single day
            tmax_day[count,:]      = np.ma.amax(data_masked[ii,:data.shape[1]], axis=0)                          ## stores T_max for each day
            tmax_day_indx[count,:] = ii[np.ma.argmax(data_masked[ii,:data.shape[1]], axis=0)]                    ## stores index of t_max for each day
            count += 1
           
        tmax = np.empty((data.shape[0],data.shape[1]))                                                           ## Space for Tmax values, as long as given records
        for i in range(data.shape[1]):                                                                           ## loop along the number of records
            xx        = jd[tmax_day_indx[:,i]]                                                                   ## defines Index of values ​​to be interpolated                           
            yy        = tmax_day[:,i]                                                                            ## defines values ​​to be interpolated 
            ii        = ~yy.mask                                                                                 ## excludes masked values for interpolation
            tmax[:,i] = np.interp(jd, xx[ii], yy[ii].compressed())                                               ## interpolates between the T_max values
            
        if swd == None:
            SFD  = (0.0119*((tmax - data_masked[:,:data.shape[1]])/data_masked[:,:data.shape[1]])**1.231)*3600.  ## converts raw data [mV] into Sap Flux density [cm3*cm-2*h-1] 
            #SFD  = SFD.filled(es)                                                                               ## filles masked elements with error sign   
        else:            
            a     = swd/2                                                                                        ## sapwood depth in cm / 2cm for needle lenght = portion of sensor in sapwood
            b     = 1-a                                                                                          ## portion of sensor in inactive xylem
            Tsw  = ((data_masked[:,:data.shape[1]] - (b*tmax))/a)                                                ## define weighted mean of ΔT in the sapwood (a) and ΔT in the inactive xylem (b)
            SFD  = (0.0119*((tmax - Tsw[:,:data.shape[1]])/Tsw[:,:data.shape[1]])**1.231)*3600                   ## converts raw data [mV] into Sap Flux density [cm3*cm-2*h-1] by using weighted mean of ΔT
            #SFD  = SFD.filled(es)                                                                               ## filles masked elements with error sign            
 
        return SFD
    
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)