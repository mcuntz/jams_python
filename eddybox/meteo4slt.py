#!/usr/bin/env python
import numpy as np
from date2dec import date2dec # ufz
from dec2date import dec2date # ufz
import os as os
import re as re

def meteo4slt(sltdir, metfile, p_t_rh, outfile,
              pat='[a-zA-Z0-9]*.slt|[a-zA-Z0-9]*.SLT', delimiter=',',
              skiprows=1, format='ascii'):
    """       
        To supply EddyFlux (Kolle & Rebmann, 2007) with meteorological
        data, it is necessary to extract for each available *.slt file the
        corresponding air pressure, air temperature and air relative humidity.
        The module looks in sltdir for available *.slt files and uses the metfile
        for syncronisation
        
        
        Definition
        ----------
        meteo4slt(sltdir, metfile, p_t_rh, outfile,
              pat='[a-zA-Z0-9]*.slt|[a-zA-Z0-9]*.SLT', delimiter=',',
              skiprows=1, format='ascii'):
        
        Input
        ----- 
        sltdir      str, path of the folder containing the *.slt files 
        metfile     str, path of the meteo file 
        p_t_rh      list, column number of Pressure, T, rH 
        outfile     str, path of the output file 
        
        
        Optional Input
        --------------
        pat         str, regular expression, describing the name pattern of
                    the *.slt files in the indir folder
        delimiter   str, column delimiter of the meteo file (default=',')
        skiprows    int, number of rows to skip at the beginning of the met file
                    e.g. header lines (default=1)
        format      str, time format of the meteo file (default='ascii')
                    
        
        Output
        ------
        outfile     file, containing Pressure, T and rH values of the meteo
                    file for each *.slt file
                    
        
        Restrictions
        ------------
        - assumes site name in slt fielname is only one character
          (e.g. W20133652300.slt AND NOT WULF20133652300.slt)
        - currently only supports format='ascii', nice would be 'eng'
        - does not work for data containning multiple years because of ugly
          doy calculation
        - works only for half hourly *.stl and half hourly meteo data
        
        
        License
        -------
        This file is part of the UFZ Python library.
    
        The UFZ Python library is free software: you can redistribute it and/or 
        modify it under the terms of the GNU Lesser General Public License as 
        published by the Free Software Foundation, either version 3 of the License,
        or (at your option) any later version.
    
        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.
    
        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not,
        see <http://www.gnu.org/licenses/>.
    
        Copyright 2014 Arndt Piayda
    
    
        History
        -------
        Written,  AP, Jul 2014
        Modified, MC, Aug 2014 - clean up and Python 3        
    """
    
    # half hour in jd
    halfhjd = 1800./86400.

    ###########################################################################
    # reading slt directory
    dirlist = os.listdir(sltdir)
    sltdates = np.array([])
    
    ###########################################################################
    # remove all files and folders from list which are not *.slt files and get size
    pat = re.compile(pat)
    for item in dirlist:
        if re.search(pat, item):
            sltdates = np.append(sltdates, item[1:-4])
    
    # replace time steps which dont fit in half hour interval
    mi = np.array([x[-2:] for x in sltdates])
    mi = np.where(mi.astype(int)<30, '00', '30')
    sltdates = np.array([sltdates[i][:-2]+mi[i] for i in range(mi.size)])
    
    ###########################################################################
    # load meteo data
    metdata = np.loadtxt(metfile, dtype='|S100', delimiter=delimiter,
                         skiprows=skiprows, usecols=[0]+p_t_rh)
    
    # get the metdate
    if format == 'ascii':
        # shift met dates one half hour back since slt time stamp marks half
        # hour start but meteo date mark half hour end
        jd = date2dec(ascii=metdata[:,0]) - halfhjd
        fulldate = dec2date(jd, fulldate=True)
        adate    = dec2date(jd, ascii=True)
        doy  = np.ceil(date2dec(ascii=adate)-
                       date2dec(yr=fulldate[0][0],
                                mo=01,dy=01,hr=00,mi=00, sc=00)).astype(int)
        doy = np.where((fulldate[3]==0) & (fulldate[4]==0), doy+1, doy)
        metdates = np.array(['%04i%03i%02i%02i'%(fulldate[0][i], doy[i],
                                                fulldate[3][i], fulldate[4][i])
                             for i in range(jd.size)])
    
    else:
        raise ValueError('meteo4slt: unknown format!')
    
    ###########################################################################
    # sync both dates
    mask = np.in1d(metdates, sltdates)
    
    ###########################################################################
    # write output
    np.savetxt(outfile, metdata[mask,1:], '%s', ',')
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
