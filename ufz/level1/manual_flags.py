#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import ufz

__all__ = ['get_manual_flags']

# --------------------------------------------------------------------

def get_manual_flags(flagfile, variable, dendro=None, julian=True):
    """
        Get start and end dates as well as flag values for a specific variable from a manual flag file.

        If start date is missing, 01.01.1900 00:00:00 will be used.
        If end date is missing, 01.01.2099 23:59:59 will be used.
        If no manual flag is found, empty lists wil be returned.

        The manual flag files is supposed to have the following structure
        var_id,start,end,flag,comment
        Lines starting with # will be ignored.
        
        
        Input
        -----
        flagfile   Filename of manual flag file
        variable   variable name in first column of flagfile


        Optional Input
        --------------
        julian     True: return julian dates (default)
                   False: return ascii strings in format DD.MM.YYYY hh:mm
        dendro     True: return also diameter at breast height and corresponding dendrometer reading
                   False: return only dates and flag


        Output
        ------
        list of start dates, list of end dates, list of flags (if dendro=False)
        list of start dates, list of end dates, list of flags, list of diameter at breast height, list of initial dendrometer readings (if dendro=False)


        Examples
        --------
        --> see __init__.py for full example of workflow

        sdate, edate, mflags = get_manual_flags('Manual_Flags_Soilnet_HH.csv"', 'Box02_Moist1')
        sdate, edate, mflags, bdh, d_ini = get_manual_flags('Manual_Flags_Soilnet_HH.csv"', 'Box02_Moist1', dendro=True)


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Aug 2015
        Modified, AW, Aug 2015, optional readout of dendrometer-data:  DBH and Dini
    """
    # Default dates
    sdef = "01.01.1900 00:00:00"
    edef = "01.01.2099 23:59:59"

    # Read manual flag file
    sdat    = ufz.sread(flagfile, comment="#", strarr=True, skip=1)
    var_id  = sdat[:,0]
    start   = sdat[:,1]
    end     = sdat[:,2]
    flag    = sdat[:,3]
    dbh     = sdat[:,4]
    dbh[dbh==''] = '0'
    d_ini   = sdat[:,5]
    d_ini[d_ini==''] = '0'
    # comment = sdat[:,6]
     
    # select variable
    ii = np.where(var_id == variable)[0]
    if ii.size > 0:
        sdate   = start[ii]
        edate   = end[ii]
        mflag   = flag[ii].astype(np.int)
        l_dbh   = dbh[ii].astype(np.float)
        l_d_ini = d_ini[ii].astype(np.float)
        # Fill default dates
        jj = np.where(sdate == "")[0]
        if jj.size > 0: sdate[jj] = sdef
        jj = np.where(edate == "")[0]
        if jj.size > 0: edate[jj] = edef
        # Julian or ascii
        if julian:
            sdate = ufz.date2dec(ascii=sdate)
            edate = ufz.date2dec(ascii=edate)
    else:
        # return empty lists if variable not in manual flag file
        sdate = list()
        edate = list()
        mflag = list()
        l_dbh   = list()
        l_d_ini = list()

    if dendro == True:
        return [sdate, edate, mflag, l_dbh, l_d_ini]  
    else:
        return [sdate, edate, mflag]   

   

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # base_folder = '/Volumes/Gruppen/tereno/CHS-Data'
    # data_folder = base_folder+'/HohesHolz'
    # location_short_name = 'HH'
    # man_file = data_folder+"/Soilnet/Manual_Flags_Soilnet_"+location_short_name+".csv"

    # v = 'Box02_Moist1'
    # print(v)
    # s, e, m = get_manual_flags(man_file, v)
    # for i in range(len(s)):
    #     print(s[i], e[i], m[i])
    # s, e, m = get_manual_flags(man_file, v, julian=False)
    # for i in range(len(s)):
    #     print(s[i], e[i], m[i])

    # v = 'Box27_Moist1'
    # print(v)
    # s, e, m = get_manual_flags(man_file, v, julian=False)
    # for i in range(len(s)):
    #     print(s[i], e[i], m[i])
