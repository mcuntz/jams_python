#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import jams

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
        This file is part of the JAMS Python package.

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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Aug 2015
        Modified, AW, Aug 2015, optional readout of dendrometer-data:  DBH and Dini
        Modified, MC, Jan 2016, igrnore leading and trailing balnk in variable ids
        Modified, BD, Oct 2016, corrected and improved error message if flag missing; introduced error messages if sdate or edate are not in the correct format
    """
    # Default dates
    sdef = "01.01.1900 00:00:00"
    edef = "01.01.2099 23:59:59"

    # Read manual flag file
    sdat    = jams.sread(flagfile, comment="#", strarr=True, skip=0)
    sdat1 = sdat[0,:]      #header   
    sdat = sdat[1:,:]      # get rid of header for the rest
    var_id  = np.array([ i.strip() for i in sdat[:,0] ])
    start   = sdat[:,1]
    end     = sdat[:,2]
    flag    = sdat[:,3]
    if dendro:
        dbh     = sdat[:,4]
        ii = np.where(dbh=='')[0]
        if ii.size > 0: dbh[ii] = '0'
        d_ini   = sdat[:,5]
        ii = np.where(d_ini=='')[0]
        if ii.size > 0: d_ini[ii] = '0'
    # comment = sdat[:,-1]
     
    # select variable
    ii = np.where(var_id == variable)[0]
    if ii.size > 0:
        sdate   = start[ii]
        edate   = end[ii]
        try:
            mflag   = flag[ii].astype(np.int)
        except:
            miss_ind = [ind for ind,fl in enumerate(flag[ii]) if fl == ''][0]
            comm_col = [m for m,com in enumerate(sdat1) if com[0:7] == 'comment']
            comm = sdat[ii,comm_col][miss_ind]
            print('get_manual_flags: forgotten flag in file '+flagfile+' - var/start/end: '+variable+'/'+str(sdate[miss_ind])+'/'+str(edate[miss_ind])+' - comment: '+comm)
            jams.encrypt.sendfail(
                'get_manual_flags: forgotten flag in file '+flagfile+' - var/start/end: '+variable+'/'+str(sdate[miss_ind])+'/'+str(edate[miss_ind])+' - comment: '+comm,
                sender='benjamin.dechant@ufz.de')
        if dendro:
            l_dbh   = dbh[ii].astype(np.float)
            l_d_ini = d_ini[ii].astype(np.float)
        # Fill default dates
        jj = np.where(sdate == "")[0]
        if jj.size > 0: sdate[jj] = sdef
        jj = np.where(edate == "")[0]
        if jj.size > 0: edate[jj] = edef
        # Julian or ascii
        if julian:
            try:
                sdate = jams.date2dec(ascii=sdate)
            except Exception as e:
                for sda_ind,sda in enumerate(sdate):
                    try:
                        sdate = jams.date2dec(ascii=sda)
                    except: 
                        pb_ind = sda_ind
                comm_col = [m for m,com in enumerate(sdat1) if com[0:7] == 'comment']
                comm = sdat[ii,comm_col][pb_ind]
                print('get_manual_flags: invalid start date in file '+flagfile+' - var/start: '+variable+'/'+str(e)[-19:-1]+' - comment: '+comm)
                jams.encrypt.sendfail(
                    'get_manual_flags: invalid start date in file '+flagfile+' - var/start: '+variable+'/'+str(e)[-19:-1]+' - comment: '+comm,
                    sender='benjamin.dechant@ufz.de')
                raise ValueError

            try:
                edate = jams.date2dec(ascii=edate)
            except Exception as e:
                for sda_ind,sda in enumerate(edate):
                    try:
                        sdate = jams.date2dec(ascii=sda)
                    except: 
                        pb_ind = sda_ind
                comm_col = [m for m,com in enumerate(sdat1) if com[0:7] == 'comment']
                comm = sdat[ii,comm_col][pb_ind]
                print('get_manual_flags: invalid end date in file '+flagfile+' - var/end: '+variable+'/'+str(e)[-19:-1]+' - comment: '+comm)
                jams.encrypt.sendfail(
                    'get_manual_flags: invalid end date in file '+flagfile+' - var/end: '+variable+'/'+str(e)[-19:-1]+' - comment: '+comm,
                    sender='benjamin.dechant@ufz.de')
                raise ValueError                        # maybe not necessary

    else:
        # return empty lists if variable not in manual flag file
        sdate = list()
        edate = list()
        mflag = list()
        if dendro:
            l_dbh   = list()
            l_d_ini = list()

    if dendro:
        return [sdate, edate, mflag, l_dbh, l_d_ini]  
    else:
        return [sdate, edate, mflag]   

   

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # base_folder = '/Volumes/Gruppen/CHS-Data'
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
