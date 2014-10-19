#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import netcdftime as nt

def dec2date(indata, calendar='standard', refdate=None, units=None,
             excelerr=True, fulldate=None, yr=False,
             mo=False, dy=False, hr=False, mi=False,
             sc=False, ascii=False, eng=False):
    """
        Converts numpy arrays with decimal date into
        numpy arrays with calendar date. Supported time formats
        are:
        standard, gregorian, julian, proleptic_gregorian,
        excel1900, excel1904, 365_day, noleap, 366_day, all_leap,
        and 360_day.

        Input is decimal date in units of days.
        Output is year, month, day, hour, minute, second
        or a combination of them. ASCII output format is possible.


        Definition
        ----------
        def dec2date(indata, calendar='standard', refdate=None, units=None,
                     excelerr=True, fulldate=None, yr=False,
                     mo=False, dy=False, hr=False, mi=False,
                     sc=False, ascii=False, eng=False):


        Input
        -----
        indata -> Input numpy array with decimal date.
                  Input date must be positive.


        Parameters
        ----------
        Calendar -> Input date format. Default value is
                   'standard'.

           'standard', 'gregorian'
                       =   Input date is standard format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until
                           05.03.1583 00:00:00 and gregorian
                           calendar from 15.03.1583 00:00:00
                           until now. Missing 10 days don't
                           exist.
           'julian'    =   Input date is julian format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until now.
           'proleptic_gregorian'
                       =   Input date is gregorian format.
                           Input is in gregorian calendar from
                           01.01.0000 00:00:00 until now.
           'excel1900' =   Input date is excel 1900 format.
                           Input date is excel date with its
                           refdate at 01.01.1900 00:00:00 until
                           now.
           'excel1904' =   Input date is excel 1904 (lotus) format.
                           Input date is excel date with its
                           refdate at 01.01.1904 00:00:00 until now.
           '365_day', 'noleap'
                       =   Input date is 365 days format. Input date
                           consists of common years only (No leap years)
                           with its refdate at 01.01.0001 00:00:00 until now.
           '366_day', 'all_leap'
                       =   Input date is 366 days format. Input date
                           consists of leap years only (No common years)
                           with its refdate at 01.01.0001 00:00:00 until now.
           '360_day'   =   Input date is 360 days format.  Input
                           date consists of years with only 360 days
                           (30 days per month) with its refdate at
                           01.01.0001 00:00:00 until now.
           'decimal'    =  Input is decimal year.
           'decimal360' =  Input is decimal year with a year of 360 days, i.e. 12 month with 30 days each.


        Optional Arguments
        ------------------
        refdate   -> Reference date for 'days since refdate' can be set by user. Input must be a
                     string in the format 'yyyy-mm-dd hh:mm:ss'.
                     Default values for different calendars are set automatically.
        units     -> Units of the the time stamp can be given. This can only be the following two:
                     'day as %Y%m%d.%f', or
                     'days since refdate', with refdate in the format 'yyyy-mm-dd hh:mm:ss'.
                     If units='day as %Y%m%d.%f' then calendar is ignored and the dates will be returned
                     assuming that %f is fractional date from 00:00:00 h.
        excelerr  -> In Excel the year 1900 is normally considered
                     as leap year, which is wrong. By default, this
                     error is taken into account (excelerr = True).
                     For excelerr = False, 1900 is considered as a
                     common year.


        Output
        ------
        fulldate -> output arrays with year, month, day, hour, minute, second
                    Default fulldate is overwritten by selection of special output yr,mn,dy,hr,mi,sc
        yr       -> output array with year
        mo       -> output array with month
        dy       -> output array with day
        hr       -> output array with hour
        mi       -> output array with minute
        sc       -> output array with second
        ascii    -> output array with strings of the format
                    'dd.mm.yyyy hh:mm:ss'
        eng      -> output array with strings of the format
                    'yyyy-mm-dd hh:mm:ss'


        Restrictions
        ------------
        Most versions of datetime do not support neagtive years,
        i.e. Julian days < 1721423.5 = 01.01.0001 00:00.

        There is an issue in netcdftime version < 0.9.5 in proleptic_gregorian for dates before year 301:
          ufz.dec2date(ufz.date2dec(ascii='01.01.0300 00:00:00', calendar='proleptic_gregorian'),
                       calendar='proleptic_gregorian')
            [300, 1, 2, 0, 0, 0]
          ufz.dec2date(ufz.date2dec(ascii='01.01.0301 00:00:00', calendar='proleptic_gregorian'),
                       calendar='proleptic_gregorian')
            [301, 1, 1, 0, 0, 0]

        Requires 'netcdftime.py' from module netcdftime available at:
        http://netcdf4-python.googlecode.com


        Examples
        --------
        # Some implementations of datetime have problems with negative years
        >>> import datetime
        >>> if datetime.MINYEAR > 0:
        ...     print('The minimum year in your datetime implementation is ', datetime.MINYEAR)
        ...     print('i.e. it does not support negative years (BC).')

        #calendar = 'standard'
        >>> year   = np.array([2000,1810,1630,1510,1271,619,1])
        >>> month  = np.array([1,4,7,9,3,8,1])
        >>> day    = np.array([5,24,15,20,18,27,1])
        >>> hour   = np.array([12,16,10,14,19,11,12])
        >>> minute = np.array([30,15,20,35,41,8,0])
        >>> second = np.array([15,10,40,50,34,37,0])
        >>> from date2dec import date2dec
        >>> decimal = date2dec(calendar='standard',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> year, month, day, hour, minute, second = dec2date(decimal, calendar= 'standard', fulldate = True)
        >>> print(year)
        [2000 1810 1630 1510 1271  619    1]
        >>> print(month)
        [1 4 7 9 3 8 1]
        >>> print(day)
        [ 5 24 15 20 18 27  1]
        >>> print(hour)
        [12 16 10 14 19 11 12]
        >>> print(minute)
        [30 15 20 35 41  8  0]
        >>> print(second)
        [15 10 40 50 34 37  0]

        # calendar = 'julian'
        >>> decimal = date2dec(calendar='julian',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> year = dec2date(decimal, calendar='julian', yr = True)
        >>> print(year)
        [2000 1810 1630 1510 1271  619    1]

        # calendar = 'proleptic_gregorian'
        >>> decimal = date2dec(calendar='proleptic_gregorian',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> ascii = dec2date(decimal, calendar='proleptic_gregorian', ascii = True)
        >>> print(ascii[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']

        # calendar = 'excel1900' WITH excelerr = True -> 1900 considered as leap year
        >>> decimal = date2dec(calendar='excel1900',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> year, day = dec2date(decimal, calendar='excel1900', yr = True, dy = True)
        >>> print(year)
        [2000 1810 1630 1510 1271  619    1]
        >>> print(day)
        [ 5 24 15 20 18 27  1]

        # calendar = 'excel1900' WITH excelerr = False -> 1900 considered as NO leap year
        # Older versions of netcdftime.py produced unnecessary output (Line 262)
        >>> decimal = date2dec(calendar='excel1900',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second,excelerr=False)
        >>> if nt.__version__ < '0.9.4':
        ...     asciidate = dec2date(decimal, calendar='excel1900', ascii = True, excelerr = False)
        ... elif nt.__version__ == '0.9.4':
        ...     asciidate = dec2date(decimal, calendar='excel1900', ascii = True, excelerr = False)
        ...     for i in range(3):
        ...         print('0 300')
        ... else:
        ...     asciidate = dec2date(decimal, calendar='excel1900', ascii = True, excelerr = False)
        ...     for i in range(7):
        ...         print('0 300')
        0 300
        0 300
        0 300
        0 300
        0 300
        0 300
        0 300

        >>> print(asciidate[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']

        #calendar = 'excel1904'
        >>> decimal = date2dec(calendar='excel1904',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> asciidate = dec2date(decimal, calendar='excel1904', ascii = True)
        >>> print(asciidate[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']
        >>> asciidate = dec2date(decimal, calendar='excel1904', ascii=True, refdate='1910-01-00 00:00:00')
        >>> print(asciidate[::4])
        ['05.01.2006 12:30:15' '18.03.1277 19:41:34']
        >>> print(dec2date(decimal[::4], calendar='excel1904', ascii=True, units='days since 1910-01-00 00:00:00'))
        ['05.01.2006 12:30:15' '18.03.1277 19:41:34']
        >>> print(dec2date(decimal[::4], calendar='excel1904', ascii=True, units='days since 1910-01-01 00:00:00'))
        ['06.01.2006 12:30:15' '19.03.1277 19:41:34']

        # check especially 1900 (no) leap year in Excel
        >>> year1   = np.array([1900,1900,1900,1900])
        >>> month1  = np.array([2,2,3,1])
        >>> day1    = np.array([28,29,1,1])
        >>> decimal = date2dec(calendar='excel1900', yr=year1, mo=month1, dy=day1)
        >>> month2, day2 = dec2date(decimal, calendar='excel1900', mo=True, dy=True)
        >>> print(month2)
        [2 2 3 1]
        >>> print(day2)
        [28 29  1  1]
        >>> decimal = date2dec(calendar='excel1900', yr=year1, mo=month1, dy=day1, excelerr=False)
        >>> month2, day2 = dec2date(decimal, calendar='excel1900', mo=True, dy=True, excelerr=False)
        >>> print(month2)
        [2 3 3 1]
        >>> print(day2)
        [28  1  1  1]
        >>> decimal = date2dec(calendar='excel1904', yr=year1, mo=month1, dy=day1)
        >>> month2, day2 = dec2date(decimal, calendar='excel1904', mo=True, dy=True)
        >>> print(month2)
        [2 3 3 1]
        >>> print(day2)
        [28  1  1  1]

        # calendar = '365_day'
        >>> decimal = date2dec(calendar='365_day',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> asciidate = dec2date(decimal, calendar='365_day', ascii = True)
        >>> print(asciidate[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']

        # calendar = '366_day'
        >>> decimal = date2dec(calendar='366_day',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> asciidate = dec2date(decimal, calendar='366_day', ascii = True)
        >>> print(asciidate[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']

        # calendar = '360_day'
        >>> decimal = date2dec(calendar='360_day',yr=year,mo=month,dy=day,hr=hour,mi=minute,sc=second)
        >>> asciidate = dec2date(decimal, calendar='360_day', ascii = True)
        >>> print(asciidate[::4])
        ['05.01.2000 12:30:15' '18.03.1271 19:41:34']

        >>> print(dec2date(719644.52101, calendar='proleptic_gregorian', ascii = True))
        28.04.1971 12:30:15
        >>> dec = date2dec(ascii='02.03.1910 03:44:55', calendar='decimal')
        >>> print(dec2date(dec, calendar='decimal', ascii=True))
        02.03.1910 03:44:55
        >>> dec = date2dec(ascii='02.03.1910 03:44:55', calendar='decimal360')
        >>> print(dec2date(dec, calendar='decimal360', ascii=True))
        02.03.1910 03:44:55
        >>> print(dec2date([dec,dec], calendar='decimal360', ascii=True))
        ['02.03.1910 03:44:55', '02.03.1910 03:44:55']
        >>> print(dec2date([[dec,dec],[dec,dec],[dec,dec]], calendar='decimal360', ascii=True)[0])
        ['02.03.1910 03:44:55', '02.03.1910 03:44:55']
        >>> print(dec2date(np.array([dec,dec]), calendar='decimal360', ascii=True))
        ['02.03.1910 03:44:55' '02.03.1910 03:44:55']
        >>> print(dec2date(np.array([[dec,dec],[dec,dec],[dec,dec]]), calendar='decimal360', ascii=True)[0:2,0])
        ['02.03.1910 03:44:55' '02.03.1910 03:44:55']

        >>> absolut = np.array([20070102.0034722, 20070102.0069444])
        >>> print(dec2date(absolut, units='day as %Y%m%d.%f', ascii=True))
        ['02.01.2007 00:05:00' '02.01.2007 00:10:00']
        >>> absolut = [20070102.0034722, 20070102.0069444]
        >>> print(dec2date(absolut, units='day as %Y%m%d.%f', ascii=True))
        ['02.01.2007 00:05:00', '02.01.2007 00:10:00']


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

        Copyright 2010-2013 Arndt Piayda, Matthias Cuntz


        History
        -------
        Written  AP, Jun 2010
        Modified MC, Feb 2012 - Input can be scalar or array
                              - Default: fulldate=True
                              - Changed checks for easier extension
                              - decimal, decimal360
                 MC, Jun 2012 - former units keyword is now called refdate
                              - units has now original meaning as in netcdftime
                              - included units='day as %Y%m%d.%f'
                 MC, Feb 2013 - solved Excel leap year problem.
                 MC, Feb 2013 - ported to Python 3
                 AP, May 2013 - solved eng output problem.
                 MC, Oct 2013 - Excel starts at 1 not at 0
                 MC, Oct 2013 - units bugs, e.g. 01.01.0001 was substracted if Julian calendar even with units
    """
    #
    # Constants
    calendars = ['standard', 'gregorian', 'julian', 'proleptic_gregorian',
                 'excel1900', 'excel1904', '365_day', 'noleap', '366_day',
                 'all_leap', '360_day', 'decimal', 'decimal360']
    #
    # Checks
    if ((nt.__version__ <= '0.9.2') & (calendar == '360_day')):
        raise ValueError("dec2date error: Your version of netcdftime.py is equal"
                         " or below 0.9.2. The 360_day calendar does not work with"
                         " arrays here. Please download a newer one.")
    calendar = calendar.lower()
    if (calendar not in calendars):
        raise ValueError("Wrong calendar! Choose: "+''.join([i+' ' for i in calendars]))
    if (refdate!= None) & (units!=None):
        raise ValueError("either refdate or units can be given.")
    #
    # Default
    if np.sum(np.array([yr,mo,dy,hr,mi,sc])) >= 1:
        ii = True
    else:
        ii = False
    if ((ascii | eng | ii) & (fulldate == None)):
        fulldate = False
    if ((not (ascii | eng | ii)) & (fulldate == None)):
        fulldate = True
    if fulldate == True:
        yr = True
        mo = True
        dy = True
        hr = True
        mi = True
        sc = True
        # Further checks
    if np.sum(np.array([ascii,fulldate,eng])) > 1:
        raise ValueError("dec2date error: Only one of ascii, fulldate, or eng can be chosen.")
    if np.sum(np.array([ascii,eng,ii])) > 1:
        raise ValueError("dec2date error: If ascii, fulldate or eng then no special selection yr,mo,dy,hr,mi,sc possible.")
    #
    # Input size and shape
    islist = type(indata) != type(np.array(indata))
    isarr = np.ndim(indata)
    if (islist & (isarr > 2)):
        raise ValueError("dec2date error: input is list > 2D; Use array input")
    if isarr == 0: indata = np.array([indata])
    else: indata = np.array(indata)
    insize  = indata.size
    inshape = indata.shape
    indata  = indata.flatten()
    #
    # depending on chosen calendar and optional set of the time refdate
    # calendar date is calculated
    if units=='day as %Y%m%d.%f':
        fdy    = indata % 1.     # day fraction
        indata = indata - fdy
        day    = indata % 100.   # day
        indata = indata - day
        tmp    = indata % 10000.
        month  = tmp / 100.     # month
        indata = indata - tmp
        year   = indata / 10000. # year
        secs   = np.rint(fdy * 86400.)
        hour   = np.floor(secs/3600.)
        secs   = secs - 3600.*hour
        minute = np.floor(secs/60.)
        second = np.rint(secs - 60.*minute)
    else:
        if (calendar == 'standard') or (calendar == 'gregorian'):
            dec0 = 0
            if units != None:
                unit = units
            elif refdate != None:
                #unit = 'days since %s' % (refdate)
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 12:00:00'
                dec0 = 1721424
            timeobj = nt.num2date(indata-dec0, unit, calendar='gregorian')
        elif calendar == 'julian':
            dec0 = 0
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 12:00:00'
                dec0 = 1721424
            timeobj = nt.num2date(indata-dec0, unit, calendar='julian')
        elif calendar == 'proleptic_gregorian':
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 00:00:00'
            timeobj = nt.num2date(indata, unit, calendar = 'proleptic_gregorian')
        elif calendar == 'excel1900':
            doerr = False
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 1900-01-00 00:00:00'
                if excelerr: doerr = True
            if doerr:
                indata1 = np.where(indata >= 61., indata-1, indata)
                timeobj = nt.num2date(indata1, unit, calendar = 'gregorian')
            else:
                timeobj = nt.num2date(indata, unit, calendar = 'gregorian')
        elif calendar == 'excel1904':
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 1904-01-00 00:00:00'
            timeobj = nt.num2date(indata, unit, calendar = 'gregorian')
        elif (calendar == '365_day') or (calendar == 'noleap'):
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 00:00:00'
            timeobj = nt.num2date(indata, unit, calendar = '365_day')
        elif (calendar == '366_day') or (calendar == 'all_leap'):
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 00:00:00'
            timeobj = nt.num2date(indata, unit, calendar = '366_day')
        elif calendar == '360_day':
            if units != None:
                unit = units
            elif refdate != None:
                unit = 'days since {0:s}'.format(refdate)
            else:
                unit = 'days since 0001-01-01 00:00:00'
            timeobj = nt.num2date(indata, unit, calendar = '360_day')
        elif calendar == 'decimal':
            fyear = np.trunc(indata)
            year  = np.array(fyear, dtype=np.int)
            leap  = (((year%4)==0) & ((year%100)!=0)) | ((year%400)==0)
            fleap = np.array(leap, dtype=np.float)
            fract_date = indata-fyear
            days_year  = 365.
            # date in hours
            fhoy  = fract_date*(days_year+fleap)*24.
            fihoy = np.trunc(fhoy)
            ihoy  = np.array(fihoy, dtype=np.int)
            # minutes
            fmoy    = (fhoy-fihoy)*60.
            fminute = np.trunc(fmoy)
            minute  = np.array(fminute, dtype=np.int)
            # seconds
            second = np.array(np.trunc((fmoy-fminute)*60.), dtype=np.int)
            # months
            fdoy  = (fihoy/24.)+1.
            fidoy = np.trunc(fdoy)
            idoy  = np.array(fidoy, dtype=np.int)
            month = np.zeros(insize, dtype=np.int)
            diy   = np.array([ [-9,0, 31, 59, 90,120,151,181,212,243,273,304,334,365],
                               [-9,0, 31, 60, 91,121,152,182,213,244,274,305,335,366] ])
            for i in range(insize):
                ii = np.squeeze(np.where(idoy[i] > np.squeeze(diy[leap[i],:])))
                month[i] = ii[-1]
            # days
            fday = np.zeros(insize, dtype=np.float)
            for i in range(insize):
                fday[i] = np.trunc(fdoy[i] - np.float(diy[leap[i],month[i]]))
            day = np.array(fday, dtype=np.int)
            # hours
            hour = ihoy % 24
        elif calendar == 'decimal360':
            fyear = np.trunc(indata)
            year  = np.array(fyear, dtype=np.int)
            fract_date = indata-fyear
            days_year  = 360.
            # date in hours
            fhoy  = fract_date*days_year*24.
            fihoy = np.trunc(fhoy)
            ihoy  = np.array(fihoy, dtype=np.int)
            # minutes
            fmoy    = (fhoy-fihoy)*60.
            fminute = np.trunc(fmoy)
            minute  = np.array(fminute, dtype=np.int)
            # seconds
            second = np.array(np.trunc((fmoy-fminute)*60.), dtype=np.int)
            # months
            fdoy  = (fihoy/24.)+1.
            fidoy = np.trunc(fdoy)
            idoy  = np.array(fidoy, dtype=np.int)
            month = np.zeros(insize, dtype=np.int)
            diy   = np.array([ -9, 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 ])
            for i in range(insize):
                ii = np.squeeze(np.where(idoy[i] > diy))
                month[i] = ii[-1]
            # days
            fday = np.zeros(insize, dtype=np.float)
            for i in range(insize):
                fday[i] = np.trunc(fdoy[i] - np.float(diy[month[i]]))
            day = np.array(fday, dtype=np.int)
            # hours
            hour = ihoy % 24
        else:
            raise ValueError("dec2date error: calendar not implemented; should have been catched before.")

        if (calendar not in ['decimal','decimal360']):
            timeobjfl = timeobj.flatten()
            year   = np.array([timeobjfl[i].year for i in range(insize)], dtype=np.int)
            month  = np.array([timeobjfl[i].month for i in range(insize)], dtype=np.int)
            day    = np.array([timeobjfl[i].day for i in range(insize)], dtype=np.int)
            hour   = np.array([timeobjfl[i].hour for i in range(insize)], dtype=np.int)
            minute = np.array([timeobjfl[i].minute for i in range(insize)], dtype=np.int)
            second = np.array([timeobjfl[i].second for i in range(insize)], dtype=np.int)
            if (calendar == 'excel1900') & excelerr:
                ii  = np.where((indata >= 60.) & (indata < 61.))[0]
                if np.size(ii) > 0:
                    month[ii] = 2
                    day[ii]   = 29
    #
    # Ascii output
    if ascii:
        output = ( ['%02d.%02d.%04d %02d:%02d:%02d' %
                      (day[i], month[i], year[i], hour[i], minute[i], second[i]) for i in range(insize)] )
        output = np.reshape(output, inshape)
        if isarr==0:
            output = output[0]
    # Ascii english output
    elif eng:
        output = ( ['%04d-%02d-%02d %02d:%02d:%02d' %
                    (year[i], month[i], day[i], hour[i], minute[i], second[i]) for i in range(insize)] )
        output = np.reshape(output, inshape)
        if isarr==0:
            output = output[0]
    else:
        # Individual output
        # if one, some or all of yr, mo, dy, hr, mi or sc is
        # choosen by the user as output, arrays for datetime
        year   = np.reshape(year, inshape)
        month  = np.reshape(month, inshape)
        day    = np.reshape(day, inshape)
        hour   = np.reshape(hour, inshape)
        minute = np.reshape(minute, inshape)
        second = np.reshape(second, inshape)
        if isarr==0:
            year   = np.int(year)
            month  = np.int(month)
            day    = np.int(day)
            hour   = np.int(hour)
            minute = np.int(minute)
            second = np.int(second)
        # filling of output list:
        output = []
        if yr: output += [year]
        if mo: output += [month]
        if dy: output += [day]
        if hr: output += [hour]
        if mi: output += [minute]
        if sc: output += [second]
        # return output arrays:
        if len(output) == 1:
            output = output[0]

    if isarr != 0:
        if islist:
            ns = np.size(inshape)
            if ns == 1:
                output = [i for i in output]
            elif ns == 2:
                loutput = [ i for i in output[:,0]]
                for i in range(np.size(output[:,0])):
                    loutput[i] = list(np.squeeze(output[i,:]))
                output = loutput
            else:
                raise ValueError("dec2date error: list output > 2D; should have been catched before.")

    return output


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
