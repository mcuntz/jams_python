#!/usr/bin/env python
import numpy as np
import netcdftime as nt

def date2dec(calendar = 'standard', units = False,
             excelerr = True, yr = 0001,
             mo = 01, dy = 01, hr = 00,
             mi = 00, sc = 00, ascii = None, eng = None):
    """
        Converts numpy arrays with calendar date into
        numpy arrays with decimal date. Supported calendar
        formats are standard, gregorian, julian, proleptic_gregorian,
        excel1900, excel1904, 365_day, noleap, 366_day, all_leap,
        360_day, decimal, or decimal360

        Input is year, month day, hour, minute,
        second or a combination of them. ASCII input
        is possible, too.

        Output is decimal date with day as unit.

        Requires 'netcdftime.py' from the module
        netcdftime available at:
        http://netcdf4-python.googlecode.com


        DEFINITION
        ----------
        def date2dec(calendar = 'standard', units = False,
             excelerr = True, yr = None,
             mo = None, dy = None, hr = None,
             mi = None, sc = None, ascii = None, eng = None):


        INPUT
        -----
        yr       -> input array with year
        mo       -> input array with month
        dy       -> input array with day
        hr       -> input array with hour
        mi       -> input array with minute
        sc       -> input array with second
        ascii    -> input array with strings of the format
                    'dd.mm.yyyy hh:mm:ss'. If seconds are
                    missing ss will be set to 00. If ascii
                    input is chosen by user, other inputs
                    will be neglected.
        eng      -> input array with strings of the format
                    'yyyy-mm-dd hh:mm:ss'. If seconds are
                    missing ss will be set to 00. If eng
                    input is chosen by user, other inputs
                    will be neglected.


        PARAMETERS
        ----------
        calendar -> Input date format. Default value is
                   'standard'.

           'standard', 'gregorian'
                       =   Input date is standard format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until
                           05.03.1583 00:00:00 and gregorian
                           calendar from 15.03.1583 00:00:00
                           until now. Missing 10 days don't
                           exsist.
           'julian'    =   Input date is julian format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until now.
           'proleptic_gregorian'
                       =   Input date is gregorian format.
                           Input is in gregorian calendar from
                           01.01.0001 00:00:00 until now.
           'excel1900' =   Input date is excel 1900 format.
                           Input date is excel date with its
                           units at 01.01.1900 00:00:00 until
                           now.
           'excel1904' =   Input date is excel 1904 (lotus) format.
                           Input date is excel date with its
                           units at 01.01.1904 00:00:00 until now.
           '365_day', 'noleap'
                       =   Input date is 365 days format. Input date
                           consists of common years only (No leap years)
                           with its units at 01.01.0001 00:00:00 until now.
           '366_day', 'all_leap'
                       =   Input date is 366 days format. Input date
                           consists of leap years only (No common years)
                           with its units at 01.01.0001 00:00:00 until now.
           '360_day'   =   Input date is 360 days format.  Input
                           date consists of years with only 360 days
                           (30 days per month)with its units at
           'decimal'    =  Output is decimal year.
           'decimal360' =  Output is decimal year with a year of 360 days, i.e. 12 month with 30 days each.


        OPTIONAL ARGUMENTS
        ------------------
        units    -> Time units can be set by user. Input must be a
                     string in the format 'yyyy-mm-dd hh:mm:ss'.
                     Default values are set automatically.
        excelerr  -> In Excel the year 1900 is normally considered
                     as leap year, which is wrong. By default, this
                     error is taken into account (excelerr = True).
                     For excelerr = False, 1900 is considered as no
                     leap year.


        OUTPUT
        ------
        output -> Output numpy array with decimal date.


        EXAMPLES
        --------
        #calendar = 'standard'
array(['05.01.2000 12:30:15', '24.04.1810 16:15:10', '15.07.1630 10:20:40',
           '20.09.1510 14:35:50', '18.03.1271 19:41:34', '27.08.0619 11:08:37',
           '23.08.-1579 20:03:41', '01.01.-4712 12:00:00'],
          dtype='|S20')
        >>> year   = np.array([2000,1810,1630,1510,1271,619,-1579,-4712])
        >>> month  = np.array([1,4,7,9,3,8,8,1])
        >>> day    = np.array([5,24,15,20,18,27,23,1])
        >>> hour   = np.array([12,16,10,14,19,11,20,12])
        >>> minute = np.array([30,15,20,35,41,8,3,0])
        >>> second = np.array([15,10,40,50,34,37,41,0])
        >>> decimal = date2dec(calendar = 'standard', \
                      yr=year, mo=month, dy=day, hr=hour, mi=minute, sc=second)
        >>> print np.round(decimal, 8)
        [ 2451549.02100694  2382262.17719907  2316600.93101852  2272848.10821759
          2185367.32053241  1947385.96431713  1144563.3358912         0.        ]
        >>> decimal = date2dec(calendar = 'standard', yr = year, mo = 06, dy = 15, hr = 12, mi = minute, sc = second)
        >>> print np.round(decimal, 8)
        [  2.45171102e+06   2.38231401e+06   2.31657101e+06   2.27275102e+06
           2.18545603e+06   1.94731301e+06   1.14449400e+06   1.66000000e+02]

        >>> a = np.array(['05.01.2000 12:30:15',\
                              '24.04.1810 16:15:10',\
                              '15.07.1630 10:20:40',\
                              '20.09.1510 14:35:50',\
                              '18.03.1271 19:41:34',\
                              '27.08. 619 11:08:37',\
                              '23.08.-1579 20:03:41',\
                              '01.01.-4712 12:00:00'])
        >>> decimal = date2dec(calendar = 'standard', ascii = a)
        >>> print np.round(decimal, 8)
        [ 2451549.02100694  2382262.17719907  2316600.93101852  2272848.10821759
          2185367.32053241  1947385.96431713  1144563.3358912         0.        ]

        #calendar = 'julian'
        >>> b = np.array(['05.01.2000 12:30:15',\
                          '24.04.1810 16:15:10',\
                          '05.07.1630 10:20:40',\
                          '20.09.1510 14:35:50',\
                          '18.03.1271 19:41:34',\
                          '27.08. 619 11:08:37',\
                          '23.08.-1579 20:03:41',\
                          '01.01.-4712 12:00:00'])
        >>> decimal = date2dec(calendar = 'julian', ascii = b)
        >>> print np.round(decimal, 8)
        [ 2451562.02100694  2382274.17719907  2316600.93101852  2272848.10821759
          2185367.32053241  1947385.96431713  1144563.3358912         0.        ]

        # calendar = 'proleptic_gregorian'
        >>> c = np.array(['28.04.1971 12:30:15',\
                          '28.02. 986 12:30:15',\
                          '31.07. 493 22:20:40',\
                          '02.01.   1 00:00:00'])
        >>> decimal = date2dec(calendar = 'proleptic_gregorian', ascii = c)
        >>> print np.round(decimal, 8)
        [  7.19644521e+05   3.59822521e+05   1.79911931e+05   1.00000000e+00]

        #calendar = 'excel1900' WITH excelerr = True -> 1900
        #considered as leap year
        >>> d = np.array(['05.01.2000 12:30:15',\
                          '27.05.1950 16:25:10',\
                          '13.08.1910 10:40:55',\
                          '01.03.1900 00:00:00',\
                          '29.02.1900 00:00:00',\
                          '28.02.1900 00:00:00',\
                          '01.01.1900 00:00:00'])
        >>> decimal = date2dec(calendar = 'excel1900', ascii = d)
        >>> print np.round(decimal, 8)
        [  3.65305210e+04   1.84106841e+04   3.87844508e+03   6.10000000e+01
           6.00000000e+01   5.90000000e+01   1.00000000e+00]

        #calendar = 'excel1900' WITH excelerr = False -> 1900
        #considered as NO leap year
        >>> e = np.array(['06.01.2000 12:30:15',\
                          '28.05.1950 16:25:10',\
                          '14.08.1910 10:40:55',\
                          '02.03.1900 00:00:00',\
                          '01.03.1900 00:00:00',\
                          '28.02.1900 00:00:00',\
                          '01.01.1900 00:00:00'])
        >>> decimal = date2dec(calendar = 'excel1900', ascii = e, excelerr = False)
        >>> print np.round(decimal, 8)
        [  3.65305210e+04   1.84106841e+04   3.87844508e+03   6.10000000e+01
           6.00000000e+01   5.90000000e+01   1.00000000e+00]

        #calendar = 'excel1904'
        >>> f = np.array(['06.01.2004 12:30:15',\
                          '28.05.1954 16:25:10',\
                          '14.08.1914 10:40:55',\
                          '02.03.1904 00:00:00',\
                          '01.03.1904 00:00:00',\
                          '29.02.1904 00:00:00',\
                          '01.01.1904 00:00:00'])
        >>> decimal = date2dec(calendar = 'excel1904', ascii = f)
        >>> print np.round(decimal, 8)
        [ 36530.52100694  18410.68414352   3878.44508102     61.             60.
             59.              0.        ]
        >>> decimal = date2dec(calendar = 'excel1904', ascii = f, units = '1890-01-01 00:00:00')
        >>> print np.round(decimal, 8)
        [ 41642.52100694  23522.68414352   8990.44508102   5173.           5172.
           5171.           5112.        ]

        #calendar = '365_day'
        >>> g = np.array(['18.08.1972 12:30:15',\
                          '25.10. 986 12:30:15',\
                          '28.11. 493 22:20:40',\
                          '01.01.   1 00:00:00'])
        >>> decimal = date2dec(calendar = '365_day', ascii = g)
        >>> print np.round(decimal, 8)
        [ 719644.52100694  359822.52100694  179911.93101852       0.        ]

        #calendar = '366_day'
        >>> h = np.array(['29.03.1967 12:30:15',\
                          '14.02. 984 12:30:15',\
                          '24.07. 492 22:20:40',\
                          '01.01.   1 00:00:00'])
        >>> decimal = date2dec(calendar = '366_day', ascii = h)
        >>> print np.round(decimal, 8)
        [ 719644.52100694  359822.52100694  179911.93101852       0.        ]

        # 350_day does not work with netcdftime.py version equal or below 0.9.2
        #calendar = '360_day'
        #>>> k = np.array(['05.01.2000 12:30:15',\
        #                  '03.07.1000 12:30:15',\
        #                  '02.10. 500 22:20:40',\
        #                  '01.01.   1 00:00:00'])
        #>>> decimal = date2dec(calendar = '360_day', ascii = k)
        #>>> print np.round(decimal, 8)
        #[ 719644.52100694  359822.52100694  179911.93101852       0.        ]

        >>> date2dec(yr=1992, mo=01, dy=26, hr=02, mi=00, sc=00, calendar='decimal')
       	1992.0685336976321
        >>> date2dec(ascii='26.01.1992 02:00', calendar='decimal360')
       	1992.0696759259258

        History
        -------
        Written  AP, Jun 2010
        Modified MC, Feb 2012 - All input can be scalar or array, also a mix
                              - Changed checks for easier extension
                              - decimal, decimal360
    """

    #
    # Checks
    calendars = ['standard', 'gregorian', 'julian', 'proleptic_gregorian',
                 'excel1900', 'excel1904', '365_day', 'noleap', '366_day',
                 'all_leap', '360_day', 'decimal', 'decimal360']
    if ((nt.__version__ <= '0.9.2') & (calendar == '360_day')):
        raise ValueError("date2dec error: Your version of netcdftime.py is equal"
                         " or below 0.9.2. The 360_day calendar does not work with"
                         " arrays here. Please download a newer one.")
    calendar = calendar.lower()
    if (calendar not in calendars):
        raise ValueError("date2dec error: Wrong calendar!"
                    " Choose: "+''.join([i+' ' for i in calendars]))
    # if ascii input is given by user, other input will be neglected
    # calculation of input size and shape
    if (ascii != None) and (eng != None):
        raise ValueError("date2dec error: 'ascii' and 'eng' mutually exclusive")
    if ascii != None:
        isarr = np.size(np.shape(ascii))
        if (isarr==0):
            ascii = np.array([ascii])
        insize   = np.size(ascii)
        outshape = np.shape(ascii)
        asciifl  = ascii.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of ascii strings to implement in datetime object. missing seconds
        # will be set to 00.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in xrange(insize):
            aa      = asciifl[i].split('.')
            dy[i]   = int(aa[0])
            mo[i]   = int(aa[1])
            tail    = aa[2].split()
            yr[i]   = int(tail[0])
            tim     = tail[1].split(':')
            hr[i]   = int(tim[0])
            mi[i]   = int(tim[1])
            if len(tim) > 2:
                sc[i] = int(tim[2])
            else:
                sc[i] = 00
            timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    if eng != None:
        isarr = np.size(np.shape(eng))
        if (isarr==0):
            eng = np.array([eng])
        insize   = np.size(eng)
        outshape = np.shape(eng)
        engfl  = eng.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of eng strings to implement in datetime object. missing seconds
        # will be set to 00.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in xrange(np.size(eng)):
            ee      = engfl[i].split('-')
            yr[i]   = int(ee[0])
            mo[i]   = int(ee[1])
            tail    = ee[2].split()
            dy[i]   = int(tail[0])
            tim     = tail[1].split(':')
            hr[i]   = int(tim[0])
            mi[i]   = int(tim[1])
            if len(tim) > 2:
                sc[i] = int(tim[2])
            else:
                sc[i] = 00
            timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    # if no ascii input, other inputs will be concidered
    # calculation of input sizes, shapes and number of axis
    if ((ascii == None) & (eng == None)):
        isarr1 = np.size(np.shape(yr))
        if (isarr1==0): yr = np.array([yr])
        isarr2 = np.size(np.shape(mo))
        if (isarr2==0): mo = np.array([mo])
        isarr3 = np.size(np.shape(dy))
        if (isarr3==0): dy = np.array([dy])
        isarr4 = np.size(np.shape(hr))
        if (isarr4==0): hr = np.array([hr])
        isarr5 = np.size(np.shape(mi))
        if (isarr5==0): mi = np.array([mi])
        isarr6 = np.size(np.shape(sc))
        if (isarr6==0): sc = np.array([sc])
        isarr = isarr1 + isarr2 + isarr3 + isarr4 + isarr5 + isarr6
        shapes = [np.shape(yr), np.shape(mo), np.shape(dy), np.shape(hr), np.shape(mi), np.shape(sc)]
        nyr    = np.size(yr)
        nmo    = np.size(mo)
        ndy    = np.size(dy)
        nhr    = np.size(hr)
        nmi    = np.size(mi)
        nsc    = np.size(sc)
        sizes  = [nyr,nmo,ndy,nhr,nmi,nsc]
        nmax   = np.amax(sizes)
        ii     = np.argmax(sizes)
        outshape = shapes[ii]
        if nyr < nmax:
            if nyr == 1: yr  = np.ones(outshape)*yr
            else: raise ValueError("date2dec error: size of yr != max input or 1.")
        if nmo < nmax:
            if nmo == 1: mo  = np.ones(outshape)*mo
            else: raise ValueError("date2dec error: size of mo != max input or 1.")
        if ndy < nmax:
            if ndy == 1: dy  = np.ones(outshape)*dy
            else: raise ValueError("date2dec error: size of dy != max input or 1.")
        if nhr < nmax:
            if nhr == 1: hr  = np.ones(outshape)*hr
            else: raise ValueError("date2dec error: size of hr != max input or 1.")
        if nmi < nmax:
            if nmi == 1: mi  = np.ones(outshape)*mi
            else: raise ValueError("date2dec error: size of mi != max input or 1.")
        if nsc < nmax:
            if nsc == 1: sc  = np.ones(outshape)*sc
            else: raise ValueError("date2dec error: size of sc != max input or 1.")
        indate  = [yr, mo, dy, hr, mi, sc]
        insize  = [np.size(i) for i in indate]
        inshape = [np.shape(i) for i in indate]
        dimnum  = [len(i) for i in inshape]
        # calculation of maximum input size and maximum number of axis for
        # reshaping the output
        indate  = [i.flatten() for i in indate]
        outsize = max(insize)
        timeobj = np.zeros(outsize, dtype=object)
        # datetime object is constructed
        for i in xrange(outsize):
            timeobj[i] = nt.datetime(indate[0][i], indate[1][i], indate[2][i], indate[3][i], indate[4][i], indate[5][i])
    # depending on chosen calendar and optional set of the time units
    # decimal date is calculated
    if (calendar == 'standard') or (calendar == 'gregorian'):
        if units == False: units = '0001-01-01 12:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='gregorian')+1721424
    elif calendar == 'julian':
        if units == False: units = '0001-01-01 12:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='julian')+1721424
    elif calendar == 'proleptic_gregorian':
        if units == False: units = '0001-01-01 00:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='proleptic_gregorian')
    elif calendar == 'excel1900':
        if units == False: units = '1900-01-00 00:00:00'
        if excelerr:
            output = nt.date2num(timeobj,'days since %s' % (units), calendar='julian')
        else:
            output = nt.date2num(timeobj,'days since %s' % (units), calendar='gregorian')
    elif calendar == 'excel1904':
        if units == False: units = '1904-01-01 00:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='gregorian')
    elif (calendar == '365_day') or (calendar == 'noleap'):
        if units == False: units = '0001-01-01 00:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='365_day')
    elif (calendar == '366_day') or (calendar == 'all_leap'):
        if units == False: units = '0001-01-01 00:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='366_day')
    elif calendar == '360_day':
        if units == False: units = '0001-01-01 00:00:00'
        output = nt.date2num(timeobj,'days since %s' % (units), calendar='360_day')
    elif calendar == 'decimal':
        ntime = np.size(yr)
        leap  = (((yr%4)==0) & ((yr%100)!=0)) | ((yr%400)==0)
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([ [-9,0, 31, 59, 90,120,151,181,212,243,273,304,334,365], \
                           [-9,0, 31, 60, 91,121,152,182,213,244,274,305,335,366] ])
        for i in xrange(ntime):
            tdy[i] = tdy[i] + np.array(diy[leap,mo[i]], dtype=np.float)
        days_year = 365.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       ((days_year+np.array(leap, dtype=np.float))*24.) )
        # for numerical stability, i.e. back and forth transforms
        import sys
        output += sys.float_info.epsilon
    elif calendar == 'decimal360':
        ntime = np.size(yr)
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([-9,  0, 30, 60, 90,120,150,180,210,240,270,300,330,360])
        for i in xrange(ntime):
            tdy[i] = tdy[i] + np.array(diy[mo[i]], dtype=np.float)
        days_year = 360.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       (days_year*24.) )
        # for numerical stability, i.e. back and forth transforms
        import sys
        output += sys.float_info.epsilon
    else:
        raise ValueError("date2dec error: calendar not implemented; should have been catched before.")


    # return of reshaped output
    output = np.reshape(output, outshape)
    if isarr == 0:
        output = np.float(output)
    return output
# END OF FUNCTION
###############################################################
# DOCTEST:
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# END
