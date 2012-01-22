#!/usr/bin/env python
import numpy as np
import netcdftime as nt

###############################################################
###############################################################
###############################################################
def dec2date(indata, calendar = 'standard', units = False,
             excelerr = True, fulldate = False, yr = False, 
             mo = False, dy = False, hr = False, mi = False,
             sc = False, ascii = False, eng = False):
    """
        Converts numpy arrays with decimal date into 
        numpy arrays with calendar date. Supported time formats
        are standard, gregorian, julian, proleptic_gregorian, 
        excel1900, excel1904, 365_day, noleap, 366_day, all_leap,
        or 360_day.
        Input is decimal date with day as unit.
        Output is year, month, day, hour, minute,
        second or a combination of them. ASCII output 
        format is possible, too. 
        
        Requires the 'numpy' package available at:
        
        http://numpy.scipy.org/
        
        Requires 'netcdftime.py' from the module 
        netcdftime available at:
        
        http://netcdf4-python.googlecode.com        
        

        DEFINITION
        ----------
        def dec2date(indata, calendar = 'standard', units = False,
             excelerr = True, fulldate = False, yr = False, 
             mo = False, dy = False, hr = False, mi = False,
             sc = False, ascii = False, eng = False):
                

        INPUT
        -----
        indata -> Input numpy array with decimal date.
                  Input date must be positive.
       

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
                           01.01.0001 00:00:00 until now.
       

        OPTIONAL ARGUMENTS
        ------------------
        units    -> Time units can be set by user. Input must be a
                     string in the format 'yyyy-mm-dd hh:mm:ss'.
                     Default values are set automatically.
        excelerr  -> In Excel the year 1900 is normally considered 
                     as leap year, which is wrong. By default, this
                     error is taken into account (excelerr = True).
                     For excelerr = False, 1900 is considered as
                     common year.
       

        OUTPUT
        ------
        fulldate -> output arrays with year, month, day, hour,
                    minute, second 
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
                    
       
        EXAMPLES
        --------
        #calendar = 'standard'         
        >>> a = np.array([2451549.02101, 2382262.17720, 2316600.93102,\
                          2272848.10822, 2185367.32053, 1947385.96432,\
                          1144563.33589, 0.0])
        >>> year, month, day, hour, minute, second \
            = dec2date(a, calendar= 'standard', fulldate = True)
        >>> print year
        [ 2000  1810  1630  1510  1271   619 -1579 -4712]
        >>> print month
        [1 4 7 9 3 8 8 1]
        >>> print day
        [ 5 24 15 20 18 27 23  1]
        >>> print hour
        [12 16 10 14 19 11 20 12]
        >>> print minute
        [30 15 20 35 41  8  3  0]
        >>> print second
        [15 10 40 50 34 37 41  0]
                         
        #calendar = 'julian'
        >>> b = np.array([2451562.02101, 2382274.17720, 2316600.93102,\
                          2272848.10822, 2185367.32053, 1947385.96432,\
                          1144563.33589, 0.0])
        >>> year = dec2date(b, calendar='julian', yr = True)
        >>> print year
        [ 2000  1810  1630  1510  1271   619 -1579 -4712]
        
        # calendar = 'proleptic_gregorian'
        >>> c = np.array([719644.52101, 359822.52101, 179911.93102, 0.0])
        >>> ascii = dec2date(c, calendar='proleptic_gregorian', ascii = True)
        >>> print ascii
        ['27.04.1971 12:30:15' '27.02.0986 12:30:15' '30.07.0493 22:20:40'
         '01.01.0001 00:00:00']
                  
        #calendar = 'excel1900' WITH excelerr = True -> 1900 
        #considered as leap year
        >>> d = np.array([36530.52101, 18410.68414, 3878.44508, 61.0,\
                          60.0, 59.0, 1.0])
        >>> year, day = dec2date(d, calendar='excel1900', yr = True, \
                              dy = True)
        >>> print year
        [2000 1950 1910 1900 1900 1900 1900]
        >>> print day
        [ 5 27 13  1 29 28  1]
        
        #calendar = 'excel1900' WITH excelerr = False -> 1900 
        #considered as NO leap year
        >>> e = np.array([36530.52101, 18410.68414, 3878.44508, 61.0,\
                          60.0, 59.0, 1.0])
        >>> asciidate = dec2date(e, calendar='excel1900', ascii = True \
                      , excelerr = False)
        0 300
        0 300
        0 300
        0 300
        
        #PREVIOUS OUTPUT IS GENERATED BY NETCDFTIME.PY IN LINE 262 UNNECESSARILY
        
        >>> print asciidate
        ['06.01.2000 12:30:15' '28.05.1950 16:25:10' '14.08.1910 10:40:55'
         '02.03.1900 00:00:00' '01.03.1900 00:00:00' '28.02.1900 00:00:00'
         '01.01.1900 00:00:00']
        
        #calendar = 'excel1904'
        >>> f = np.array([36530.52101, 18410.68414, 3878.44508, 61.0,\
                          60.0, 59.0, 0.0])
        >>> asciidate = dec2date(f, calendar='excel1904', \
                                 ascii = True)
        >>> print asciidate
        ['06.01.2004 12:30:15' '28.05.1954 16:25:10' '14.08.1914 10:40:55'
         '02.03.1904 00:00:00' '01.03.1904 00:00:00' '29.02.1904 00:00:00'
         '01.01.1904 00:00:00']
        >>> asciidate = dec2date(f, calendar='excel1904', \
                                 ascii = True, units = '1910-01-01 00:00:00')
        >>> print asciidate
        ['06.01.2010 12:30:15' '28.05.1960 16:25:10' '14.08.1920 10:40:55'
         '03.03.1910 00:00:00' '02.03.1910 00:00:00' '01.03.1910 00:00:00'
         '01.01.1910 00:00:00']
         
        #calendar = '365_day'
        >>> g = np.array([719644.52101, 359822.52101, 179911.93102, 0.0])
        >>> asciidate = dec2date(g, calendar='365_day', ascii = True)
        >>> print asciidate
        ['18.08.1972 12:30:15' '25.10.0986 12:30:15' '28.11.0493 22:20:40'
         '01.01.0001 00:00:00']
         
        #calendar = '366_day'
        >>> h = np.array([ 719644.52101, 359822.52101, 179911.93102, 0.0])
        >>> asciidate = dec2date(h, calendar='366_day', ascii = True)
        >>> print asciidate
        ['29.03.1967 12:30:15' '14.02.0984 12:30:15' '24.07.0492 22:20:40'
         '01.01.0001 00:00:00']
         
        #calendar = '360_day'
        #>>> k = np.array([ 719644.52101, 359822.52101, 179911.93102, 0.0])
        #>>> asciidate = dec2date(k, calendar='360_day', ascii= True)
        #>>> print asciidate        
        #['05.01.2000 12:30:15' '03.07.1000 12:30:15' '02.10.0500 22:20:40'
        # '01.01.0001 00:00:00']


        History
        -------
        Written AP, Jun 2010
    """
    
    if (int(nt.__version__[0]) <= 0) and\
       (int(nt.__version__[2]) <= 9) and\
       (int(nt.__version__[4]) <= 2) and\
       (calendar == '360_day'):
        raise ValueError("date2dec error: Your version of netcdftime.py is equal"
                         " or below 0.9.2. The 360_day calendar does not work with"
                         " arrays here. Please download a newer one.")
    
    calendar = calendar.lower()
    # calculation of input size and shape
    insize   = np.size(indata)
    inshape  = np.shape(indata)

    # if user input of calendar is undefined:
    if (calendar != 'standard') and (calendar != 'gregorian') \
        and (calendar != 'julian') and\
       (calendar != 'proleptic_gregorian') and\
       (calendar != 'excel1900') and (calendar != 'excel1904') and\
       (calendar != '365_day') and (calendar != 'noleap') and (calendar != '366_day')\
        and (calendar != 'all_leap') and (calendar != '360_day'):
        raise ValueError("dec2date error: Wrong calendar!"
                    " Choose 'standard', 'gregorian', 'julian',"
                    " 'proleptic_gregorian',"
                    " 'excel1900', 'excel1904'," 
                    " '365_day', 'noleap', '366_day', 'all_leap' or '360_day'")

    # depending on chosen calendar and optional set of the time units
    # calendar date is calculated
    if   (calendar == 'standard') or (calendar == 'gregorian'):
        if units == False:
            units = '0001-01-01 12:00:00'
        timeobj = nt.num2date(indata-1721424, 
                  'days since %s' % (units),
                  calendar='gregorian')
    elif calendar == 'julian':
        if units == False:
            units = '0001-01-01 12:00:00'
        timeobj = nt.num2date(indata-1721424, 
                  'days since %s' % (units),
                  calendar='julian')
    elif calendar == 'proleptic_gregorian':
        if units == False:
            units = '0001-01-00 00:00:00'
        timeobj = nt.num2date(indata, 
                  'days since %s' % (units),
                  calendar = 'proleptic_gregorian')
    elif calendar == 'excel1900':
        if units == False:
            units = '1900-01-00 00:00:00'
        if excelerr:
            timeobj = nt.num2date(indata, 
                      'days since %s' % (units), 
                      calendar = 'julian')
        else:
            timeobj = nt.num2date(indata, 
                      'days since %s' % (units), 
                      calendar = 'gregorian')
    elif calendar == 'excel1904':
        if units == False:
            units = '1904-01-01 00:00:00'
        timeobj = nt.num2date(indata, 
                  'days since %s' % (units),
                  calendar = 'gregorian')
    elif (calendar == '365_day') or (calendar == 'noleap'):
        if units == False:
            units = '0001-01-01 00:00:00'
        timeobj = nt.num2date(indata, 
                  'days since %s' % (units),
                  calendar = '365_day')
    elif (calendar == '366_day') or (calendar == 'all_leap'):
        if units == False:
            units = '0001-01-01 00:00:00'
        timeobj = nt.num2date(indata, 
                  'days since %s' % (units),
                  calendar = '366_day')
    elif calendar == '360_day':
        if units == False:
            units = '0001-01-01 00:00:00'
        timeobj = nt.num2date(indata, 
                  'days since %s' % (units),
                  calendar = '360_day')        

    timeobjfl = timeobj.flatten()
    # generating ascii output
    if ascii :  
        asciiout = [' ']*insize
        for i in xrange(insize):
            asciiout[i] = '%02d.%02d.%04d %02d:%02d:%02d' %(
                          timeobjfl[i].day, 
                          timeobjfl[i].month, 
                          timeobjfl[i].year, 
                          timeobjfl[i].hour, 
                          timeobjfl[i].minute, 
                          timeobjfl[i].second)
        asciiout = np.reshape(asciiout, inshape)
    # if user choose 'ascii' output and another output 
    # variable simultaneously, output will be 'ascii' ONLY.
        if ascii & (yr or mo or dy or hr or mi or sc or fulldate or eng):
                print ("dec2date warning: If you choose"
                       " ascii, output will be ascii ONLY")
        return asciiout 
    if eng :  
        engout = [' ']*insize
        for i in xrange(insize):
            engout[i] = '%04d-%02d-%02d %02d:%02d:%02d' %(
                          timeobjfl[i].year,
                          timeobjfl[i].month,
                          timeobjfl[i].day,
                          timeobjfl[i].hour, 
                          timeobjfl[i].minute, 
                          timeobjfl[i].second)
        engout = np.reshape(engout, inshape)
    # if user choose 'eng' output and another output 
    # variable simultaneously, output will be 'ascii' ONLY.
        if eng & (yr or mo or dy or hr or mi or sc or fulldate):
                print ("dec2date warning: If you choose"
                       " eng, output will be eng ONLY")
        return engout  
    # if fulldate is selected for output, the entire calendar
    # date is returned.     
    if fulldate:
        yr = True
        mo = True
        dy = True
        hr = True
        mi = True
        sc = True
    # if user missed to choose output:
    if (yr == False) and (mo == False) and \
       (dy == False) and (hr == False) and \
       (mi == False) and (sc == False) and \
       (fulldate == False):    
        raise ValueError("dec2date error: Missing output "
                         "specification! Choose fulldate,"
                         " yr, mo, dy, hr, mi, sc, "
                         "a combination of them or ascii")
    # if one, some or all of yr, mo, dy, hr, mi or sc is 
    # choosen by the user as output, arrays for datetime 
    # objects are initialised:            
    if yr:
        year   = np.zeros(insize, dtype=np.int)
    if mo:
        month  = np.zeros(insize, dtype=np.int)
    if dy:
        day    = np.zeros(insize, dtype=np.int)
    if hr:
        hour   = np.zeros(insize, dtype=np.int)
    if min:
        minute = np.zeros(insize, dtype=np.int)
    if sc:
        second = np.zeros(insize, dtype=np.int)
    # output arrays are filled with datetime objects.
    for i in xrange(insize):       
        if yr:
            year[i]   = int(timeobjfl[i].year)
        if mo:
            month[i]  = int(timeobjfl[i].month)
        if dy:
            day[i]    = int(timeobjfl[i].day)
        if hr:
            hour[i]   = int(timeobjfl[i].hour)
        if mi:
            minute[i] = int(timeobjfl[i].minute)
        if sc:
            second[i] = int(timeobjfl[i].second)
    # reshaping of arrays
    if yr:
        year   = np.reshape(year, inshape)
    if mo:
        month  = np.reshape(month, inshape)
    if dy:
        day    = np.reshape(day, inshape)
    if hr:
        hour   = np.reshape(hour, inshape)
    if mi:
        minute = np.reshape(minute, inshape)
    if sc:
        second = np.reshape(second, inshape)
    # filling of output array:
    output = []
    if yr:
        output += [year]
    if mo:
        output += [month]
    if dy:
        output += [day]
    if hr:
        output += [hour]
    if mi:
        output += [minute]
    if sc:
        output += [second]
    # return output arrays:
    if len(output) == 1:
        return output[0]
    else:
        return output
    
# END OF FUNCTION
###############################################################
###############################################################
###############################################################
# DOCTEST:
if __name__ == '__main__':
    import doctest
    doctest.testmod()
# END
