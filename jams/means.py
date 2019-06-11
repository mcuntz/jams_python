#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.date2dec import date2dec
from jams.dec2date import dec2date

def means(date, dat, year=False, month=False, day=False, hour=False, half_hour=False, minute=False,
          meanday=False, meanmonth=False, seasonal=False,
          sum=False, max=False, min=False, onlydat=False):
    """
        Calculate daily, monthly, yearly, etc. means of data depending on date stamp.

        date is Julian date. dat will be averaged over columns, i.e. first dimension.
        If no option is given, the mean will be over the whole first column.

        Returns centred dates, averages.
            Yearly   dates are centred at June 15, 12:00h.
            Monthly  dates are centred at 15th, 12:00h.
            Daily    dates are centred at 12:00h.
            Hourly   dates are centred at 30 min.
            Half hourly dates are centred at 15 and 45 min.
            Minutely dates are centred at 30 sec.
            Mean daily dates centred on 30 min of 01. January of first year.
            Mean monthly dates centred on 15th of month, 12:00h of first year.
            Seasonal dates are centred at 12:00h of each day of first (leap) year.


        Definition
        ----------
        def means(date, dat, year=False, month=False, day=False, hour=False, half_hour=False, minute=False,
                  meanday=False, meanmonth=False, seasonal=False,
                  sum=False, max=False, min=False, onlydat=False):


        Input
        -----
        date      1D array with Julian dates
        dat       ND (masked-)array


        Optional Input
        --------------
        year      if True, yearly   means.
        month     if True, monthly  means.
        day       if True, daily    means.
        hour      if True, hourly   means.
        half_hour if True, half-hourly means.
        minute    if True, minutely means.
        meanday   if True, mean daily cycle.
        meanmonth if True, mean monthly cycle.
        seasonal  if True, seasonal cycle.
        sum       if True, calculate sums instead of means.
        max       if True, calculate maxima instead of means.
        min       if True, calculate minima instead of means.
        onlydat   if True, return only meandat, else return [outdate, meandat]


        Output
        ------
        outdate    centred date on year, month, day, etc.
        meandat    averaged data in 1st dimension.

        If meanday==True, then all hours will be written; as a masked-array if hours are missing.
        If meanmonth==True, then all months will be written; as a masked-array if months are missing.
        If seasonal==True: input data has to be spaced <= days, otherwise consider meanmonth.
        If seasonal==True, then all days will be written; as a masked-array if days are missing.


        Note
        ----
        If input date signifies the end of the time step, the user should remove half a time step
        before using the routine. Otherwise the routine would have to guess the time step, which
        is error prone.
        For example, remove 15 minutes in decimal days from time steps of half an hour:
        date = jams.date2dec(ascii=adate) - 15./(24.*60.)


        Examples
        --------
        >>> from autostring import astr
        >>> from date2dec import date2dec
        >>> dates  = ['01.01.1990 12:00:00', '01.02.1990 12:00:00', '01.03.1990 12:00:00',
        ...           '01.01.1990 12:10:00', '01.01.1990 13:00:00', '01.01.1990 14:00:00']
        >>> jdates = date2dec(ascii=dates)
        >>> x      = np.arange(len(jdates))+1.
        >>> odates, xout = means(jdates, x)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 3.500

        >>> odates, xout = means(jdates, x, year=True)
        >>> print(astr(xout, 3, pp=True))
        ['3.500']

        >>> odates, xout = means(jdates, x, month=True)
        >>> print(astr(xout, 3, pp=True))
        ['4.000' '2.000' '3.000']

        >>> odates, xout = means(jdates, x, day=True)
        >>> print(astr(xout, 3, pp=True))
        ['4.000' '2.000' '3.000']

        >>> odates, xout = means(jdates, x, hour=True)
        >>> print(astr(xout, 3, pp=True))
        ['2.500' '5.000' '6.000' '2.000' '3.000']

        >>> odates, xout = means(jdates, x, half_hour=True)
        >>> print(astr(xout, 3, pp=True))
        ['2.500' '5.000' '6.000' '2.000' '3.000']

        >>> odates, xout = means(jdates, x, meanday=True)
        >>> print(astr(xout[10:16], 3, pp=True))
        ['--   ' '--   ' '2.500' '5.000' '6.000' '--   ']

        >>> odates, xout = means(jdates, x, meanmonth=True)
        >>> print(astr(xout[0:5], 3, pp=True))
        ['4.000' '2.000' '3.000' '--   ' '--   ']

        >>> odates, xout = means(jdates, x, seasonal=True)
        >>> print(astr(xout[0:5], 3, pp=True))
        ['4.000' '--   ' '--   ' '--   ' '--   ']

        >>> print(astr(means(jdates, x, month=True, onlydat=True), 3, pp=True))
        ['4.000' '2.000' '3.000']

        # Masked arrays
        >>> x         = np.ma.array(x, mask=np.zeros(x.size, dtype=np.bool))
        >>> x.mask[0] = True
        >>> odates, xout = means(jdates, x)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 4.000

        >>> odates, xout = means(jdates, x, year=True)
        >>> print(astr(xout, 3, pp=True))
        ['4.000']

        >>> odates, xout = means(jdates, x, month=True)
        >>> print(astr(xout, 3, pp=True))
        ['5.000' '2.000' '3.000']

        >>> odates, xout = means(jdates, x, day=True)
        >>> print(astr(xout, 3, pp=True))
        ['5.000' '2.000' '3.000']

        # sum
        >>> odates, xout = means(jdates, x, sum=True)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 20.000
        >>> odates, xout = means(jdates, x, year=True, sum=True)
        >>> print(astr(xout, 3, pp=True))
        ['20.000']
        >>> odates, xout = means(jdates, x, month=True, sum=True)
        >>> print(astr(xout, 3, pp=True))
        ['15.000' ' 2.000' ' 3.000']
        >>> odates, xout = means(jdates, x, day=True, sum=True)
        >>> print(astr(xout, 3, pp=True))
        ['15.000' ' 2.000' ' 3.000']

        # max
        >>> odates, xout = means(jdates, x, max=True)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 6.000
        >>> odates, xout = means(jdates, x, year=True, max=True)
        >>> print(astr(xout, 3, pp=True))
        ['6.000']
        >>> odates, xout = means(jdates, x, month=True, max=True)
        >>> print(astr(xout, 3, pp=True))
        ['6.000' '2.000' '3.000']
        >>> odates, xout = means(jdates, x, day=True, max=True)
        >>> print(astr(xout, 3, pp=True))
        ['6.000' '2.000' '3.000']

        # min
        >>> odates, xout = means(jdates, x, min=True)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 2.000
        >>> odates, xout = means(jdates, x, year=True, min=True)
        >>> print(astr(xout, 3, pp=True))
        ['2.000']
        >>> odates, xout = means(jdates, x, month=True, min=True)
        >>> print(astr(xout, 3, pp=True))
        ['4.000' '2.000' '3.000']
        >>> odates, xout = means(jdates, x, day=True, min=True)
        >>> print(astr(xout, 3, pp=True))
        ['4.000' '2.000' '3.000']

        # 2D and masked arrays
        >>> x  = np.repeat(x,2).reshape((x.size,2))
        >>> odates, xout = means(jdates, x)
        >>> print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
        2.448e+06 ['4.000' '4.000']

        >>> odates, xout = means(jdates, x, year=True)
        >>> print(astr(xout, 3, pp=True))
        [['4.000' '4.000']]

        >>> odates, xout = means(jdates, x, month=True)
        >>> print(astr(xout, 3, pp=True))
        [['5.000' '5.000']
         ['2.000' '2.000']
         ['3.000' '3.000']]

        >>> odates, xout = means(jdates, x, day=True)
        >>> print(astr(xout, 3, pp=True))
        [['5.000' '5.000']
         ['2.000' '2.000']
         ['3.000' '3.000']]

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2013-2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jul 2013
        Modified, MC, Jul 2013 - meanday
                  MC, Apr 2014 - max, min
                  MC, Jun 2015 - onlydat, meanmonth
                  MC, Oct 2018 - seasonal
                  MC, Jun 2019 - bug in minute averaging: compred unique minutes to hours instead of minutes
                               - half_hour inspired by Robin Leucht for level1 data
                               - did not take keyword 'retrospective' from Robin Leucht but added a Note above for this case
    """
    # Constants
    myundef  = 9e33
    ismasked = type(dat) == np.ma.core.MaskedArray

    # Assure array
    if not ismasked: dat = np.array(dat)

    # Assure ND-array
    isone = False
    if dat.ndim == 1:
        isone = True
        dat = dat[:,np.newaxis]

    # Check options
    allopts = day+month+year+hour+half_hour+minute+meanday+meanmonth+seasonal
    assert allopts <= 1, "only one averaging option day, month, year, etc. possible"

    # Check aggregation
    allaggs = sum+max+min
    assert allaggs <= 1, "only one aggregation option sum, min, max possible"

    # average 1st dim
    if allopts == 0:
        dout = 0.5*(date[-1]+date[0])
        if sum:
            out = np.ma.sum(dat, 0)
        elif max:
            out = np.ma.amax(dat, 0)
        elif min:
            out = np.ma.amin(dat, 0)
        else:
            out = np.ma.mean(dat, 0)
        if isone:
            if onlydat:
                return out[0]
            else:
                return dout, out[0]            
        else:
            if onlydat:
                return out
            else:
                return dout, out
    else:
        yr, mo, dy, hr, mn, sc = dec2date(date)

        # year
        if year:
            yrs  = np.unique(yr)
            nout = yrs.size
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                ii = np.where(yr==i)[0]
                if np.size(ii) > 0:
                    dout[zahl] = date2dec(yr=i, mo=6, dy=15, hr=12)
                    if sum:
                        out[zahl,:] = np.ma.sum(dat[ii,:],0)
                    elif max:
                        out[zahl,:] = np.ma.amax(dat[ii,:],0)
                    elif min:
                        out[zahl,:] = np.ma.amin(dat[ii,:],0)
                    else:
                        out[zahl,:] = np.ma.mean(dat[ii,:],0)
                    zahl += 1

        # month
        if month:
            yrs  = np.unique(yr)
            mos  = np.unique(mo)
            nout = yrs.size*mos.size
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                for j in mos:
                    ii = np.where((yr==i) & (mo==j))[0]
                    if np.size(ii) > 0:
                        dout[zahl]  = date2dec(yr=i, mo=j, dy=15, hr=12)
                        if sum:
                            out[zahl,:] = np.ma.sum(dat[ii,:],0)
                        elif max:
                            out[zahl,:] = np.ma.amax(dat[ii,:],0)
                        elif min:
                            out[zahl,:] = np.ma.amin(dat[ii,:],0)
                        else:
                            out[zahl,:] = np.ma.mean(dat[ii,:],0)
                        zahl += 1

        # day
        if day:
            yrs  = np.unique(yr)
            mos  = np.unique(mo)
            dys  = np.unique(dy)
            nout = yrs.size*mos.size*dys.size
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                for j in mos:
                    for k in dys:
                        ii = np.where((yr==i) & (mo==j) & (dy==k))[0]
                        if np.size(ii) > 0:
                            dout[zahl]  = date2dec(yr=i, mo=j, dy=k, hr=12)
                            if sum:
                                out[zahl,:] = np.ma.sum(dat[ii,:],0)
                            elif max:
                                out[zahl,:] = np.ma.amax(dat[ii,:],0)
                            elif min:
                                out[zahl,:] = np.ma.amin(dat[ii,:],0)
                            else:
                                out[zahl,:] = np.ma.mean(dat[ii,:],0)
                            zahl += 1

        # hour
        if hour:
            yrs  = np.unique(yr)
            mos  = np.unique(mo)
            dys  = np.unique(dy)
            hrs  = np.unique(hr)
            nout = yrs.size*mos.size*dys.size*hrs.size
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                for j in mos:
                    for k in dys:
                        for l in hrs:
                            ii = np.where((yr==i) & (mo==j) & (dy==k) & (hr==l))[0]
                            if np.size(ii) > 0:
                                dout[zahl]  = date2dec(yr=i, mo=j, dy=k, hr=l, mi=30)
                                if sum:
                                    out[zahl,:] = np.ma.sum(dat[ii,:],0)
                                elif max:
                                    out[zahl,:] = np.ma.amax(dat[ii,:],0)
                                elif min:
                                    out[zahl,:] = np.ma.amin(dat[ii,:],0)
                                else:
                                    out[zahl,:] = np.ma.mean(dat[ii,:],0)
                                zahl += 1
        # half_hour
        if half_hour:
            yrs  = np.unique(yr)
            mos  = np.unique(mo)
            dys  = np.unique(dy)
            hrs  = np.unique(hr)
            nout = yrs.size*mos.size*dys.size*hrs.size*2
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                for j in mos:
                    for k in dys:
                        for l in hrs:
                            for m in range(2):
                                ii = np.where((yr==i) & (mo==j) & (dy==k) & (hr==l) & ((mn//30)==m))[0]
                                if np.size(ii) > 0:
                                    dout[zahl]  = date2dec(yr=i, mo=j, dy=k, hr=l, mi=15+m*30)
                                    if sum:
                                        out[zahl,:] = np.ma.sum(dat[ii,:],0)
                                    elif max:
                                        out[zahl,:] = np.ma.amax(dat[ii,:],0)
                                    elif min:
                                        out[zahl,:] = np.ma.amin(dat[ii,:],0)
                                    else:
                                        out[zahl,:] = np.ma.mean(dat[ii,:],0)
                                    zahl += 1

        # minute
        if minute:
            yrs  = np.unique(yr)
            mos  = np.unique(mo)
            dys  = np.unique(dy)
            hrs  = np.unique(hr)
            mns  = np.unique(mn)
            nout = yrs.size*mos.size*dys.size*hrs.size*mns.size
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in yrs:
                for j in mos:
                    for k in dys:
                        for l in hrs:
                            for m in mns:
                                ii = np.where((yr==i) & (mo==j) & (dy==k) & (hr==l) & (mn==m))[0]
                                if np.size(ii) > 0:
                                    dout[zahl]  = date2dec(yr=i, mo=j, dy=k, hr=l, mi=m, sc=30)
                                    if sum:
                                        out[zahl,:] = np.ma.sum(dat[ii,:],0)
                                    elif max:
                                        out[zahl,:] = np.ma.amax(dat[ii,:],0)
                                    elif min:
                                        out[zahl,:] = np.ma.amin(dat[ii,:],0)
                                    else:
                                        out[zahl,:] = np.ma.mean(dat[ii,:],0)
                                    zahl += 1

        # Mean daily
        if meanday:
            nout = 24
            hrs  = range(nout)
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in hrs:
                dout[zahl] = date2dec(yr=yr[0], mo=1, dy=1, hr=i, mi=30)
                ii = np.where(hr==i)[0]
                if np.size(ii) > 0:
                    if sum:
                        out[zahl,:] = np.ma.sum(dat[ii,:],0)
                    elif max:
                        out[zahl,:] = np.ma.amax(dat[ii,:],0)
                    elif min:
                        out[zahl,:] = np.ma.amin(dat[ii,:],0)
                    else:
                        out[zahl,:] = np.ma.mean(dat[ii,:],0)
                zahl += 1
            if np.any(out==myundef):
                out = np.ma.array(out, mask=(out==myundef), keep_mask=True)

        # Mean monthly
        if meanmonth:
            nout = 12
            mos  = range(1,nout+1)
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in mos:
                dout[zahl] = date2dec(yr=yr[0], mo=i, dy=15, hr=12, mi=0)
                ii = np.where(mo==i)[0]
                if np.size(ii) > 0:
                    if sum:
                        out[zahl,:] = np.ma.sum(dat[ii,:],0)
                    elif max:
                        out[zahl,:] = np.ma.amax(dat[ii,:],0)
                    elif min:
                        out[zahl,:] = np.ma.amin(dat[ii,:],0)
                    else:
                        out[zahl,:] = np.ma.mean(dat[ii,:],0)
                zahl += 1
            if np.any(out==myundef):
                out = np.ma.array(out, mask=(out==myundef), keep_mask=True)

        # Seasonal
        if seasonal:
            dim = np.array([[-9,31,28,31,30,31,30,31,31,30,31,30,31],
                            [-9,31,29,31,30,31,30,31,31,30,31,30,31]])
            # leaps = np.array((((yr%4)==0) & ((yr%100)!=0)) | ((yr%400)==0)).astype(np.int)
            leaps = (((yr%4)==0) & ((yr%100)!=0)) | ((yr%400)==0)
            leap  = np.any(leaps)
            ileap = int(leap)
            if leap:
                iileap = np.argmax(leaps)
                nout = 366
            else:
                iileap = 0
                nout = 365
            dys = range(1,nout+1)
            dout = np.ones(nout)*myundef
            if ismasked:
                out  = np.ma.ones([nout]+list(dat.shape[1:]))*myundef
            else:
                out  = np.ones([nout]+list(dat.shape[1:]))*myundef
            zahl = 0
            for i in range(1,13):
                for j in range(1,dim[ileap,i]+1):
                    dout[zahl] = date2dec(yr=yr[iileap], mo=i, dy=j, hr=12, mi=0)
                    ii = np.where((mo==i) & (dy==j))[0]
                    if np.size(ii) > 0:
                        if sum:
                            out[zahl,:] = np.ma.sum(dat[ii,:],0)
                        elif max:
                            out[zahl,:] = np.ma.amax(dat[ii,:],0)
                        elif min:
                            out[zahl,:] = np.ma.amin(dat[ii,:],0)
                        else:
                            out[zahl,:] = np.ma.mean(dat[ii,:],0)
                    zahl += 1
            if np.any(out==myundef):
                out = np.ma.array(out, mask=(out==myundef), keep_mask=True)

    ii = np.where(dout != myundef)[0]
    if len(ii) > 0:
        dout = dout[ii]
        out  = out[ii,:]
    else:
        raise ValueError("all values undefined after mean.")

    if isone:
        if onlydat:
            return out[:,0]
        else:
            return dout, out[:,0]
    else:
        if onlydat:
            return out
        else:
            return dout, out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # from date2dec import date2dec
    # dates  = ['01.01.1990 12:00:00', '01.02.1990 12:00:00', '01.03.1990 12:00:00',
    #          '01.01.1990 12:10:00', '01.01.1990 13:00:00', '01.01.1990 14:00:00']
    # jdates = date2dec(ascii=dates)
    # x      = np.arange(len(jdates))+1.
    # print(dates)
    # #['01.01.1990 12:00:00', '01.02.1990 12:00:00', '01.03.1990 12:00:00', '01.01.1990 12:10:00', '01.01.1990 13:00:00', '01.01.1990 14:00:00']
    # print(x)
    # #[ 1.  2.  3.  4.  5.  6.]
    # odates, xout = means(jdates, x)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 3.500
    # odates, xout = means(jdates, x, year=True)
    # print(astr(xout, 3, pp=True))
    # #['3.500']
    # odates, xout = means(jdates, x, month=True)
    # print(astr(xout, 3, pp=True))
    # #['4.000' '2.000' '3.000']
    # odates, xout = means(jdates, x, day=True)
    # print(astr(xout, 3, pp=True))
    # #['4.000' '2.000' '3.000']
    # odates, xout = means(jdates, x, meanday=True)
    # print(astr(xout[10:16], 3, pp=True))
    # #['--   ' '--   ' '2.500' '5.000' '6.000' '--   ']

    # # mask
    # x         = np.ma.array(x, mask=np.zeros(x.size, dtype=np.bool))
    # x.mask[0] = True
    # odates, xout = means(jdates, x)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 4.000
    # odates, xout = means(jdates, x, year=True)
    # print(astr(xout, 3, pp=True))
    # #['4.000']
    # odates, xout = means(jdates, x, month=True)
    # print(astr(xout, 3, pp=True))
    # #['5.000' '2.000' '3.000']
    # odates, xout = means(jdates, x, day=True)
    # print(astr(xout, 3, pp=True))
    # #['5.000' '2.000' '3.000']

    # # sum
    # odates, xout = means(jdates, x, sum=True)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 20.000
    # odates, xout = means(jdates, x, year=True, sum=True)
    # print(astr(xout, 3, pp=True))
    # #['20.000']
    # odates, xout = means(jdates, x, month=True, sum=True)
    # print(astr(xout, 3, pp=True))
    # #['15.000' ' 2.000' ' 3.000']
    # odates, xout = means(jdates, x, day=True, sum=True)
    # print(astr(xout, 3, pp=True))
    # #['15.000' ' 2.000' ' 3.000']

    # # max
    # odates, xout = means(jdates, x, max=True)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 6.000
    # odates, xout = means(jdates, x, year=True, max=True)
    # print(astr(xout, 3, pp=True))
    # #['6.000']
    # odates, xout = means(jdates, x, month=True, max=True)
    # print(astr(xout, 3, pp=True))
    # #['6.000' ' 2.000' ' 3.000']
    # odates, xout = means(jdates, x, day=True, max=True)
    # print(astr(xout, 3, pp=True))
    # #['6.000' ' 2.000' ' 3.000']

    # # min
    # odates, xout = means(jdates, x, min=True)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 2.000
    # odates, xout = means(jdates, x, year=True, min=True)
    # print(astr(xout, 3, pp=True))
    # #['2.000']
    # odates, xout = means(jdates, x, month=True, min=True)
    # print(astr(xout, 3, pp=True))
    # #['4.000' ' 2.000' ' 3.000']
    # odates, xout = means(jdates, x, day=True, min=True)
    # print(astr(xout, 3, pp=True))
    # #['4.000' ' 2.000' ' 3.000']

    # # 2D & mask
    # x  = np.repeat(x,2).reshape((x.size,2))
    # odates, xout = means(jdates, x)
    # print(astr(odates, 3, pp=True), astr(xout, 3, pp=True))
    # #2.448e+06 ['4.000' '4.000']
    # odates, xout = means(jdates, x, year=True)
    # print(astr(xout, 3, pp=True))
    # #[['4.000' '4.000']]
    # odates, xout = means(jdates, x, month=True)
    # print(astr(xout, 3, pp=True))
    # #[['5.000' '5.000']
    # # ['2.000' '2.000']
    # # ['3.000' '3.000']]
    # odates, xout = means(jdates, x, day=True)
    # print(astr(xout, 3, pp=True))
    # #[['5.000' '5.000']
    # # ['2.000' '2.000']
    # # ['3.000' '3.000']]
