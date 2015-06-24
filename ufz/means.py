#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.date2dec import date2dec
from ufz.dec2date import dec2date

def means(date, dat, year=False, month=False, day=False, hour=False, minute=False,
          meanday=False, meanmonth=False, sum=False, max=False, min=False, onlydat=False):
    """
        Calculate daily, monthly, yearly, etc. means of data depending on date stamp.

        date is Julian date. dat will be averaged over columns, i.e. first dimension.
        If no option is given, the mean will be over the whole first column.

        Returns centred dates, averages.
            Yearly   dates are centred at June 15, 12:00h.
            Monthly  dates are centred at 15th, 12:00h.
            Daily    dates are centred at 12:00h.
            Hourly   dates are centred at 30 min.
            Minutely dates are centred at 30 sec.
            Mean daily dates centred on 30 min of 01. January of first year.
            Mean monthly dates centred on 15th of month, 12:00h of first year.


        Definition
        ----------
        def means(date, dat, year=False, month=False, day=False, hour=False, minute=False,
                  meanday=False, meanmonth=False, sum=False, max=False, min=False, onlydat=False):


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
        minute    if True, minutely means.
        meanday   if True, mean daily cycle.
        meanmonth if True, mean monthly cycle.
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

        >>> odates, xout = means(jdates, x, meanday=True)
        >>> print(astr(xout[10:16], 3, pp=True))
        ['--   ' '--   ' '2.500' '5.000' '6.000' '--   ']

        >>> odates, xout = means(jdates, x, meanmonth=True)
        >>> print(astr(xout[0:5], 3, pp=True))
        ['4.000' '2.000' '3.000' '--   ' '--   ']

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

        Copyright 2013-2014 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2013
        Modified, MC, Jul 2013 - meanday
                  MC, Apr 2014 - max, min
                  MC, Jun 2015 - onlydat, meanmonth
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
    allopts = day+month+year+hour+minute+meanday+meanmonth
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
                                ii = np.where((yr==i) & (mo==j) & (dy==k) & (hr==l) & (hr==m))[0]
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
            hrs  = np.arange(24)
            nout = hrs.size
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
            mos  = np.arange(12)+1
            nout = mos.size
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
