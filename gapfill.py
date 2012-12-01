#!/usr/bin/env python
import numpy as np
from dec2date import *

def gapfill(datein, datain, rgin, tairin, vpdin,
            data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
            rg_dev=50., tair_dev=2.5, vpd_dev=5.,
            longestmarginalgap=60, undef=np.nan, ddof=1,
            err=False, shape=False):
    """
        Fills gaps of flux data from Eddy covariance measurements according to
        Reichstein et al. (2005).
        If there is a gap in the data, look for similar meteorological conditions
        (defined as maximum possible deviations) in a certain time window and fill
        with the average of these 'similar' values.

        The routine can also do the same search for similar meteorological conditions
        for every data point and calculate its standard deviation as a measure of uncertainty.


        Definition
        ----------
        def gapfill(datein, datain, rgin, tairin, vpdin,
                    data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
                    rg_dev=50., tair_dev=2.5, vpd_dev=5,
                    longestmarginalgap=60, undef=np.nan, ddof=1,
                    err=False, shape=False):


        Input
        -----
        datein        julian days
        datain        fluxes to fill
        rgin          global radiation [W m-2]
        tairin        air Temperature [deg C]
        vpdin         vapour pressure deficit [hPa]


        Optional Input
        -------------
        data_flag     flags of fluxes: False=good quality; True=bad data (default: False)
        rg_flag       flags of global radiation: False=good quality; True=bad data (default: False)
        tair_flag     flags of air temperature: False=good quality; True=bad data (default: False)
        vpd_flag      flags of vapour pressure deficit: False=good quality; True=bad data (default: False)


        Parameters
        ----------
        rg_dev               threshold for maximum deviation of global radiation (default: 50)
        tair_dev             threshold for maximum deviation of air Temperature (default: 2.5)
        vpd_dev              threshold for maximum deviation of vpd (default: 5)
        longestmarginalgap   avoid extraploation into a gap longer than longestmarginalgap days (default: 60)
        undef                undefined values in data  (default: np.nan)
        ddof                 Degrees of freedom tu use in calculation of standard deviation
                             for error estimate (default: 1)
        err                  if True, fill every data point, i.e. used for error generation (default: False)
        shape                if False then outputs are 1D arrays;
                             if True, output have the same shape as datain
                             if a shape tuple is given, then this tuple is used to reshape


        Ouput
        -----
        if err=False:
            filled_data, quality_class
        else:
            err_data


        Restrictions
        ------------
        If err=True, there is no error estimate if there are no meteorological
          conditions in the vicinity of the data point.


        Literature
        ----------
        Reichstein et al. (2005)
            On the separation of net ecosystem exchange into assimilation and ecosystem
            respiration: review and improved algorithm.
            Global Change Biology 11, 1424-1439


        Examples
        --------
        >>> import numpy as np
        >>> from fread import *
        >>> from date2dec import *
        >>> ifile = 'nee2gpp_test.csv'
        >>> undef = -9999.
        >>> dat   = fread(ifile, skip=2, transpose=True)
        >>> ndat  = dat.shape[1]
        >>> head  = fread(ifile, skip=2, header=True)
        >>> head1 = head[0]
        >>> ihead = dict(zip(head1, range(len(head1))))
        >>> for ii in xrange(len(head1)): exec(head1[ii].lower() + ' = ' + 'dat[ihead["'+head1[ii]+'"],:]')
        >>> year  = np.ones(day.shape, dtype=day.dtype) * 1998.
        >>> hh    = hour.astype(np.int)
        >>> mn    = np.round((hour-hh)*60.)
        >>> y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
        >>> jdate = y0 + day

        >>> nee_f, nee_qc = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
        >>> print nee_qc[11000:11020]
        [ 1.  1.  1.  1.  1.  1.  1.  1.  1.  2.  2.  2.  2.  2.  2.  2.  2.  2.
          2.  2.]
        >>> print nee_f[11000:11020]
        [  2.14         2.76333333   1.87666667  -2.765       -6.59        -8.63454545
         -18.6775     -15.63333333 -19.61       -15.536      -12.4025     -15.32857143
         -14.255      -14.52333333 -13.95       -14.155      -12.90666667
         -12.90666667 -12.725      -16.3       ]

        >>> nee_std = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True, err=True)
        >>> print nee_std[11000:11020]
        [  1.64647097e+00   1.31720664e+00   1.86345110e+00   4.14701901e+00
           2.23898638e+00   3.29813997e+00   5.37231406e+00   1.31184234e+01
           6.47709812e+00  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03
          -9.99900000e+03  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03
          -9.99900000e+03  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03]

        >>> nee_err     = np.ones(nee_std.shape, dtype=np.int)*(-1)
        >>> kk          = np.where((nee_std!=undef) & (nee_f!=0.))[0]
        >>> nee_err[kk] = np.abs(nee_std[kk]/nee_f[kk]*100.).astype(np.int)
        >>> print nee_err[11000:11020]
        [ 76  47  99 149  33  38  28  83  33  -1  -1  -1  -1  -1  -1  -1  -1  -1
          -1  -1]


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

        Copyright 2012 Matthias Cuntz


        History
        -------
        Written  MC, Mar 2012 - modified gap_filling.py
    """

    # -------------------------------------------------------------
    # Checks
    # remember shape is any
    inshape = datain.shape
    date = np.squeeze(datein)
    data = np.squeeze(datain)
    rg   = np.squeeze(rgin)
    tair = np.squeeze(tairin)
    vpd  = np.squeeze(vpdin)
    if np.ndim(date) != 1: raise ValueError('squeezed dates must be 1D array.')
    if np.ndim(data) != 1: raise ValueError('squeezed data must be 1D array.')
    if np.ndim(rg)   != 1: raise ValueError('squeezed rg must be 1D array.')
    if np.ndim(tair) != 1: raise ValueError('squeezed tair must be 1D array.')
    if np.ndim(vpd)  != 1: raise ValueError('squeezed vpd must be 1D array.')

    # check flags
    ndata = data.size
    if (data_flag != None):
        data_flg = np.squeeze(data_flag)
    else:
        data_flg = np.where(np.squeeze(data) == undef, True, False)
    if (rg_flag != None):
        rg_flg = np.squeeze(rg_flag)
    else:
        rg_flg = np.where(np.squeeze(rg) == undef, True, False)
    if (tair_flag != None):
        tair_flg = np.squeeze(tair_flag)
    else:
        tair_flg = np.where(np.squeeze(tair) == undef, True, False)
    if (vpd_flag != None):
        vpd_flg = np.squeeze(vpd_flag)
    else:
        vpd_flg = np.where(np.squeeze(vpd) == undef, True, False)

    if ((date.size != ndata) | (rg.size != ndata) | (tair.size != ndata) |
        (vpd.size != ndata) | (data_flg.size != ndata) | (rg_flg.size != ndata) |
        (tair_flg.size != ndata) | (vpd_flg.size != ndata)):
        raise ValueError('inputs must have the same size.')
        
    # -------------------------------------------------------------
    # Parameters

    # number of data points per week; basic factor of the time
    # window
    ddate    = np.diff(date)
    if np.any((ddate-ddate[0]) > 1e-4): # 1e-4 are ca 2.secs
        raise ValueError('dates must be equally spaced.')
    week    = np.int(np.around(7./ddate[0]))
    nperday = week / 7
    #hour    = (np.array(np.floor((date-np.trunc(date))*24.), dtype=np.int) + 12) % 24
    hour, mi = dec2date(date, hr=True, mi=True)
    hour     = hour + mi/60.

    ndata = data.size
    if err:
        # error estimate
        data_std  = np.ones(ndata)*undef
    else:
        # array for filled values
        data_fill = np.where(~data_flg, data, undef)
        # gap quality classes
        quality   = np.zeros(ndata)

    #------------------------------------------------------------------
    # Gap filling

    # flag for all meteorological conditions
    meteo_flg  = (~tair_flg) & (~vpd_flg) & (~rg_flg)
    # flag for all meteorological conditions and data
    total_flag = meteo_flg & (~data_flg)

    # Check for large margins at beginning
    largemargin = np.zeros(ndata)
    firstvalid  = np.amin(np.squeeze(np.where(~data_flg)))
    lastvalid   = np.amax(np.squeeze(np.where(~data_flg)))
    nn          = nperday*longestmarginalgap
    if firstvalid > nn:        largemargin[0:(firstvalid-nn)] = 1
    if lastvalid < (ndata-nn): largemargin[(lastvalid+nn):]   = 1

    # Fill loop over all data points
    for j in xrange(ndata):
        if not err:
            # no reason to go further, no gap -> continue
            if (~(data_flg[j])) | (largemargin[j] == 1): continue
        # 3 Methods
        #   1. tair, vpd and global radiation;
        #   2. just global radiation;
        #   3. none meteorolgical conditions: take the mean of +- hour

        # for better overview: dynamic calculation of radiation threshold
        # minimum 20; maximum 50 [Wm-2] according to private correspondence
        # with Markus Reichstein
        rg_devmax = np.maximum(20.,np.minimum(rg[j],rg_dev))

        # Method 1: all met conditions
        if meteo_flg[j]:
            # search for values around the met-conditions in a window of time
            # (one week in the first iteration and odd weeks in the next)
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                           (np.abs(tair[win]-tair[j]) < tair_dev) &
                           (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    # assign also quality category of gap filling
                    quality[j]   = 1
                continue
            else: # --> extend time window to two weeks
                j1  = j - np.arange(1,2*week+1) + 1
                j2  = j + np.arange(1,2*week)
                jj  = np.append(j1,j2)
                win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
                conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                               (np.abs(tair[win]-tair[j]) < tair_dev) &
                               (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        quality[j]   = 1
                    continue

        # if you come here, no error estimate
        if err: continue

        # If nothing is found under similar meteo within two weeks,
        # look for global radiation within one week ->

        # Method 2: just global radiation available
        if (~(rg_flg[j])):
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(rg[win]-rg[j]) < rg_devmax) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j]   = 1
                continue

        # If still nothing is found under similar rg within one week,
        # take the same hour within 1-7 days

        # Method 3: same hour
        enough = False
        for i in xrange(2):
            t_win = (nperday * (2*i+1))/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (~(data_flg[win]))
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    if i == 0:
                        quality[j] = 1
                    else:
                        quality[j] = 2
                break

        # sanity check
        if err:
            if data_std[j]  != undef: continue
        else:
            if data_fill[j] != undef: continue

        # If still nothing is found, start a new cycle with increased window size
        # Method 4: same as 1 but for 3-12 weeks
        if meteo_flg[j]:
            for multi in xrange(3,12):
                j1  = j - np.arange(1,multi*week+1) + 1
                j2  = j + np.arange(1,multi*week)
                jj  = np.append(j1,j2)
                win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
                conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                               (np.abs(tair[win]-tair[j]) < tair_dev) &
                               (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        # assign also quality category of gap filling
                        if multi <= 2:
                            quality[j] = 1
                        elif multi > 4:
                            quality[j] = 3
                        else:
                            quality[j] = 2
                    break
            # Check because continue does not support to jump out of two loops
            if err:
                if data_std[j]  != undef: continue
            else:
                if data_fill[j] != undef: continue

        # Method 5: same as 2 but for 2-12 weeks
        if (~(rg_flg[j])):
            for multi in xrange(2,12):
                j1  = j - np.arange(1,multi*week+1) + 1
                j2  = j + np.arange(1,multi*week)
                jj  = np.append(j1,j2)
                win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
                # get boolean array where meteo-conditions are in a given width
                conditions = ( (np.abs(rg[win]  -rg[j]) < rg_devmax) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        if multi ==0:
                            quality[j] = 1
                        elif multi <= 2:
                            quality[j] = 2
                        else:
                            quality[j] = 3
                    break
            if err:
                if data_std[j]  != undef: continue
            else:
                if data_fill[j] != undef: continue

        # Method 6: same as 3 but for 3-120 days
        for i in xrange(3,120):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (~(data_flg[win]))
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j] = 3
                break

    if shape != False:
        if shape != True:
            if err:
                return np.reshape(data_std,shape)
            else:
                return np.reshape(data_fill,shape), np.reshape(quality,shape)
        else:
            if err:
                return np.reshape(data_std,inshape)
            else:
                return np.reshape(data_fill,inshape), np.reshape(quality,inshape)
    else:
        if err:
            return data_std
        else:
            return data_fill, quality


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # print 'Read data'
    # from fread import *
    # from date2dec import *
    # ifile = 'nee2gpp_test.csv'
    # undef = -9999.
    # dat  = fread(ifile, skip=2, transpose=True)
    # ndat = dat.shape[1]
    # # Create variables with names from the header
    # head  = fread(ifile, skip=2, header=True)
    # head1 = head[0]
    # ihead = dict(zip(head1, range(len(head1))))
    # for ii in xrange(len(head1)):
    #     exec(head1[ii].lower() + ' = ' + 'dat[ihead["'+head1[ii]+'"],:]')
    # # Date
    # year  = np.ones(day.shape, dtype=day.dtype) * 1998.
    # hh    = hour.astype(np.int)
    # mn    = np.round((hour-hh)*60.)
    # y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
    # jdate = y0 + day
    # print 'Fill data'
    # nee_f, nee_qc = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
    # print nee_qc[11000:11020]
    # #[ 1.  1.  1.  1.  1.  1.  1.  1.  1.  2.  2.  2.  2.  2.  2.  2.  2.  2.
    # #  2.  2.]
    # print nee_f[11000:11020]
    # #[  2.14         2.76333333   1.87666667  -2.765       -6.59        -8.63454545
    # # -18.6775     -15.63333333 -19.61       -15.536      -12.4025     -15.32857143
    # # -14.255      -14.52333333 -13.95       -14.155      -12.90666667
    # # -12.90666667 -12.725      -16.3       ]

    # nee_std = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True, err=True)
    # print nee_std[11000:11020]
    # #[  1.64647097e+00   1.31720664e+00   1.86345110e+00   4.14701901e+00
    # #   2.23898638e+00   3.29813997e+00   5.37231406e+00   1.31184234e+01
    # #   6.47709812e+00  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03
    # #  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03
    # #  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03  -9.99900000e+03]

    # nee_err     = np.ones(nee_std.shape, dtype=np.int)*(-1)
    # kk          = np.where((nee_std!=undef) & (nee_f!=0.))[0]
    # nee_err[kk] = np.abs(nee_std[kk]/nee_f[kk]*100.).astype(np.int)
    # print nee_err[11000:11020]
    # #[ 76  47  99 149  33  38  28  83  33  -1  -1  -1  -1  -1  -1  -1  -1  -1
    # #  -1  -1]
