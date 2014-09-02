#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from dec2date import dec2date

def gapfill(date, data, rg, tair, vpd,
            data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
            rg_dev=50., tair_dev=2.5, vpd_dev=5.,
            longgap=60, fullday=False, undef=-9999, ddof=1,
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
        def gapfill(date, data, rg, tair, vpd,
                    data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
                    rg_dev=50., tair_dev=2.5, vpd_dev=5.,
                    longgap=60, fullday=False, undef=np.nan, ddof=1,
                    err=False, shape=False):


        Input
        -----
        date        1D-array of julian days
        data        ND-array of fluxes to fill
        rg          1D-array of global radiation [W m-2]
        tair        1D-array of air Temperature [deg C]
        vpd         1D-array of vapour pressure deficit [hPa]


        Optional Input
        -------------
        data_flag     1D-array of flags of fluxes: False=good quality; True=bad data (default: False)
        rg_flag       1D-array of flags of global radiation: False=good quality; True=bad data (default: False)
        tair_flag     1D-array of flags of air temperature: False=good quality; True=bad data (default: False)
        vpd_flag      1D-array of flags of vapour pressure deficit: False=good quality; True=bad data (default: False)


        Parameters
        ----------
        rg_dev               threshold for maximum deviation of global radiation (default: 50)
        tair_dev             threshold for maximum deviation of air Temperature (default: 2.5)
        vpd_dev              threshold for maximum deviation of vpd (default: 5)
        longgap   avoid extraploation into a gap longer than longgap days (default: 60)
        fullday              If True, move beginning of large gap to start of next day
                                      and move end of large gap to end of last day (default: False)
        undef                undefined values in data  (default: -9999, np.nan not allowed)
        ddof                 Degrees of freedom to use in calculation of standard deviation
                             for error estimate (default: 1)
        err                  if True, fill every data point, i.e. used for error generation (default: False)
        shape                if False then outputs are 1D arrays;
                             if True, output have the same shape as input data;
                             if a tuple is given, then this tuple is used to reshape. (default: False)


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
        Works not with undef=np.nan


        Literature
        ----------
        Reichstein et al. (2005)
            On the separation of net ecosystem exchange into assimilation and ecosystem
            respiration: review and improved algorithm.
            Global Change Biology 11, 1424-1439


        Examples
        --------
        >>> import numpy as np
        >>> from fread import fread
        >>> from date2dec import date2dec
        >>> from autostring import astr
        >>> ifile = 'test_gapfill.csv' # Tharandt 1998 = Online tool example file
        >>> undef = -9999.
        >>> dat   = fread(ifile, skip=2, transpose=True)
        >>> ndat  = dat.shape[1]
        >>> head  = fread(ifile, skip=2, header=True)
        >>> head1 = head[0]
        >>> ihead = dict(list(zip(head1, list(range(len(head1))))))
        >>> for ii in range(len(head1)): exec(head1[ii].lower() + ' = ' + 'dat[ihead["'+head1[ii]+'"],:]')
        >>> year  = np.ones(day.shape, dtype=day.dtype) * 1998.
        >>> hh    = hour.astype(np.int)
        >>> mn    = np.round((hour-hh)*60.)
        >>> y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
        >>> jdate = y0 + day

        # 1D fill
        >>> nee_f, nee_qc = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
        >>> print(astr(nee_qc[11006:11012],0,pp=True))
        ['1' '1' '1' '2' '2' '2']
        >>> print(astr(nee_f[11006:11012],3,pp=True))
        ['-18.678' '-15.633' '-19.610' '-15.536' '-12.402' '-15.329']

        # 1D err
        >>> nee_std = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True, err=True)
        >>> print(astr(nee_std[11006:11012],3,pp=True))
        ['    5.372' '   13.118' '    6.477' '-9999.000' '-9999.000' '-9999.000']

        >>> nee_err     = np.ones(nee_std.shape, dtype=np.int)*(-1)
        >>> kk          = np.where((nee_std!=undef) & (nee_f!=0.))[0]
        >>> nee_err[kk] = np.abs(nee_std[kk]/nee_f[kk]*100.).astype(np.int)
        >>> print(astr(nee_err[11006:11012],pp=True))
        [' 28' ' 83' ' 33' ' -1' ' -1' ' -1']
    
        # 2D fill
        >>> nee2 = np.dstack([nee,nee,nee])
        >>> print(nee2.shape)
        (1, 17520, 3)
        >>> nee_f2, nee_qc2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
        >>> print(nee_f2.shape, nee_qc2.shape)
        (1, 17520, 3) (17520,)
        >>> nee_f, nee_qc   = gapfill(jdate, nee,  rg, tair, vpd, data_flag=(qcnee>1), undef=undef)
        >>> nee_f2, nee_qc2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef)
        >>> print(nee_f2.shape, nee_qc2.shape)
        (17520, 3) (17520,)
        >>> print(np.all(nee_f==nee_f2[:,0]), np.all(nee_f==nee_f2[:,1]), np.all(nee_f==nee_f2[:,2]))
        True True True

        # 2D err
        >>> nee_std  = gapfill(jdate, nee,  rg, tair, vpd, data_flag=(qcnee>1), undef=undef, err=True)
        >>> nee_std2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, err=True)
        >>> print(np.all(nee_std==nee_std2[:,0]), np.all(nee_std==nee_std2[:,1]), np.all(nee_std==nee_std2[:,2]))
        True True True


        License
        -------
        This file is part of the UFZ Python library.

        It is NOT released under the GNU Lesser General Public License, yet.

        If you use this routine, please contact Matthias Cuntz.

        Copyright 2012-2014 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2012 - modified gap_filling.py
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
                               - data ND-array
                               - longestmarginalgap was only working at beginning and end of time series
                                 renamed to longgap
                               - fullday
    """

    # -------------------------------------------------------------
    # Checks

    # check input types
    tma = np.ma.core.MaskedArray
    assert type(date) != tma, 'dates cannot be masked array.'
    if type(data) == tma: data = data.filled(undef)
    if type(rg)   == tma: rg   = rg.filled(undef)
    if type(tair) == tma: tair = tair.filled(undef)
    if type(vpd)  == tma: vpd  = vpd.filled(undef)

    # remember shape if any
    inshape = data.shape
    date = date.squeeze()
    data = data.squeeze()
    rg   = rg.squeeze()
    tair = tair.squeeze()
    vpd  = vpd.squeeze()
    assert date.ndim == 1, 'squeezed dates must be 1D array.'
    assert rg.ndim   == 1, 'squeezed rg must be 1D array.'
    assert tair.ndim == 1, 'squeezed tair must be 1D array.'
    assert vpd.ndim  == 1, 'squeezed vpd must be 1D array.'
    assert data.ndim <= 2, 'squeezed data must be 1D or 2D array.'
    
    # make 2D array
    if data.ndim == 1:
        is2d = False
        data = data[:,np.newaxis]
    else:
        is2d = True

    # check flags
    ndata  = data.shape[0]
    ndata2 = data.shape[1]
    if (data_flag != None):
        data_flg = data_flag.squeeze()
    else:
        data_flg = np.where(data[:,0] == undef, True, False)
    if (rg_flag != None):
        rg_flg = rg_flag.squeeze()
    else:
        rg_flg = np.where(rg == undef, True, False)
    if (tair_flag != None):
        tair_flg = tair_flag.squeeze()
    else:
        tair_flg = np.where(tair == undef, True, False)
    if (vpd_flag != None):
        vpd_flg = vpd_flag.squeeze()
    else:
        vpd_flg = np.where(vpd == undef, True, False)

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
    nperday = week // 7
    #hour    = (np.array(np.floor((date-np.trunc(date))*24.), dtype=np.int) + 12) % 24
    hour, mi = dec2date(date, hr=True, mi=True)
    hour    = hour + mi/60.
    day     = np.floor(date-0.5).astype(np.int)

    if err:
        # error estimate
        data_std  = np.ones(data.shape)*undef
    else:
        # array for filled values
        data_fill = data.copy()
        ii = np.where(data_flg)[0]
        if ii.size > 0: data_fill[ii,:] = undef
        # gap quality classes
        quality = np.zeros(ndata)

    #------------------------------------------------------------------
    # Large margins

    # Check for large margins at beginning
    largegap    = np.zeros(ndata, dtype=np.bool)
    firstvalid  = np.amin(np.where(~data_flg)[0])
    lastvalid   = np.amax(np.where(~data_flg)[0])
    nn          = np.int(nperday*longgap)
    if firstvalid > nn:        largegap[0:(firstvalid-nn)] = True
    if lastvalid < (ndata-nn): largegap[(lastvalid+nn):]   = True

    # search largegap - code from maskgroup.py
    index  = []
    length = []
    count  = 0
    for i in range(ndata):
        if i==0:
            if data_flg[i]:
                index += [i]
                count  = 1
        if i>0:
            if data_flg[i] and not data_flg[i-1]:
                index += [i]
                count  = 1
            elif data_flg[i]:
                count += 1
            elif not data_flg[i] and data_flg[i-1]:
                length += [count]
                count = 0
            else:
                pass
    if count>0:
        length += [count]
    # set largegap
    for i in range(len(index)):
        if length[i] > nn:
            largegap[index[i]:index[i]+length[i]] = True

    # set or unset rest of days in large gaps
    if fullday:
        for i in range(ndata-1):
            if largegap[i] and not largegap[i+1]:   # end of large margin
                largegap[np.where(day == day[i])[0]] = False
            elif not largegap[i] and largegap[i+1]: # beginning of large margin
                largegap[np.where(day == day[i])[0]] = False
            else:
                continue

    #------------------------------------------------------------------
    # Gap filling

    # flag for all meteorological conditions
    meteo_flg  = (~tair_flg) & (~vpd_flg) & (~rg_flg)
    # flag for all meteorological conditions and data
    total_flag = meteo_flg & (~data_flg)

    # Fill loop over all data points
    for j in range(ndata):
        if not err:
            # no reason to go further, no gap -> continue
            if (~(data_flg[j])) | largegap[j]: continue
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
                nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                dat = np.ma.array(data[win,:], mask=nconditions)
                if err:
                    data_std[j,:]  = np.ma.std(dat, axis=0, ddof=ddof)
                else:
                    data_fill[j,:] = np.ma.mean(dat, axis=0)
                    # assign also quality category of gap filling
                    quality[j]     = 1
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
                    nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                    dat = np.ma.array(data[win,:], mask=nconditions)
                    if err:
                        data_std[j,:]  = np.ma.std(dat, axis=0, ddof=ddof)
                    else:
                        data_fill[j,:] = np.ma.mean(dat, axis=0)
                        quality[j]     = 1
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
                nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                dat = np.ma.array(data[win,:], mask=nconditions)
                if err:
                    data_std[j,:]  = np.ma.std(dat, axis=0, ddof=ddof)
                else:
                    data_fill[j,:] = np.ma.mean(dat, axis=0)
                    quality[j]     = 1
                continue

        # If still nothing is found under similar rg within one week,
        # take the same hour within 1-7 days

        # Method 3: same hour
        enough = False
        for i in range(2):
            t_win = (nperday * (2*i+1))//2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (~(data_flg[win]))
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                dat = np.ma.array(data[win,:], mask=nconditions)
                if err:
                    data_std[j,:]  = np.ma.std(dat, axis=0, ddof=ddof)
                else:
                    data_fill[j,:] = np.ma.mean(dat, axis=0)
                    if i == 0:
                        quality[j] = 1
                    else:
                        quality[j] = 2
                break

        # sanity check
        if err:
            if data_std[j,0]  != undef: continue
        else:
            if data_fill[j,0] != undef: continue

        # If still nothing is found, start a new cycle with increased window size
        # Method 4: same as 1 but for 3-12 weeks
        if meteo_flg[j]:
            for multi in range(3,12):
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
                    nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                    dat = np.ma.array(data[win,:], mask=nconditions)
                    if err:
                        data_std[j,:]  = np.ma.std(dat, axis=0, ddof=ddof)
                    else:
                        data_fill[j,:] = np.ma.mean(dat, axis=0)
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
                if data_std[j,0]  != undef: continue
            else:
                if data_fill[j,0] != undef: continue

        # Method 5: same as 2 but for 2-12 weeks
        if (~(rg_flg[j])):
            for multi in range(2,12):
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
                    nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                    dat = np.ma.array(data[win,:], mask=nconditions)
                    if err:
                        data_std[j,:] = np.ma.std(dat, axis=0, ddof=ddof)
                    else:
                        data_fill[j,:] = np.ma.mean(dat, axis=0)
                        if multi ==0:
                            quality[j] = 1
                        elif multi <= 2:
                            quality[j] = 2
                        else:
                            quality[j] = 3
                    break
            if err:
                if data_std[j,0]  != undef: continue
            else:
                if data_fill[j,0] != undef: continue

        # Method 6: same as 3 but for 3-120 days
        for i in range(3,120):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.unique(np.sort(np.clip(jj,0,ndata-1)))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (~(data_flg[win]))
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                nconditions = np.repeat(~conditions,ndata2).reshape(conditions.size,ndata2)
                dat = np.ma.array(data[win,:], mask=nconditions)
                if err:
                    data_std[j,:] = np.ma.std(dat, axis=0, ddof=ddof)
                else:
                    data_fill[j,:] = np.ma.mean(dat, axis=0)
                    quality[j] = 3
                break

    if shape != False:
        if shape != True:
            ishape = shape
        else:
            ishape = inshape
        if err:
            return np.reshape(data_std,ishape)
        else:
            if np.prod(ishape) == quality.size: # 1D input with dimensions of 1
                return np.reshape(data_fill,ishape), np.reshape(quality,ishape)
            else:                               # 2D data but only 1D fill quality flags
                return np.reshape(data_fill,ishape), quality
    else:
        if err:
            if not is2d: data_std = data_std.squeeze()
            return data_std
        else:
            if not is2d: data_fill = data_fill.squeeze()
            return data_fill, quality


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

#     print('Read data')
#     import numpy as np
#     from fread import fread
#     from date2dec import date2dec
#     from autostring import astr
#     ifile = 'test_gapfill.csv' # Tharandt 1998 = Online tool example file
#     undef = -9999.
#     dat   = fread(ifile, skip=2, transpose=True)
#     ndat  = dat.shape[1]
#     head  = fread(ifile, skip=2, header=True)
#     head1 = head[0]
#     ihead = dict(list(zip(head1, list(range(len(head1))))))
#     for ii in range(len(head1)): exec(head1[ii].lower() + ' = ' + 'dat[ihead["'+head1[ii]+'"],:]')
#     year  = np.ones(day.shape, dtype=day.dtype) * 1998.
#     hh    = hour.astype(np.int)
#     mn    = np.round((hour-hh)*60.)
#     y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
#     jdate = y0 + day
# 
#     # 1D fill
#     nee_f, nee_qc = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
#     print(astr(nee_qc[11006:11012],0,pp=True))
#     # ['1' '1' '1' '2' '2' '2']
#     print(astr(nee_f[11006:11012],3,pp=True))
#     # ['-18.678' '-15.633' '-19.610' '-15.536' '-12.402' '-15.329']
# 
#     # 1D err
#     nee_std = gapfill(jdate, nee, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True, err=True)
#     print(astr(nee_std[11006:11012],3,pp=True))
#     # ['    5.372' '   13.118' '    6.477' '-9999.000' '-9999.000' '-9999.000']
# 
#     nee_err     = np.ones(nee_std.shape, dtype=np.int)*(-1)
#     kk          = np.where((nee_std!=undef) & (nee_f!=0.))[0]
#     nee_err[kk] = np.abs(nee_std[kk]/nee_f[kk]*100.).astype(np.int)
#     print(astr(nee_err[11006:11012],pp=True))
#     # [' 28' ' 83' ' 33' ' -1' ' -1' ' -1']
#     
#     # 2D fill
#     nee2 = np.dstack([nee,nee,nee])
#     print(nee2.shape)
#     # (1, 17520, 3)
#     nee_f2, nee_qc2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, shape=True)
#     print(nee_f2.shape, nee_qc2.shape)
#     # (1, 17520, 3) (17520,)
#     nee_f, nee_qc   = gapfill(jdate, nee,  rg, tair, vpd, data_flag=(qcnee>1), undef=undef)
#     nee_f2, nee_qc2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef)
#     print(nee_f2.shape, nee_qc2.shape)
#     # (17520, 3) (17520,)
#     print(np.all(nee_f==nee_f2[:,0]), np.all(nee_f==nee_f2[:,1]), np.all(nee_f==nee_f2[:,2]))
#     # True True True
# 
#     # 2D err
#     nee_std  = gapfill(jdate, nee,  rg, tair, vpd, data_flag=(qcnee>1), undef=undef, err=True)
#     nee_std2 = gapfill(jdate, nee2, rg, tair, vpd, data_flag=(qcnee>1), undef=undef, err=True)
#     print(np.all(nee_std==nee_std2[:,0]), np.all(nee_std==nee_std2[:,1]), np.all(nee_std==nee_std2[:,2]))
#     # True True True
