#!/usr/bin/env python

#from sys import float_info
import numpy as np
#import pdb
from dec2date import * # from ufz

###############################################################
###############################################################
###############################################################
def gapfill(date, data, rg, tair, vpd,
            data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
            rg_dev=50., tair_dev=2.5, vpd_dev=5,
            longestmarginalgap=60, undef=-9999.,
            err=False):
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
                    rg_dev=50., tair_dev=2.5, vpd_dev=5,
                    longestmarginalgap=60, undef=-9999.,
                    err=False):


        Input
        -----
        date          julian days
        data          fluxes to fill
        rg            global radiation [W m-2]
        tair          air Temperature [deg C]
        vpd           vapour pressure deficit [hPa]


        Optional Input
        -------------
        data_flag     flags of fluxes: 0=good quality; >0=bad data (default: 0)
        rg_flag       flags of global radiation: 0=good quality; >0=bad data (default: 0)
        tair_flag     flags of air temperature: 0=good quality; >0=bad data (default: 0)
        vpd_flag      flags of vapour pressure deficit: 0=good quality; >0=bad data (default: 0)


        Parameters
        ----------
        rg_dev               threshold for maximum deviation of global radiation (default: 50)
        tair_dev             threshold for maximum deviation of air Temperature (default: 2.5)
        vpd_dev              threshold for maximum deviation of vpd (default: 5)
        longestmarginalgap   avoid extraploation into a gap longer than longestmarginalgap days (default: 60)
        undef                undefined values in data  (default: -9999.)
        err                  if True, fill every data point, i.e. used for error generation (default: False)


        Ouput
        -----
        if err=False:
            filled_data, quality_class
        else:
            err_data


        Restrictions
        ------------
        if err=True, there is no error estimate if there are no meteorological
        conditions in the vicinity of the data point.


        Literature
        ----------
        Reichstein et al. (2005) On the separation of net ecosystem exchange into
        assimilation and ecosystem respiration: review and improved algorithm.
	Global Change Biology,11,9, p. 1424-1439.


        Examples
        --------
        >>> print 'To be done!'


        History
        -------
        Written  MC, Mar 2012 - modified gap_filling.py
    """

    # -------------------------------------------------------------
    # Checks

    if np.size(np.shape(date)) != 1:
        raise ValueError('Error gapfill: dates must be 1D array.')
    if np.size(np.shape(data)) != 1:
        raise ValueError('Error gapfill: data must be 1D array.')
    if np.size(np.shape(rg)) != 1:
        raise ValueError('Error gapfill: rg must be 1D array.')
    if np.size(np.shape(tair)) != 1:
        raise ValueError('Error gapfill: tair must be 1D array.')
    if np.size(np.shape(vpd)) != 1:
        raise ValueError('Error gapfill: vpd must be 1D array.')

    # check flags
    if (data_flag != None):
        data_flg = data_flag
    else:
        data_flg = np.where(data == undef, 1, 0)
    if (rg_flag != None):
        rg_flg = rg_flag
    else:
        rg_flg = np.zeros(ndata)
    if (tair_flag != None):
        tair_flg = tair_flag
    else:
        tair_flg = np.zeros(ndata)
    if (vpd_flag != None):
        vpd_flg = vpd_flag
    else:
        vpd_flg = np.zeros(ndata)

    # -------------------------------------------------------------
    # Parameters

    tperday = 48
    
    # epsilon: Machine variable for floating point calculations
    eps = 1e-9

    # number of data points per week; basic factor of the time
    # window
    week = tperday * 7
    ddate = date-np.roll(date,1)
    ddate[0] = ddate[1]
    if np.any((ddate-ddate[0]) > eps):
        raise ValueError('Error gapfill: dates must be equally spaced.')
    print week
    week = np.int(np.around(7./ddate[0]))
    print week
    print tperday
    nperday = week / 7
    print nperday

    (years,months,days,hours,minutes,sc) = dec2date(date,\
            fulldate=True)
    hour = hours + minutes/60
    print hour
    hour = (date-np.trunc(date))*24.
    print hour

    ndata = np.size(data)
    if err:
        # array for filled values
        data_fill = np.ones(ndata)*undef
        # gap quality classes
        quality   = np.zeros(ndata)
    else:
        # error estimate
        data_std  = np.ones(ndata)*undef

    #------------------------------------------------------------------
    # Gap filling

    # flag for all meteorological conditions
    meteo_flg  = (tair_flg == 0) & (vpd_flg == 0) & (rg_flg == 0)
    # flag for all meteorological conditions and data
    total_flag = meteo_flg & (data_flg == 0)

    # Check for large margins at beginning
    largemargin = np.zeros(ndata)
    firstvalid  = np.amin(np.squeeze(np.where(data_flg==0)))
    lastvalid   = np.amax(np.squeeze(np.where(data_flg==0)))
    nn          = nperday*longestmarginalgap
    if firstvalid > nn:
        largemargin[0:(firstvalid-nn)] = 1
    if lastvalid < (nadata-nn):
        largemargin[(lastvalid+nn):]   = 1

    # Fill loop over all data points
    for j in xrange(nadata):
        if not err:
            # no reason to go further, no gap -> continue
            if data_flg[j] == 0 or largemargin[j] == 1: continue
        # 3 Methods
        #   1. tair, vpd and global radiation;
        #   2. just global radiation;
        #   3. none meteorolgical conditions: take the mean of +- hour

        # for better overview: dynamic calculation of radiation threshold
        # minimum 20; maximum 50 [Wm-2] according to private correspondence
        # with Markus Reichstein
        rg_devmax = np.maximum(20,np.minimum(rg[j],rg_dev))

        # Method 1: all met conditions
        if meteo_flg[j]:
            # search for values around the met-conditions in a window of time
            # (one week in the first iteration and odd weeks in the next)
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,nadata-1))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(np.abs(rg[win]  -rg[j])  -rg_devmax) < eps) &
                           (np.abs(np.abs(tair[win]-tair[j])-tair_dev)  < eps) &
                           (np.abs(np.abs(vpd[win] -vpd[j]) -vpd_dev)   < eps) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    # assign also quality category of gap filling
                    quality[j]   = 1
                continue
            else: # --> extend time window to two weeks
                j1  = j - np.arange(1,2*week+1) + 1
                j2  = j + np.arange(1,2*week)
                jj  = np.append(j1,j2)
                win = np.sort(np.clip(jj,0,nadata-1))
                conditions = ( (np.abs(np.abs(rg[win]  -rg[j])  -rg_devmax) < eps) &
                               (np.abs(np.abs(tair[win]-tair[j])-tair_dev)  < eps) &
                               (np.abs(np.abs(vpd[win] -vpd[j]) -vpd_dev)   < eps) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        quality[j]   = 2
                    continue

        # if you come here, no error estimate
        if err: continue

        # If nothing is found under similar meteo within two weeks,
        # look for global radiation within one week ->

        # Method 2: just global radiation available
        if rg_flg[j] == 0:
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,nadata-1))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(np.abs(rg[win]  -rg[j])  -rg_devmax) < eps) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j]   = 3
                continue

        # If still nothing is found under similar rg within one week,
        # take the same hour within 1-7 days

        # Method 3: same hour
        for i in xrange(3):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,nadata-1))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (data_flg[win] == 0)
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j]   = 4
                continue

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
                win = np.sort(np.clip(jj,0,nadata-1))
                conditions = ( (np.abs(np.abs(rg[win]  -rg[j])  -rg_devmax) < eps) &
                               (np.abs(np.abs(tair[win]-tair[j])-tair_dev)  < eps) &
                               (np.abs(np.abs(vpd[win] -vpd[j]) -vpd_dev)   < eps) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat)
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
        if rg_flg[j] == 0:
            for multi in xrange(2,12):
                j1  = j - np.arange(1,multi*week+1) + 1
                j2  = j + np.arange(1,multi*week)
                jj  = np.append(j1,j2)
                win = np.sort(np.clip(jj,0,nadata-1))
                # get boolean array where meteo-conditions are in a given width
                conditions = ( (np.abs(np.abs(rg[win]  -rg[j])  -rg_devmax) < eps) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat)
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

        # Method 6: same as 3 but for 3-120 weeks
        for i in range(3,120):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,nadata-1))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (data_flg[win] == 0)
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j] = 3
                break

    if err:
        return data_std
    else:
        return data_fill, quality


if __name__ == '__main__':
    import doctest
    doctest.testmod()
