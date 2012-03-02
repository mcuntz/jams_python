#!/usr/bin/env python
import numpy as np
import scipy.optimize as opt #curve_fit, fmin

#def partition(nee, t, isday):
if True:
    import ufz
    dat   = ufz.fread('partition_test.csv', skip=2, transpose=True)
    dates = ufz.date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:]) 
    NEE   = np.squeeze(dat[7,:])
    rg    = np.squeeze(dat[10,:])
    tair  = np.squeeze(dat[11,:])
    tsoil = np.squeeze(dat[12,:])
    undef = -9999.
    ii    = np.squeeze(np.where((NEE != undef) & (tair != undef) & (rg != undef)))
    nee   = NEE[ii]
    t     = tair[ii]
    isday = np.where(rg[ii] > 10., 1, 0)
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
    # Fit functions
    def lloyd(p,t):
        tref = 227.13
        t0   = 283.15
        return p[0]*np.exp(p[1]*(1./(tref-t)-1./(t-t0)))
    
    def lloyd_fmin(p,t,nee):
        return np.sum((nee-lloyd(p,t))**2)
    
    def lloyd_amin(p,t,nee):
        return np.sum(np.abs(nee-lloyd(p,t)))


    # -------------------------------------------------------------
    # Checks

    if np.size(np.shape(nee)) != 1:
        raise ValueError('Error partition: nee must be 1D array.')
    if np.size(np.shape(t)) != 1:
        raise ValueError('Error partition: t must be 1D array.')
    if np.size(np.shape(isday)) != 1:
        raise ValueError('Error partition: isday must be 1D array.')

    ndata = np.size(nee)
    if ((np.size(t) != ndata) | (np.size(isday) != ndata)):
        raise ValueError('Error partition: inputs must have the same size.')

    # -------------------------------------------------------------
    # Partition

    ii = np.squeeze(np.where(np.array(isday,dtype=np.int) == 0))
    xx = t[ii]
    yy = nee[ii]

    # Minpack**2
    plf = opt.fmin(lloyd_fmin, np.array([2.,200.]), args=(xx,yy))
    print 'Lloyd fmin**2: ', plf

    # Minpack-abs
    pla = opt.fmin(lloyd_amin, np.array([2.,200.]), args=(xx,yy))
    print 'Lloyd abs: ', pla



if __name__ == '__main__':
    import doctest
    doctest.testmod()
