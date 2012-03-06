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
    undef = -9999.
    ii    = np.squeeze(np.where((NEE != undef) & (tair != undef) & (rg != undef)))
    jul   = dates[ii]
    nee   = NEE[ii]
    t     = tair[ii] + 273.15
    rad   = rg[ii]
    isday = np.where(rad > 10., 1, 0)
    local = True
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
    def lloyd(t, p0, p1):
        tref = 283.15
        t0   = 227.13
        return p0*np.exp(p1*(1./(tref-t0)-1./(t-t0)))

    def rref(et, p0):
        return p0*et

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

    ii  = np.squeeze(np.where(np.array(isday,dtype=np.int) == 0))
    nii = np.size(ii)
    jj  = jul[ii] # night time, undef removed julian day
    tt  = t[ii]   #      "    ,        "      temperature
    nn  = nee[ii] #      "    ,        "      nee

    if not local: # Global relationship for comparison
        p, c = opt.curve_fit(lloyd, tt, nn, p0=[2.,200.])      # global parameter, global cov matrix
        #nii  = np.size(tt)
        #res  = np.sum((nn-lloyd(tt,p[0],p[1]))**2)            # global residuals
        #s    = np.sqrt(np.diag(c * res)/np.float(nii-2))      # global param error
        #r    = c[0,1]/np.sqrt(c[0,0]*c[1,1])                  # corr between parameters
        #print 'Global curve fit: ', p, ' +- ', s, 'corr=', r
        ndates   = np.size(dates)
        Reco     = np.ones(ndates)*undef
        ii       = np.squeeze(np.where((tair != undef)))
        t        = tair[ii] + 273.15
        Reco[ii] = lloyd(t, p[0], p[1])

    else: # Local relationships
        locp = [] # local param
        locs = [] # local err
        #locc = [] # local cov
        dmin = np.int(np.floor(np.amin(jj))) # be aware that julian days starts at noon, i.e. 1.0 is 12h
        dmax = np.int(np.ceil(np.amax(jj)))  # so the search will be from noon to noon and therefore includes total nights
        for i in xrange(dmin,dmax,5):
            iii  = np.squeeze(np.where((jj>=i) & (jj<(i+14))))
            niii = np.size(iii)
            if niii > 6:
                if (np.amax(tt[iii])-np.amin(tt[iii])) >= 5.:
                    p, c = opt.curve_fit(lloyd, tt[iii], nn[iii], p0=[2.,200.])
                    res  = np.sum((nn[iii]-lloyd(tt[iii],p[0],p[1]))**2)
                    s    = np.sqrt(np.diag(c * res)/np.float(niii-2))
                    locp = locp + [p]
                    locs = locs + [s]
                    #locc = locc + [c]
        if len(locp) == 0:
            raise ValueError('No local relationship found.')
        locp   = np.squeeze(np.array(locp))
        locs   = np.squeeze(np.array(locs))
        #locc   = np.squeeze(np.array(locc))
        # best local relationship (avg of first 3 best)
        iii  = np.squeeze(np.where((locp[:,1] > 0.) & (locp[:,1] < 450.) & (np.abs(locs[:,1]/locp[:,1]) < 0.5)))
        niii = np.size(iii)
        if niii==0:
            raise ValueError('No good local relationship found.')
        elif niii==1:
            bestp = locp[iii,:]
            bests = locs[iii,:]
            #cc    = locc[iii,:]
        elif niii==2:
            bestp = np.mean(locp[iii,:],axis=0)
            bests = np.mean(locs[iii,:],axis=0)
            #ls    = locs[iii,:]
            #iis   = np.argsort(ls[:,1])
            #cc    = locc[iii[iis[0]],:]
        else:
            lp  = locp[iii,:]
            ls  = locs[iii,:]
            iis = np.argsort(ls[:,1])
            bestp = np.mean(lp[iis[0:3],:],axis=0)
            bests = np.mean(ls[iis[0:3],:],axis=0)
            #cc = locc[iii[iis[0]],:]
        #corr = cc[0,1]/np.sqrt(cc[0,0]*cc[1,1])
        #print 'Local best curve fit: ', bestp, ' +- ', bests, 'corr=', corr
                        
        refp  = [] # Rref param
        refii = [] # mean index of data points
        E0 = bestp[1]
        et = lloyd(tt, 1., E0)
        for i in xrange(dmin,dmax,5):
            iii  = np.squeeze(np.where((jj>=i) & (jj<(i+4))))
            niii = np.size(iii)
            if niii > 3:
                p = np.sum(nn[iii])/np.sum(et[iii])
                refp  = refp  + [p]
                refii = refii + [np.int((iii[0]+iii[-1])/2)]
        if len(refp) == 0:
            raise ValueError('No ref relationship found.')
        refp   = np.squeeze(np.array(refp))
        refii  = np.squeeze(np.array(refii))

        Rref   = np.interp(dates, jj[refii], refp)
        ndates = np.size(dates)
        Reco   = np.ones(ndates)*undef
        ii    = np.squeeze(np.where((tair != undef)))
        t     = tair[ii] + 273.15
        Reco[ii]  = lloyd(t, Rref[ii], E0)

    ii   = np.squeeze(np.where((NEE != undef) & (tair != undef)))
    GPP  = np.ones(ndates)*undef
    t    = tair[ii] + 273.15
    GPP[ii]  = Reco[ii] - NEE[ii]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
