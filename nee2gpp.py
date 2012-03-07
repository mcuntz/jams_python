#!/usr/bin/env python
import numpy as np
import scipy.optimize as opt #curve_fit, fmin

def nee2gpp(dates, nee, t, isday, undef=np.nan, method='local', shape=False, masked=False):
# if True:
#     import ufz
#     dat   = ufz.fread('nee2gpp_test.csv', skip=2, transpose=True)
#     dates = ufz.date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
#     nee   = np.squeeze(dat[5,:])
#     rg    = np.squeeze(dat[6,:])
#     tair  = np.squeeze(dat[7,:])
#     undef = -9999.
#     method = 'local'
#     isday = np.where(rg > 10., True, False)
#     t    = np.where(tair == undef, undef, tair+273.15)
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data.

        It uses either a fit of Reco vs. temperature to all nighttime data,
        or several fits over the season of Reco vs. temperature as in Reichstein et al. (2005),
        or the daytime method of Lasslop et al. (2010) (not implemented yet),
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp(dates, nee, t, isday, undef=np.nan, method='local'):


        Input
        -----
        Inputs are 1D arrays that can be masked or not.
        dates         julian days
        nee           net ecosystem exchange (uptake is <0) [umol m-2 s-1]
        t             temperature [K]


        Parameters
        ----------
        undef        undefined values in data  (default: np.nan)
                     Input arrays will be masked at undef, keeping the original mask
        method       if 'global':  fit of Reco vs. temperature to all nighttime data
                     if 'local':   method of Reichstein et al. (2005)
                     if 'lasslop': method of Lasslop et al. (2010) (not implemented yet)
        shape        if False then outputs are 1D arrays;
                     if True, output have the same shape as datain
                     if a shape tuple is given, then this tuple is used to reshape
        masked       if False: outputs are undef where nee and t are masked or undef
                     if True:  return masked arrays where outputs would be undef


        Ouput
        -----
        GPP, Reco    photosynthesis, ecosystem respiration


        Restrictions
        ------------
        None.


        Literature
        ----------
        Reichstein et al. (2005)
            On the separation of net ecosystem exchange into assimilation and ecosystem
            respiration: review and improved algorithm.            
            Global Change Biology 11, 1424-1439

        Lasslop et al. (2010)
            Separation of net ecosystem exchange into assimilation and respiration using
            a light response curve approach: critical issues and global evaluation
            Global Change Biology 16, 187-208


        Examples
        --------
        >>> import ufz
        >>> dat   = ufz.fread('nee2gpp_test.csv', skip=2, transpose=True)
        >>> dates = ufz.date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.27224371e+00
           8.15630940e+00   1.04248029e+01   8.25499825e+00   1.09568661e+01]
        >>> print Reco[1120:1128]
        [ 1.60969988  1.72185913  1.87729631  2.03726371  2.2244794   2.44427293
          2.66330825  2.90459609]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.27224371e+00
           8.15630940e+00   1.04248029e+01   8.25499825e+00   1.09568661e+01]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.62210448e+00
           8.46924520e+00   1.06904080e+01   8.46784095e+00   1.11070112e+01]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
        >>> print GPP[1120:1128]
        [-- -- -- 4.27224370609 8.1563093998 10.4248029293 8.25499825486
         10.9568660884]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
        >>> print GPP[1120:1128]
        [[ -9.99900000e+03]
         [ -9.99900000e+03]
         [ -9.99900000e+03]
         [  4.27224371e+00]
         [  8.15630940e+00]
         [  1.04248029e+01]
         [  8.25499825e+00]
         [  1.09568661e+01]]


        History
        -------
        Written  MC, Mar 2012
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

    # remember shape is any
    inshape = np.shape(nee)
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    # Check squeezed shape
    if np.size(np.shape(dates)) != 1:
        raise ValueError('Error nee2gpp: squeezed dates must be 1D array.')
    if np.size(np.shape(nee)) != 1:
        raise ValueError('Error nee2gpp: squeezed nee must be 1D array.')
    if np.size(np.shape(t)) != 1:
        raise ValueError('Error nee2gpp: squeezed t must be 1D array.')
    if np.size(np.shape(isday)) != 1:
        raise ValueError('Error nee2gpp: squeezed isday must be 1D array.')
    ndata = np.size(dates)
    if ((np.size(nee) != ndata) | (np.size(t) != ndata) | (np.size(isday) != ndata)):
        raise ValueError('Error nee2gpp: inputs must have the same size.')

    # -------------------------------------------------------------
    # Transform to masked array with 1D mask
    matype = type(np.ma.array([1]))
    if type(nee) != matype:
        nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if np.size(nee.mask) == 1:
            nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    if type(t) != matype:
        t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if np.size(t.mask) == 1:
            t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    if type(isday) != matype:
        isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if np.size(isday.mask) == 1:
            isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    # mask also undef
    if np.ma.any(nee   == undef): nee   = np.ma.array(nee,   mask=(nee==undef),   keep_mask=True)
    if np.ma.any(t     == undef): t     = np.ma.array(t,     mask=(t==undef),     keep_mask=True)
    if np.ma.any(isday == undef): isday = np.ma.array(isday, mask=(isday==undef), keep_mask=True)


    # -------------------------------------------------------------
    # Partition

    # Select valid nighttime
    mask = isday | nee.mask | t.mask | isday.mask
    ii   = np.squeeze(np.where(~mask))
    nii  = np.size(ii)
    jul  = dates[ii]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])

    #
    # Global relationship
    if method.lower() == 'global':
        p, c = opt.curve_fit(lloyd, tt, net, p0=[2.,200.])      # global parameter, global cov matrix
        Reco     = np.ones(ndata)*undef
        ii       = np.squeeze(np.where(~t.mask))
        Reco[ii] = lloyd(t[ii], p[0], p[1])

    #
    # Local relationship
    elif method.lower() == 'local':
        # 1. each 5 days, in 15 day period, fit if range of T > 5
        locp = [] # local param
        locs = [] # local err
        #locc = [] # local cov
        dmin = np.int(np.floor(np.amin(jul))) # be aware that julian days starts at noon, i.e. 1.0 is 12h
        dmax = np.int(np.ceil(np.amax(jul)))  # so the search will be from noon to noon and thus includes all nights
        for i in xrange(dmin,dmax,5):
            iii  = np.squeeze(np.where((jul>=i) & (jul<(i+14))))
            niii = np.size(iii)
            if niii > 6:
                if (np.amax(tt[iii])-np.amin(tt[iii])) >= 5.:
                    p, c = opt.curve_fit(lloyd, tt[iii], net[iii], p0=[2.,200.]) # params, covariance
                    res  = np.sum((net[iii]-lloyd(tt[iii],p[0],p[1]))**2)        # residuals
                    s    = np.sqrt(np.diag(c * res)/np.float(niii-2))            # std err
                    #s    = np.sqrt(np.diag(c * res))                             # std dev
                    locp = locp + [p]
                    locs = locs + [s]
                    #locc = locc + [c]
        if len(locp) == 0:
            raise ValueError('No local relationship found.')
        locp   = np.squeeze(np.array(locp))
        locs   = np.squeeze(np.array(locs))
        #locc   = np.squeeze(np.array(locc))

        # 2. E0 = avg of best 3
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

        # 3. Refit Rref with fixed E0, each 4 days
        refp  = [] # Rref param
        refii = [] # mean index of data points
        E0 = bestp[1]
        et = lloyd(tt, 1., E0)
        for i in xrange(dmin,dmax,4):
            iii  = np.squeeze(np.where((jul>=i) & (jul<(i+4))))
            niii = np.size(iii)
            if niii > 3:
                p = np.sum(net[iii]*et[iii])/np.sum(et[iii]**2)
                refp  = refp  + [p]
                refii = refii + [np.int((iii[0]+iii[-1])/2)]
        if len(refp) == 0:
            raise ValueError('No ref relationship found.')
        refp   = np.squeeze(np.array(refp))
        refii  = np.squeeze(np.array(refii))

        # 4. Interpol Rref
        Rref   = np.interp(dates, jul[refii], refp)

        # 5. Calc Reco
        Reco   = np.ones(ndata)*undef
        ii     = np.squeeze(np.where(~t.mask))
        Reco[ii]  = lloyd(t[ii], Rref[ii], E0)

    #
    # Include new methods here
    else:
        raise ValueError('Error nee2gpp: method not implemented yet.')

    # GPP
    GPP  = np.ones(ndata)*undef
    ii   = np.squeeze(np.where(~(t.mask | nee.mask)))
    GPP[ii] = Reco[ii] - nee[ii]

    if masked:
        GPP  = np.ma.array(GPP,  mask=(GPP == undef))
        Reco = np.ma.array(Reco, mask=(Reco == undef))

    if shape != False:
        if shape != True:
            return np.reshape(GPP,shape), np.reshape(Reco,shape)
        else:
            return np.reshape(GPP,inshape), np.reshape(Reco,inshape)
    else:
        return GPP, Reco

if __name__ == '__main__':
    import doctest
    doctest.testmod()
