#!/usr/bin/env python
import numpy as np
import scipy.optimize as opt # curve_fit, fmin, fmin_tnc

# ----------------------------------------------------------------------
def nee2gpp(dates, nee, t, isday, rg=False, vpd=False, undef=np.nan,
            method='local', shape=False, masked=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data.

        It uses either a fit of Reco vs. temperature to all nighttime data,
        or several fits over the season of Reco vs. temperature as in Reichstein et al. (2005),
        or the daytime method of Lasslop et al. (2010),
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp(dates, nee, t, isday, rg=False, vpd=False, undef=np.nan,
                    method='local', shape=False, masked=False):


        Input
        -----
        Inputs are 1D arrays that can be masked or not.
        dates         julian days
        nee           net ecosystem exchange (uptake is <0) [umol m-2 s-1]
        t             temperature [K]


        Optional Input
        --------------
        If method = 'day' | 'lasslop', extra inputs are
        rg            global radiation, i.e. shortwave down [W m-2]
        vpd           vapour pressure deficit [Pa]


        Parameters
        ----------
        undef        undefined values in data  (default: np.nan)
                     Input arrays will be masked at undef, keeping the original mask
        method       if 'global':                fit of Reco vs. temperature to all nighttime data
                     if 'local' | 'reichstein':  method of Reichstein et al. (2005)
                     if 'day'   | 'lasslop':     method of Lasslop et al. (2010)
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
        >>> from fread import * # from ufz
        >>> from date2dec import * # from ufz
        >>> dat   = fread('nee2gpp_test.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
           8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
        >>> print Reco[1120:1128]
        [ 1.78172209  1.90616886  2.07856924  2.2560362   2.46373274  2.70757535
          2.95064665  3.2184422 ]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
           8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.33166157e+00
           8.18228013e+00   1.04092252e+01   8.19395317e+00   1.08427448e+01]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
        >>> print GPP[1120:1128]
        [-- -- -- 4.49101619818 8.39556273706 10.6881053462 8.54233664766
         11.2707121977]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
        >>> print GPP[1120:1128]
        [[ -9.99900000e+03]
         [ -9.99900000e+03]
         [ -9.99900000e+03]
         [  4.49101620e+00]
         [  8.39556274e+00]
         [  1.06881053e+01]
         [  8.54233665e+00]
         [  1.12707122e+01]]
        >>> VPD = np.squeeze(dat[8,:])
        >>> vpd = np.where(VPD == undef, undef, VPD*100.)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, rg, vpd, undef=undef, method='day')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   2.89693618e+00
           6.77103400e+00   9.06351370e+00   6.95696901e+00   9.77798943e+00]
        >>> print Reco[1120:1128]
        [ 0.35174817  0.42088838  0.53124809  0.66195618  0.839204    1.0829837
          1.36527901  1.72571943]


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

        Copyright 2012 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
        Modified MC, Nov 2012 - wrapper for individual routines nee2gpp_reichstein etc.
    """

    # Global relationship in Reichstein et al. (2005)
    if method.lower() == 'global':
        return nee2gpp_global(dates, nee, t, isday, undef=undef, shape=shape, masked=masked)
    # Local relationship = Reichstein et al. (2005)
    elif ((method.lower() == 'local') | (method.lower() == 'reichstein')):
        return nee2gpp_reichstein(dates, nee, t, isday, undef=undef, shape=shape, masked=masked)
    # Lasslop et al. (2010) method
    elif ((method.lower() == 'day') | (method.lower() == 'lasslop')):
        return nee2gpp_lasslop(dates, nee, t, isday, rg, vpd, undef=undef, shape=shape, masked=masked)
    # Include new methods here
    else:
        raise ValueError('Error nee2gpp: method not implemented yet.')



# ----------------------------------------------------------------------
def nee2gpp_global(dates, nee, t, isday, undef=np.nan,
            shape=False, masked=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data, using a fit of Reco vs. temperature to all nighttime data,
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp_global(dates, nee, t, isday, undef=np.nan, shape=False, masked=False):


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


        Examples
        --------
        >>> from fread import * # from ufz
        >>> from date2dec import * # from ufz
        >>> dat   = fread('nee2gpp_test.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.33166157e+00
           8.18228013e+00   1.04092252e+01   8.19395317e+00   1.08427448e+01]


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

        Copyright 2012 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
        Modified MC, Nov 2012 - individual routine
    """

    # Checks

    # remember shape if any
    inshape = nee.shape
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    # Check squeezed shape
    if dates.ndim != 1: raise ValueError('Error nee2gpp_global: squeezed dates must be 1D array.')
    if nee.ndim   != 1: raise ValueError('Error nee2gpp_global: squeezed nee must be 1D array.')
    if t.ndim     != 1: raise ValueError('Error nee2gpp_global: squeezed t must be 1D array.')
    if isday.ndim != 1: raise ValueError('Error nee2gpp_global: squeezed isday must be 1D array.')
    ndata = dates.size
    if ((nee.size != ndata) | (t.size != ndata) | (isday.size != ndata)):
        raise ValueError('Error nee2gpp_global: inputs must have the same size.')

    # Transform to masked array with 1D mask
    matype = type(np.ma.array([1]))
    if type(nee) != matype:
        nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if nee.mask.size == 1:
            nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    if type(t) != matype:
        t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if t.mask.size == 1:
            t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    if type(isday) != matype:
        isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if isday.mask.size == 1:
            isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee))  : nee   = np.ma.array(nee,   mask=np.isnan(nee),   keep_mask=True)
        if np.ma.any(np.isnan(t))    : t     = np.ma.array(t,     mask=np.isnan(t),     keep_mask=True)
        if np.ma.any(np.isnan(isday)): isday = np.ma.array(isday, mask=np.isnan(isday), keep_mask=True)
    else:
    	if np.ma.any(nee   == undef): nee   = np.ma.array(nee,   mask=(nee==undef),   keep_mask=True)
    	if np.ma.any(t     == undef): t     = np.ma.array(t,     mask=(t==undef),     keep_mask=True)
    	if np.ma.any(isday == undef): isday = np.ma.array(isday, mask=(isday==undef), keep_mask=True)

    # Partition - Global relationship as in Reichstein et al. (2005)

    # Select valid nighttime
    mask = isday | nee.mask | t.mask | isday.mask
    ii   = np.squeeze(np.where(~mask))
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    # p, c     = opt.curve_fit(lloyd, tt, net, p0=[2.,200.]) # global parameter, global cov matrix
    p        = opt.fmin(cost_lloyd, [2.,200.], args=(tt, net), disp=False)
    Reco     = np.ones(ndata)*undef
    ii       = np.squeeze(np.where(~t.mask))
    Reco[ii] = lloyd(t[ii], p[0], p[1])

    # GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.squeeze(np.where(~(t.mask | nee.mask)))
    GPP[ii] = Reco[ii] - nee[ii]

    # Return
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



# ----------------------------------------------------------------------
def nee2gpp_reichstein(dates, nee, t, isday, rg=False, vpd=False, undef=np.nan,
            shape=False, masked=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data, using several fits of Reco vs. temperature of nighttime data
        over the season, as in Reichstein et al. (2005), in order to calculate Reco
        and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp_reichstein(dates, nee, t, isday, undef=np.nan, shape=False, masked=False):


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


        Examples
        --------
        >>> from fread import * # from ufz
        >>> from date2dec import * # from ufz
        >>> dat   = fread('nee2gpp_test.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
           8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
        >>> print Reco[1120:1128]
        [ 1.78172209  1.90616886  2.07856924  2.2560362   2.46373274  2.70757535
          2.95064665  3.2184422 ]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='Reichstein')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
           8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
        >>> print GPP[1120:1128]
        [-- -- -- 4.49101619818 8.39556273706 10.6881053462 8.54233664766
         11.2707121977]
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
        >>> print GPP[1120:1128]
        [[ -9.99900000e+03]
         [ -9.99900000e+03]
         [ -9.99900000e+03]
         [  4.49101620e+00]
         [  8.39556274e+00]
         [  1.06881053e+01]
         [  8.54233665e+00]
         [  1.12707122e+01]]


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

        Copyright 2012 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
        Modified MC, Nov 2012 - individual routine
    """

    # Checks

    # remember shape if any
    inshape = nee.shape
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    # Check squeezed shape
    if dates.ndim != 1: raise ValueError('Error nee2gpp_reichstein: squeezed dates must be 1D array.')
    if nee.ndim   != 1: raise ValueError('Error nee2gpp_reichstein: squeezed nee must be 1D array.')
    if t.ndim     != 1: raise ValueError('Error nee2gpp_reichstein: squeezed t must be 1D array.')
    if isday.ndim != 1: raise ValueError('Error nee2gpp_reichstein: squeezed isday must be 1D array.')
    ndata = dates.size
    if ((nee.size != ndata) | (t.size != ndata) | (isday.size != ndata)):
        raise ValueError('Error nee2gpp_reichstein: inputs must have the same size.')

    # Transform to masked array with 1D mask
    matype = type(np.ma.array([1]))
    if type(nee) != matype:
        nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if nee.mask.size == 1:
            nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    if type(t) != matype:
        t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if t.mask.size == 1:
            t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    if type(isday) != matype:
        isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if isday.mask.size == 1:
            isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee))  : nee   = np.ma.array(nee,   mask=np.isnan(nee),   keep_mask=True)
        if np.ma.any(np.isnan(t))    : t     = np.ma.array(t,     mask=np.isnan(t),     keep_mask=True)
        if np.ma.any(np.isnan(isday)): isday = np.ma.array(isday, mask=np.isnan(isday), keep_mask=True)
    else:
    	if np.ma.any(nee   == undef): nee   = np.ma.array(nee,   mask=(nee==undef),   keep_mask=True)
    	if np.ma.any(t     == undef): t     = np.ma.array(t,     mask=(t==undef),     keep_mask=True)
    	if np.ma.any(isday == undef): isday = np.ma.array(isday, mask=(isday==undef), keep_mask=True)


    # Partition - Local relationship = Reichstein et al. (2005)

    # Select valid nighttime
    mask = isday | nee.mask | t.mask | isday.mask
    ii   = np.squeeze(np.where(~mask))
    jul  = dates[ii]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    # 1. each 5 days, in 15 day period, fit if range of T > 5
    locp = [] # local param
    locs = [] # local err
    #locc = [] # local cov
    dmin = np.int(np.floor(np.amin(jul))) # be aware that julian days starts at noon, i.e. 1.0 is 12h
    dmax = np.int(np.ceil(np.amax(jul)))  # so the search will be from noon to noon and thus includes all nights
    for i in xrange(dmin,dmax,5):
        iii  = np.squeeze(np.where((jul>=i) & (jul<(i+14))))
        niii = iii.size
        if niii > 6:
            if (np.amax(tt[iii])-np.amin(tt[iii])) >= 5.:
                p, c = opt.curve_fit(lloyd, tt[iii], net[iii], p0=[2.,200.]) # params, covariance
                # p = opt.fmin(cost_lloyd, [2.,200.], args=(tt[iii], net[iii]), disp=False)
                res  = np.sum((net[iii]-lloyd(tt[iii],p[0],p[1]))**2)        # residuals
                s    = np.sqrt(np.diag(c * res)/np.float(niii-2))            # std err
                # s    = np.sqrt(np.diag(c * res))                           # std dev
                locp = locp + [p]
                locs = locs + [s]
                # locc = locc + [c]
    if len(locp) == 0:
        raise ValueError('Error nee2gpp_reichstein: No local relationship found.')
    locp   = np.squeeze(np.array(locp))
    locs   = np.squeeze(np.array(locs))
    # locc   = np.squeeze(np.array(locc))

    # 2. E0 = avg of best 3
    iii  = np.squeeze(np.where((locp[:,1] > 0.) & (locp[:,1] < 450.) & (np.abs(locs[:,1]/locp[:,1]) < 0.5)))
    niii = iii.size
    if niii==0:
        raise ValueError('Error nee2gpp_reichstein: No good local relationship found.')
    elif niii==1:
        bestp = locp[iii,:]
        bests = locs[iii,:]
        # cc    = locc[iii,:]
    elif niii==2:
        bestp = np.mean(locp[iii,:],axis=0)
        bests = np.mean(locs[iii,:],axis=0)
        # ls    = locs[iii,:]
        # iis   = np.argsort(ls[:,1])
        # cc    = locc[iii[iis[0]],:]
    else:
        lp    = locp[iii,:]
        ls    = locs[iii,:]
        iis   = np.argsort(ls[:,1])
        bestp = np.mean(lp[iis[0:3],:],axis=0)
        bests = np.mean(ls[iis[0:3],:],axis=0)
        # cc = locc[iii[iis[0]],:]
    # corr = cc[0,1]/np.sqrt(cc[0,0]*cc[1,1])

    # 3. Refit Rref with fixed E0, each 4 days
    refp  = [] # Rref param
    refii = [] # mean index of data points
    E0    = bestp[1]
    et    = lloyd(tt, 1., E0)
    for i in xrange(dmin,dmax,4):
        iii  = np.squeeze(np.where((jul>=i) & (jul<(i+4))))
        niii = iii.size
        if niii > 3:
            # Calc directly minisation of (nee-p*et)**2
            # p = np.sum(net[iii]*et[iii])/np.sum(et[iii]**2)
            # p, c = opt.curve_fit(lloyd_rref, et[iii], net[iii], p0=[2.])
            p     = opt.fmin(cost_lloyd_rref, [2.], args=(et[iii], net[iii]), disp=False)
            refp  = refp  + [p]
            refii = refii + [np.int((iii[0]+iii[-1])/2)]
    if len(refp) == 0:
        raise ValueError('Error nee2gpp_reichstein: No ref relationship found.')
    refp  = np.squeeze(np.array(refp))
    refii = np.squeeze(np.array(refii))

    # 4. Interpol Rref
    Rref = np.interp(dates, jul[refii], refp)

    # 5. Calc Reco
    Reco     = np.ones(ndata)*undef
    ii       = np.squeeze(np.where(~t.mask))
    Reco[ii] = lloyd(t[ii], Rref[ii], E0)

    # GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.squeeze(np.where(~(t.mask | nee.mask)))
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



# ----------------------------------------------------------------------
def nee2gpp_lasslop(dates, nee, t, isday, rg, vpd, undef=np.nan,
            shape=False, masked=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data, using the daytime method of Lasslop et al. (2010),
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp_lasslop(dates, nee, t, isday, rg, vpd, undef=np.nan,
                    shape=False, masked=False):


        Input
        -----
        Inputs are 1D arrays that can be masked or not.
        dates         julian days
        nee           net ecosystem exchange (uptake is <0) [umol m-2 s-1]
        t             temperature [K]
        rg            global radiation, i.e. shortwave down [W m-2]
        vpd           vapour pressure deficit [Pa]


        Parameters
        ----------
        undef        undefined values in data  (default: np.nan)
                     Input arrays will be masked at undef, keeping the original mask
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
        Lasslop et al. (2010)
            Separation of net ecosystem exchange into assimilation and respiration using
            a light response curve approach: critical issues and global evaluation
            Global Change Biology 16, 187-208


        Examples
        --------
        >>> from fread import * # from ufz
        >>> from date2dec import * # from ufz
        >>> dat   = fread('nee2gpp_test.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> VPD = np.squeeze(dat[8,:])
        >>> vpd = np.where(VPD == undef, undef, VPD*100.)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, rg, vpd, undef=undef, method='day')
        >>> print GPP[1120:1128]
        [ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   2.89693618e+00
           6.77103400e+00   9.06351370e+00   6.95696901e+00   9.77798943e+00]
        >>> print Reco[1120:1128]
        [ 0.35174817  0.42088838  0.53124809  0.66195618  0.839204    1.0829837
          1.36527901  1.72571943]


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

        Copyright 2012 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
        Modified MC, Nov 2012 - individual routine
    """

    # Checks

    # remember shape if any
    inshape = nee.shape
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    # Check squeezed shape
    if dates.ndim != 1: raise ValueError('Error nee2gpp_lasslop: squeezed dates must be 1D array.')
    if nee.ndim   != 1: raise ValueError('Error nee2gpp_lasslop: squeezed nee must be 1D array.')
    if t.ndim     != 1: raise ValueError('Error nee2gpp_lasslop: squeezed t must be 1D array.')
    if isday.ndim != 1: raise ValueError('Error nee2gpp_lasslop: squeezed isday must be 1D array.')
    ndata = dates.size
    if ((nee.size != ndata) | (t.size != ndata) | (isday.size != ndata)):
        raise ValueError('Error nee2gpp_lasslop: inputs must have the same size.')
    if rg.ndim  != 1: raise ValueError('Error nee2gpp_lasslop: squeezed rg must be 1D array.')
    if vpd.ndim != 1: raise ValueError('Error nee2gpp_lasslop: squeezed vpd must be 1D array.')
    if ((rg.size != ndata) | (vpd.size != ndata)):
        raise ValueError('Error nee2gpp_lasslop: lasslop inputs must have the same size as other inputs.')

    # Transform to masked array with 1D mask
    matype = type(np.ma.array([1]))
    if type(nee) != matype:
        nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if nee.mask.size == 1:
            nee = np.ma.array(nee, mask=np.zeros(ndata,dtype=np.bool))
    if type(t) != matype:
        t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if t.mask.size == 1:
            t = np.ma.array(t, mask=np.zeros(ndata,dtype=np.bool))
    if type(isday) != matype:
        isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if isday.mask.size == 1:
            isday = np.ma.array(isday, mask=np.zeros(ndata,dtype=np.bool))
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee))  : nee   = np.ma.array(nee,   mask=np.isnan(nee),   keep_mask=True)
        if np.ma.any(np.isnan(t))    : t     = np.ma.array(t,     mask=np.isnan(t),     keep_mask=True)
        if np.ma.any(np.isnan(isday)): isday = np.ma.array(isday, mask=np.isnan(isday), keep_mask=True)
    else:
    	if np.ma.any(nee   == undef): nee   = np.ma.array(nee,   mask=(nee==undef),   keep_mask=True)
    	if np.ma.any(t     == undef): t     = np.ma.array(t,     mask=(t==undef),     keep_mask=True)
    	if np.ma.any(isday == undef): isday = np.ma.array(isday, mask=(isday==undef), keep_mask=True)

    if type(rg) != matype:
        rg = np.ma.array(rg, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if rg.mask.size == 1:
            rg = np.ma.array(rg, mask=np.zeros(ndata,dtype=np.bool))
    if type(vpd) != matype:
        vpd = np.ma.array(vpd, mask=np.zeros(ndata,dtype=np.bool))
    else:
        if vpd.mask.size == 1:
            vpd = np.ma.array(vpd, mask=np.zeros(ndata,dtype=np.bool))
    # mask also undef
    if np.ma.any(rg   == undef): rg  = np.ma.array(rg,  mask=(rg==undef),  keep_mask=True)
    if np.ma.any(vpd == undef):  vpd = np.ma.array(vpd, mask=(vpd==undef), keep_mask=True)


    # Partition - Lasslop et al. (2010) method
    do_lgpp = False
    mask  = nee.mask | t.mask | isday.mask | rg.mask | vpd.mask
    # night
    nmask = isday | mask
    nii   = np.squeeze(np.where(~nmask))
    njul  = dates[nii]
    ntt   = np.ma.compressed(t[nii])
    nnet  = np.ma.compressed(nee[nii])
    aRref = np.mean(nnet)
    # day
    dmask = (~isday) | mask
    dii   = np.squeeze(np.where(~dmask))
    djul  = dates[dii]
    dtt   = np.ma.compressed(t[dii])
    dnet  = np.ma.compressed(nee[dii])
    drg   = np.ma.compressed(rg[dii])
    dvpd  = np.ma.compressed(vpd[dii])
    # starting values for optim
    aalpha = 0.01
    qnet   = np.sort(dnet)
    nqnet  = qnet.size
    abeta0 = np.abs(qnet[np.floor(0.97*nqnet)]-qnet[np.ceil(0.03*nqnet)])
    ak     = 0.
    # out
    lE0    = []
    lalpha = []
    if do_lgpp:
        lbeta0 = []
        lk     = []
    lRref  = []
    lii    = []
    dmin = np.int(np.floor(np.amin(dates)))
    dmax = np.int(np.ceil(np.amax(dates)))
    zaehl = -1
    for i in xrange(dmin,dmax,2):
        good = True
        # 1. Estimate E0 from nighttime data
        iii  = np.squeeze(np.where((njul>=i) & (njul<(i+12))))
        niii = iii.size
        if niii > 3:
            # p, c = opt.curve_fit(lloyd, ntt[iii], nnet[iii], p0=[aRref,100.])
            p  = opt.fmin(cost_lloyd, [aRref,100.], args=(ntt[iii], nnet[iii]), disp=False)
            E0 = np.maximum(p[1], 50.)
        else:
            if zaehl >= 0:
                E0 = lE0[zaehl]
            else:
                # large gap at beginning of data set, i.e. skip the period
                good = False
                continue
        # 2. Estimate alpha, k, beta0, Rref from daytime data
        iii  = np.squeeze(np.where((djul>=i) & (djul<(i+4))))
        niii = iii.size
        if niii > 3:
            et     = lloyd(dtt[iii], 1., E0)
            again  = True
            ialpha = aalpha
            ibeta0 = abeta0
            ik     = ak
            iRref  = aRref
            bounds = [[None,None],[None,None],[None,None],[None,None]]
            while again:
                again = False
                p, nfeval, rc  = opt.fmin_tnc(cost_lasslop, [ialpha,ibeta0,ik,iRref], bounds=bounds,
                                              args=(drg[iii], et, dvpd[iii], dnet[iii]),
                                              approx_grad=True, disp=False)
                # if parameters beyond some bounds, set params and redo the optim or skip
                if ((p[0] < 0.) | (p[0] > 0.22)): # alpha
                    again = True
                    if zaehl >= 0:
                        bounds[0] = [lalpha[zaehl],lalpha[zaehl]]
                        ialpha    = lalpha[zaehl]
                    else:
                        bounds[0] = [0.,0.]
                        ialpha    = 0.
                if p[1] < 0.:                    # beta0
                    bounds[1] = [0.,0.]
                    ibeta0    = 0.
                    again = True
                if p[1] > 250.:
                    good = False
                    continue
                if p[2] < 0.:                    # k
                    bounds[2] = [0.,0.]
                    ik        = 0.
                    again = True
                if p[3] < 0:                     # Rref
                    good = False
                    continue
            if good:
                lalpha = lalpha + [p[0]]
                if do_lgpp:
                    lbeta0 = lbeta0 + [p[1]]
                    lk     = lk     + [p[2]]
                lRref  = lRref  + [p[3]]
                lii    = lii    + [np.int((iii[0]+iii[-1])/2)]
            else:
                continue
        else:
            continue
        lE0    = lE0 + [E0]
        zaehl += 1
    if len(lE0) == 0:
        raise ValueError('Error nee2gpp_lasslop: No day relationship found.')
    lE0    = np.squeeze(np.array(lE0))
    if do_lgpp:
        lalpha = np.squeeze(np.array(lalpha))
        lbeta0 = np.squeeze(np.array(lbeta0))
        lk     = np.squeeze(np.array(lk))
    lRref  = np.squeeze(np.array(lRref))
    lii    = np.squeeze(np.array(lii))

    # 3. Interpol E0 and Rref
    E0   = np.interp(dates, djul[lii], lE0)
    Rref = np.interp(dates, djul[lii], lRref)

    # 4. Calc Reco
    Reco     = np.ones(ndata)*undef
    ii       = np.squeeze(np.where(~t.mask))
    Reco[ii] = lloyd(t[ii], Rref[ii], E0[ii])

    # 5. Calc GPP from light response for check
    if do_lgpp:
        alpha    = np.interp(dates, djul[lii], lE0)
        beta0    = np.interp(dates, djul[lii], lbeta0)
        k        = np.interp(dates, djul[lii], lk)
        et       = lloyd(t, 1., E0)
        lmask    = t.mask | isday.mask | rg.mask | vpd.mask
        ii       = np.squeeze(np.where(~lmask))
        lgpp     = np.zeros(ndata)
        lgpp[ii] = lasslop(rg[ii], et[ii], vpd[ii], alpha[ii], beta0[ii], k[ii], Rref[ii]) - Reco[ii]

    # GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.squeeze(np.where(~(t.mask | nee.mask)))
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



# -------------------------------------------------------------
# Fit functions
def lloyd(T, Rref, E0):
    """ Lloyd & Taylor (1994)
        T       Temperature [k]
        Rref    Respiration at Tref=10 degC [umol(C) m-2 s-1]
        E0      Activation energy [K]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return Rref*np.exp(E0*(1./(Tref-T0)-1./(T-T0)))

def cost_lloyd(p, T, NEE):
    """ Cost function for Lloyd """
    return np.sum(np.abs(NEE-lloyd(T, p[0], p[1])))



def lloyd_rref(et, Rref):
    """ If E0 is know in Lloyd & Taylor (1994) then one can calc
        the exponential term outside the routine and the fitting
        becomes linear """
    return Rref*et

def cost_lloyd_rref(p, et, NEE):
    """ Cost function for rref """
    return np.sum(np.abs(NEE-lloyd_rref(et, p[0])))



def lasslop(Rg, et, VPD, alpha, beta0, k, Rref):
    """ Lasslop et al. (2010) is basically the rectangular, hyperbolic
        light-response of NEE as by Falge et al. (2001), where the
        respiration is calculated with Lloyd & Taylor (1994), and the
        maximum canopy uptake rate at light saturation decreases
        exponentially with VPD as in Koerner (1995)
        Rg      Global radiation [W m-2]
        et      Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
        VPD     Vapour Pressure Deficit [Pa]
        alpha   Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
        beta0   Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
        k       e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
        Rref    Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = Rref*et
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.maximum(np.minimum(-k*(VPD-VPD0), 600.), -600.)
    beta  = np.where(VPD > VPD0, beta0*np.exp(kk), beta0)
    return -alpha*beta*Rg/(alpha*Rg+beta) + gamma

def cost_lasslop(p, Rg, et, VPD, NEE):
    """ Cost function for Lasslop """
    return np.sum(np.abs(NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3])))



# -------------------------------------------------------------
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # from fread import * # from ufz
    # from date2dec import * # from ufz
    # dat   = fread('nee2gpp_test.csv', skip=2, transpose=True)
    # dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    # NEE   = np.squeeze(dat[5,:])
    # rg    = np.squeeze(dat[6,:])
    # tair  = np.squeeze(dat[7,:])
    # undef = -9999.
    # isday = np.where(rg > 10., True, False)
    # tt    = np.where(tair == undef, undef, tair+273.15)
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
    # print GPP[1120:1128]
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
    # #  8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
    # print Reco[1120:1128]
    # #[ 1.78172209  1.90616886  2.07856924  2.2560362   2.46373274  2.70757535
    # #  2.95064665  3.2184422 ]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
    # print GPP[1120:1128]
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
    # #   8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
    # print GPP[1120:1128]
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.33166157e+00
    # #   8.18228013e+00   1.04092252e+01   8.19395317e+00   1.08427448e+01]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
    # print GPP[1120:1128]
    # #[-- -- -- 4.49101619818 8.39556273706 10.6881053462 8.54233664766
    # # 11.2707121977]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
    # print GPP[1120:1128]
    # #[[ -9.99900000e+03]
    # # [ -9.99900000e+03]
    # # [ -9.99900000e+03]
    # # [  4.49101620e+00]
    # # [  8.39556274e+00]
    # # [  1.06881053e+01]
    # # [  8.54233665e+00]
    # # [  1.12707122e+01]]
    # VPD = np.squeeze(dat[8,:])
    # vpd = np.where(VPD == undef, undef, VPD*100.)
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, rg, vpd, undef=undef, method='day')
    # print GPP[1120:1128]
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   2.89693618e+00
    # #   6.77103400e+00   9.06351370e+00   6.95696901e+00   9.77798943e+00]
    # print Reco[1120:1128]
    # #[ 0.35174817  0.42088838  0.53124809  0.66195618  0.839204    1.0829837
    # #  1.36527901  1.72571943]
