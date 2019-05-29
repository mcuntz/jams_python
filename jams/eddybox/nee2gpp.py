#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.optimize as opt # curve_fit, fmin, fmin_tnc
import jams.functions as functions # from jams
from jams.mad import mad # from jams
import warnings
# import pdb

# ----------------------------------------------------------------------
def nee2gpp(dates, nee, t, isday, rg=False, vpd=False, undef=np.nan,
            method='reichstein', shape=False, masked=False, nogppnight=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data.

        It uses either
          1. a fit of Reco vs. temperature to all nighttime data, or
          2. several fits over the season of Reco vs. temperature as in Reichstein et al. (2005), or
          3. the daytime method of Lasslop et al. (2010),
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp(dates, nee, t, isday, rg=False, vpd=False, undef=np.nan,
                    method='reichstein', shape=False, masked=False):


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
        method       if 'global' | 'falge':      fit of Reco vs. temperature to all nighttime data
                     if 'local'  | 'reichstein': method of Reichstein et al. (2005)
                     if 'day'    | 'lasslop':    method of Lasslop et al. (2010)
        shape        if False then outputs are 1D arrays;
                     if True, output have the same shape as datain
                     if a shape tuple is given, then this tuple is used to reshape
        masked       if False: outputs are undef where nee and t are masked or undef
                     if True:  return masked arrays where outputs would be undef

        If method = 'night' | 'reichstein', extra parameters are
        nogppnight   if True:  Resp=NEE, GPP=0 at night, GPP always positive
                     if False: Resp=lloyd_taylor, GPP=Resp-NEE at night (default)


        Ouput
        -----
        GPP, Reco    photosynthesis, ecosystem respiration


        Restrictions
        ------------
        Negative respiration possible at night when gpp is forced to 0 with nogppnight=True


        Literature
        ----------
        Falge et al. (2001)
            Gap filling strategies for defensible annual sums of net ecosystem exchange
            Acricultural and Forest Meteorology 107, 43-69

        Lasslop et al. (2010)
            Separation of net ecosystem exchange into assimilation and respiration using
            a light response curve approach: critical issues and global evaluation
            Global Change Biology 16, 187-208

        Reichstein et al. (2005)
            On the separation of net ecosystem exchange into assimilation and ecosystem
            respiration: review and improved algorithm.
            Global Change Biology 11, 1424-1439


        Examples
        --------
        >>> from jams.fread import fread # from jams
        >>> from jams.date2dec import date2dec # from jams
        >>> dat   = fread('test_nee2gpp.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> from jams.autostring import astr
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.491' '    8.396' '   10.688' '    8.542' '   11.271']
        >>> print(astr(Reco[1120:1128],3,pp=True))
        ['1.782' '1.906' '2.079' '2.256' '2.464' '2.708' '2.951' '3.218']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.491' '    8.396' '   10.688' '    8.542' '   11.271']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='reichstein')
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.491' '    8.396' '   10.688' '    8.542' '   11.271']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['--    ' '--    ' '--    ' ' 4.491' ' 8.396' '10.688' ' 8.542' '11.271']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
        >>> print(astr(GPP[1120:1128],3,pp=True))
        [['-9999.000']
         ['-9999.000']
         ['-9999.000']
         ['    4.491']
         ['    8.396']
         ['   10.688']
         ['    8.542']
         ['   11.271']]
        >>> VPD = np.squeeze(dat[8,:])
        >>> vpd = np.where(VPD == undef, undef, VPD*100.)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, rg, vpd, undef=undef, method='day')
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    2.897' '    6.771' '    9.064' '    6.957' '    9.778']
        >>> print(astr(Reco[1120:1128],3,pp=True))
        ['0.352' '0.421' '0.531' '0.662' '0.839' '1.083' '1.365' '1.726']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2014 Matthias Cuntz, Arndt Piayda - mc (at) macu (dot) de

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

        Copyright 2012-2013 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
                 MC, Nov 2012 - wrapper for individual routines nee2gpp_reichstein etc.
                 MC, Feb 2013 - ported to Python 3
                 MC, May 2013 - replaced cost functions by generel cost function cost_abs if possible
                 AP, Aug 2014 - replaced fmin with fmin_tnc to permit params<0,
                                permit gpp<0 at any time if nogppnight=True 
    """

    # Global relationship in Reichstein et al. (2005)
    if ((method.lower() == 'global') | (method.lower() == 'falge')):
        return nee2gpp_falge(dates, nee, t, isday, undef=undef, shape=shape, masked=masked)
    # Local relationship = Reichstein et al. (2005)
    elif ((method.lower() == 'local') | (method.lower() == 'reichstein')):
        return nee2gpp_reichstein(dates, nee, t, isday, undef=undef, shape=shape, masked=masked, nogppnight=nogppnight)
    # Lasslop et al. (2010) method
    elif ((method.lower() == 'day') | (method.lower() == 'lasslop')):
        return nee2gpp_lasslop(dates, nee, t, isday, rg, vpd, undef=undef, shape=shape, masked=masked, nogppnight=nogppnight)
    # Include new methods here
    else:
        raise ValueError('Error nee2gpp: method not implemented yet.')



# ----------------------------------------------------------------------
def nee2gpp_falge(dates, nee, t, isday, undef=np.nan,
            shape=False, masked=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data, using a fit of Reco vs. temperature to all nighttime data,
        in order to calculate Reco and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp_falge(dates, nee, t, isday, undef=np.nan, shape=False, masked=False):


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
        Falge et al. (2001)
            Gap filling strategies for defensible annual sums of net ecosystem exchange
            Acricultural and Forest Meteorology 107, 43-69


        Examples
        --------
        >>> from jams.fread import fread # from jams
        >>> from jams.date2dec import date2dec # from jams
        >>> dat   = fread('test_nee2gpp.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
        >>> from jams.autostring import astr
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.332' '    8.182' '   10.409' '    8.194' '   10.843']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2013 Matthias Cuntz, Arndt Piayda - mc (at) macu (dot) de

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

        Copyright 2012-2013 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
                 MC, Nov 2012 - individual routine
                 MC, Feb 2013 - ported to Python 3
    """

    # Checks

    # remember shape if any
    inshape = nee.shape
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    # Check squeezed shape
    if dates.ndim != 1: raise Error('Error nee2gpp_falge: squeezed dates must be 1D array.')
    if nee.ndim   != 1: raise Error('Error nee2gpp_falge: squeezed nee must be 1D array.')
    if t.ndim     != 1: raise Error('Error nee2gpp_falge: squeezed t must be 1D array.')
    if isday.ndim != 1: raise Error('Error nee2gpp_falge: squeezed isday must be 1D array.')
    ndata = dates.size
    if ((nee.size != ndata) | (t.size != ndata) | (isday.size != ndata)):
        raise Error('Error nee2gpp_falge: inputs must have the same size.')

    # Transform to masked array with 1D mask
    nee   = np.ma.array(nee, mask=False)
    t     = np.ma.array(t, mask=False)
    isday = np.ma.array(isday, mask=False)
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee)):   nee[np.isnan(nee)]     = np.ma.masked
        if np.ma.any(np.isnan(t)):     t[np.isnan(t)]         = np.ma.masked
        if np.ma.any(np.isnan(isday)): isday[np.isnan(isday)] = np.ma.masked
    else:
    	if np.ma.any(nee==undef):   nee[nee==undef]     = np.ma.masked
    	if np.ma.any(t==undef):     t[t==undef]         = np.ma.masked
    	if np.ma.any(isday==undef): isday[isday==undef] = np.ma.masked

    # Partition - Global relationship as in Falge et al. (2001)

    # Select valid nighttime
    mask = isday | nee.mask | t.mask | isday.mask
    ii   = np.where(~mask)[0]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    # p, c     = opt.curve_fit(functions.lloyd_fix, tt, net, p0=[2.,200.]) # global parameter, global cov matrix
    #p        = opt.fmin(functions.cost_lloyd_fix, [2.,200.], args=(tt, net), disp=False)
    p        = opt.fmin(functions.cost_abs, [2.,200.], args=(functions.lloyd_fix_p, tt, net), disp=False)
    Reco     = np.ones(ndata)*undef
    ii       = np.where(~t.mask)[0]
    Reco[ii] = functions.lloyd_fix(t[ii], p[0], p[1])

    # GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.where(~(t.mask | nee.mask))[0]
    GPP[ii] = Reco[ii] - nee[ii]

    # Return
    if masked:
        if np.isnan(undef):
            GPP  = np.ma.array(GPP,  mask=np.isnan(GPP))
            Reco = np.ma.array(Reco, mask=np.isnan(Reco))
        else:
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
            shape=False, masked=False, nogppnight=False):
    """
        Calculate photosynthesis (GPP) and ecosystem respiration (Reco) from original
        Eddy flux data, using several fits of Reco vs. temperature of nighttime data
        over the season, as in Reichstein et al. (2005), in order to calculate Reco
        and then GPP = Reco - NEE.


        Definition
        ----------
        def nee2gpp_reichstein(dates, nee, t, isday, undef=np.nan, shape=None, masked=False):


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
        shape        if False then outputs are 1D arrays (default)
                     if True, output have the same shape as datain
                     if a shape tuple is given, then this tuple is used to reshape
        masked       if False: outputs are undef where nee and t are masked or undef (default)
                     if True:  return masked arrays where outputs would be undef
        nogppnight   if True:  Resp=NEE, GPP=0 at night
                     if False: Resp=lloyd_taylor, GPP=Resp-NEE at night (default)


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
        >>> from jams.fread import fread # from jams
        >>> from jams.date2dec import date2dec # from jams
        >>> dat   = fread('test_nee2gpp.csv', skip=2, transpose=True)
        >>> dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
        >>> NEE   = np.squeeze(dat[5,:])
        >>> rg    = np.squeeze(dat[6,:])
        >>> tair  = np.squeeze(dat[7,:])
        >>> undef = -9999.
        >>> isday = np.where(rg > 10., True, False)
        >>> tt    = np.where(tair == undef, undef, tair+273.15)
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
        >>> from jams.autostring import astr
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.491' '    8.396' '   10.688' '    8.542' '   11.271']
        >>> print(astr(Reco[1120:1128],3,pp=True))
        ['1.782' '1.906' '2.079' '2.256' '2.464' '2.708' '2.951' '3.218']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='Reichstein')
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    4.491' '    8.396' '   10.688' '    8.542' '   11.271']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['--    ' '--    ' '--    ' ' 4.491' ' 8.396' '10.688' ' 8.542' '11.271']
        >>> GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
        >>> print(astr(GPP[1120:1128],3,pp=True))
        [['-9999.000']
         ['-9999.000']
         ['-9999.000']
         ['    4.491']
         ['    8.396']
         ['   10.688']
         ['    8.542']
         ['   11.271']]


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2013 Matthias Cuntz, Arndt Piayda - mc (at) macu (dot) de

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

        Copyright 2012-2013 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
                 MC, Nov 2012 - individual routine
                 MC, Feb 2013 - ported to Python 3
    """
    # Checks

    # remember shape if any
    if shape != False:
        if shape != True:
            inshape = shape
        else:
            inshape = nee.shape
    dates   = np.squeeze(dates)
    nee     = np.squeeze(nee)
    t       = np.squeeze(t)
    isday   = np.squeeze(isday)
    if shape == False: inshape = nee.shape
    # Check squeezed shape
    if dates.ndim != 1: raise ValueError('Error nee2gpp_reichstein: squeezed dates must be 1D array.')
    if nee.ndim   != 1: raise ValueError('Error nee2gpp_reichstein: squeezed nee must be 1D array.')
    if t.ndim     != 1: raise ValueError('Error nee2gpp_reichstein: squeezed t must be 1D array.')
    if isday.ndim != 1: raise ValueError('Error nee2gpp_reichstein: squeezed isday must be 1D array.')
    ndata = dates.size
    if ((nee.size != ndata) | (t.size != ndata) | (isday.size != ndata)):
        raise ValueError('Error nee2gpp_reichstein: inputs must have the same size.')

    # Transform to masked array with 1D mask
    nee   = np.ma.array(nee, mask=False)
    t     = np.ma.array(t, mask=False)
    isday = np.ma.array(isday, mask=False)
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee)):   nee[np.isnan(nee)]     = np.ma.masked
        if np.ma.any(np.isnan(t)):     t[np.isnan(t)]         = np.ma.masked
        if np.ma.any(np.isnan(isday)): isday[np.isnan(isday)] = np.ma.masked
    else:
    	if np.ma.any(nee==undef):   nee[nee==undef]     = np.ma.masked
    	if np.ma.any(t==undef):     t[t==undef]         = np.ma.masked
    	if np.ma.any(isday==undef): isday[isday==undef] = np.ma.masked

    # Partition - Local relationship = Reichstein et al. (2005)

    # Select valid nighttime
    mask = isday | nee.mask | t.mask | isday.mask
    ii   = np.where(~mask)[0]
    if (ii.size==0):
        print('Warning nee2gpp_reichstein: no valid nighttime data.')
        if masked:
            GPP  = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
            Reco = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
        else:
            GPP  = np.ones(np.reshape(nee,inshape))*undef
            Reco = np.ones(np.reshape(nee,inshape))*undef
        return GPP, Reco
    jul  = dates[ii]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    # 1. each 5 days, in 15 day period, fit if range of T > 5
    locp = [] # local param
    locs = [] # local err
    dmin = np.int(np.floor(np.amin(jul))) # be aware that julian days starts at noon, i.e. 1.0 is 12h
    dmax = np.int(np.ceil(np.amax(jul)))  # so the search will be from noon to noon and thus includes all nights
    for i in range(dmin,dmax,5):
        iii  = np.where((jul>=i) & (jul<(i+14)))[0]
        niii = iii.size
        if niii > 6:
            tt1  = tt[iii]
            net1 = net[iii]
            mm   = ~mad(net1, z=4.5) # make fit more robust by removing outliers
            if (np.ptp(tt[iii]) >= 5.) & (np.sum(mm) > 6):
                # print(i)
                #p     = opt.fmin(functions.cost_lloyd_fix, [2.,200.], args=(tt1[mm], net1[mm]), disp=False) # robust params
                
                p, temp1, temp2 = opt.fmin_tnc(functions.cost_lloyd_fix, [2.,200.], bounds=[[0.,None],[0.,None]],
                                              args=(tt1[mm], net1[mm]),
                                              approx_grad=True, disp=False)
                
                try:
                    p1, c = opt.curve_fit(functions.lloyd_fix, tt1[mm], net1[mm], p0=p, maxfev=10000) # params, covariance
                    if np.all(np.isfinite(c)): # possible return of curvefit: c=inf
                        s = np.sqrt(np.diag(c))
                    else:
                        s = 10.*np.abs(p)
                except:
                    s = 10.*np.abs(p)
                locp += [p]
                locs += [s]
                # if ((s[1]/p[1])<0.5) & (p[1] > 0.): pdb.set_trace()
    if len(locp) == 0:
        raise ValueError('Error nee2gpp_reichstein: No local relationship found.')
        print('Warning nee2gpp_reichstein: No local relationship found.')
        if masked:
            GPP  = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
            Reco = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
        else:
            GPP  = np.ones(np.reshape(nee,inshape))*undef
            Reco = np.ones(np.reshape(nee,inshape))*undef
        return GPP, Reco
    locp   = np.squeeze(np.array(locp).astype(np.float))
    locs   = np.squeeze(np.array(locs).astype(np.float))
    # 2. E0 = avg of best 3
    # Reichstein et al. (2005), p. 1430, 1st paragraph.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        iii  = np.where((locp[:,1] > 0.) & (locp[:,1] < 450.) & (np.abs(locs[:,1]/locp[:,1]) < 0.5))[0]
    niii = iii.size
    if niii==0:
        # raise ValueError('Error nee2gpp_reichstein: No good local relationship found.')
        # loosen the criteria: take the best three estimates anyway
        iii   = np.where((locp[:,1] > 0.))[0]
        niii = iii.size
        if niii<1:
            raise ValueError('Error nee2gpp_reichstein: No E0>0 found.')
            print('Warning nee2gpp_reichstein: No E0>0 found.')
            if masked:
                GPP  = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
                Reco = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
            else:
                GPP  = np.ones(np.reshape(nee,inshape))*undef
                Reco = np.ones(np.reshape(nee,inshape))*undef
            return GPP, Reco
        lp    = locp[iii,:]
        ls    = locs[iii,:]
        iis   = np.argsort(ls[:,1])
        bestp = np.mean(lp[iis[0:np.minimum(3,niii)],:],axis=0)
        bests = np.mean(ls[iis[0:np.minimum(3,niii)],:],axis=0)
    elif niii==1:
        bestp = np.squeeze(locp[iii,:])
        bests = np.squeeze(locs[iii,:])
    elif niii==2:
        bestp = np.mean(locp[iii,:],axis=0)
        bests = np.mean(locs[iii,:],axis=0)
        # ls    = locs[iii,:]
        # iis   = np.argsort(ls[:,1])
    else:
        lp    = locp[iii,:]
        ls    = locs[iii,:]
        iis   = np.argsort(ls[:,1])
        bestp = np.mean(lp[iis[0:3],:],axis=0)
        bests = np.mean(ls[iis[0:3],:],axis=0)

    # 3. Refit Rref with fixed E0, each 4 days
    refp  = [] # Rref param
    refii = [] # mean index of data points
    E0    = bestp[1]
    et    = functions.lloyd_fix(tt, 1., E0)
    for i in range(dmin,dmax,4):
        iii  = np.where((jul>=i) & (jul<(i+4)))[0]
        niii = iii.size
        if niii > 3:
            # Calc directly minisation of (nee-p*et)**2
            # p = np.sum(net[iii]*et[iii])/np.sum(et[iii]**2)
            # p, c = opt.curve_fit(functions.lloyd_only_rref, et[iii], net[iii], p0=[2.])
            #p      = opt.fmin(functions.cost_lloyd_only_rref, [2.], args=(et[iii], net[iii]), disp=False)
            #p = opt.fmin(functions.cost_abs, [2.], args=(functions.lloyd_only_rref_p, et[iii], net[iii]), disp=False)
            
            p, temp1, temp2 = opt.fmin_tnc(functions.cost_abs, [2.], bounds=[[0.,None]],
                                              args=(functions.lloyd_only_rref_p, et[iii], net[iii]),
                                              approx_grad=True, disp=False)
                
            refp  += [p]
            refii += [np.int((iii[0]+iii[-1])//2)]
    if len(refp) == 0:
        raise ValueError('Error nee2gpp_reichstein: No ref relationship found.')
        print('Warning nee2gpp_reichstein: No ref relationship found.')
        if masked:
            GPP  = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
            Reco = np.ma.array(np.reshape(nee,inshape), mask=np.ones(inshape, dtype=np.bool))
        else:
            GPP  = np.ones(np.reshape(nee,inshape))*undef
            Reco = np.ones(np.reshape(nee,inshape))*undef
        return GPP, Reco
    refp  = np.squeeze(np.array(refp))
    refii = np.squeeze(np.array(refii))

    # 4. Interpol Rref
    Rref = np.interp(dates, jul[refii], refp)

    # 5. Calc Reco
    Reco     = np.ones(ndata)*undef
    ii       = np.where(~t.mask)[0]
    Reco[ii] = functions.lloyd_fix(t[ii], Rref[ii], E0)

    # 6. Calc GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.where(~(t.mask | nee.mask))[0]
    GPP[ii] = Reco[ii] - nee[ii]

    # 7. Set GPP=0 at night, if wanted
    if nogppnight:
        mask = isday | nee.mask | t.mask | isday.mask # night
        ii   = np.where(~mask)[0]
        Reco[ii] = nee[ii]
        GPP[ii]  = 0.
        # and prohibit negative gpp at any time
        mask = nee.mask | t.mask | (GPP>0.)
        ii   = np.where(~mask)[0]
        Reco[ii] -= GPP[ii]
        GPP[ii]  = 0.

    if masked:
        if np.isnan(undef):
            GPP  = np.ma.array(GPP,  mask=np.isnan(GPP))
            Reco = np.ma.array(Reco, mask=np.isnan(Reco))
        else:
            GPP  = np.ma.array(GPP,  mask=(GPP==undef))
            Reco = np.ma.array(Reco, mask=(Reco==undef))


    return GPP.reshape(inshape), Reco.reshape(inshape)

# ----------------------------------------------------------------------
def nee2gpp_lasslop(dates, nee, t, isday, rg, vpd, undef=np.nan,
                    shape=False, masked=False, nogppnight=False):
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
        nogppnight   if True:  Resp=NEE, GPP=0 at night
                     if False: Resp=lloyd_taylor, GPP=Resp-NEE at night (default)


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
        >>> from jams.fread import fread # from jams
        >>> from jams.date2dec import date2dec # from jams
        >>> dat   = fread('test_nee2gpp.csv', skip=2, transpose=True)
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
        >>> from jams.autostring import astr
        >>> print(astr(GPP[1120:1128],3,pp=True))
        ['-9999.000' '-9999.000' '-9999.000' '    2.897' '    6.771' '    9.064' '    6.957' '    9.778']
        >>> print(astr(Reco[1120:1128],3,pp=True))
        ['0.352' '0.421' '0.531' '0.662' '0.839' '1.083' '1.365' '1.726']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2013 Matthias Cuntz, Arndt Piayda - mc (at) macu (dot) de

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

        Copyright 2012-2013 Matthias Cuntz, Arndt Piayda


        History
        -------
        Written  MC, Mar 2012
        Modified AP, Mar 2012 - undef=np.nan
                 MC, Nov 2012 - individual routine
                 MC, Feb 2013 - ported to Python 3
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
    nee   = np.ma.array(nee, mask=False)
    t     = np.ma.array(t, mask=False)
    isday = np.ma.array(isday, mask=False)
    rg    = np.ma.array(rg, mask=False)
    vpd   = np.ma.array(vpd, mask=False)
    # mask also undef
    if np.isnan(undef):
        if np.ma.any(np.isnan(nee)):   nee[np.isnan(nee)]     = np.ma.masked
        if np.ma.any(np.isnan(t)):     t[np.isnan(t)]         = np.ma.masked
        if np.ma.any(np.isnan(isday)): isday[np.isnan(isday)] = np.ma.masked
        if np.ma.any(np.isnan(rg)):    rg[np.isnan(rg)]       = np.ma.masked
        if np.ma.any(np.isnan(vpd)):   vpd[np.isnan(vpd)]     = np.ma.masked
    else:
    	if np.ma.any(nee==undef):   nee[nee==undef]     = np.ma.masked
    	if np.ma.any(t==undef):     t[t==undef]         = np.ma.masked
    	if np.ma.any(isday==undef): isday[isday==undef] = np.ma.masked
    	if np.ma.any(rg==undef):    rg[rg==undef]       = np.ma.masked
    	if np.ma.any(vpd==undef):   vpd[vpd==undef]     = np.ma.masked

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
    for i in range(dmin,dmax,2):
        good = True
        # 1. Estimate E0 from nighttime data
        iii  = np.squeeze(np.where((njul>=i) & (njul<(i+12))))
        niii = iii.size
        if niii > 3:
            # p, c = opt.curve_fit(functions.lloyd_fix, ntt[iii], nnet[iii], p0=[aRref,100.])
            #p  = opt.fmin(functions.cost_lloyd_fix, [aRref,100.], args=(ntt[iii], nnet[iii]), disp=False)
            #p  = opt.fmin(functions.cost_abs, [aRref,100.], args=(functions.lloyd_fix_p, ntt[iii], nnet[iii]), disp=False)
            p, temp1, temp2 = opt.fmin_tnc(functions.cost_abs, [aRref,100.], bounds=[[0.,None],[0.,None]],
                                              args=(functions.lloyd_fix_p, ntt[iii], nnet[iii]),
                                              approx_grad=True, disp=False)
            
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
            et     = functions.lloyd_fix(dtt[iii], 1., E0)
            again  = True
            ialpha = aalpha
            ibeta0 = abeta0
            ik     = ak
            iRref  = aRref
            bounds = [[None,None],[None,None],[None,None],[None,None]]
            while again:
                again = False
                p, nfeval, rc  = opt.fmin_tnc(functions.cost_lasslop, [ialpha,ibeta0,ik,iRref], bounds=bounds,
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
    Reco[ii] = functions.lloyd_fix(t[ii], Rref[ii], E0[ii])

    # 5. Calc GPP from light response for check
    if do_lgpp:
        alpha    = np.interp(dates, djul[lii], lE0)
        beta0    = np.interp(dates, djul[lii], lbeta0)
        k        = np.interp(dates, djul[lii], lk)
        et       = functions.lloyd_fix(t, 1., E0)
        lmask    = t.mask | isday.mask | rg.mask | vpd.mask
        ii       = np.squeeze(np.where(~lmask))
        lgpp     = np.zeros(ndata)
        lgpp[ii] = functions.lasslop(rg[ii], et[ii], vpd[ii], alpha[ii], beta0[ii], k[ii], Rref[ii]) - Reco[ii]

    # 6. GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.squeeze(np.where(~(t.mask | nee.mask)))
    GPP[ii] = Reco[ii] - nee[ii]

    # 7. Set GPP=0 at night, if wanted
    if nogppnight:
        mask = isday | nee.mask | t.mask | isday.mask # night
        ii   = np.where(~mask)[0]
        Reco[ii] = nee[ii]
        GPP[ii]  = 0.
        # and prohibit negative gpp at any time
        mask = nee.mask | t.mask | (GPP>0.)
        ii   = np.where(~mask)[0]
        Reco[ii] -= GPP[ii]
        GPP[ii]  = 0.


    if masked:
        if np.isnan(undef):
            GPP  = np.ma.array(GPP,  mask=np.isnan(GPP))
            Reco = np.ma.array(Reco, mask=np.isnan(Reco))
        else:
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
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from fread import fread # from jams
    # from date2dec import date2dec # from jams
    # dat   = fread('test_nee2gpp.csv', skip=2, transpose=True)
    # dates = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    # NEE   = np.squeeze(dat[5,:])
    # rg    = np.squeeze(dat[6,:])
    # tair  = np.squeeze(dat[7,:])
    # undef = -9999.
    # isday = np.where(rg > 10., True, False)
    # tt    = np.where(tair == undef, undef, tair+273.15)
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
    # print(GPP[1120:1128])
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
    # #  8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
    # print(Reco[1120:1128])
    # #[ 1.78172209  1.90616886  2.07856924  2.2560362   2.46373274  2.70757535
    # #  2.95064665  3.2184422 ]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local')
    # print(GPP[1120:1128])
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.49101620e+00
    # #   8.39556274e+00   1.06881053e+01   8.54233665e+00   1.12707122e+01]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='global')
    # print(GPP[1120:1128])
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   4.33166157e+00
    # #   8.18228013e+00   1.04092252e+01   8.19395317e+00   1.08427448e+01]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', masked=True)
    # print(GPP[1120:1128])
    # #[-- -- -- 4.49101619818 8.39556273706 10.6881053462 8.54233664766
    # # 11.2707121977]
    # GPP, Reco = nee2gpp(dates, NEE, tt, isday, undef=undef, method='local', shape=(np.size(NEE),1))
    # print(GPP[1120:1128])
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
    # print(GPP[1120:1128])
    # #[ -9.99900000e+03  -9.99900000e+03  -9.99900000e+03   2.89693618e+00
    # #   6.77103400e+00   9.06351370e+00   6.95696901e+00   9.77798943e+00]
    # print(Reco[1120:1128])
    # #[ 0.35174817  0.42088838  0.53124809  0.66195618  0.839204    1.0829837
    # #  1.36527901  1.72571943]
