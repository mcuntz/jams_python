#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None,
                si=True, sti=True,
                mean=False, wmean=False,
                saltelli=False):
    """
        Calculates the first-order Si and total STi variance-based sensitivity indices
        after Saltelli et al. (2008) with improvements of Saltelli (2002)
        and Mai et al. (2014).

        Definition
        ----------
        def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None,
                        si=True, sti=True,
                        mean=False, wmean=False,
                        saltelli=False):


        Optional Input
        --------------
        Either
          s        ns*(k+2) model outputs
          ns       base sample number
        or
          ya       ns model outputs
          yb       second ns model outputs
          yc       (ns,k) model outputs
        with k: # or parameters. See Saltelli et al. (2008, p. 164ff) for definitions.
        
        si         if True (default): output Si
        sti        if True (default): output STi
        mean       if True: output mean Si and/or STi (if 2D/3D input)
        wmean      if True: output variance weighted mean Si and/or STi (if 2D/3D input)
        saltelli   if True: original Saltelli without changes of Mai et al. (2014)

        If 2D input (3D for yc), then Si/STi will be calculated for each element in first dimension.
        Assumes that first dimension is time, i.e. calculates indices per time step.
        mean and wmean give then the mean and variance weighted mean of the time steps, resp.


        Output
        ------
        First-order sensitivity indices Si and Total sensitivity indices STi.


        Restrictions
        ------------
        The right order and size is assumed if s and ns are given.
        If three positional arguments are give, i.e. would be s, ns and ya, then it is assumed
        that the user gave actually ya, yb and yc.


        References
        ----------
        Mai et a. (2014) ...
        Saltelli, A. (2002) Making best use of model evaluations to compute sensitivity indices.
            Computer Physics Communications, 145(2), 280-297.
        Saltelli, A. et al. (2008) Global sensitivity analysis. The primer.
            John Wiley & Sons Inc., NJ, USA, ISBN 978-0-470-05997-5 (pp. 1-292)


        Examples
        --------

        >>> import numpy as np
        >>> import scipy.stats as stats
        >>> import ufz
        >>> from sobol import i4_sobol_generate

        # ------------------------------------------------------------------------
        #                                1D
        # ------------------------------------------------------------------------

        # seed for reproducible results in doctest
        >>> np.random.seed(1)
        
        # number of parameter sets
        >>> nsets = 1000

        # number of parameter
        >>> npara = 3

        # randomly sample nsets with Latin-Hypercube
        >>> dist = [stats.uniform, stats.uniform, stats.uniform]
        >>> pars = [(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi)]
        >>> setsa  = ufz.lhs(dist, pars, nsets)
        >>> setsb  = ufz.lhs(dist, pars, nsets)

        # Saltelli has two constants in the model
        >>> a = 0.5
        >>> b = 2.0

        # generate model output yA and yB
        >>> iya = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
        >>> iyb = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])

        # generate model output yCi
        >>> iyc = np.empty((npara,nsets))*np.nan
        >>> for i in range(npara):
        ...     tmpset = np.copy(setsb)
        ...     tmpset[i,:] = setsa[i,:]
        ...     iyc[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])

        # theoretical results: has to converge for n->Infinity to this values
        >>> vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
        >>> theo_si  = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])
        >>> theo_sti = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy, (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])

        # Theoretical results
        >>> print('T :: si  =',ufz.astr(theo_si,3,pp=True))
        T :: si  = ['0.372' '0.000' '0.000']
        >>> print('T :: sti =',ufz.astr(theo_sti,3,pp=True))
        T :: sti = ['1.000' '0.000' '0.628']
    
        # Give 3 positional arguments
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc)
        >>> from autostring import astr
        >>> print('S :: si  =',astr(isi1,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']
        >>> print('S :: sti =',astr(isti1,3,pp=True))
        S :: sti = ['1.037' '0.000' '0.647']

        # Give 3 optional arguments
        >>> isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc)
        >>> print('S :: si  =',astr(isi2,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']

        # Give 2 positional arguments
        >>> s = np.concatenate((iya, iyb, np.ravel(iyc)))
        >>> ns = iya.size
        >>> isi3, isti3 = sobol_index(s, ns)
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']

        # Give 2 optional arguments
        >>> isi4, isti4 = sobol_index(s=s, ns=ns)
        >>> print('S :: si  =',astr(isi4,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']

        # Give 2 positional arguments and original Saltelli
        >>> isi3, isti3 = sobol_index(s, ns, saltelli=True)
        >>> print('Sal :: si  =',astr(isi3,3,pp=True))
        Sal :: si  = [' 0.361' ' 0.000' ' 0.005']
        >>> print('Sal :: sti =',astr(isti3,3,pp=True))
        Sal :: sti = [' 1.050' '-0.020' ' 0.648']

        # 2 optional arguments and no STi output
        >>> isi5, = sobol_index(s=s, ns=ns, si=True, sti=False)
        >>> print('S :: si  =',astr(isi5,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']

        # 2 optional arguments and no Si output
        >>> isti5, = sobol_index(s=s, ns=ns, si=False, sti=True)
        >>> print('S :: sti =',astr(isti5,3,pp=True))
        S :: sti = ['1.037' '0.000' '0.647']

        # ------------------------------------------------------------------------
        #                                2D
        # ------------------------------------------------------------------------
        # number of timepoints
        >>> ntime = 100
        >>> iya = np.empty((ntime,nsets))*np.nan
        >>> iyb = np.empty((ntime,nsets))*np.nan
        >>> iyc = np.empty((ntime,npara,nsets))*np.nan
        >>> s   = np.empty((ntime,(npara+2)*nsets))*np.nan
        >>> ns = nsets
        >>> theo_si  = []
        >>> theo_sti = []
        >>> for t in range(ntime):
        ...     a = 1.0*t + 50.0
        ...     b = 2.0
        ...     iya[t,:] = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
        ...     iyb[t,:] = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
        ...     for i in range(npara):
        ...         tmpset = np.copy(setsb)
        ...         tmpset[i,:] = setsa[i,:]
        ...         iyc[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])
        ...     s[t,:] = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc[t,:,:])))
        ...     vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
        ...     theo_si  = theo_si + [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])]
        ...     theo_sti = theo_sti + [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy, (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])]

        # -------------------------------------------------------------------------------------------
        # 1D model output --> ya, yb are 2D and yc is 3D
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc)

        #
        # SI :: 1st time point
        >>> print('S :: si[t0]  =',astr(isi1[0,:],3,pp=True))
        S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
        >>> print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
        T :: si[t0]  = ['0.325' '0.127' '0.000']

        #
        # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi1[2,:],3,pp=True))
        S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']

        #
        # STI :: 1st time point
        >>> print('S :: sti[t0] =',astr(isti1[0,:],3,pp=True))
        S :: sti[t0] = ['0.876' '0.138' '0.560']
        >>> print('T :: sti[t0] =',astr(theo_sti[0],3,pp=True))
        T :: sti[t0] = ['0.873' '0.127' '0.548']

        #
        # STI :: 3rd time point
        >>> print('S :: sti[t2] =',astr(isti1[2,:],3,pp=True))
        S :: sti[t2] = ['0.866' '0.147' '0.554']
        >>> print('T :: sti[t2] =',astr(theo_sti[2],3,pp=True))
        T :: sti[t2] = ['0.864' '0.136' '0.543']

        # -------------------------------------------------------------------------------------------
        >>> isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc, mean=True)

        #
        # SI :: 1st time point
        >>> print('S :: si[t0]  =',astr(isi2[0,:],3,pp=True))
        S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
        >>> print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
        T :: si[t0]  = ['0.325' '0.127' '0.000']

        #
        # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi2[2,:],3,pp=True))
        S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']

        #
        # SI :: mean over all timepoints
        >>> print('S :: si_m    =',astr(msi2,3,pp=True))
        S :: si_m    = [' 0.266' ' 0.334' '-0.005']

        #
        # STI :: mean over all timepoints
        >>> print('S :: sti_m   =',astr(msti2,3,pp=True))
        S :: sti_m   = ['0.631' '0.361' '0.412']

        # -------------------------------------------------------------------------------------------
        >>> isi3, isti3 = sobol_index(s=s, ns=ns)

        #
        # SI :: 1st time point
        >>> print('S :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
        >>> print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
        T :: si[t0]  = ['0.325' '0.127' '0.000']

        #
        # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']

        # -------------------------------------------------------------------------------------------
        >>> isi3, isti3 = sobol_index(s, ns, saltelli=True)

        # SI :: 1st time point
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.347' ' 0.120' '-0.002']
        >>> print('Tal :: si[t0]  =',astr(theo_si[0],3,pp=True))
        Tal :: si[t0]  = ['0.325' '0.127' '0.000']

        # SI :: 3rd time point
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.344' ' 0.129' '-0.002']
        >>> print('Tal :: si[t2]  =',astr(theo_si[2],3,pp=True))
        Tal :: si[t2]  = ['0.321' '0.136' '0.000']
    
        # -------------------------------------------------------------------------------------------
        >>> isi4, isti4, wsi2, wsti2 = sobol_index(s, ns, wmean=True)

        #
        # SI :: 1st time point
        >>> print('S :: si[t0]  =',astr(isi4[0,:],3,pp=True))
        S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
        >>> print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
        T :: si[t0]  = ['0.325' '0.127' '0.000']

        #
        # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi4[2,:],3,pp=True))
        S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']

        #
        # SI :: weighted mean over all timepoints
        >>> print('S :: si_w    =',astr(wsi2,3,pp=True))
        S :: si_w    = [' 0.257' ' 0.360' '-0.005']

        #
        # STI :: weighted mean over all timepoints
        >>> print('S :: sti_w   =',astr(wsti2,3,pp=True))
        S :: sti_w   = ['0.602' '0.387' '0.395']

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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Aug 2013 - Si per time step
                  JM, Sep 2013 - changes of Mai et al. (2014), theoretical example in docstring
                  MC, Sep 2013 - saltelli
    """
    # Check input
    if (si==False) & (sti==False):
        raise ValueError('No output chosen: si=False and sti=False.')
    # s, ns or ya, yb, yc
    isone = False
    if ((s != None) & (ns != None) & (ya == None)):
        if s.ndim == 1:
            isone = True
            s = s[np.newaxis,:]
        nsa   = ns
        ntime = s.shape[0]
        ss    = s.shape[1]
        nn    = (ss//nsa)-2
        iyA   = s[:,0:nsa]
        iyB   = s[:,nsa:2*nsa]
        iyC   = np.reshape(s[:,2*nsa:],(ntime,nn,nsa))
    else:
        if not (((s != None) & (ns != None) & (ya != None)) | ((ya != None) & (yb !=None) & (yc !=None))):
            raise ValueError('Either s and ns has to be given or ya, yb and yc.')
        if ((s != None) & (ns != None) & (ya != None)):
            if s.ndim == 1:
                isone = True
                s  = s[np.newaxis,:]
                ns = ns[np.newaxis,:]
                ya = ya[np.newaxis,:,:]
            iyA     = s
            iyB     = ns
            iyC     = ya
        else:
            if ya.ndim == 1:
                isone = True
                ya = ya[np.newaxis,:]
                yb = yb[np.newaxis,:]
                yc = yc[np.newaxis,:,:]
            iyA     = ya
            iyB     = yb
            iyC     = yc
        ntime = iyA.shape[0]
        nsa   = iyA.shape[1]
        if not ((nsa == iyB.shape[1]) & (nsa == iyC.shape[2])):
            raise ValueError('ya and yb must have same size as yc[1].')
        nn = iyC.shape[1]
    fnsa  = 1. / np.float(nsa)
    fnsa1 = 1. / np.float(nsa-1)
    f02   = fnsa*np.sum(iyA*iyB, axis=1)
    if saltelli:
        iyAiyA = np.mean(iyA**2, axis=1)
        varA   = iyAiyA - f02
    else:
        f0B    = fnsa*np.sum(iyB,    axis=1)
        iyBiyB = np.mean(iyB**2,     axis=1)
        yab    = np.append(iyA, iyB, axis=1)
        varAB  = np.var(yab, ddof=1, axis=1)
        varB   = iyBiyB - f0B**2
    isi  = np.empty((ntime,nn))
    isti = np.empty((ntime,nn))
    for i in range(nn):
        iyCi = iyC[:,i,:]
        if saltelli:
            iyAiyC    = fnsa1*np.sum(iyA*iyCi, axis=1)
            isi[:,i]  = (iyAiyC - f02) / varA
            iyBiyC    = fnsa1*np.sum(iyB*iyCi, axis=1)
            isti[:,i] = 1. - (iyBiyC - f02) / varA
        else:
            iyAiyC    = fnsa*np.sum(iyA*iyCi, axis=1)
            isi[:,i]  = (iyAiyC - f02) / varAB
            iyBiyC    = fnsa*np.sum(iyB*iyCi, axis=1)
            isti[:,i] = 1. - (iyBiyC - f0B**2) / varB

    if not isone:
        # simple mean
        if mean:
            msi  = np.mean(isi,  axis=0)
            msti = np.mean(isti, axis=0)

        # weighted mean
        if wmean:
            if saltelli:
                denom = 1./np.sum(varA, axis=0)
                varA  = varA[:,np.newaxis]
                wsi   = np.sum(isi*varA,  axis=0) * denom
                wsti  = np.sum(isti*varA, axis=0) * denom
            else:
                denomAB = 1./np.sum(varAB, axis=0)
                denomB  = 1./np.sum(varB,  axis=0)
                varB    = varB[:, np.newaxis]
                varAB   = varAB[:,np.newaxis]
                wsi     = np.sum(isi*varAB, axis=0) * denomAB
                wsti    = np.sum(isti*varB, axis=0) * denomB

    # Output
    if isone:
        isi  = isi[0,:]
        isti = isti[0,:]
        
    out = []
    if si:  out += [isi]
    if sti: out += [isti]
    if not isone:
        if mean:
            if si:  out += [msi]
            if sti: out += [msti]
        if wmean:
            if si:  out += [wsi]
            if sti: out += [wsti]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # import numpy as np
    # import scipy.stats as stats
    # import ufz
    # from sobol import i4_sobol_generate

    # # ------------------------------------------------------------------------
    # #                                1D
    # # ------------------------------------------------------------------------

    # # seed for reproducible results in doctest
    # np.random.seed(1)
    # # number of parameter sets
    # nsets = 1000
    # # number of parameter
    # npara = 3
    # # randomly sample nsets with Latin-Hypercube
    # dist = [stats.uniform, stats.uniform, stats.uniform]
    # pars = [(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi)]
    # setsa  = ufz.lhs(dist, pars, nsets)
    # setsb  = ufz.lhs(dist, pars, nsets)
    # # Saltelli has two constants in the model
    # a = 0.5
    # b = 2.0
    # # generate model output yA and yB
    # iya = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
    # iyb = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
    # # generate model output yCi
    # iyc = np.empty((npara,nsets))*np.nan
    # for i in range(npara):
    #     tmpset = np.copy(setsb)
    #     tmpset[i,:] = setsa[i,:]
    #     iyc[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])
    # # theoretical results: has to converge for n->Infinity to this values
    # vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
    # theo_si  = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])
    # theo_sti = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy, (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])

    # # Theoretical results
    # print('T :: si  =',ufz.astr(theo_si,3,pp=True))
    # #T :: si  = ['0.372' '0.000' '0.000']
    # print('T :: sti =',ufz.astr(theo_sti,3,pp=True))
    # #T :: sti = ['1.000' '0.000' '0.628']
    
    # # Give 3 positional arguments
    # isi1, isti1 = sobol_index(iya,iyb,iyc)
    # from autostring import astr
    # print('S :: si  =',astr(isi1,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']
    # print('S :: sti =',astr(isti1,3,pp=True))
    # #S :: sti = ['1.037' '0.000' '0.647']

    # # Give 3 optional arguments
    # isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc)
    # print('S :: si  =',astr(isi2,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']

    # # Give 2 positional arguments
    # s = np.concatenate((iya, iyb, np.ravel(iyc)))
    # ns = iya.size
    # isi3, isti3 = sobol_index(s, ns)
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']

    # # Give 2 optional arguments
    # isi4, isti4 = sobol_index(s=s, ns=ns)
    # print('S :: si  =',astr(isi4,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']

    # # 2 optional arguments and no STi output
    # isi5, = sobol_index(s=s, ns=ns, si=True, sti=False)
    # print('S :: si  =',astr(isi5,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']

    # # 2 optional arguments and no Si output
    # isti5, = sobol_index(s=s, ns=ns, si=False, sti=True)
    # print('S :: sti =',astr(isti5,3,pp=True))
    # #S :: sti = ['1.037' '0.000' '0.647']

    # # Give 2 positional arguments and original Saltelli
    # isi3, isti3 = sobol_index(s, ns, saltelli=True)
    # print('Sal :: si  =',astr(isi3,3,pp=True))
    # print('Sal :: sti =',astr(isti3,3,pp=True))
    # #Sal :: si  = [' 0.361' ' 0.000' ' 0.005']
    # #Sal :: sti = [' 1.050' '-0.020' ' 0.648']

    # # ------------------------------------------------------------------------
    # #                                2D
    # # ------------------------------------------------------------------------
    # # number of timepoints
    # ntime = 100
    # iya = np.empty((ntime,nsets))*np.nan
    # iyb = np.empty((ntime,nsets))*np.nan
    # iyc = np.empty((ntime,npara,nsets))*np.nan
    # s   = np.empty((ntime,(npara+2)*nsets))*np.nan
    # ns = nsets
    # theo_si  = []
    # theo_sti = []
    # for t in range(ntime):
    #     a = 1.0*t + 50.0
    #     b = 2.0
    #     # generate model output yA and yB
    #     iya[t,:] = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
    #     iyb[t,:] = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
    #     # generate model output yCi
    #     for i in range(npara):
    #         tmpset = np.copy(setsb)
    #         tmpset[i,:] = setsa[i,:]
    #         iyc[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])
    #     # for 2 optional argument call
    #     s[t,:] = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc[t,:,:])))
    #     # theoretical results: has to converge for n->Infinity to this values
    #     vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
    #     theo_si  = theo_si + [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])]
    #     theo_sti = theo_sti + [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy, (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])]
    # # print('Theoretical Si:  ',ufz.astr(theo_si,3,pp=True))
    # # print('Theoretical STi: ',ufz.astr(theo_sti,3,pp=True))

    # # -------------------------------------------------------------------------------------------
    # # 1D model output --> ya, yb are 2D and yc is 3D
    # isi1, isti1 = sobol_index(iya,iyb,iyc)
    # #
    # # SI :: 1st time point
    # print('S :: si[t0]  =',astr(isi1[0,:],3,pp=True))
    # #S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
    # print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #T :: si[t0]  = ['0.325' '0.127' '0.000']
    # #
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi1[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']
    # #
    # # STI :: 1st time point
    # print('S :: sti[t0] =',astr(isti1[0,:],3,pp=True))
    # #S :: sti[t0] = ['0.876' '0.138' '0.560']
    # print('T :: sti[t0] =',astr(theo_sti[0],3,pp=True))
    # #T :: sti[t0] = ['0.873' '0.127' '0.548']
    # #
    # # STI :: 3rd time point
    # print('S :: sti[t2] =',astr(isti1[2,:],3,pp=True))
    # #S :: sti[t2] = ['0.866' '0.147' '0.554']
    # print('T :: sti[t2] =',astr(theo_sti[2],3,pp=True))
    # #T :: sti[t2] = ['0.864' '0.136' '0.543']

    # # -------------------------------------------------------------------------------------------
    # isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc, mean=True)
    # #
    # # SI :: 1st time point
    # print('S :: si[t0]  =',astr(isi2[0,:],3,pp=True))
    # #S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
    # print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #T :: si[t0]  = ['0.325' '0.127' '0.000']
    # #
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi2[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']
    # #
    # # SI :: mean over all timepoints
    # print('S :: si_m    =',astr(msi2,3,pp=True))
    # #S :: si_m    = [' 0.266' ' 0.334' '-0.005']
    # #
    # # STI :: mean over all timepoints
    # print('S :: sti_m   =',astr(msti2,3,pp=True))
    # #S :: sti_m   = ['0.631' '0.361' '0.412']

    # # -------------------------------------------------------------------------------------------
    # isi3, isti3 = sobol_index(s=s, ns=ns)
    # #
    # # SI :: 1st time point
    # print('S :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
    # print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #T :: si[t0]  = ['0.325' '0.127' '0.000']
    # #
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']

    # # -------------------------------------------------------------------------------------------
    # isi3, isti3 = sobol_index(s, ns, saltelli=True)
    # #
    # # SI :: 1st time point
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.347' ' 0.120' '-0.002']
    # print('Tal :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #Tal :: si[t0]  = ['0.325' '0.127' '0.000']
    # #
    # # SI :: 3rd time point
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.344' ' 0.129' '-0.002']
    # print('Tal :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #Tal :: si[t2]  = ['0.321' '0.136' '0.000']
    
    # # -------------------------------------------------------------------------------------------
    # isi4, isti4, wsi2, wsti2 = sobol_index(s, ns, wmean=True)
    # #
    # # SI :: 1st time point
    # print('S :: si[t0]  =',astr(isi4[0,:],3,pp=True))
    # #S :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
    # print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #T :: si[t0]  = ['0.325' '0.127' '0.000']
    # #
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi4[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.334' ' 0.125' '-0.002']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']
    # #
    # # SI :: weighted mean over all timepoints
    # print('S :: si_w    =',astr(wsi2,3,pp=True))
    # #S :: si_w    = [' 0.257' ' 0.360' '-0.005']
    # #
    # # STI :: weighted mean over all timepoints
    # print('S :: sti_w   =',astr(wsti2,3,pp=True))
    # #S :: sti_w   = ['0.602' '0.387' '0.395']
