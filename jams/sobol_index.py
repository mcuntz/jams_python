#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None,
                si=True, sti=True,
                mean=False, wmean=False,
                method='Mai1999'):
    """
        Calculates the first-order Si and total STi variance-based sensitivity indices
        summarised in Saltelli et al. (2010) with improvements of Mai et al. (2014).


        Definition
        ----------
        def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None,
                        si=True, sti=True,
                        mean=False, wmean=False,
                        method='Mai1999'):


        Optional Input
        --------------
        Either
          s        ns*(k+2) model outputs
          ns       base sample number
        or
          ya       ns model outputs f(A)
          yb       second ns model outputs f(B)
          yc       (ns,k) model outputs f(A_B) or f(B_A)
        with k: # or parameters. See Saltelli et al. (2008, 2010) for definitions.
        yc can be either f(B_A), i.e. yb with one column of ya
        or vice versa f(A_B), i.e. ya with one column of yb. See Saltelli et al. (2010).

        If 2D input (3D for yc), then Si/STi will be calculated for each element in first dimension.
        Assumes that first dimension is time, i.e. calculates indices per time step.
        mean and wmean give then the mean and variance weighted mean of the time steps, resp.

        si         if True (default): output Si
        sti        if True (default): output STi
        mean       if True: output mean Si and/or STi (if 2D/3D input)
        wmean      if True: output variance weighted mean Si and/or STi (if 2D/3D input)
        method     string, case-insensitive (default: 'Mai1999')
                   'Saltelli2008' - The formulation presented in 'The Primer'. (yc=f(B_A))
                                    Si  = (1/n*sum_j(f(A)_j*f(B_A^i)_j) - mean(f(A))^2)/var(f(A))
                                    STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(A))^2))/var(f(A))
                   'Homma1996'    - Si takes mixed A, B term for squared mean in nominator. (yc=f(B_A))
                                    Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var(f(A))
                                    STi as Saltelli2008
                   'Sobol2007'    - STi takes A, A_B instead of B, B_A (needs f(A_B) and f(B_A), not avail.)
                                    Si as Homma1996
                                    STi = 1/n*sum_j(f(A)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
                   'Saltelli2010' - Si takes B, A_B instead of A, B_A (yc=f(A_B))
                                    Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
                                    STi as Sobol2007
                   'Jansen1999'   - Calculate Si and STi by expectation(variance) instead of variance(expectation) (yc=f(A_B))
                                    Si  = (var(f(A)) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var(f(A))
                                    STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
                   'Mai2012'      - Si of Homma 1996 but denominator is variance of A and B (yc=f(B_A)) and
                                    STi of Homma 1996 but mean and variance of f(A) replaced by f(B)
                                    Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var([f(A),f(B)])
                                    STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(B)^2)))/var(f(B))
                   'Mai2013'      - Si of Saltelli2010 but denominator is variance of A and B (yc=f(A_B))
                                    Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
                                    STi as Sobol2007/Saltelli2010
                   'Mai2014'      - SI of Jansen1999  but denominator is variance of A and B (yc=f(A_B))
                                    Si  = (var([f(A),f(B)]) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var([f(A),f(B)])
                                    STi as Jansen1999
                   'Mai1999'      - SI of Mai2013 (Saltelli2010 with var([f(A),f(B)])) and STi of Jansen1999 (yc=f(A_B))
                                    Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
                                    STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))


        Output
        ------
        First-order sensitivity indices Si and Total sensitivity indices STi.


        Restrictions
        ------------
        1. The right order and size is assumed if s and ns are given.
        2. If three positional arguments are given, i.e. would be s, ns and ya, then it is assumed
           that the user provided ya, yb and yc.
        3. yc is either f(A_B) or f(B_A) depending on the method.
        4. method='Sobol2007' would need f(A_B) and f(B_A) and is therefore not available here.


        References
        ----------
        Mai et al. (2014) ...
        Saltelli, A. et al. (2008) Global sensitivity analysis. The primer.
            John Wiley & Sons Inc., NJ, USA, ISBN 978-0-470-05997-5 (pp. 1-292)
        Saltelli, A. et al. (2010), Variance based sensitivity analysis of model output. Design and estimator
            for the total sensitivity index. Computer Physics Communications 181, 259-270.


        Examples
        --------

        >>> import numpy as np
        >>> import scipy.stats as stats
        >>> from lhs import lhs
        >>> from autostring import astr

        >>> # ------------------------------------------------------------------------
        >>> #                                1D
        >>> #
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)
        >>> # number of parameter sets
        >>> nsets = 1000
        >>> # number of parameters
        >>> npara = 3
        >>> # randomly sample nsets with Latin-Hypercube
        >>> dist = [stats.uniform, stats.uniform, stats.uniform]
        >>> pars = [(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi)]
        >>> setsa  = lhs(dist, pars, nsets)
        >>> setsb  = lhs(dist, pars, nsets)
        >>> # Saltelli has two constants in the model
        >>> a = 0.5
        >>> b = 2.0
        >>> # generate model output yA and yB
        >>> iya = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
        >>> iyb = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
        >>> # generate model output yCi = f(B_A)
        >>> iyc = np.empty((npara,nsets))*np.nan
        >>> for i in range(npara):
        ...     tmpset = np.copy(setsb)
        ...     tmpset[i,:] = setsa[i,:]
        ...     iyc[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])
        >>> # generate model output yCi = f(A_B)
        >>> iyc2 = np.empty((npara,nsets))*np.nan
        >>> for i in range(npara):
        ...     tmpset = np.copy(setsa)
        ...     tmpset[i,:] = setsb[i,:]
        ...     iyc2[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])

        >>> # Theoretical results: has to converge for n->Infinity to this values
        >>> vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
        >>> theo_si  = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])
        >>> theo_sti = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy,
        ...                      (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])
        >>> # Theoretical results
        >>> print('T :: si  =',astr(theo_si,3,pp=True))
        T :: si  = ['0.372' '0.000' '0.000']
        >>> print('T :: sti =',astr(theo_sti,3,pp=True))
        T :: sti = ['1.000' '0.000' '0.628']

        >>> # 3 positional arguments
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc2)
        >>> print('S :: si  =',astr(isi1,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']
        >>> print('S :: sti =',astr(isti1,3,pp=True))
        S :: sti = ['1.014' '0.000' '0.654']

        >>> # 3 optional arguments
        >>> isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc2)
        >>> print('S :: si  =',astr(isi2,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']

        >>> # 2 positional arguments
        >>> s2 = np.concatenate((iya, iyb, np.ravel(iyc2)))
        >>> ns = iya.size
        >>> isi3, isti3 = sobol_index(s2, ns)
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']

        >>> # 2 optional arguments
        >>> isi4, isti4 = sobol_index(s=s2, ns=ns)
        >>> print('S :: si  =',astr(isi4,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']

        >>> # 2 optional arguments and no STi output
        >>> isi5 = sobol_index(s=s2, ns=ns, si=True, sti=False)
        >>> print('S :: si  =',astr(isi5,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']

        >>> # 2 optional arguments and no Si output
        >>> isti5 = sobol_index(s=s2, ns=ns, si=False, sti=True)
        >>> print('S :: sti =',astr(isti5,3,pp=True))
        S :: sti = ['1.014' '0.000' '0.654']

        >>> # Method=Saltelli2008
        >>> s = np.concatenate((iya, iyb, np.ravel(iyc)))
        >>> isi3, isti3 = sobol_index(s, ns, method='Saltelli2008')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = ['0.367' '0.011' '0.016']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = [' 1.038' '-0.019' ' 0.641']

        >>> # Method=Homma1996
        >>> isi3, isti3 = sobol_index(s, ns, method='homma1996')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.356' ' 0.000' ' 0.005']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = [' 1.038' '-0.019' ' 0.641']

        >>> # Method=Mai2012
        >>> isi3, isti3 = sobol_index(s, ns, method='mai2012')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.353' ' 0.000' ' 0.005']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = ['1.037' '0.000' '0.647']

        >>> # Method=saltelli2010
        >>> isi3, isti3 = sobol_index(s2, ns, method='saltelli2010')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.348' ' 0.000' '-0.049']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = [' 0.984' ' 0.000' ' 0.633']

        >>> # Method=jansen1999
        >>> isi3, isti3 = sobol_index(s2, ns, method='jansen1999')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.319' ' 0.002' '-0.068']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = ['1.014' '0.000' '0.654']

        >>> # Method=mai2013
        >>> isi3, isti3 = sobol_index(s2, ns, method='mai2013')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = [' 0.984' ' 0.000' ' 0.633']

        >>> # Method=mai2014
        >>> isi3, isti3 = sobol_index(s2, ns, method='mai2014')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.325' ' 0.011' '-0.059']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = ['1.014' '0.000' '0.654']

        >>> # Method=mai
        >>> isi3, isti3 = sobol_index(s2, ns, method='mai1999')
        >>> print('S :: si  =',astr(isi3,3,pp=True))
        S :: si  = [' 0.345' ' 0.000' '-0.049']
        >>> print('S :: sti =',astr(isti3,3,pp=True))
        S :: sti = ['1.014' '0.000' '0.654']

        >>> # ------------------------------------------------------------------------
        >>> #                                2D
        >>> #
        >>> # number of timepoints
        >>> ntime = 100
        >>> iya   = np.empty((ntime,nsets))*np.nan
        >>> iyb   = np.empty((ntime,nsets))*np.nan
        >>> iyc   = np.empty((ntime,npara,nsets))*np.nan
        >>> iyc2  = np.empty((ntime,npara,nsets))*np.nan
        >>> s     = np.empty((ntime,(npara+2)*nsets))*np.nan
        >>> s2    = np.empty((ntime,(npara+2)*nsets))*np.nan
        >>> ns    = nsets
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
        ...         iyc[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + (
        ...                      b*tmpset[2,:]**4  * np.sin(tmpset[0,:]))
        ...     for i in range(npara):
        ...         tmpset = np.copy(setsa)
        ...         tmpset[i,:] = setsb[i,:]
        ...         iyc2[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + (
        ...                       b*tmpset[2,:]**4  * np.sin(tmpset[0,:]))
        ...     s[t,:]  = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc[t,:,:])))
        ...     s2[t,:] = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc2[t,:,:])))
        ...     vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
        ...     theo_si  += [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])]
        ...     theo_sti += [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy,
        ...                            (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])]

        >>> # 1D model output --> ya, yb are 2D and yc is 3D
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc2)

        >>> # SI :: 1st time point
        >>> print('S :: si[t0]  =',astr(isi1[0,:],3,pp=True))
        S :: si[t0]  = [' 0.290' ' 0.140' '-0.031']
        >>> print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
        T :: si[t0]  = ['0.325' '0.127' '0.000']
        >>> # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi1[2,:],3,pp=True))
        S :: si[t2]  = [' 0.286' ' 0.149' '-0.030']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']
        >>> # STI :: 1st time point
        >>> print('S :: sti[t0] =',astr(isti1[0,:],3,pp=True))
        S :: sti[t0] = ['0.881' '0.126' '0.568']
        >>> print('T :: sti[t0] =',astr(theo_sti[0],3,pp=True))
        T :: sti[t0] = ['0.873' '0.127' '0.548']
        >>> # STI :: 3rd time point
        >>> print('S :: sti[t2] =',astr(isti1[2,:],3,pp=True))
        S :: sti[t2] = ['0.871' '0.134' '0.562']
        >>> print('T :: sti[t2] =',astr(theo_sti[2],3,pp=True))
        T :: sti[t2] = ['0.864' '0.136' '0.543']

        >>> # Mean
        >>> isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc2, mean=True)
        >>> print('S :: si_m    =',astr(msi2,3,pp=True))
        S :: si_m    = [' 0.207' ' 0.367' '-0.015']
        >>> print('S :: sti_m   =',astr(msti2,3,pp=True))
        S :: sti_m   = ['0.643' '0.351' '0.415']

        >>> # Weighted mean
        >>> isi4, isti4, wsi2, wsti2 = sobol_index(s2, ns, wmean=True)
        >>> print('S :: si_w    =',astr(wsi2,3,pp=True))
        S :: si_w    = [' 0.197' ' 0.394' '-0.014']
        >>> print('S :: sti_w   =',astr(wsti2,3,pp=True))
        S :: sti_w   = ['0.616' '0.377' '0.397']

        >>> # 2 positional arguments
        >>> isi3, isti3 = sobol_index(s=s2, ns=ns)
        >>> # SI :: 3rd time point
        >>> print('S :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        S :: si[t2]  = [' 0.286' ' 0.149' '-0.030']
        >>> print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
        T :: si[t2]  = ['0.321' '0.136' '0.000']

        >>> # Method=Saltelli2008
        >>> isi3, isti3 = sobol_index(s, ns, method='saltelli2008')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = ['0.354' '0.130' '0.011']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = ['0.352' '0.139' '0.010']

        >>> # Method=Homma1996
        >>> isi3, isti3 = sobol_index(s, ns, method='homma1996')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.342' ' 0.118' '-0.002']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.340' ' 0.127' '-0.002']

        >>> # Method=Mai2012
        >>> isi3, isti3 = sobol_index(s, ns, method='Mai2012')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.334' ' 0.125' '-0.002']

        >>> # Method=Saltelli2010
        >>> isi3, isti3 = sobol_index(s2, ns, method='saltelli2010')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.294' ' 0.142' '-0.031']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.291' ' 0.152' '-0.031']

        >>> # Method=Jansen1999
        >>> isi3, isti3 = sobol_index(s2, ns, method='jansen1999')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.284' ' 0.133' '-0.072']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.281' ' 0.143' '-0.072']

        >>> # Method=Mai2013
        >>> isi3, isti3 = sobol_index(s2, ns, method='Mai2013')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.290' ' 0.140' '-0.031']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.286' ' 0.149' '-0.030']

        >>> # Method=Mai2014
        >>> isi3, isti3 = sobol_index(s2, ns, method='mai2014')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.295' ' 0.147' '-0.056']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.293' ' 0.156' '-0.055']

        >>> # Method=Mai
        >>> isi3, isti3 = sobol_index(s2, ns, method='mai1999')
        >>> print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
        Sal :: si[t0]  = [' 0.290' ' 0.140' '-0.031']
        >>> print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
        Sal :: si[t2]  = [' 0.286' ' 0.149' '-0.030']

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Matthias Cuntz - mc (at) macu (dot) de

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


        History
        -------
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Aug 2013 - Si per time step
                  JM, Sep 2013 - changes of Mai et al. (2014), theoretical example in docstring
                  MC, Sep 2013 - saltelli
                  MC, Sep 2013 - method, removed saltelli
                  MC, Apr 2014 - assert
    """
    # Check input
    assert (si+sti) > 0, 'No output chosen: si=False and sti=False.'
    # s, ns or ya, yb, yc
    isone = False
    if ((s is not None) and (ns is not None) and (ya is None)):
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
        assert (((s is not None) and (ns is not None) and (ya is not None)) or ((ya is not None) and (yb is not None) and (yc is not None))), 'Either s and ns has to be given or ya, yb and yc.'
        if ((s is not None) and (ns is not None) and (ya is not None)):
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

    mm = method.lower()
    if mm == 'saltelli2008':
        meanA = np.mean(iyA, axis=1)
        varA  = np.var(iyA,  axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = (np.mean(iyA*iyCi, axis=1) - meanA**2) / varA
            isti[:,i] = 1. - (np.mean(iyB*iyCi, axis=1) - meanA**2) / varA
    elif mm == 'homma1996':
        meanA = np.mean(iyA, axis=1)
        varA  = np.var(iyA,  axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = np.mean(iyA*(iyCi - iyB), axis=1) / varA
            isti[:,i] = 1. - (np.mean(iyB*iyCi, axis=1) - meanA**2) / varA
            #old code isti[:,i] = 1. - (np.mean(iyB*(iyCi-iyA), axis=1)) / varA
    elif mm == 'sobol2007':
        raise ValueError('Sobol2007 would need f(A_B) and f(B_A). It is thus not implemented here.')
    elif mm == 'saltelli2010':
        varA  = np.var(iyA,  axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = np.mean(iyB*(iyCi - iyA), axis=1) / varA
            isti[:,i] = np.mean(iyA*(iyA - iyCi), axis=1) / varA
    elif mm == 'jansen1999':
        varA  = np.var(iyA,  axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = 1. - np.sum((iyB -iyCi)**2, axis=1) / (2.*nsa * varA)
            isti[:,i] = np.sum((iyA -iyCi)**2, axis=1) / (2.*nsa * varA)
    elif mm == 'mai2012':
        meanB = np.mean(iyB, axis=1)
        varB  = np.var(iyB,  axis=1)
        varAB = np.var(np.append(iyA, iyB, axis=1), axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = np.mean(iyA*(iyCi - iyB), axis=1) / varAB
            isti[:,i] = 1. - (np.mean(iyB*iyCi, axis=1) - meanB**2) / varB
    elif mm == 'mai2013':
        varA  = np.var(iyA,  axis=1)
        varAB = np.var(np.append(iyA, iyB, axis=1), axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = np.mean(iyB*(iyCi - iyA), axis=1) / varAB
            isti[:,i] = np.mean(iyA*(iyA - iyCi), axis=1) / varA
    elif mm == 'mai2014':
        varA  = np.var(iyA,  axis=1)
        varAB = np.var(np.append(iyA, iyB, axis=1), axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = 1. - np.sum((iyB -iyCi)**2, axis=1) / (2.*nsa * varAB)
            isti[:,i] = np.sum((iyA -iyCi)**2, axis=1) / (2.*nsa * varA)
    elif mm == 'mai1999':
        varA  = np.var(iyA,  axis=1)
        varAB = np.var(np.append(iyA, iyB, axis=1), axis=1)
        isi  = np.empty((ntime,nn))
        isti = np.empty((ntime,nn))
        for i in range(nn):
            iyCi      = iyC[:,i,:]
            isi[:,i]  = np.mean(iyB*(iyCi - iyA), axis=1) / varAB
            isti[:,i] = np.sum((iyA -iyCi)**2, axis=1) / (2.*nsa * varA)
    else:
        raise ValueError('method unknown: {0}.'.format(method))

    if not isone:
        # simple mean
        if mean:
            msi  = np.mean(isi,  axis=0)
            msti = np.mean(isti, axis=0)

        # weighted mean
        if wmean:
            if (mm == 'mai2012'):
                denomAB = 1./np.sum(varAB)
                denomB  = 1./np.sum(varB)
                varB    = varB[:, np.newaxis]
                varAB   = varAB[:,np.newaxis]
                wsi     = np.sum(isi*varAB, axis=0) * denomAB
                wsti    = np.sum(isti*varB, axis=0) * denomB
            elif (mm == 'mai1999') | (mm == 'mai2013') | (mm == 'mai2014'):
                denomAB = 1./np.sum(varAB)
                denomA  = 1./np.sum(varA)
                varA    = varA[:, np.newaxis]
                varAB   = varAB[:,np.newaxis]
                wsi     = np.sum(isi*varAB, axis=0) * denomAB
                wsti    = np.sum(isti*varA, axis=0) * denomA
            else:
                denom = 1./np.sum(varA)
                varA  = varA[:,np.newaxis]
                wsi   = np.sum(isi*varA,  axis=0) * denom
                wsti  = np.sum(isti*varA, axis=0) * denom

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

    if len(out) == 1:
        return out[0]
    else:
        return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # import scipy.stats as stats
    # from lhs import lhs
    # from autostring import astr

    # # ------------------------------------------------------------------------
    # #                                1D
    # #
    # # seed for reproducible results in doctest
    # np.random.seed(1)
    # # number of parameter sets
    # nsets = 1000
    # # number of parameters
    # npara = 3
    # # randomly sample nsets with Latin-Hypercube
    # dist = [stats.uniform, stats.uniform, stats.uniform]
    # pars = [(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi),(-np.pi,2.0*np.pi)]
    # setsa  = lhs(dist, pars, nsets)
    # setsb  = lhs(dist, pars, nsets)
    # # Saltelli has two constants in the model
    # a = 0.5
    # b = 2.0
    # # generate model output yA and yB
    # iya = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
    # iyb = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
    # # generate model output yCi = f(B_A)
    # iyc = np.empty((npara,nsets))*np.nan
    # for i in range(npara):
    #     tmpset = np.copy(setsb)
    #     tmpset[i,:] = setsa[i,:]
    #     iyc[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])
    # # generate model output yCi = f(A_B)
    # iyc2 = np.empty((npara,nsets))*np.nan
    # for i in range(npara):
    #     tmpset = np.copy(setsa)
    #     tmpset[i,:] = setsb[i,:]
    #     iyc2[i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + b*tmpset[2,:]**4  * np.sin(tmpset[0,:])

    # # Theoretical results: has to converge for n->Infinity to this values
    # vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
    # theo_si  = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])
    # theo_sti = np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy,
    #                      (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])
    # # Theoretical results
    # print('T :: si  =',astr(theo_si,3,pp=True))
    # #T :: si  = ['0.372' '0.000' '0.000']
    # print('T :: sti =',astr(theo_sti,3,pp=True))
    # #T :: sti = ['1.000' '0.000' '0.628']
    # #
    # # 3 positional arguments
    # isi1, isti1 = sobol_index(iya,iyb,iyc2)
    # print('S :: si  =',astr(isi1,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']
    # print('S :: sti =',astr(isti1,3,pp=True))
    # #S :: sti = ['1.014' '0.000' '0.654']

    # # 3 optional arguments
    # isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc2)
    # print('S :: si  =',astr(isi2,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']

    # # 2 positional arguments
    # s2 = np.concatenate((iya, iyb, np.ravel(iyc2)))
    # ns = iya.size
    # isi3, isti3 = sobol_index(s2, ns)
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']

    # # 2 optional arguments
    # isi4, isti4 = sobol_index(s=s2, ns=ns)
    # print('S :: si  =',astr(isi4,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']

    # # 2 optional arguments and no STi output
    # isi5 = sobol_index(s=s2, ns=ns, si=True, sti=False)
    # print('S :: si  =',astr(isi5,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']

    # # 2 optional arguments and no Si output
    # isti5 = sobol_index(s=s2, ns=ns, si=False, sti=True)
    # print('S :: sti =',astr(isti5,3,pp=True))
    # #S :: sti = ['1.014' '0.000' '0.654']

    # # Method=Saltelli2008
    # s = np.concatenate((iya, iyb, np.ravel(iyc)))
    # isi3, isti3 = sobol_index(s, ns, method='Saltelli2008')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = ['0.367' '0.011' '0.016']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = [' 1.038' '-0.019' ' 0.641']

    # # Method=Homma1996
    # isi3, isti3 = sobol_index(s, ns, method='homma1996')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.356' ' 0.000' ' 0.005']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = [' 1.038' '-0.019' ' 0.641']

    # # Method=Mai2012
    # isi3, isti3 = sobol_index(s, ns, method='mai2012')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.353' ' 0.000' ' 0.005']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = ['1.037' '0.000' '0.647']

    # # Method=saltelli2010
    # isi3, isti3 = sobol_index(s2, ns, method='saltelli2010')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.348' ' 0.000' '-0.049']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = [' 0.984' ' 0.000' ' 0.633']

    # # Method=jansen1999
    # isi3, isti3 = sobol_index(s2, ns, method='jansen1999')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.319' ' 0.002' '-0.068']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = ['1.014' '0.000' '0.654']

    # # Method=mai2013
    # isi3, isti3 = sobol_index(s2, ns, method='mai2013')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = [' 0.984' ' 0.000' ' 0.633']

    # # Method=mai2014
    # isi3, isti3 = sobol_index(s2, ns, method='mai2014')
    # print('S :: si  =',astr(isi3,3,pp=True))
    # #S :: si  = [' 0.345' ' 0.000' '-0.049']
    # print('S :: sti =',astr(isti3,3,pp=True))
    # #S :: sti = ['1.014' '0.000' '0.654']

    # # ------------------------------------------------------------------------
    # #                                2D
    # #
    # # number of timepoints
    # ntime = 100
    # iya   = np.empty((ntime,nsets))*np.nan
    # iyb   = np.empty((ntime,nsets))*np.nan
    # iyc   = np.empty((ntime,npara,nsets))*np.nan
    # iyc2  = np.empty((ntime,npara,nsets))*np.nan
    # s     = np.empty((ntime,(npara+2)*nsets))*np.nan
    # s2    = np.empty((ntime,(npara+2)*nsets))*np.nan
    # ns    = nsets
    # theo_si  = []
    # theo_sti = []
    # for t in range(ntime):
    #     a = 1.0*t + 50.0
    #     b = 2.0
    #     iya[t,:] = np.sin(setsa[0,:])  + a*(np.sin(setsa[1,:]))**2  + b*setsa[2,:]**4  * np.sin(setsa[0,:])
    #     iyb[t,:] = np.sin(setsb[0,:])  + a*(np.sin(setsb[1,:]))**2  + b*setsb[2,:]**4  * np.sin(setsb[0,:])
    #     for i in range(npara):
    #         tmpset = np.copy(setsb)
    #         tmpset[i,:] = setsa[i,:]
    #         iyc[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + (
    #                      b*tmpset[2,:]**4  * np.sin(tmpset[0,:]))
    #     for i in range(npara):
    #         tmpset = np.copy(setsa)
    #         tmpset[i,:] = setsb[i,:]
    #         iyc2[t,i,:] = np.sin(tmpset[0,:])  + a*(np.sin(tmpset[1,:]))**2  + (
    #                       b*tmpset[2,:]**4  * np.sin(tmpset[0,:]))
    #     s[t,:]  = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc[t,:,:])))
    #     s2[t,:] = np.concatenate((iya[t,:], iyb[t,:], np.ravel(iyc2[t,:,:])))
    #     vy = 0.5 + a**2/8.0 + b*np.pi**4/5.0+b**2*np.pi**8/18.0
    #     theo_si  += [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0)/vy, (a**2/8.0)/vy, 0.0])]
    #     theo_sti += [np.array([(0.5+b*np.pi**4/5.0+b**2*np.pi**8/50.0+8.0/225.0*b**2*np.pi**8)/vy,
    #                            (a**2/8.0)/vy, (8.0/225.0*b**2*np.pi**8)/vy])]

    # # 1D model output --> ya, yb are 2D and yc is 3D
    # isi1, isti1 = sobol_index(iya,iyb,iyc2)

    # # SI :: 1st time point
    # print('S :: si[t0]  =',astr(isi1[0,:],3,pp=True))
    # #S :: si[t0]  = [' 0.295' ' 0.147' '-0.056']
    # print('T :: si[t0]  =',astr(theo_si[0],3,pp=True))
    # #T :: si[t0]  = ['0.325' '0.127' '0.000']
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi1[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.293' ' 0.156' '-0.055']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']
    # # STI :: 1st time point
    # print('S :: sti[t0] =',astr(isti1[0,:],3,pp=True))
    # #S :: sti[t0] = ['0.881' '0.126' '0.568']
    # print('T :: sti[t0] =',astr(theo_sti[0],3,pp=True))
    # #T :: sti[t0] = ['0.873' '0.127' '0.548']
    # # STI :: 3rd time point
    # print('S :: sti[t2] =',astr(isti1[2,:],3,pp=True))
    # #S :: sti[t2] = ['0.871' '0.134' '0.562']
    # print('T :: sti[t2] =',astr(theo_sti[2],3,pp=True))
    # #T :: sti[t2] = ['0.864' '0.136' '0.543']

    # # Mean
    # isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc2, mean=True)
    # print('S :: si_m    =',astr(msi2,3,pp=True))
    # #S :: si_m    = [' 0.231' ' 0.377' '-0.037']
    # print('S :: sti_m   =',astr(msti2,3,pp=True))
    # #S :: sti_m   = ['0.643' '0.351' '0.415']
    # #
    # # Weighted mean
    # isi4, isti4, wsi2, wsti2 = sobol_index(s2, ns, wmean=True)
    # print('S :: si_w    =',astr(wsi2,3,pp=True))
    # #S :: si_w    = [' 0.224' ' 0.404' '-0.035']
    # print('S :: sti_w   =',astr(wsti2,3,pp=True))
    # #S :: sti_w   = ['0.616' '0.377' '0.397']

    # # 2 positional arguments
    # isi3, isti3 = sobol_index(s=s2, ns=ns)
    # # SI :: 3rd time point
    # print('S :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #S :: si[t2]  = [' 0.293' ' 0.156' '-0.055']
    # print('T :: si[t2]  =',astr(theo_si[2],3,pp=True))
    # #T :: si[t2]  = ['0.321' '0.136' '0.000']

    # # Method=Saltelli2008
    # isi3, isti3 = sobol_index(s, ns, method='saltelli2008')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = ['0.354' '0.130' '0.011']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = ['0.352' '0.139' '0.010']

    # # Method=Homma1996
    # isi3, isti3 = sobol_index(s, ns, method='homma1996')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.342' ' 0.118' '-0.002']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.340' ' 0.127' '-0.002']

    # # Method=Mai2012
    # isi3, isti3 = sobol_index(s, ns, method='Mai2012')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.336' ' 0.116' '-0.002']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.334' ' 0.125' '-0.002']

    # # Method=Saltelli2010
    # isi3, isti3 = sobol_index(s2, ns, method='saltelli2010')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.294' ' 0.142' '-0.031']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.291' ' 0.152' '-0.031']

    # # Method=Jansen1999
    # isi3, isti3 = sobol_index(s2, ns, method='jansen1999')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.284' ' 0.133' '-0.072']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.281' ' 0.143' '-0.072']

    # # Method=Mai2013
    # isi3, isti3 = sobol_index(s2, ns, method='Mai2013')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.290' ' 0.140' '-0.031']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.286' ' 0.149' '-0.030']

    # # Method=Mai2014
    # isi3, isti3 = sobol_index(s2, ns, method='mai2014')
    # print('Sal :: si[t0]  =',astr(isi3[0,:],3,pp=True))
    # #Sal :: si[t0]  = [' 0.295' ' 0.147' '-0.056']
    # print('Sal :: si[t2]  =',astr(isi3[2,:],3,pp=True))
    # #Sal :: si[t2]  = [' 0.293' ' 0.156' '-0.055']
