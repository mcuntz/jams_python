#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.const import eps

def zacharias(h, clay, sand, db, params=None, thetar=False, thetas=False, lnalpha=False, n=False):
    """
        Soil water content with the van Genuchten equation and
        the pedotransfer functions of Zacharias et al. (2007).


        Definition
        ----------
        def zacharias(h, clay, sand, db, params=None, thetar=False, thetas=False, lnalpha=False, n=False):


        Input
        -----
        h          pressure head, scalar or array [cm], 0=saturation, 15000=wilting point
        clay       clay content, scalar or array [%, i.e. 0-100]
        sand       sand content, scalar or array [%, i.e. 0-100]
        db         bulk density, scalar or array [g/cm3], quartz=2.65


        Optional Input
        --------------
        params     Parameter for Zacharias et al. (2007) pedotransfer functions
                   If None, values from Zacharias et al. will be taken that are different
                   between sandy and non-sandy soil (<66.5% sand)


        Options
        -------
        thetar     If True, outputs residual water content thetar as well [m3 m-3]
        thetas     If True, outputs saturation water content thetas as well [m3 m-3]
        lnalpha    If True, outpus logarithm of shape parameter alpha as well [1/cm]
        n          If True, output exponent n as well


        Output
        ------
        Soil water content theta [m^3 m^-3]


        Restrictions
        ------------
        Does not check the validity of the parameter set, i.e. negative soil moistures
        can occur, for example.
        Use zacharias_check to check the parameter set first.


        Examples
        --------
        >>> h = np.array([0.0000000, 0.0000000, 10.000000, 31.622777,
        ...               100.00000, 199.52623, 199.52623,
        ...               501.18723, 2511.8864, 15848.932])
        >>> sand = np.array([12.800000, 61.600000, 17.200000, 85.800000,
        ...                  16.500000, 12.800000, 61.600000,
        ...                  17.200000, 85.800000, 16.500000])
        >>> clay = np.array([30.500000, 17.200000, 25.500000, 8.9000000,
        ...                  28.100000, 30.500000, 17.200000,
        ...                  25.500000, 8.9000000, 28.100000])
        >>> rho = np.array([1.2100000, 1.3400000, 1.4600000, 1.6300000,
        ...                 1.3000000, 1.2100000, 1.3400000,
        ...                 1.4600000, 1.6300000, 1.3000000])
        >>> from autostring import astr
        >>> print(astr(zacharias(h, clay, sand, rho),3,pp=True))
        ['0.500' '0.453' '0.421' '0.245' '0.393' '0.381' '0.285' '0.313' '0.039' '0.221']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jun 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Nov 2016 - const.tiny -> const.eps
    """
    #
    # Check input
    ih = np.where(h==0., eps, h)
    if np.any(ih < 0.) | np.any(ih > 1e6):
        raise ValueError('h must be >=0 and <= 1e6 (=pf6)')
    iclay = np.where(clay==0., eps, clay)
    if np.any(iclay < 0.) | np.any(iclay > 100.):
        raise ValueError('clay must be >=0 and <= 100.')
    isand = np.where(sand==0., eps, sand)
    if np.any(isand < 0.) | np.any(isand > 100.):
        raise ValueError('sand must be >=0 and <= 100.')
    idb = np.array(db)
    if np.any(idb < 0.) | np.any(idb > 2.65):
        raise ValueError('db must be >=0 and <= 2.65.')
    nn = np.size(isand)
    if (np.size(iclay) != nn) | (np.size(idb) != nn) | (np.size(ih) != nn):
        raise ValueError('h, sand, clay, and db must have the same sizes.')
    if params is not None:
        if np.size(params) != 15:
            raise ValueError('size(params) must be 15.')
    # save output shape
    ns    = np.shape(isand)
    iclay = np.ravel(iclay)
    isand = np.ravel(isand)
    idb   = np.ravel(idb)
    # Take right params
    par0  = np.empty(nn)
    par1  = np.empty(nn)
    par2  = np.empty(nn)
    par3  = np.empty(nn)
    par4  = np.empty(nn)
    par5  = np.empty(nn)
    par6  = np.empty(nn)
    par7  = np.empty(nn)
    par8  = np.empty(nn)
    par9  = np.empty(nn)
    par10 = np.empty(nn)
    par11 = np.empty(nn)
    par12 = np.empty(nn)
    par13 = np.empty(nn)
    par14 = np.empty(nn)
    if params is not None:
        # Either params given
        par0[:]  = params[0]
        par1[:]  = params[1]
        par2[:]  = params[2]
        par3[:]  = params[3]
        par4[:]  = params[4]
        par5[:]  = params[5]
        par6[:]  = params[6]
        par7[:]  = params[7]
        par8[:]  = params[8]
        par9[:]  = params[9]
        par10[:] = params[10]
        par11[:] = params[11]
        par12[:] = params[12]
        par13[:] = params[13]
        par14[:] = params[14]
    else:
        # or take Zacharias
        parclay = np.array([ 0.,     0.,     0.,
                             0.788,  0.001, -0.263,
                            -0.648,  0.044, -3.168,  0.023,
                             1.392,  1.212, -0.704, -0.418, -0.024])
        parsand = np.array([ 0.,     0.,     0.,
                             0.890, -0.001, -0.322,
                            -4.197,  0.076, -0.276, 0.013,
                            -2.562,  3.750, -0.016, 7e-9,  4.004])
        for i in range(nn):
            if isand[i] < 66.5:
                par = parclay
            else:
                par = parsand
            par0[i]  = par[0]
            par1[i]  = par[1]
            par2[i]  = par[2]
            par3[i]  = par[3]
            par4[i]  = par[4]
            par5[i]  = par[5]
            par6[i]  = par[6]
            par7[i]  = par[7]
            par8[i]  = par[8]
            par9[i]  = par[9]
            par10[i] = par[10]
            par11[i] = par[11]
            par12[i] = par[12]
            par13[i] = par[13]
            par14[i] = par[14]

    # Zacharias pedotransfer
    ithetar = par0  + par1*iclay                        + par2*idb
    ithetas = par3  + par4*iclay                        + par5*idb
    ilna    = par6  + par7*iclay                        + par8*idb + par9*isand
    inn     = par10 + par11*np.exp(par12*np.log(iclay))            + par13*np.exp(par14*np.log(isand))
    imm     = 1. - 1./np.where(inn != 0., inn, 1e-3)
    # van Genuchten sign
    # limit exp to 600 and log to eps so that no over- and underflows occur
    expmax  = 600.
    lnah    = np.log(np.maximum(np.exp(ilna)*ih, eps))
    ahn     = np.exp(np.minimum(inn*lnah, expmax))
    denom   = np.maximum(np.exp(np.minimum(imm*np.log(1.+ahn), expmax)), eps)
    itheta  = np.where(ih <= eps, ithetas, ithetar + (ithetas-ithetar)/denom)
    # Output
    itheta = np.reshape(itheta, ns)
    if nn==1: itheta = np.float(itheta)
    if (thetar==False) & (thetas==False) & (lnalpha==False) & (n==False):
        return itheta
    else:
        out = [itheta]
        if thetar==True:
            ithetar = np.reshape(ithetar, ns)
            if nn==1: ithetar = np.float(ithetar)
            out = out + [ithetar]
        if thetas==True:
            if nn==1: ithetas = np.float(ithetas)
            ithetas = np.reshape(ithetas, ns)
            out = out + [ithetas]
        if lnalpha==True:
            if nn==1: ilna = np.float(ilna)
            ilna = np.reshape(ilna, ns)
            out = out + [ilna]
        if n==True:
            if nn==1: inn = np.float(inn)
            inn = np.reshape(inn, ns)
            out = out + [inn]
        return out

def zacharias_check(params, sand=None, clay=None):
    """
        Checks if a given parameter set is valid for all possible soils with
        the van Genuchten equation and the pedotransfer functions of Zacharias et al. (2007).

        Definition
        ----------
        def zacharias_check(params, sand=None, clay=None):


        Input
        -----
        params     array[15] with parameters a1, b1, ..., d4 of Zacharias et al. (2007)


        Optional Input
        --------------
        clay       If given: 1. < clay < 99. then calc sand content < clay
        sand       If given: 1. < sand < 99. then calc sand content > sand


        Output
        ------
        Boolean    True:  valid parameter set for all soils
                   False: not valid in at least one extreme case


        Examples
        --------
        >>> parclay = np.array([ 0.,     0.,     0.,
        ...                      0.788,  0.001, -0.263,
        ...                     -0.648,  0.044, -3.168,  0.023,
        ...                      1.392,  1.212, -0.704, -0.418, -0.024])
        >>> print(zacharias_check(parclay))
        True
        >>> print(zacharias_check(parclay, clay=66))
        True
        >>> parsand = np.array([ 0.,     0.,     0.,
        ...                      0.890, -0.001, -0.322,
        ...                     -4.197,  0.076, -0.276, 0.013,
        ...                     -2.562,  3.750, -0.016, 7e-9,  4.004])
        >>> print(zacharias_check(parsand))
        False
        >>> print(zacharias_check(parsand, sand=66))
        True


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2016 Matthias Cuntz - mc (at) macu (dot) de

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
    """
    #
    # Check input
    if params is not None:
        if np.size(params) != 15: raise ValueError('size(params) must be 15.')
    # Check ranges
    isgood = True
    if (params[0]  <   0.)   |  (params[0]  >  1.):   isgood = False
    if (params[1]  <  -0.01) |  (params[1]  >  0.01): isgood = False
    if (params[2]  <  -0.5)  |  (params[2]  >  0.5):  isgood = False
    if (params[3]  <   0.)   |  (params[3]  >  1.):   isgood = False
    if (params[4]  <  -0.01) |  (params[4]  >  0.01): isgood = False
    if (params[5]  <  -0.5)  |  (params[5]  >  0.5):  isgood = False
    if (params[6]  < -10.)   |  (params[6]  > 10.):   isgood = False
    if (params[7]  <  -1.)   |  (params[7]  >  1.):   isgood = False
    if (params[8]  <  -5.)   |  (params[8]  >  5.):   isgood = False
    if (params[9]  <  -1.)   |  (params[9]  >  1.):   isgood = False
    if (params[10] < -10.)   |  (params[10] > 10.):   isgood = False
    if (params[11] < -10.)   |  (params[11] > 10.):   isgood = False
    if (params[12] <  -5.)   |  (params[12] >  5.):   isgood = False
    if (params[13] < -10.)   |  (params[13] > 10.):   isgood = False
    if (params[14] <  -5.)   |  (params[14] >  5.):   isgood = False
    # Soil ranges
    h    = [1., 1e6]
    if clay is not None:
        rclay = [98., 1.,   99.-clay]
        rsand = [1.,  clay, clay]
    elif sand is not None:
        rclay = [1.,  1.,   99.-sand]
        rsand = [98., sand, sand]
    else:
        rclay = [1.,  1., 98.]
        rsand = [98., 1.,  1.]
    db   = [0.5, 2.3]
    nn   = 2*3*2
    ih    = np.empty(nn)
    iclay = np.empty(nn)
    isand = np.empty(nn)
    idb   = np.empty(nn)
    zaehl = 0
    for i in [0,1]:
        for j in [0,1,2]:
            for k in [0,1]:
                ih[zaehl]    = h[i]
                iclay[zaehl] = rclay[j]
                isand[zaehl] = rsand[j]
                idb[zaehl]   = db[k]
                zaehl += 1
    # Calc Zacharias
    itheta, ithetar, ithetas, ilna, inn = zacharias(ih, iclay, isand, idb, params=params,
                                                    thetar=True, thetas=True, lnalpha=True, n=True)
    #ia  = np.exp(ilna)
    #imm = 1. - 1./inn
    if np.any(ithetar < 0.) | np.any(ithetar > 1.): isgood = False
    if np.any(ithetas < 0.) | np.any(ithetas > 1.) | np.any(ithetas < ithetar): isgood = False
    #if np.any(ilna > 0.): isgood = False
    if np.any(inn < 1.) | np.any(inn > 10.): isgood = False
    #if np.any(imm < 0.): isgood = False
    #if np.any(ia > 1.): isgood = False
    if np.any(itheta < ithetar) | np.any(itheta > ithetas): isgood = False

    return isgood


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # parclay = np.array([ 0.,     0.,     0.,
    #                      0.788,  0.001, -0.263,
    #                      -0.648,  0.044, -3.168,  0.023,
    #                      1.392,  1.212, -0.704, -0.418, -0.024])
    # print zacharias_check(parclay)
    # print zacharias_check(parclay, clay=66)
    # parsand = np.array([ 0.,     0.,     0.,
    #                      0.890, -0.001, -0.322,
    #                      -4.197,  0.076, -0.276, 0.013,
    #                      -2.562,  3.750, -0.016, 7e-9,  4.004])
    # print zacharias_check(parsand)
    # print zacharias_check(parsand, sand=66)

    # h = np.array([0.0000000, 0.0000000, 10.000000, 31.622777,
    #               100.00000, 199.52623, 199.52623,
    #               501.18723, 2511.8864, 15848.932])
    # sand = np.array([12.800000, 61.600000, 17.200000, 85.800000,
    #                  16.500000, 12.800000, 61.600000,
    #                  17.200000, 85.800000, 16.500000])
    # clay = np.array([30.500000, 17.200000, 25.500000, 8.9000000,
    #                  28.100000, 30.500000, 17.200000,
    #                  25.500000, 8.9000000, 28.100000])
    # rho = np.array([1.2100000, 1.3400000, 1.4600000, 1.6300000,
    #                 1.3000000, 1.2100000, 1.3400000,
    #                 1.4600000, 1.6300000, 1.3000000])
    # # # From IDL code
    # # zthetafit = np.array([0.50026998, 0.45277997, 0.42107488, 0.24487569,
    # #                       0.39254300, 0.38144154, 0.28500143,
    # #                       0.31320843, 0.038888339, 0.22115273])
    # print zacharias(h, clay, sand, rho)
    # # [ 0.50027     0.45278     0.42107491  0.24487575  0.39254299  0.38144151
    # #   0.28500141  0.31320842  0.03888839  0.22115267]

