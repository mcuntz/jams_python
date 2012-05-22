#!/usr/bin/env python
import numpy as np

def zacharias(h, clay, sand, db, params=None, thetar=False, thetas=False, lnalpha=False, n=False, check=True):
    """
        Soil water content with the van Genuchten equation and
        the pedotransfer functions of Zacharias et al. (2007).

        Definition
        ----------
        def zacharias(h, clay, sand, db, params=None,
                      thetar=False, thetas=False, lnalpha=False, n=False, check=True):


        Input
        -----
        h          number array


        Optional Input
        --------------
        powten     Power of ten array
                   If missing, simple round (ceil, floor) is taken.
        ceil       ceil instead of round to the nearest power of ten
        floor      floor instead of round to the nearest power of ten


        Output
        ------
        Soil water content: theta [m^3 m^-3].


        Restrictions
        ------------
        Powten is exactly opposite of decimal keyword of numpy.around.

        From numpy.around documentation:
        'For values exactly halfway between rounded decimal values,
        Numpy rounds to the nearest even value. Thus 1.5 and 2.5 round to 2.0,
        -0.5 and 0.5 round to 0.0, etc. Results may also be surprising due to the
        inexact representation of decimal fractions in the IEEE floating point
        standard and errors introduced when scaling by powers of ten.'


        Examples
        --------
        >>> h = np.array([0.0000000, 0.0000000, 10.000000, 31.622777, \
                          100.00000, 199.52623, 199.52623, \
                          501.18723, 2511.8864, 15848.932])
        >>> sand = np.array([12.800000, 61.600000, 17.200000, 85.800000, \
                             16.500000, 12.800000, 61.600000, \
                             17.200000, 85.800000, 16.500000])
        >>> clay = np.array([30.500000, 17.200000, 25.500000, 8.9000000, \
                             28.100000, 30.500000, 17.200000, \
                             25.500000, 8.9000000, 28.100000])
        >>> rho = np.array([1.2100000, 1.3400000, 1.4600000, 1.6300000, \
                            1.3000000, 1.2100000, 1.3400000, \
                            1.4600000, 1.6300000, 1.3000000])
        >>> print zacharias(h, clay, sand, rho)
        [ 0.50027     0.45278     0.42107491  0.24487575  0.39254299  0.38144151
          0.28500141  0.31320842  0.03888839  0.22115267]


        History
        -------
        Written, MC, Jun 2011
    """
    #
    # Check input
    tiny = np.finfo(np.float).eps
    ih = np.where(h==0., tiny, h)
    if np.any(ih < 0.) | np.any(ih > 1e6):
        raise ValueError('h must be >=0 and <= 1e6 (=pf6)')
    iclay = np.where(clay==0., tiny, clay)
    if np.any(iclay < 0.) | np.any(iclay > 100.):
        raise ValueError('clay must be >=0 and <= 100.')
    isand = np.where(sand==0., tiny, sand)
    if np.any(isand < 0.) | np.any(isand > 100.):
        raise ValueError('sand must be >=0 and <= 100.')
    idb = np.array(db)
    if np.any(idb < 0.) | np.any(idb > 2.65):
        raise ValueError('db must be >=0 and <= 2.65.')
    nn = np.size(isand)
    if (np.size(iclay) != nn) | (np.size(idb) != nn) | (np.size(ih) != nn):
        raise ValueError('h, sand, clay, and db must have the same sizes.')
    if params != None:
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
    if params != None:
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
        for i in xrange(nn):
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
    imm     = 1. - 1./inn
    # van Genuchten
    itheta  = np.where(ih <= tiny, ithetas, 
                       ithetar + (ithetas-ithetar)/np.exp(imm*np.log(1.+np.exp(inn*np.log(np.exp(ilna)*ih)))))    
    # Output
    itheta = np.reshape(itheta, ns)
    if nn==1: itheta = itheta[0]
    if (thetar==False) & (thetas==False) & (lnalpha==False) & (n==False):
        return itheta
    else:
        out = [itheta]
        if thetar==True:
            ithetar = np.reshape(ithetar, ns)
            if nn==1: ithetar = ithetar[0]
            out = out + [ithetar]
        if thetas==True:
            if nn==1: ithetas = ithetas[0]
            ithetas = np.reshape(ithetas, ns)
            out = out + [ithetas]
        if lnalpha==True:
            if nn==1: ilna = ilna[0]
            ilna = np.reshape(ilna, ns)
            out = out + [ilna]
        if n==True:
            if nn==1: inn = inn[0]
            inn = np.reshape(inn, ns)
            out = out + [inn]
        return out


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # h = np.array([0.0000000, 0.0000000, 10.000000, 31.622777, \
    #               100.00000, 199.52623, 199.52623, \
    #               501.18723, 2511.8864, 15848.932])
    # sand = np.array([12.800000, 61.600000, 17.200000, 85.800000, \
    #                  16.500000, 12.800000, 61.600000, \
    #                  17.200000, 85.800000, 16.500000])
    # clay = np.array([30.500000, 17.200000, 25.500000, 8.9000000, \
    #                  28.100000, 30.500000, 17.200000, \
    #                  25.500000, 8.9000000, 28.100000])
    # rho = np.array([1.2100000, 1.3400000, 1.4600000, 1.6300000, \
    #                 1.3000000, 1.2100000, 1.3400000, \
    #                 1.4600000, 1.6300000, 1.3000000])
    # # # From IDL code
    # # zthetafit = np.array([0.50026998, 0.45277997, 0.42107488, 0.24487569, \
    # #                       0.39254300, 0.38144154, 0.28500143, \
    # #                       0.31320843, 0.038888339, 0.22115273])
    # print zacharias(h, clay, sand, rho)
    # # [ 0.50027     0.45278     0.42107491  0.24487575  0.39254299  0.38144151
    # #   0.28500141  0.31320842  0.03888839  0.22115267]
