#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from jams.date2dec   import date2dec
from jams.dec2date   import dec2date
from jams.autostring import astr
from jams.eddybox    import itc
from jams.eddybox    import spikeflag
from jams.eddybox    import ustarflag

def fluxflag(fluxfile, metfile, outdir, swdr, T, lat, delimiter=[',',','],
             skiprows=[1,1], format=['ascii','ascii'], limit=0.3,
             analyzer='LI7000', spike_f=True, itc_f=True, ustar_f=True,
             novalue=-9999, plot=False):    
    '''
    Quality flag calculation for Eddy Covariance data originating from EddyFlux.
    Including stationarity, exceeding limits, exceeding change rates, exceeding
    variances, integral turbulence characteristics, spike detection and ustar
    thresholds. Missing values are flagged with 2. Ustar thresholding works
    ONLY for a data set of ONE FULL year.  
    
    
    Definition
    ----------
    fluxflag(fluxfile, metfile, outdir, Rg, T, lat, delimiter=[',',','],
             skiprows=[1,1], format=['ascii','ascii'], limit=0.3,
             analyzer='LI7000', spike_f=True, itc_f=True, ustar_f=True,
             novalue=-9999, plot=False):
    
    
    Input
    ----- 
    fluxfile    str, path and file name of EddyFlux output file
    metfile     str, path and file name of the meteorology file (must be in
                sync with fluxfile)
    outdir      str, path of the output folder
    swdr        int, column number of short wave incoming radiation in metfile,
                column number starts with 0 which is first data column.
                Used for swdr>0=isday
    T           int, column number of air temperature in metfile, column number
                starts with 0 which is first data column. Used for ustarflag
    lat         latitude of tower position in decimal degrees. Used for itcflag
    
                        
    Optional Input
    --------------
    delimiter   list of str, delimiters of fluxfile and metfile
                (default: [',',','])
    skiprows    list of int, lines to skip at the beginning of fluxfile and
                metfile, e.g. header lines (default: [1,1])
    format      list of str, time formats of fluxfile and metfile, 'ascii' and 
                'eng' possible (default: ['ascii','ascii'])
    limit       float, relative deviation limit from the itc model above which
                values are flagged with itcflag. e.g. 0.3 means 30% above and
                30% below the model (default: 0.3)
    analyzer    str, analyzer type. 'LI7000' and 'LI7500' possible. Used for
                exceeding variance limits. (default: 'LI7000')
    spike_f     bool, if True spike detection and flagging is performed
    itc_f       bool, if True itc calculation and flagging is performed
    ustar_f     bool, if True ustar threshold calculation and flagging is
                performed
    novalue     int/float, missing value of fluxfile and metfile
                (default: -9999)
    plot        bool, if True subroutines (itc, spike and ustar) perform
                plotting (default: False)
    
    
    Output
    ------
    flags.csv     file containing all calculated flags
    fluxflags.csv file containing fluxes with the sum of all respective flags
                  for each flux. Where flag>2, flag=2.
    flags.log     file containing statistics for each flag
    
    
    Restrictions
    ------------
    - ustar flagging works ONLY for a data set of ONE FULL year
    - ustar flagging works ONLY for half hourly time steps
    
    
    License
    -------
    This file is part of the JAMS Python package.
    
    It is NOT released under the GNU Lesser General Public License, yet.
    
    If you use this routine, please contact Arndt Piayda.
    
    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Aug 2014
    '''
    ###########################################################################
    # reading input files
    d = np.loadtxt(fluxfile, dtype='|S100', delimiter=delimiter[0], skiprows=skiprows[0])
    m = np.loadtxt(metfile,  dtype='|S100', delimiter=delimiter[1], skiprows=skiprows[1])
    
    if format[0]=='ascii':
        datev   = date2dec(ascii=d[:,0])
    elif format[0]=='eng':
        datev   = date2dec(eng=d[:,0])
    else:
        raise ValueError('fluxflag: unknown format')    
    if format[1]=='ascii':
        datem   = date2dec(ascii=m[:,0])
    elif format[1]=='eng':
        datem   = date2dec(eng=m[:,0])
    else:
        raise ValueError('fluxflag: unknown format')
    
    val = np.where(d[:,1:]=='', str(novalue), d[:,1:]).astype(np.float)
    val = np.where(val==novalue, np.NaN, val)
    met = np.where(m[:,1:]=='', str(novalue), m[:,1:]).astype(np.float)
    met = np.where(met==novalue, np.NaN, met)
    print('LOADING FILES COMPLETED')
    
    ###########################################################################
    # calculate isday
    isday = met[:,swdr]>0.
    fluxes = val[:,[28,29,30,31,5]] # corr. H, corr. LE, corr. E, corr. C, tau
    zeta  = val[:,26]
    varu  = val[:,41]
    ustar = val[:,4]
    varw  = val[:,8]
    vart  = val[:,9]
    rho   = val[:,62]
    H     = val[:,28]
    T     = met[:,T]
        
    ###########################################################################
    # calculate standard flags
    sttest2  = np.where(np.isnan(val[:,37:38]), np.NaN,
                        np.unpackbits(val[:,37:38].astype('uint8'), axis=1))
    h2oexli  = np.where(np.isnan(val[:,56]), np.NaN,
                        np.where(val[:,56] > 100, 2, 0))
    c02exli  = np.where(np.isnan(val[:,55]), np.NaN,
                        np.where(val[:,55] > 100, 2, 0))
    h2oexcr  = np.where(np.isnan(val[:,50]), np.NaN,
                        np.where(val[:,50] > 50, 2, 0))
    co2excr  = np.where(np.isnan(val[:,49]), np.NaN,
                        np.where(val[:,49] > 50, 2, 0))
    texcr    = np.where(np.isnan(val[:,48]), np.NaN,
                        np.where(val[:,48] > 50, 2, 0))
    tauexcr  = np.where(np.isnan(val[:,45]) ^ np.isnan(val[:,47]),
                        np.NaN, np.where((val[:,45]>50) ^ (val[:,47]>50), 2, 0))
    zetaexcr = np.where(np.isnan(val[:,45]) ^ np.isnan(val[:,47]) ^ np.isnan(val[:,48]),
                        np.NaN, np.where((val[:,45] > 50) ^ (val[:,47] > 50) ^ (val[:,48] > 50), 2, 0))
    wexvar   = np.where(np.isnan(val[:,8]), np.NaN,
                        np.where(val[:,8] > 3, 1, 0))
    texvar   = np.where(np.isnan(val[:,9]), np.NaN,
                        np.where(val[:,9] > 3, 1, 0))
    if analyzer=='LI7000':
        h2oexvar3  = np.where(np.isnan(val[:,11]), np.NaN,
                              np.where(val[:,11] > 3, 1, 0))
        h2oexvar10 = np.where(np.isnan(val[:,11]), np.NaN,
                              np.where(val[:,11] > 10, 1, 0))
        co2exvar3  = np.where(np.isnan(val[:,10]), np.NaN,
                              np.where(val[:,10] > 3, 1, 0))
    elif analyzer=='LI7500':
        h2oexvar3  = np.where(np.isnan(val[:,11]), np.NaN,
                              np.where(val[:,11] > 5000, 1, 0))
        h2oexvar10 = np.where(np.isnan(val[:,11]), np.NaN,
                              np.where(val[:,11] > 17300, 1, 0))
        co2exvar3  = np.where(np.isnan(val[:,10]), np.NaN,
                              np.where(val[:,10] > 3.5, 1, 0))
    else:
        raise ValueError('fluxflag: unknown analyzer')
    wvar0    = np.where(np.isnan(val[:,8]), np.NaN,
                        np.where(val[:,8] < 0.00001, 2, 0))
    h2ovar0  = np.where(np.isnan(val[:,11]), np.NaN,
                        np.where(val[:,11] < 0.00001, 2, 0))
    co2var0  = np.where(np.isnan(val[:,10]), np.NaN,
                        np.where(val[:,10] < 0.00001, 2, 0))
    print('CALCULATE STANDARD FLAGS COMPLETE')
    
    ###########################################################################
    # itc flag calculation
    if itc_f:
        itcu, itcw, itct = itc.itc(H, zeta, ustar, varu, varw, vart, rho, lat,
                               limit, '%s/itc'%outdir, plot=plot)
        print('ITC CALCUATION COMPLETED')
    else:
        itcu = np.zeros_like(H, dtype=np.int)
        itcw = np.zeros_like(H, dtype=np.int)
        itct = np.zeros_like(H, dtype=np.int)
    
    ###########################################################################
    # summing flags for each flux
    Hflag   = np.nansum(np.vstack((sttest2[:,6], itcw, texcr, wexvar, texvar,
                                   wvar0)).transpose(), 1).astype(int)
    Hflag[Hflag>1]=2
    LEflag  = np.nansum(np.vstack((sttest2[:,5], itcw, h2oexli, h2oexcr, wexvar,
                                   h2oexvar3, wvar0, h2ovar0)).transpose(), 1).astype(int)
    LEflag[LEflag>1]=2
    Eflag   = np.nansum(np.vstack((sttest2[:,5], itcw, h2oexli, h2oexcr, wexvar,
                                   h2oexvar3, wvar0, h2ovar0)).transpose(), 1).astype(int)
    Eflag[Eflag>1]=2
    Cflag   = np.nansum(np.vstack((sttest2[:,4], itcw, c02exli, co2excr, wexvar,
                                   h2oexvar10, wvar0, co2var0)).transpose(), 1).astype(int) #, co2exvar3
    Cflag[Cflag>1]=2
    Tauflag = np.nansum(np.vstack((sttest2[:,3], itcu, itcw, tauexcr, wexvar,
                                   wvar0)).transpose(), 1).astype(int)
    Tauflag[Tauflag>1]=2
        
    ###########################################################################
    # spike detection
    #inflag = np.zeros_like(fluxes, dtype=np.int)
    inflag = np.vstack((Hflag, LEflag, Eflag, Cflag, Tauflag)).transpose()
    if spike_f:
        spikef = spikeflag.spikeflag(datev, fluxes, inflag, isday, '%s/spike'%outdir, z=7,
                           deriv=1, udef=np.NaN, spike_v=2, plot=plot)
        print('SPIKE DETECTION COMPLETED')
    else:
        spikef = np.zeros_like(inflag, dtype=np.int)

    inflag += spikef
    
    ###########################################################################
    # ustar flagging
    if ustar_f:
        C_ustar_T             = np.vstack((fluxes[:,3], ustar, T)).transpose()
        C_ustar_T_flags       = np.isnan(C_ustar_T).astype(np.int)
        C_ustar_T_flags[:,0] += inflag[:,3]
        ustarf = ustarflag.ustarflag(datev, C_ustar_T, C_ustar_T_flags, isday,
                           '%s/ustar'%outdir, ustar_v=2, plot=plot)
        print('USTAR FLAGGING COMPLETED')
    else:
        ustarf = np.zeros_like(datev, dtype=np.int)

    ###########################################################################
    # the big gathering :-D
    header = np.array(['sttest2_H', ' sttest2_LE/E', ' sttest2_C', ' sttest2_tau',
                       ' h2oexli', ' c02exli', ' h2oexcr', ' co2excr', ' texcr',
                       ' tauexcr', ' zetaexcr', ' wexvar', ' texvar', ' h2oexvar3',
                       ' h2oexvar10', ' co2exvar3', ' wvar0', ' h2ovar0', ' co2var0',
                       ' itcu', ' itcw', ' itct', ' spike_H', ' spike_LE', ' spike_E',
                       ' spike_C', ' spike_tau', ' ustar'])
    maxlen = [len(x) for x in header]
    flags  = np.vstack((sttest2[:,6], sttest2[:,5], sttest2[:,4], sttest2[:,3],
                        h2oexli, c02exli, h2oexcr, co2excr, texcr, tauexcr,
                        zetaexcr, wexvar, texvar, h2oexvar3, h2oexvar10,
                        co2exvar3, wvar0, h2ovar0, co2var0, itcu, itcw, itct,
                        spikef[:,0], spikef[:,1], spikef[:,2], spikef[:,3],
                        spikef[:,4], ustarf)).transpose()
    flags = np.where(np.isnan(flags), 2, flags).astype(np.int)
    flags_str = np.zeros_like(flags, dtype='|S100')
    for i, item in enumerate(maxlen):
        flags_str[:,i] = astr(flags[:,i], prec=item)
    
    # save flag file
    output = np.hstack((d[:,0:1], flags_str))
    np.savetxt('%s/flags.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header))[np.newaxis,:],
                          output)), '%s', delimiter=',')

    ###########################################################################
    # summing flags for each flux
    Hflag   = np.nansum(np.vstack((sttest2[:,6], itcw, texcr, wexvar, texvar,
                                   wvar0, spikef[:,0])).transpose(), 1).astype(int)
    Hflag[Hflag>1]=2
    LEflag  = np.nansum(np.vstack((sttest2[:,5], itcw, h2oexli, h2oexcr, wexvar,
                                   h2oexvar3, wvar0, h2ovar0, spikef[:,1])).transpose(), 1).astype(int)
    LEflag[LEflag>1]=2
    Eflag   = np.nansum(np.vstack((sttest2[:,5], itcw, h2oexli, h2oexcr, wexvar,
                                   h2oexvar3, wvar0, h2ovar0, spikef[:,2])).transpose(), 1).astype(int)
    Eflag[Eflag>1]=2
    Cflag   = np.nansum(np.vstack((sttest2[:,4], itcw, c02exli, co2excr, wexvar,
                                   h2oexvar10, co2exvar3, wvar0, co2var0, spikef[:,3],
                                   ustarf)).transpose(), 1).astype(int) ### why h20exvar10
    Cflag[Cflag>1]=2
    Tauflag = np.nansum(np.vstack((sttest2[:,3], itcu, itcw, tauexcr, wexvar,
                                   wvar0, spikef[:,4])).transpose(), 1).astype(int)
    Tauflag[Tauflag>1]=2
    header2     = np.array(['          H', '   Hflag', '         LE', '  LEflag', '          E', '   Eflag', '          C',
                            '   Cflag', '        TAU', ' TAUflag'])
    fluxes[np.isnan(fluxes)] = novalue
    fluxes_str  = np.array([['%11.5f'%x for x in y] for y in fluxes])
    fluxes_str  = np.where(fluxes_str=='%11.5f'%novalue, ' '*(11-len(str(novalue)))+str(novalue), fluxes_str)
    Hflag_str   = astr(Hflag, prec=8)
    LEflag_str  = astr(LEflag, prec=8)
    Eflag_str   = astr(Eflag, prec=8)
    Cflag_str   = astr(Cflag, prec=8)
    Tauflag_str = astr(Tauflag, prec=8)
    
    output = np.hstack((d[:,0:1], fluxes_str.repeat(2, axis=1)))
    output[:,2::2] = np.vstack((Hflag_str, LEflag_str, Eflag_str, Cflag_str,
                               Tauflag_str)).transpose()
    np.savetxt('%s/fluxflags.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header2))[np.newaxis,:],
                          output)), '%s', delimiter=',')
    
    ###########################################################################
    # write log file with stats
    log = open('%s/flags.log'%outdir, 'w')
    log.write('flag, 0, 1, 2\n')
    for i, item in enumerate(header):
        hist = tuple(np.histogram(flags[:,i], [-0.5,0.5,1.5,2.5])[0])
        log.write(item.strip()+', %i, %i, %i\n'%hist)    
    log.close()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
