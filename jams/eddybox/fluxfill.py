#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.eddybox import gapfill
from jams.esat import esat
from jams.date2dec import date2dec
from jams.autostring import astr

def fluxfill(fluxfile, metfile, outdir, rg, tair, rh, delimiter=[',',','],
             skiprows=[1,1], format=['ascii','ascii'], undef=-9999, plot=False):
    '''
    Wrapper function for gapfill with file management and plotting.
    
    
    Definition
    ----------
    fluxfill(fluxfile, metfile, outdir, rg, tair, rh, delimiter=[',',','],
             skiprows=[1,1], format=['ascii','ascii'], undef=-9999, plot=False):
    
    
    Input
    ----- 
    fluxfile    str, path and file name of fluxflag output file containing
                fluxes and flags
    metfile     str, path and file name of the meteorology file (must be in
                sync with fluxfile)
    outdir      str, path of the output folder
    rg          int, column number of global radiation [W m-2] in metfile,
                column number starts with 0 which is first data column.
    tair        int, column number of air temperature [deg C] in metfile, column
                number starts with 0 which is first data column.
    rh          int, column number of relative humidity [%] in metfile, column
                number starts with 0 which is first data column.
                
                        
    Optional Input
    --------------
    delimiter   list of str, delimiters of fluxfile and metfile
                (default: [',',','])
    skiprows    list of int, lines to skip at the beginning of fluxfile and
                metfile, e.g. header lines (default: [1,1])
    format      list of str, time formats of fluxfile and metfile, 'ascii' and 
                'eng' possible (default: ['ascii','ascii'])
    undef       int/float, missing value of fluxfile and metfile
                (default: -9999, np.nan not possible)
    plot        bool, if True performs plotting (default: False)
    
    
    Output
    ------
    fluxfilled.csv file containing fluxes with original flags and quality flags
                   of the gap filling. Fluxes are filled where flag>1. 
    fluxfilled.pdf  plot of each gap filled flux.
        
    
    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2014 Arndt Piayda

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
    Written,  AP, Aug 2014
    '''       

    ###########################################################################
    # reading input files
    d = np.loadtxt(fluxfile, dtype='|S100', delimiter=delimiter[0], skiprows=skiprows[0])
    m = np.loadtxt(metfile,  dtype='|S100', delimiter=delimiter[1], skiprows=skiprows[1])

    assert (d.shape[1]==11) | (d.shape[1]==19), 'fluxfill: fluxfile must be from fluxflag or profile2storage and have 11 or 19 cols'
    
    if format[0]=='ascii':
        datev   = date2dec(ascii=d[:,0])
    elif format[0]=='eng':
        datev   = date2dec(eng=d[:,0])
    else:
        raise ValueError('fluxfill: unknown format')    
    if format[1]=='ascii':
        datem   = date2dec(ascii=m[:,0])
    elif format[1]=='eng':
        datem   = date2dec(eng=m[:,0])
    else:
        raise ValueError('fluxfill: unknown format')
    
    val = np.where(d[:,1:]=='', str(undef), d[:,1:]).astype(np.float)
    met = np.where(m[:,1:]=='', str(undef), m[:,1:]).astype(np.float)
    
    ###########################################################################
    # assign variables
    if (d.shape[1]==11):
        data      = val[:,0::2] # corr. H, corr. LE, corr. E, corr. C, tau
        data_flag = (val[:,1::2]>1) | (val[:,0::2]==undef)
    else:
        data      = val[:,[0,1, 3,4, 6,7, 9,10, 12]] # corr. H, corr. H+s, corr. LE, corr. LE+s, corr. E, corr. E+s, corr. C, corr. C+s, tau
        data_flag = np.repeat((val[:,[2,5,8,11,13]]>1) | (val[:,[0,1, 3,4, 6,7, 9,10, 12]]==undef), 2, axis=1)[:,:-1]
    rg        = met[:,rg]
    rg_flag   = rg==undef
    tair      = met[:,tair]
    tair_flag = tair==undef
    rh        = met[:,rh]
    rh_flag   = rh==undef
    vpd       = np.empty_like(tair)
    vpd[tair_flag | rh_flag]    = undef
    vpd[~(tair_flag | rh_flag)] = (1.-rh[~(tair_flag | rh_flag)]/100.)*esat(tair[~(tair_flag | rh_flag)]+273.14)
    vpd_flag  = vpd==undef
    
    flux = np.zeros_like(data)
    flag = np.zeros_like(data_flag, dtype=np.int)
    
    ###########################################################################
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pp1 = pdf.PdfPages(outdir+'/fluxfilled.pdf')
        
    ###########################################################################
    # do gapfill
    for i in range(data.shape[1]):
        flux[:,i], flag[:,i] = gapfill.gapfill(datev, data[:,i], rg, tair, vpd,
                                       data_flag[:,i], rg_flag, tair_flag,
                                       vpd_flag, rg_dev=50., tair_dev=2.5,
                                       vpd_dev=5., longgap=60, fullday=False,
                                       undef=undef, ddof=1, err=False,
                                       shape=False)
    
        #######################################################################
        # plot
        if plot:
            majticks = mpl.dates.MonthLocator(bymonthday=1)
            format_str='%d %m %Y %H:%M'
            date01 = date2dec(yr=1, mo=1, dy=2, hr=0, mi=0, sc=0)
            
            fig1 = plt.figure(1)
            sub1 = fig1.add_subplot(111)
            l1 =sub1.plot(datev-date01, flux[:,i], '-g')
            l2 =sub1.plot(datev-date01, np.ma.array(data[:,i],mask=data_flag[:,i]), '-b')
            
            sub1.set_xlim(datev[0]-date01,datev[-1]-date01)
            sub1.xaxis.set_major_locator(majticks)
            sub1.xaxis.set_major_formatter(mpl.dates.DateFormatter(format_str))
            fig1.autofmt_xdate()
            
            plt.show()
            fig1.savefig(pp1, format='pdf')
    
    if plot:
        pp1.close()
    
    ###########################################################################
    # prepare output and save file
    flux_str       = np.array([['%11.5f'%x for x in y] for y in flux])

    if (d.shape[1]==11):
        header         = np.array(['          H', '   Hflag', '     Hgf',
                                   '         LE', '  LEflag', '    LEgf',
                                   '          E', '   Eflag', '     Egf',
                                   '          C', '   Cflag', '     Cgf',
                                   '        TAU', ' TAUflag', '   TAUgf'])
        output         = np.hstack((d[:,0:1], flux_str.repeat(3, axis=1)))
        output[:,2::3] = astr(val[:,1::2].astype(int), prec=8)
        output[:,3::3] = astr(flag, prec=8)
    else:
        header         = np.array(['          H', '       H+sT', '   Hflag', '     Hgf',
                                   '         LE', '     LE+sLE', '  LEflag', '    LEgf',
                                   '          E', '       E+sE', '   Eflag', '     Egf',
                                   '          C', '       C+sC', '   Cflag', '     Cgf',
                                   '        TAU',    ' TAUflag',             '   TAUgf',
                                   '         sT', '        sLE', '         sE', '         sC'])    
        output                   = np.hstack((d[:,0:1], np.insert(flux_str, [2,2,4,4,6,6,8,8], flux_str[:,:-1], axis=1),
                                              flux_str[:,-1:].repeat(2, axis=1), d[:,-4:]))
        output[:,[3,7,11,15,18]] = astr(val[:,[2,5,8,11,13]].astype(int), prec=8)
        output[:,[4,8,12,16,19]] = astr(flag[:,[0,2,4,6,8]], prec=8)
    
    np.savetxt('%s/fluxfilled.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header))[np.newaxis,:],
                          output)), '%s', delimiter=',')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
