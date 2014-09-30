import numpy as np
from gapfill import gapfill
from esat import esat
from date2dec import date2dec
from autostring import astr

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
    This file is part of the UFZ Python library.

    The UFZ Python library is free software: you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    The UFZ Python library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with The UFZ Python library.  If not,
    see <http://www.gnu.org/licenses/>.

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Aug 2014
    '''       

    ###########################################################################
    # reading input files
    d = np.loadtxt(fluxfile, dtype='|S100', delimiter=delimiter[0], skiprows=skiprows[0])
    m = np.loadtxt(metfile,  dtype='|S100', delimiter=delimiter[1], skiprows=skiprows[1])
    
    if d.shape[1]!=11:
        raise ValueError('fluxfill: fluxfile must be from fluxflag and have 11 cols')
    
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
    data      = val[:,0::2] # corr. H, corr. LE, corr. E, corr. C, tau
    data_flag = val[:,1::2]>1
    rg        = met[:,rg]
    rg_flag   = rg==undef
    tair      = met[:,tair]
    tair_flag = tair==undef
    rh        = met[:,rh]
    vpd       = (1.-rh/100.)*esat(tair+273.14)    
    vpd_flag  = vpd==undef
    
    flux = np.zeros_like(data)
    flag = np.zeros_like(data_flag, dtype=np.int)
    
    ###########################################################################
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pp1 = pdf.PdfPages(outdir+'/fluxfilled.pdf')
        
    ###########################################################################
    # do gapfill
    for i in range(data.shape[1]):
        flux[:,i], flag[:,i] = gapfill(datev, data[:,i], rg, tair, vpd,
                                       data_flag[:,i], rg_flag, tair_flag,
                                       vpd_flag, rg_dev=50., tair_dev=2.5,
                                       vpd_dev=5., longgap=60, fullday=False,
                                       undef=undef, ddof=1, err=False,
                                       shape=False)
    
        #######################################################################
        # plot
        if plot:
            fig1 = plt.figure(1)
            sub1 = fig1.add_subplot(111)
            l1 =sub1.plot(datev, flux[:,i], '-g')
            l2 =sub1.plot(datev, np.ma.array(data[:,i],mask=data_flag[:,i]), '-b')
            plt.show()
            fig1.savefig(pp1, format='pdf')
    
    if plot:
        pp1.close()
    
    ###########################################################################
    # prepare output and save file
    header         = np.array(['          H', '   Hflag', '     Hgf',
                               '         LE', '  LEflag', '    LEfg',
                               '          E', '   Eflag', '     Efg',
                               '          C', '   Cflag', '     Cfg',
                               '        TAU', ' TAUflag', '   TAUfg'])
    flux_str       = np.array([['%11.5f'%x for x in y] for y in flux])
    output         = np.hstack((d[:,0:1], flux_str.repeat(3, axis=1)))
    output[:,2::3] = astr(val[:,1::2].astype(int), prec=8)
    output[:,3::3] = astr(flag, prec=8)    
    
    np.savetxt('%s/fluxfilled.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header))[np.newaxis,:],
                          output)), '%s', delimiter=',')

if __name__ == '__main__':
    import doctest
    doctest.testmod()