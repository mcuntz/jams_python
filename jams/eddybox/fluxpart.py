#!/usr/bin/env python
import numpy as np
from jams.eddybox import nee2gpp
from jams.esat import esat
from jams.date2dec import date2dec

def fluxpart(fluxfile, metfile, outdir, swdr, tair, rh, method='local',
             nogppnight=False, delimiter=[',',','], skiprows=[1,1],
             format=['ascii','ascii'], undef=-9999, plot=False):
    '''
    Wrapper function for nee2gpp with file management and plotting.
    
    
    Definition
    ----------
    fluxpart(fluxfile, metfile, outdir, rg, tair, rh, method='local',
             delimiter=[',',','], skiprows=[1,1], format=['ascii','ascii'],
             undef=-9999, plot=False):
    
    
    Input
    ----- 
    fluxfile    str, path and file name of fluxflag or fluxfill output file
                containing fluxes and flags
    metfile     str, path and file name of the meteorology file (must be in
                sync with fluxfile)
    outdir      str, path of the output folder
    swdr        int, column number of short wave downward radiation [W m-2] in
                metfile, column number starts with 0 which is first data column.
                swdr is used for lasslopp and for swdr>0=isday
    tair        int, column number of air temperature [deg C] in metfile, column
                number starts with 0 which is first data column.
    rh          int, column number of relative humidity [%] in metfile, column
                number starts with 0 which is first data column.
                

    Optional Input
    --------------
    method      str, if 'global', fit of Reco vs. temperature to all nighttime data
                     if 'local' | 'reichstein',  method of Reichstein et al. (2005)
                     if 'day'   | 'lasslop',     method of Lasslop et al. (2010)
                     (default: 'local')
    nogppnight  if True:  Resp=NEE, GPP=0 at night, GPP always positive
                if False: Resp=lloyd_taylor, GPP=Resp-NEE at night (default)
    delimiter   list of str, delimiters of fluxfile and metfile
                (default: [',',','])
    skiprows    list of int, lines to skip at the beginning of fluxfile and
                metfile, e.g. header lines (default: [1,1])
    format      list of str, time formats of fluxfile and metfile, 'ascii' and 
                'eng' possible (default: ['ascii','ascii'])
    undef       int/float, missing value of fluxfile and metfile
                (default: -9999, np.nan is not possible)
    plot        bool, if True performs plotting (default: False)
    
    
    Output
    ------
    fluxpart.csv file containing fluxes with original flags and quality flags
                 of the gap filling plus gpp and reco fluxes. gpp and reco get
                 flags of c flux. 
    fluxpart.pdf plot of c flux, gpp and reco
        
    
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
    
    assert (d.shape[1]==16) | (d.shape[1]==24), 'fluxfill: fluxfile must be from fluxfill and have 16 or 19 cols'
    
    if format[0]=='ascii':
        datev   = date2dec(ascii=d[:,0])
    elif format[0]=='eng':
        datev   = date2dec(eng=d[:,0])
    else:
        raise ValueError('fluxpart: unknown format')    
    if format[1]=='ascii':
        datem   = date2dec(ascii=m[:,0])
    elif format[1]=='eng':
        datem   = date2dec(eng=m[:,0])
    else:
        raise ValueError('fluxpart: unknown format')
    
    val = np.where(d[:,1:]=='', str(undef), d[:,1:]).astype(np.float)
    met = np.where(m[:,1:]=='', str(undef), m[:,1:]).astype(np.float)
    
    ###########################################################################
    # assign variables
    if (d.shape[1]==11):
        nee       = val[:,9] #corr. C
    else:
        nee       = val[:,12] # oder 13 #corr. C
        nee_stor  = val[:,13] # oder 13 #corr. C
    swdr      = met[:,swdr]
    isday     = swdr>0.
    tair      = met[:,tair]+273.15
    rh        = met[:,rh]
    vpd       = (1.-rh/100.)*esat(tair)    
    
    ###########################################################################
    # do partitioning
    gpp, reco           = nee2gpp(datev, nee,      tair, isday, rg=swdr,
                                  vpd=vpd, undef=undef, method=method,
                                  shape=False, masked=False, nogppnight=nogppnight)
    gpp_stor, reco_stor = nee2gpp(datev, nee_stor, tair, isday, rg=swdr,
                                  vpd=vpd, undef=undef, method=method,
                                  shape=False, masked=False, nogppnight=nogppnight)
    
    #######################################################################
    # plot
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pp1 = pdf.PdfPages(outdir+'/fluxpart.pdf')
        
        majticks = mpl.dates.MonthLocator(bymonthday=1)
        format_str='%d %m %Y %H:%M'
        date01 = date2dec(yr=1, mo=1, dy=2, hr=0, mi=0, sc=0)
        
        fig1 = plt.figure(1)
        if (d.shape[1]==11):
            sub1 = fig1.add_subplot(111)
        else:
            sub1 = fig1.add_subplot(211)
            sub2 = fig1.add_subplot(212)
        l1 =sub1.plot(datev-date01, nee, '-k', label='nee')
        l2 =sub1.plot(datev-date01, gpp, '-g', label='gpp')
        l3 =sub1.plot(datev-date01, reco, '-r', label='reco')
        if (d.shape[1]!=11):
            l4 =sub2.plot(datev-date01, nee_stor, '-k', label='nee')
            l5 =sub2.plot(datev-date01, gpp_stor, '-g', label='gpp')
            l6 =sub2.plot(datev-date01, reco_stor, '-r', label='reco')
        
        sub1.set_xlim(datev[0]-date01,datev[-1]-date01)
        sub1.xaxis.set_major_locator(majticks)
        sub1.xaxis.set_major_formatter(mpl.dates.DateFormatter(format_str))
        if (d.shape[1]!=11):
            sub2.set_xlim(datev[0]-date01,datev[-1]-date01)
            sub2.xaxis.set_major_locator(majticks)
            sub2.xaxis.set_major_formatter(mpl.dates.DateFormatter(format_str))
        fig1.autofmt_xdate()
         
        plt.legend(loc='best')
        plt.show()
        fig1.savefig(pp1, format='pdf')
        pp1.close()
    
    ###########################################################################
    # prepare output and save file
    if (d.shape[1]==11):
        header         = np.array(['          H', '   Hflag', '     Hgf',
                                   '         LE', '  LEflag', '    LEgf',
                                   '          E', '   Eflag', '     Egf',
                                   '          C', '   Cflag', '     Cgf',
                                   '        TAU', ' TAUflag', '   TAUgf',
                                   '        GPP', ' GPPflag', '   GPPgf',
                                   '       Reco', 'Recoflag', '  Recogf'])
        flux     = np.vstack((gpp,reco)).transpose()
        flux_str = np.array([['%11.5f'%x for x in y] for y in flux]).repeat(3, axis=1)
        flux_str[:,[1,2,4,5]] = np.tile(d[:,11:13],2)
        output   = np.hstack((d[:,:], flux_str))
    else:
        header         = np.array(['          H', '       H+sT', '   Hflag', '     Hgf',
                                   '         LE', '     LE+sLE', '  LEflag', '    LEgf',
                                   '          E', '       E+sE', '   Eflag', '     Egf',
                                   '          C', '       C+sC', '   Cflag', '     Cgf',
                                   '        TAU',    ' TAUflag',             '   TAUgf',
                                   '         sT', '        sLE', '         sE', '         sC',
                                   '        GPP', '     GPP+sC', ' GPPflag', '   GPPfg',
                                   '       Reco', '    Reco+sC', 'Recoflag', '  Recofg'])
        flux     = np.vstack((gpp,gpp_stor,reco,reco_stor)).transpose()
        flux_str = np.array([['%11.5f'%x for x in y] for y in flux])
        flux_str = np.insert(flux_str, [2,2,4,4], flux_str, axis=1)
        flux_str[:,[2,3,6,7]] = np.tile(d[:,11:13],2)
        output   = np.hstack((d[:,:], flux_str))
    
    np.savetxt('%s/fluxpart.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header))[np.newaxis,:],
                          output)), '%s', delimiter=',')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
