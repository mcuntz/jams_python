import numpy as np
from nee2gpp import nee2gpp
from esat import esat
from date2dec import date2dec

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
    This file is part of the UFZ Python library.
    
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
    
    if d.shape[1]!=16:
        raise ValueError('fluxpart: fluxfile must be from fluxfill and have 16 cols')
    
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
    nee       = val[:,9] #corr. C
    swdr      = met[:,swdr]
    isday     = swdr>0.
    tair      = met[:,tair]+273.15
    rh        = met[:,rh]
    vpd       = (1.-rh/100.)*esat(tair)    
    
    ###########################################################################
    # do partitioning
    gpp, reco = nee2gpp(datev, nee, tair, isday, rg=swdr, vpd=vpd, undef=undef,
                        method=method, shape=False, masked=False,
                        nogppnight=nogppnight)
    
    #######################################################################
    # plot
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pp1 = pdf.PdfPages(outdir+'/fluxpart.pdf')

        fig1 = plt.figure(1)
        sub1 = fig1.add_subplot(111)
        l1 =sub1.plot(datev, nee, '-k', label='nee')
        l2 =sub1.plot(datev, gpp, '-g', label='gpp')
        l3 =sub1.plot(datev, reco, '-r', label='reco')
        plt.legend(loc='best')
        plt.show()
        fig1.savefig(pp1, format='pdf')
        pp1.close()
    
    ###########################################################################
    # prepare output and save file
    header         = np.array(['          H', '   Hflag', '     Hgf',
                               '         LE', '  LEflag', '    LEfg',
                               '          E', '   Eflag', '     Efg',
                               '          C', '   Cflag', '     Cfg',
                               '        TAU', ' TAUflag', '   TAUfg',
                               '        GPP', ' GPPflag', '   GPPfg',
                               '       Reco', 'Recoflag', '  Recofg'])
    flux = np.vstack((gpp,reco)).transpose()
    flux_str       = np.array([['%11.5f'%x for x in y] for y in flux]).repeat(3, axis=1)
    flux_str[:,[1,2,4,5]] = np.tile(d[:,11:13],2)
    output         = np.hstack((d[:,:], flux_str))
    np.savetxt('%s/fluxpart.csv'%outdir,
               np.vstack((np.concatenate((['            date'], header))[np.newaxis,:],
                          output)), '%s', delimiter=',')

if __name__ == '__main__':
    import doctest
    doctest.testmod()