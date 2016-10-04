#!/usr/bin/env python
import numpy as np
from jams.date2dec import date2dec
from jams.const import cheat_quartz, cheat_water, cheat_air, density_quartz
from scipy.interpolate import splrep, splint
import scipy.optimize as opt

def energyclosure(fluxfile, metfile, outdir, Rn, G, swdr, Ts=None, theta=None,
                  depths=None, por=None, method='year', force0=True,
                  delimiter=[',',','], skiprows=[1,1], format=['ascii','ascii'],
                  undef=-9999, plot=False):
    '''
    Calculation of energy balance closure and correction of Eddy Covariance data
    with different methods. Possible calculation of soil heat flux from soil
    temperature and moisture if not given.
    
    
    Definition
    ----------
    energyclosure(fluxfile, metfile, outdir, Rn, G, swdr, Ts=None, theta=None,
                  depths=None, por=None, method='year', force0=True,
                  delimiter=[',',','], skiprows=[1,1], format=['ascii','ascii'],
                  undef=-9999, plot=False):
    
    
    Input
    ----- 
    fluxfile    str, path and file name of fluxflag or fluxfill or fluxpart
                output file containing fluxes and flags
    metfile     str, path and file name of the meteorology file (must be in
                sync with fluxfile)
    outdir      str, path of the output folder
    Rn          int, column number of net radiation [W m-2] in metfile, column
                number starts with 0 which is first data column.
    G           int, column number of soil heat flux [W m-2] in metfile, column
                number starts with 0 which is first data column. If soil heat
                flux not available, set to False
    swdr        int, column number of short wave downward radiation [W m-2] in
                metfile, column number starts with 0 which is first data column.
                swdr is used for swdr>0=isday
                

    Optional Input
    --------------
    Ts          list(M) of int, if G is not given, column numbers of soil
                temperatures in certain depths in metfile, column number starts
                with 0 which is first data column.
    theta       list(M) of int, if G is not given, column numbers of soil
                moistures in certain depths in metfile, column number starts
                with 0 which is first data column.
    depths      list(M) or array(M) of int, if G is not given, depths of Ts and
                theta measurements [m], positively increasing with depths
                e.g. [0.05, 0.30, 0.60]
    por         float, if G is not given, porosity of the soil [-], must be
                bigger or equal than np.amax(theta)     
    method      str, method of how energy balance closure is calculated and
                applied.
                    if 'year', fit of whole year daytime flag=0 data
                    Rn-G vs. H+Le and application of the fraction to all H, Le
                    and E values
                    implement new methods here
                    implement new methods here
                    implement new methods here
                    (default: 'year')
    force0      bool, if method='year', True forces fit through origin, if False
                fit is allowed to have intercept (default: True)
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
    fluxclosed.csv file containing fluxes and flags where depending on the
                        method chosen fluxes are corrected for energy balance
                        closure
    fluxclosed.pdf if method='year', plot of whole year daytime flag=0 data
                   Rn-G vs. H+Le
        
    fluxclosed_stor.pdf if method='year', plot of whole year daytime flag=0 data
                   Rn-G vs. H+Le including storages
        
    
    Restrictions
    ------------
    Works only with half hourly time steps


    License
    -------
    This file is part of the JAMS Python package.
    
    It is NOT released under the GNU Lesser General Public License, yet.
    
    If you use this routine, please contact Arndt Piayda.
    
    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Sep 2014
    '''       
    ###########################################################################
    # reading input files
    d = np.loadtxt(fluxfile, dtype='|S100', delimiter=delimiter[0])
    m = np.loadtxt(metfile,  dtype='|S100', delimiter=delimiter[1])
    
    assert (d.shape[1]==22) | (d.shape[1]==32), 'energyclosure: fluxfile must be from fluxpart and have 22 or 32 cols'
        
    if format[0]=='ascii':
        datev   = date2dec(ascii=d[skiprows[0]:,0])
    elif format[0]=='eng':
        datev   = date2dec(eng=d[skiprows[0]:,0])
    else:
        raise ValueError('energyclosure: unknown format')    
    if format[1]=='ascii':
        datem   = date2dec(ascii=m[skiprows[1]:,0])
    elif format[1]=='eng':
        datem   = date2dec(eng=m[skiprows[1]:,0])
    else:
        raise ValueError('energyclosure: unknown format')
    
    val = np.where(d[skiprows[0]:,1:]=='', str(undef), d[skiprows[0]:,1:]).astype(np.float)
    met = np.where(m[skiprows[1]:,1:]=='', str(undef), m[skiprows[1]:,1:]).astype(np.float)
    
    ###########################################################################
    # assign variables
    if (d.shape[1]==22):
        H         = val[:,0]
        Hflag     = val[:,1]
        Le        = val[:,3]
        Leflag    = val[:,4]
        E         = val[:,6]
        Eflag     = val[:,7]
    else:
        H         = val[:,0]
        Hflag     = val[:,2]
        Le        = val[:,4]
        Leflag    = val[:,6]
        E         = val[:,8]
        Eflag     = val[:,10]    
        H_stor      = val[:,1]
        Hflag_stor  = val[:,2]
        Le_stor     = val[:,5]
        Leflag_stor = val[:,6]
        E_stor      = val[:,9]
        Eflag_stor  = val[:,10]    
    Rn        = met[:,Rn]
    Ts        = met[:,Ts]
    swdr      = met[:,swdr]
    isday     = swdr>0.

    ###########################################################################
    # check if soil heat flux is given or needs to be calculated
    if not G:
        if Ts==None or theta==None or depths==None or por==None:
            raise ValueError('energyclosure: if G is not given, Ts, theta, depths and rhos are needed to calculate G')
        else:
            G=soilheatflux(Ts, theta, depths, rhos)
    else:
        G = met[:,G]
        
    ###########################################################################
    # apply energy balance closure methods
    # yearly approach
    if (method.lower() == 'year'):
        closure = energyclosure_year(H, Hflag, Le, Leflag, Rn, G, isday, outdir,
                                     force0=force0, undef=undef, plot='energyclosure.pdf')
        if force0:
            H_closed, Le_closed, E_closed = H/closure, Le/closure, E/closure
        else:
            H_closed, Le_closed, E_closed = H/closure[0], Le/closure[0], E/closure[0]

        if (d.shape[1]==32):
            closure_stor = energyclosure_year(H_stor, Hflag_stor, Le_stor, Leflag_stor,#
                                              Rn, G, isday, outdir, force0=force0,
                                              undef=undef, plot='energyclosure_stor.pdf')
            if force0:
                H_stor_closed, Le_stor_closed, E_stor_closed = H_stor/closure, Le_stor/closure, E_stor/closure
            else:
                H_stor_closed, Le_stor_closed, E_stor_closed = H_stor/closure[0], Le_stor/closure[0], E_stor/closure[0]
            
    # Include new methods here
    # Include new methods here
    # Include new methods here
    else:
        raise ValueError('energyclosure: method not implemented yet.')

    H_closed[H==undef]   = undef
    Le_closed[Le==undef] = undef
    E_closed[E==undef]   = undef
    if (d.shape[1]==32):
        H_stor_closed[H_stor==undef]   = undef
        Le_stor_closed[Le_stor==undef] = undef
        E_stor_closed[E_stor==undef]   = undef
    
    ###########################################################################
    # prepare and write output
    if (d.shape[1]==22):    
        flux = np.vstack((H_closed,Le_closed,E_closed)).transpose()
    else:
        flux = np.vstack((H_closed,H_stor_closed,Le_closed,Le_stor_closed,E_closed,E_stor_closed)).transpose()        
    
    flux_str = np.array([['%11.5f'%x for x in y] for y in flux])
    flux_str = np.where(flux_str=='%11.5f'%undef, ' '*(11-len(str(undef)))+str(undef), flux_str)
    
    if (d.shape[1]==22):    
        d[skiprows[0]:,[1,4,7]] = flux_str
    else:
        d[skiprows[0]:,[1,2,5,6,9,10]] = flux_str
    
    np.savetxt('%s/fluxclosed.csv'%outdir, d, '%s', delimiter=',')

def soilheatflux(Ts, theta, depths, por, undef=-9999):
    '''
    Calculates soil heat flux from soil temperatures and soil moistures at
    different depths. It's best do include surface temperatures at 0 m.
    
    
    Definition
    ----------
    soilheatflux(Ts, theta, depths, por, undef=-9999):
    
    
    Input
    ----- 
    Ts          array(N,M), soil temperatures [K or degC]
    theta       array(N,M), soil moisture [-]
    depths      list(M) or array(M), depths of the Ts and theta measurements [m]
                positively increasing with depths e.g. [0.05, 0.30, 0.60]
    por         float, porosity of the soil [-] must be bigger or equal
                than np.amax(theta)
               

    Optional Input
    --------------
    undef       int/float, missing value of Ts and theta
                (default: -9999, np.nan is not possible)
    
    
    Output
    ------
    G           array(N), soil heat flux [W m-2]
    
    
    Restrictions
    ------------
    Works only with half hourly time steps
    
    
    License
    -------
    This file is part of the JAMS Python package.

    The JAMS Python package is free software: you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    The JAMS Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with The JAMS Python package.  If not,
    see <http://www.gnu.org/licenses/>.

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Sep 2014
    '''       
    if (np.size(depths)!=Ts.shape[1]) or (Ts.shape!=theta.shape):
        raise ValueError('soilheatflux: dimensions mismatch')
    
    Ts    = np.ma.array(Ts, mask=Ts==undef)
    theta = np.ma.array(theta, mask=theta==undef)
    
    # calc soil density and soil heat capacity
    rhos  = (1.-por)*density_quartz
    csoil = cheat_quartz*(1.-por) + cheat_water*theta + cheat_air*(1-theta)*por
    
    # calculate heat flux    
    T_diff = (Ts[:-1,:] - Ts[1:,:])/1800. * rhos * csoil[:-1,:]
    G      = np.ma.masked_all_like(Ts[:,0])
    for i, item in enumerate(T_diff):
        if not item.mask.any():
            tck  = splrep(depths, item, k=1)
            G[i+1] = splint(depths[0],depths[-1],tck)
    G[0]=G[1]
    
    return np.ma.filled(G, fill_value=undef)

def energyclosure_year(H, Hflag, Le, Leflag, Rn, G, isday, outdir, force0=True,
                       undef=-9999, plot=False):
    '''
    Calculates energy balance closure with the 'year' approach. A straight line
    is fitted to Rn-g vs. H+Le and the fraction is returned.
    
    
    Definition
    ----------
    energyclosure_year(H, Le, Rn, G, isday, outdir, force0=force0, undef=undef,
                       plot=plot):
    
    
    Input
    ----- 
    H           array(N), sensible heat flux [W m-2]
    Hflag       array(N), quality flag of sensible heat flux, 0 where good quality
    Le          array(N), latent heat flux [W m-2]
    Leflag      array(N), quality flag of latent heat flux, 0 where good quality
    Rn          array(N), net radiation [W m-2]
    G           array(N), soil heat flux [W m-2]
    isday       array(N), bool, True if day, False if night, fitting is done
                only to day data
    outdir      str, path of the output folder
    

    Optional Input
    --------------
    force0      bool, if True fitting line is forced to origin and output is
                only the fraction, if False fitted line can have a intercept and
                output is then (fraction, intercept), (default: True)
    undef       int/float, missing value of input (default: -9999, np.nan is
                not possible)
    plot        bool, if False, no fitting plot is made and saved to outdir
                (default: False), otherwise give file name like 'energyclosure.pdf'
    
    
    Output
    ------
    fraction    float, energy balance closure gap, e.g 0.8 means 80% of
                available energy is observed, 20% are missing. list(2) of floats
                with energy balance closure gap and intercept if force0=False
    
    
    License
    -------
    This file is part of the JAMS Python package.

    The JAMS Python package is free software: you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    The JAMS Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with The JAMS Python package.  If not,
    see <http://www.gnu.org/licenses/>.

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Sep 2014
    '''       
    ###########################################################################
    # variables
    H  = np.ma.array(H,  mask=(H==undef)  | (~isday) |(Hflag>0))
    Le = np.ma.array(Le, mask=(Le==undef) | (~isday) |(Leflag>0))
    Rn = np.ma.array(Rn, mask=(Rn==undef) | (~isday))
    G  = np.ma.array(G,  mask=(G==undef)  | (~isday))
    
    x = Rn-G
    y = H+Le
    
    ###########################################################################
    # fit
    if force0:
        p_guess = [1]
        func    = lin2
    else:
        p_guess = [0,1]
        func    = lin

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        x_mod = np.arange(np.ma.min(x),np.ma.max(x),(np.ma.max(x)-np.ma.min(x))/50.)
        fig1=plt.figure('pre-fit')
        fig1.add_subplot(111, aspect='equal')
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')                   
        plt.plot(x, y, 'bo', label='n: %i'%np.ma.count(x+y))         
        plt.plot(x_mod, func(x_mod,p_guess), 'r-')
        plt.xlabel('Rn-G')
        plt.ylabel('H+Le')
        plt.legend(loc='best')
        plt.show()
    
    p_opt = opt.fmin(absdif, p_guess, args=(x, y, func), disp=0)
    
    if plot:
        fig1=plt.figure('post-fit')                                        
        fig1.add_subplot(111, aspect='equal')
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')                   
        plt.plot(x, y, 'bo', label='n: %i'%np.ma.count(x+y))         
        if force0:
            plt.plot(x_mod, func(x_mod,p_opt), 'r-', label='slope: %0.2f'%p_opt)
        else:
            plt.plot(x_mod, func(x_mod,p_opt), 'r-', label='slope: %0.2f'%p_opt[0])
        plt.xlabel('Rn-G')
        plt.ylabel('H+Le')
        plt.legend(loc='best')
        plt.show()
        
        pp1 = pdf.PdfPages('%s/'%(outdir)+plot)
        fig1.savefig(pp1, format='pdf')
        pp1.close()

    return p_opt

def absdif(p,x,y,func):
    '''
    energyclosure: objective function
    '''
    return np.ma.sum(np.ma.abs(y-func(x, p)))

def lin(x, p):
    '''
    energyclosure: linear function
    '''    
    return p[0] * x + p[1]

def lin2(x, p):
    '''
    energyclosure: linear function with no intercept
    '''
    return p[0] * x

if __name__ == '__main__':
    import doctest
    doctest.testmod()
