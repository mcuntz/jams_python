#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
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
