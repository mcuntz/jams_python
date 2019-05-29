#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.sread import sread
from jams.fread import fread
from scipy import interpolate
import csv
import sys
import os
import time as t
import math as m
import shutil as sh
from scipy.optimize import fmin

def eddycorr(indir, sltdir, cfile, hfile, meteofile, outfile, novalue=-9999,
             histstep=10, attach=True, plot=False):
    '''
    Moves EddyCorr files (cfile=35_corr.csv and hfile=36_corr.csv) from sltdir
    to indir after they have been created by EddyCorr (Kolle & Rebmann, 2007).
    Time lags between wind and concentrations are plotted and user selects
    thresholds. Water lag is correlated against rH from meteofile for missing
    values. Plots and lags are saved to outfile in indir and attached to the
    meteofile for being read by EddyFlux.


    Definition
    ----------
    eddycorr(indir, sltdir, cfile, hfile, meteofile, outfile, novalue=-9999,
             histstep=10, attach=True):


    Input
    -----
    indir       str, path of the folder where results will be saved
    sltdir      str, path of the folder containing the *.slt files and EddyCorr
                files
    cfile       str, name of the carbon lag file, e.g. 35_corr.csv (EddyCorr)
    hfile       str, name of the water lag file, e.g. 36_corr.csv (EddyCorr)
    meteofile   str, path of the meteorological file for EddyFlux (in sync with
                available *.slt files, e.g. made by meteo4slt) where the lags
                will be attached
    outfile     str, name of the output file


    Optional Input
    --------------
    novalue     int, novalue in meteofile (default=-9999)
    histstep    int, histogram steps for plotting (default=10)
    attach      bool, if True, lags will be attached to meteofile.


    Output
    ------
    lags_X.csv     file containing the lags
    c_lags_X.pdf   plot of the final carbon lags
    h_lags_X.pdf   plot of the final water lags
    laglog_X_X.log log file


    License
    -------
    This file is part of the JAMS Python package.

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

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Jul 2014
    Modified, AP, Aug 2014 - major bug fix
    '''

    ############################################################################
    # move correlation files from sltdir to indir
    sh.copy(os.path.join(sltdir, cfile), indir)
    sh.copy(os.path.join(sltdir, hfile), indir)

    ############################################################################
    # reading input file
    header = sread(os.path.join(indir,cfile))[0]
    doys   = np.array(sread(os.path.join(indir,cfile), nc=1, skip=1), dtype='|S16')
    #doys   = np.array([x[5:12] for x in doys.flatten()], dtype = '|S7')
    day    = np.array([x[5:8] for x in doys.flatten()], dtype = '|S7').astype(float)
    hour   = np.array([x[8:10] for x in doys.flatten()], dtype = '|S2').astype(float)
    min    = np.array([x[10:12] for x in doys.flatten()], dtype = '|S2').astype(float)
    doysfloat   = day + (hour + min/60.)/24.
    c      = np.array(fread(os.path.join(indir,cfile), skip=1, cskip=1))
    h      = np.array(fread(os.path.join(indir,hfile), skip=1, cskip=1))
    try:
        m      = np.array(fread(meteofile, cskip=2, nc=1))
        m      = np.where(m.astype(int) == novalue, np.NaN, m)
    except TypeError:
        print('EddyCorrWarning: No meteorology file found, no fitting with rH can be done!')
        m      = False
    cout   = np.empty((np.shape(c)[0],1))
    hout   = np.empty((np.shape(h)[0],1))

    ############################################################################
    # checking for breakpoints
    print("\nBREAKPOINT CHECK:\n")
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        # plot lag c
        plt.figure(1)
        plt.plot(c[:,0], c[:,2], 'bo', label='minlag (sam)')
        plt.plot(c[:,0], c[:,5], 'ro', label='maxlag (sam)')
        plt.xlabel('DOY')
        plt.ylabel('Lags [Samples]')
        plt.title('c')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()
        plt.show()

    breaks = [0]
    inp = True
    while inp:
        inp = raw_input("Breakpoints? [ddd.ddd or n]: ")
        if inp.lower() == 'n':
            inp = False
        else:
            try:
                breaks += [np.min(np.where(np.abs(doysfloat-float(inp))<0.02083))]
            except (ValueError, IndexError):
                print('EddyCorrWarning: type in floats(with . not ,), integers or n, nothing else!')
            inp = True
    breaks += [np.shape(c)[0]]

    ############################################################################
    # calling the calculation function
    if np.all(m) != False:
        for i in range(len(breaks)-1):
            cout[breaks[i]:breaks[i+1],0], hout[breaks[i]:breaks[i+1],0] = calc(c[breaks[i]:breaks[i+1]],
                                                                            h[breaks[i]:breaks[i+1]],
                                                                            m[breaks[i]:breaks[i+1]],
                                                                            doys[breaks[i]:breaks[i+1]],
                                                                            histstep, indir, plot=plot)
    else:
        for i in range(len(breaks)-1):
            cout[breaks[i]:breaks[i+1],0], hout[breaks[i]:breaks[i+1],0] = calc(c[breaks[i]:breaks[i+1]],
                                                                            h[breaks[i]:breaks[i+1]],
                                                                            m,
                                                                            doys[breaks[i]:breaks[i+1]],
                                                                            histstep, indir, plot=plot)

    ################################################################################
    # writing output file
    output = csv.writer(open(os.path.join(indir,outfile), 'w'))

    for i in xrange(np.shape(cout)[0]):
        output.writerow([np.abs(cout[i][0].astype(int)), np.abs(hout[i][0].astype(int))])

    ################################################################################
    # attaching lags to meteo file
    if attach:
        meteo = np.loadtxt(meteofile, dtype='|S21', delimiter=',')
        meteo = np.concatenate((meteo, np.abs(cout.astype(int)), np.abs(hout.astype(int))), 1)
        np.savetxt(meteofile, meteo, fmt='%s', delimiter=',')

    ############################################################################
    # plot c lag output

    print("\nFINAL LAG OUTPUT:\n")
    if plot:
        fig1 = plt.figure(7)
        sub = fig1.add_subplot(111)
        sub.plot(c[:,0], cout.astype(int), 'ro')
        plt.xlabel('doy')
        plt.ylabel('lag (sam)')
        sub.set_title('c lag ouput')
        sub.axis('auto')
        sub.grid('on')

        # plot h lag output
        fig2 = plt.figure(8)
        sub = fig2.add_subplot(111)
        sub.plot(h[:,0], hout.astype(int), 'ro')
        plt.xlabel('doy')
        plt.ylabel('h lag (sam)')
        sub.set_title('h lag ouput')
        sub.axis('auto')
        sub.grid('on')
        plt.show()

        ############################################################################
        # save figures
        pp1 = pdf.PdfPages(os.path.join(indir,'c_lags%s.pdf'%(outfile[4:-4])))
        pp2 = pdf.PdfPages(os.path.join(indir,'h_lags%s.pdf' %(outfile[4:-4])))
        fig1.savefig(pp1, format='pdf')
        fig2.savefig(pp2, format='pdf')
        pp1.close()
        pp2.close()

############################################################################
# calculations
def calc(c, h, m, doys, histstep, indir, plot=False):
    '''
    part of eddycorr: plotting and fitting of hlag to rH and witing log file
    '''
    ############################################################################
    # histograms of lag distribution c
    histcmin, bin_edgescmin = np.histogram(c[:,2], bins=np.arange(np.min(c[:,2]),
                                           np.max(c[:,2])+histstep,histstep),
                                           range=None, normed=False, weights=None)
    histcmax, bin_edgescmax = np.histogram(c[:,5], bins=np.arange(np.min(c[:,5]),
                                           np.max(c[:,5])+histstep,histstep),
                                           range=None, normed=False, weights=None)

    # histograms of lag distribution h
    histhmin, bin_edgeshmin = np.histogram(h[:,2], bins=np.arange(np.min(h[:,2]),
                                           np.max(h[:,2])+histstep,histstep),
                                           range=None, normed=False, weights=None)
    histhmax, bin_edgeshmax = np.histogram(h[:,5], bins=np.arange(np.min(h[:,5]),
                                           np.max(h[:,5])+histstep,histstep),
                                           range=None, normed=False, weights=None)

    ############################################################################
    # plot lag c
    print("\nCO2 LAGS:\n")
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        plt.figure(1)
        plt.plot(c[:,0], c[:,2], 'bo', label='minlag (sam)')
        plt.plot(c[:,0], c[:,5], 'ro', label='maxlag (sam)')
        plt.xlabel('DOY')
        plt.ylabel('Lags [Samples]')
        plt.title('c')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()

        # plot histogram lag c
        plt.figure(2)
        plt.bar(bin_edgescmin[0:-1:1], histcmin, width=histstep/2-1, color='b', label='minlag (sam)')
        plt.bar(bin_edgescmax[0:-1:1]+histstep/2, histcmax, width=histstep/2-1, color='r', label='maxlag (sam)')
        plt.xlabel('Classes of lags [samples]')
        plt.ylabel('count')
        plt.title('c')
        plt.grid('on')
        plt.legend()
        plt.show()

    ctop = float(raw_input("Top of CO2 lag range: "))
    cbottom = float(raw_input("Bottom of CO2 lag range: "))

    ############################################################################
    # plot lag h
    print("\nH2O LAGS:\n")
    if plot:
        plt.figure(3)
        plt.plot(h[:,0], h[:,2], 'bo', label='minlag (sam)')
        plt.plot(h[:,0], h[:,5], 'ro', label='maxlag (sam)')
        plt.xlabel('DOY')
        plt.ylabel('Lags [Samples]')
        plt.title('h')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()

        # plot histogram lag h
        plt.figure(4)
        plt.bar(bin_edgeshmin[0:-1:1], histhmin, width=histstep/2-1, color='b', label='minlag (sam)')
        plt.bar(bin_edgeshmax[0:-1:1]+histstep/2, histhmax, width=histstep/2-1, color='r', label='maxlag (sam)')
        plt.xlabel('Classes of lags [samples]')
        plt.ylabel('count')
        plt.title('h')
        plt.grid('on')
        plt.legend()
        plt.show()

    print("! Only maxlag (sam) can be used for H2O !")
    htop    = float(raw_input("Top of H2O lag range: "))
    hmedian = float(raw_input("Median of  H2O lag range: "))
    hbottom = float(raw_input("Bottom of H2O lag range: "))

    ############################################################################
    # preparing data for regression

    print("\nCHECKING lag-rH FIT:\n")

    hsub = h[np.where((h[:,5]<=htop)&(h[:,5]>=hbottom))[0],5]
    msub = m[np.where((h[:,5]<=htop)&(h[:,5]>=hbottom))[0],0]
    if plot:
        # plot max lag vs. rH
        fig3 = plt.figure(6)
        sub = fig3.add_subplot(111)
        sub.plot(msub, hsub, 'ro', label='measurements')
        plt.xlabel('rH [%]')
        plt.ylabel('maxlag (sam)')
        sub.set_title('h maxlag (sam) - rH correlation')
        sub.axis('auto')
        sub.grid('on')
        plt.legend()
        plt.show()

    pbt = 't'
    while pbt=='t':
        mmax = float(raw_input("Top of rel. hum: "))
        hmin = float(raw_input("Maximum of  H2O lag: "))
        valid=(msub<mmax) & (hsub<hmin)
        p_guess1 = float(raw_input("Initial guess - Lag offset: "))*-1.
        p_guess2 = float(raw_input("Initial guess - rH multiplier: "))

        p, ff = fit(-hsub[valid], msub[valid], func=f, p_guess=[p_guess1,p_guess2],plot=plot)
        print('Offset=%f'%(p[0]*-1.), 'Multiplier=%f'%(p[1]))
        pbt = raw_input("(g)ood fit - proceed, (b)ad fit - proceed, (t)ry again: ")
        if pbt not in ['g', 'b', 't']:
            pbt = 't'

    ############################################################################
    # lags for c for output
    cout = np.where((c[:,2]<=ctop) & (c[:,2]>=cbottom), c[:,2],
                    np.where((c[:,5]<=ctop) & (c[:,5]>=cbottom), c[:,5],
                             (ctop+cbottom)/2))

    # lags for h for output
    if pbt=="g":
        hout = np.where((h[:,5]>htop)^(h[:,5]<hbottom), -f_inv(m[:,0], p), h[:,5])
        hout = np.where((hout>htop)^(hout<hbottom), hmedian, hout)
        hout = np.where(np.isnan(hout), hmedian, hout)
    else:
        hout = np.where((h[:,5]>htop)^(h[:,5]<hbottom), hmedian, h[:,5])
        hout = np.where(np.isnan(hout), hmedian, hout)

    ################################################################################
    # writing log file
    log = open(os.path.join(indir,'lag_%i_%02i_%02i_%02i_%02i_%02i.log' %(
                                                    t.localtime()[0], t.localtime()[1],
                                                    t.localtime()[2], t.localtime()[3],
                                                    t.localtime()[4], t.localtime()[5],))
                                                    , 'w')
    log.write('From doy %s %s:%s to doy %s %s:%s: \n' %(doys.flatten()[0][5:8], doys.flatten()[0][8:10], doys.flatten()[0][10:12], doys.flatten()[-1][5:8], doys.flatten()[-1][8:10], doys.flatten()[-1][10:12]))
    log.write('Top of CO2 lag range: %i\n' %(int(ctop)))
    log.write('Bottom of CO2 lag range: %i\n' %(int(cbottom)))
    log.write('Top of H2O lag range: %i\n' %(int(htop)))
    log.write('Bottom of H2O lag range: %i\n' %(int(hbottom)))
    log.write('Median of  H2O lag range: %i\n' %(int(hmedian)))
    if pbt=="g":
        log.write('Moisture dependency:\n')
        log.write('Top of rH used for fitting: %f\n'%mmax)
        log.write('Maximum of  H2O lag used for fitting: %f\n'%hmin)
        log.write('Initial guess - Lag offset: %f\n'%(p_guess1*-1.))
        log.write('Initial guess - rH multiplier: %f\n'%p_guess2)
        log.write('Moisture dependency: -maxlag_h[sam] = %f + exp(rH/%f)'%(p[0], p[1]))
    else:
        log.write('Median used instead of moisture dependency.')

    return cout, hout

############################################################################
# fitting
def fit(x,y,func,p_guess,plot=False):
    '''
    eddycorr: fitting procedure
    '''
    xfit=np.ma.array(x, mask=np.isnan(x))
    yfit=np.ma.array(y, mask=np.isnan(y))

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        x_mod = np.arange(np.ma.min(xfit),np.ma.max(xfit),(np.ma.max(xfit)-np.ma.min(xfit))/50.)
        plt.figure('pre-fit')
        plt.plot(xfit, yfit, 'bo')
        plt.plot(x_mod, func(x_mod,p_guess), 'k-')
        plt.show()

    # fit
    p_opt = fmin(absdif, p_guess, args=(xfit, yfit, func), disp=0)

    if plot:
        plt.figure('post-fit')
        plt.plot(xfit, yfit, 'bo')
        plt.plot(x_mod, func(x_mod,p_opt), 'k-')
        plt.show()

    return p_opt, func

############################################################################
# h lag model
def f(x, p):
    '''
    eddycorr: model for relating relative humidity and -maxlag_h
    '''
    return p[1]*np.log(x-p[0])

def f_inv(x, p):
    '''
    eddycorr: model for relating relative humidity and -maxlag_h
    '''
    return p[0]+np.exp(x/p[1])

def absdif(p,x,y,func):
    '''
    eddycorr: objective function
    '''
    return np.sum(np.abs(y-func(x, p)))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
