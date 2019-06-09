#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.fread import fread
from jams.sread import sread
from jams.date2dec import date2dec
import time as t
import shutil as sh
import re
import os as os

################################################################################
def eddyspec(indir, cfile, hfile, rawfile, sltdir, tspfile='34_specmean.csv',
             cspfile='35_specmean.csv', hspfile='36_specmean.csv',
             novalue=-9999, plot=False):
    '''
    Provides plots with carbon and water fluxes from the rawfile together
    with the respective time lags in the cfile and hfile for the user to select
    'golden days' with sufficient flux and constant lags (about 6 half hour
    values during midday) to perform spectrum analysis with EddySpec and
    SpecMean (Kolle & Rebmann, 2007). When user selected desired time period,
    the respective *.slt files are printed to the console. The user can than
    use EddySpec and SpecMean with these files and the script continues after
    user is finished. Inductances for water and carbon are fitted and saved
    together with spectrum plots. It is recommended to repeat the spectrum
    analysis for different times within the year.


    Definition
    ----------
    eddyspec(indir, cfile, hfile, rawfile, sltdir, tspfile='34_specmean.csv',
             cspfile='35_specmean.csv', hspfile='36_specmean.csv', novalue=-9999):


    Input
    -----
    indir       str, path of the folder where results will be saved
    cfile       str, path of the carbon lag file
    hfile       str, path of the water lag file
    rawfile     str, path of the raw flux file
    sltdir      str, path of the folder containing the *.slt files (will contain
                EddySpec and SpecMean files after use)


    Optional Input
    --------------
    tspfile     str, name of the average temperature spectrum file
                (default: '34_specmean.csv')
    cspfile     str, name of the average carbon spectrum file
                (default: '35_specmean.csv')
    hspfile     str, name of the average water spectrum file
                (default: '36_specmean.csv')
    novalue     int, novalue in rawfile (default=-9999)


    Output
    ------
    (for each run in the year, a new folder is created in indir and EddySpec and
    SpecMean files are moved here from sltdir)
    c_specw.pdf    plot of the carbon spectrum
    h_specw.pdf    plot of the water spectrum
    spec_X_X.log   log file containing used *.slt files and inductances


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2014 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  AP, Jul 2014
    Modified, MC, Aug 2014 - clean up and Python 3
    '''
    ############################################################################
    # reading input files
    #lags
    lagdate   = np.array(sread('%s' %(cfile), nc=1, skip=1), dtype='|S16')
    day       = np.array([x[5:8] for x in lagdate.flatten()], dtype = '|S7').astype(float)
    hour      = np.array([x[8:10] for x in lagdate.flatten()], dtype = '|S2').astype(float)
    min       = np.array([x[10:12] for x in lagdate.flatten()], dtype = '|S2').astype(float)
    doysfloat = day + (hour + min/60.)/24.
    cl        = np.array(fread('%s' %(cfile), skip=1, cskip=1))
    hl        = np.array(fread('%s' %(hfile), skip=1, cskip=1))
    # fluxes
    fluxdate = date2dec(ascii = np.array(sread('%s' %(rawfile), nc=1, skip=1), dtype='|S16'))
    year     = np.array(sread('%s' %(rawfile), nc=1, skip=1), dtype='|S16')
    year     = np.array([x[6:10] for x in year.flatten()], dtype = '|S4').astype(int)
    fluxdate = fluxdate.flatten() - date2dec(yr = year, mo = 1, dy = 1, hr = 0, mi = 0, sc = 0)
    cf       = np.array(fread('%s' %(rawfile), skip=1, cskip=4, nc=1))
    cf       = np.where(cf == novalue, np.nan, cf).flatten()
    hf       = np.array(fread('%s' %(rawfile), skip=1, cskip=1, nc=1))
    hf       = np.where(hf == novalue, np.nan, hf).flatten()

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        ############################################################################
        # plot lag c
        plt.figure(1)
        plt.plot(cl[:,0], cl[:,2], 'bo', label='minlag (sam)')
        plt.plot(cl[:,0], cl[:,5], 'ro', label='maxlag (sam)')
        plt.xlabel('DOY')
        plt.ylabel('Lags [Samples]')
        plt.title('c lags')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()

        ############################################################################
        # plot flux c
        plt.figure(2)
        plt.plot(fluxdate, cf, 'r-', label='cflux')
        plt.xlabel('DOY')
        plt.ylabel('cflux [?]')
        plt.title('c fluxes')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()

        ############################################################################
        # plot lag h
        plt.figure(3)
        plt.plot(hl[:,0], hl[:,2], 'bo', label='minlag (sam)')
        plt.plot(hl[:,0], hl[:,5], 'ro', label='maxlag (sam)')
        plt.xlabel('DOY')
        plt.ylabel('Lags [Samples]')
        plt.title('h lags')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()

        ############################################################################
        # plot flux h
        plt.figure(4)
        plt.plot(fluxdate, hf, 'b-', label='hflux')
        plt.xlabel('DOY')
        plt.ylabel('hflux [?]')
        plt.title('h fluxes')
        plt.axis('auto')
        plt.grid('on')
        plt.legend()
        plt.show()

    interval = []
    inp = True
    while inp:
        inpfrom = raw_input("From doy [ddd.ddd]: ")
        inpto = raw_input("To doy [ddd.ddd]: ")
        try:
            interval += [np.min(np.where(np.abs(doysfloat-float(inpfrom))<0.02083))]
            interval += [np.min(np.where(np.abs(doysfloat-float(inpto))<0.02083))]
            inp = False
        except (ValueError, IndexError):
            print('EddySpecWarning: type in floats(with . not ,), nothing else!')
            inp = True
    inpclag = raw_input("C lag: ")
    inphlag = raw_input("H lag: ")

    print(lagdate[interval[0]:interval[1]+1])
    print('\nDO EDDYSPEC AND SPECMEAN NOW!\n')
    inp = raw_input("Mean spectrums ready? [y/n]: ").lower()
    if inp == "y":
        pass
    else:
        import sys
        sys.exit()

    ############################################################################
    # move spectrum files from slt folder to spec folder
    if len(os.listdir(indir))==0:
        os.mkdir(indir+'/1')
        indir = indir+'/1'
    else:
        os.mkdir(indir+'/%i'%(int(os.listdir(indir)[-1])+1))
        indir = indir+'/'+os.listdir(indir)[-1]
    specpat  = re.compile('[0-9]*_specmean.csv')
    specpat2  = re.compile('[a-zA-Z0-9]*_[0-9]*_spec.csv')
    sltdirlist = os.listdir(sltdir)
    for file in sltdirlist:
        if bool(re.search(specpat, file)) | bool(re.search(specpat2, file)):
            sh.move('%s/%s' %(sltdir, file), indir)

    ############################################################################
    # conductance fitting
    # reading mean spec files
    tsp = np.array(fread('%s/%s' %(indir,tspfile), skip=1, nc=2))
    csp = np.array(fread('%s/%s' %(indir,cspfile), skip=1, nc=2))
    hsp = np.array(fread('%s/%s' %(indir,hspfile), skip=1, nc=2))

    # filter input specs for log10(input) != NaN
    global fcoc, fcoh
    fcoc = csp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(csp[:,1])))))[0],1]
    fcoh = hsp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(hsp[:,1])))))[0],1]
    xc = tsp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(csp[:,1])))))[0],0]
    yc = tsp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(csp[:,1])))))[0],1]
    xh = tsp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(hsp[:,1])))))[0],0]
    yh = tsp[np.where((np.invert(np.isnan(np.log10(tsp[:,1])))) & (np.invert(np.isnan(np.log10(hsp[:,1])))))[0],1]

    # calculate deviations from t spec for different inductences and select the one with the smallest residual
    con = np.arange(0,5,0.01)
    devc = np.empty_like(con)
    devh = np.empty_like(con)

    for i in xrange(np.shape(con)[0]):
        devc[i] = np.nansum(np.log10(yc)-modc(xc,con[i]))
        devh[i] = np.nansum(np.log10(yh)-modh(xh,con[i]))

    conc = con[np.argmin(np.abs(devc))]
    conh = con[np.argmin(np.abs(devh))]

    if plot:
        ############################################################################
        # plot c spec
        fig1 = plt.figure(5)
        sub = fig1.add_subplot(111)
        sub.plot(np.log10(tsp[:,0]), np.log10(tsp[:,1]), 'ro-', label='T Spec')
        sub.plot(np.log10(csp[:,0]), np.log10(csp[:,1]), 'bo-', label='C Spec')
        sub.plot(np.log10(xc), modc(xc,conc),'go-', label='C Spec with %.2f Induc.' %(conc))
        plt.xlabel('F[Hz]')
        plt.ylabel('F[Hz]*Co(w,X)')
        sub.set_title('c Spectrum')
        sub.axis('auto')
        sub.grid('on')
        sub.legend(loc='best')

        ############################################################################
        # plot h spec
        fig2 = plt.figure(6)
        sub = fig2.add_subplot(111)
        sub.plot(np.log10(tsp[:,0]), np.log10(tsp[:,1]), 'ro-', label='T Spec')
        sub.plot(np.log10(hsp[:,0]), np.log10(hsp[:,1]), 'bo-', label='H Spec')
        sub.plot(np.log10(xh), modh(xh,conh),'go-', label='H Spec %.2f Induc.' %(conh))
        plt.xlabel('F[Hz]')
        plt.ylabel('F[Hz]*Co(w,X)')
        sub.set_title('h Spectrum')
        sub.axis('auto')
        sub.grid('on')
        sub.legend(loc='best')

        plt.show()

        ############################################################################
        # save figures
        pp1 = pdf.PdfPages('%s/c_spec%s.pdf' %(indir,rawfile[6:-4]))
        pp2 = pdf.PdfPages('%s/h_spec%s.pdf' %(indir,rawfile[6:-4]))
        fig1.savefig(pp1, format='pdf')
        fig2.savefig(pp2, format='pdf')
        pp1.close()
        pp2.close()

    ################################################################################
    # writing log file
    log = open('%s/spec_%i_%02i_%02i_%02i_%02i_%02i.log' %(indir,
                                                    t.localtime()[0], t.localtime()[1],
                                                    t.localtime()[2], t.localtime()[3],
                                                    t.localtime()[4], t.localtime()[5],)
                                                    , 'w')
    log.write('Inductence determination for the mean of:\n')
    for item in lagdate[interval[0]:interval[1]+1]:
        log.write('%s\n'%(item[0]))
    log.write('\nC lag: %s\n' %(inpclag))
    log.write('H lag: %s\n' %(inphlag))
    log.write('\nC Inductence: %.2f\n' %(conc))
    log.write('H Inductence: %.2f\n' %(conh))

# inductance fit functions for c
def modc(fc, conc):
    '''
    part of eddyspec: inductance model for c
    '''
    y = 1./(1. + 4. * np.pi**2 * fc**2 * conc**2)
    return np.log10(fcoc/y)

# inductance fit functions for h
def modh(fh, conh):
    '''
    part of eddyspec: inductance model for h
    '''
    y = 1./(1. + 4. * np.pi**2 * fh**2 * conh**2)
    return np.log10(fcoh/y)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
