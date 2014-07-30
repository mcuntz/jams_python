import numpy as np
import sread, fread
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import pylab as pl
from scipy import interpolate
import csv
import sys
import time as t
import math as m
import shutil as sh

global missing_package
missing_package = False
try:
    from scipy.optimize import curve_fit # requires at least scipy 0.8.0
except ImportError:
    missing_package = True
    
def eddycorr(indir, sltdir, cfile, hfile, meteofile, outfile, novalue=-9999,
             histstep=10, attach=True):
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
                available *.slt files) where the lags will be attached
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
    Written,  AP, Jul 2014
    '''
        
    ############################################################################
    # missing package warning
    if missing_package == True:
        raise ValueError('eddycorr: scipy.optimize.curve_fit can not be found!')
    
    ############################################################################
    # move correlation files from sltdir to indir
    sh.move('%s\%s' %(sltdir, cfile), indir)
    sh.move('%s\%s' %(sltdir, hfile), indir)
    
    ############################################################################
    # reading input file
    header = sread.sread('%s/%s' %(indir,cfile))[0]
    doys   = np.array(sread.sread('%s/%s' %(indir,cfile), nc=1, skip=1), dtype='|S16')
    #doys   = np.array([x[5:12] for x in doys.flatten()], dtype = '|S7')
    day    = np.array([x[5:8] for x in doys.flatten()], dtype = '|S7').astype(float)
    hour   = np.array([x[8:10] for x in doys.flatten()], dtype = '|S2').astype(float)
    min    = np.array([x[10:12] for x in doys.flatten()], dtype = '|S2').astype(float)
    doysfloat   = day + (hour + min/60.)/24.
    c      = np.array(fread.fread('%s/%s' %(indir,cfile), skip=1, cskip=1))
    h      = np.array(fread.fread('%s/%s' %(indir,hfile), skip=1, cskip=1))
    try:
        m      = np.array(fread.fread(meteofile, cskip=1, nc=1))
        m      = np.where(m.astype(int) == novalue, pl.NaN, m)
    except TypeError: 
        print 'EddyCorrWarning: No meteorology file found, no fitting with rH can be done!'
        m      = False
    cout   = np.empty((np.shape(c)[0],1))
    hout   = np.empty((np.shape(h)[0],1))
    
    ############################################################################
    # checking for breakpoints
    print("\nBREAKPOINT CHECK:\n")
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
                print 'EddyCorrWarning: type in floats(with . not ,), integers or n, nothing else!'
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
                                                                            histstep, indir)
    else:
        for i in range(len(breaks)-1):
            cout[breaks[i]:breaks[i+1],0], hout[breaks[i]:breaks[i+1],0] = calc(c[breaks[i]:breaks[i+1]],
                                                                            h[breaks[i]:breaks[i+1]],
                                                                            m,
                                                                            doys[breaks[i]:breaks[i+1]],
                                                                            histstep, indir)
    
    ################################################################################
    # writing output file
    output = csv.writer(open('%s/%s' %(indir,outfile), 'wb')) 
                                                            
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
    pp1 = pdf.PdfPages('%s/c_lags%s.pdf' %(indir,outfile[4:-4]))
    pp2 = pdf.PdfPages('%s/h_lags%s.pdf' %(indir,outfile[4:-4]))
    fig1.savefig(pp1, format='pdf')
    fig2.savefig(pp2, format='pdf')
    pp1.close()
    pp2.close()

############################################################################
# h lag model
def mod(h,ms,hs,ho):
    '''
    relates water time lag with relative humidity
    '''
    return hs*np.exp(-h*ho)+ms

############################################################################
# calculations
def calc(c, h, m, doys, histstep, indir):
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
    
    print "! Only maxlag (sam) can be used for H2O !"
    htop    = float(raw_input("Top of H2O lag range: "))
    hmedian = float(raw_input("Median of  H2O lag range: "))
    hbottom = float(raw_input("Bottom of H2O lag range: "))          

    ############################################################################
    # preparing data for regression
    
    print("\nCHECKING lag-rH FIT:\n")        
    
    try:
        hsub = np.transpose(h[np.where((h[:,5]<=htop)&(h[:,5]>=hbottom)),5])
        hsubreg = -hsub.flatten()+np.max(hsub)
        msub = np.transpose(m[np.where((h[:,5]<=htop)&(h[:,5]>=hbottom)),0])
        msubreg = msub.flatten()-np.min(msub)
        
        # curve fitting 
        popti, pcovi = curve_fit(mod, hsubreg, msubreg, maxfev = 1000)
        
        # coefficient of determination
        sstot=np.sum((msubreg-np.mean(msubreg))**2)
        sserr=np.sum((msubreg-mod(hsubreg,popti[0],popti[1],popti[2]))**2)
        rsquare=1-sserr/sstot
    
        ############################################################################
        # plot rH
        plt.figure(5)                                        
        plt.plot(h[:,0], m[:,0], 'bo')                      
        plt.xlabel('DOY')
        plt.ylabel('rH [%]')                    
        plt.title('rH')
        plt.axis('auto')
        plt.grid('on')                  
              
        # plot max lag vs. rH
        fig3 = plt.figure(6)
        sub = fig3.add_subplot(111)                                       
        sub.plot(hsub, msub, 'ro', label='measurements')
        sub.plot(-np.arange(0, np.max(hsubreg)+histstep, histstep)+np.max(hsub),
                 mod(np.arange(0, np.max(hsubreg)+histstep, histstep),popti[0],popti[1],popti[2])+np.min(msub)
                 ,'b-', label='exponential fit \n r^2 = %s' %(str(round(rsquare,2))))         
        plt.xlabel('maxlag (sam)')                               
        plt.ylabel('rH [%]')                    
        sub.set_title('h maxlag (sam) - rH correlation')
        sub.axis('auto')
        sub.grid('on')                            
        plt.legend()
        plt.show()
        
        tm = "n"
        decision = raw_input("Fit sufficient? (Y/N): ").lower()
        if decision == "y":
            pp3 = pdf.PdfPages('%s/h_fit%s.pdf' %(indir,outfile[4:-4]))
            fig3.savefig(pp3, format='pdf')
            pp3.close()
        else:
            tm = raw_input("Use median instead? (Y/N): ").lower()
            if tm == "y":
                pass
            else:
                sys.exit()
                
    except (RuntimeError, TypeError):
            print 'EddyCorrWarning: No sufficient fit can be done!'   
            tm = raw_input("Use median instead? (Y/N): ").lower()
            if tm == "y":
                pass
            else:
                sys.exit()
    
    ############################################################################
    # lags for c for output
    cout = np.where((c[:,2]<=ctop) & (c[:,2]>=cbottom), c[:,2],
                    np.where((c[:,5]<=ctop) & (c[:,5]>=cbottom), c[:,5],
                             (ctop+cbottom)/2))
    
    # lags for h for output
    if tm =="n":
        hout = np.where((h[:,5]>htop)^(h[:,5]<hbottom), (np.log(( m[:,0].flatten()-np.min(msub) -popti[0])/popti[1])/popti[2])+np.max(hsub),
                                                        h[:,5])
        hout = np.where((np.isnan(hout)) & (m[:,0]>mod(-htop+np.max(hsub),popti[0],popti[1],popti[2])+np.min(msub)), htop, hout)
        hout = np.where((np.isnan(hout)) & (m[:,0]<mod(-hbottom+np.max(hsub),popti[0],popti[1],popti[2])+np.min(msub)), hbottom, hout)
        hout = np.where(np.isnan(hout), hmedian, hout)
        hout = np.where(hout<hbottom, hbottom, hout)
    else:
        hout = np.where((h[:,5]>htop)^(h[:,5]<hbottom), hmedian, h[:,5])
        hout = np.where(np.isnan(hout), hmedian, hout)
    
    ################################################################################
    # writing log file
    log = open('%s/lag_%i_%02i_%02i_%02i_%02i_%02i.log' %(indir, 
                                                    t.localtime()[0], t.localtime()[1], 
                                                    t.localtime()[2], t.localtime()[3],
                                                    t.localtime()[4], t.localtime()[5],)
                                                    , 'w')
    log.write('From doy %s %s:%s to doy %s %s:%s: \n' %(doys.flatten()[0][5:8], doys.flatten()[0][8:10], doys.flatten()[0][10:12], doys.flatten()[-1][5:8], doys.flatten()[-1][8:10], doys.flatten()[-1][10:12]))
    log.write('Top of CO2 lag range: %i\n' %(int(ctop)))    
    log.write('Bottom of CO2 lag range: %i\n' %(int(cbottom)))    
    log.write('Top of H2O lag range: %i\n' %(int(htop)))
    log.write('Bottom of H2O lag range: %i\n' %(int(hbottom)))    
    log.write('Median of  H2O lag range: %i\n' %(int(hmedian)))
    if tm != "y":
        log.write('Fit equation: rH[%] = rH_scale * exp(-maxlag_h[sam] * maxlag_h_scale) + rH_offset\n')
        log.write('rH_scale : %f\nmaxlag_h_scale : %f\nrH_offset : %f\n' %(popti[1], popti[2], popti[0]))
        log.write('r^2 : %s' %(str(round(rsquare,2))))
    else:        
        log.write('No fit between rH and maxlag_h[sam] could be performed.')
    
    return cout, hout
