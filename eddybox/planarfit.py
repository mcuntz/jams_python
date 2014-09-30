import numpy as np
import sread, fread
import csv
import os as os
import re
from math import pi
from scipy.stats import skew
import sys

def planarfit(indirpf, rawfile, outfile, pfmat='pfitmatrix.csv',
              pf0file='pfitdata0.csv', pf1file='pfitdata1.csv',
              pf2file='pfitdata2.csv', histsteps=50, plot=False):
    '''
    Extracts raw wind speeds from the raw flux file of EddyFlux as input for
    EDDYPFit. When EDDYFit is finished, the script loads the results and
    do plots. If user is satisfied, results are saved. 
    
    
    Definition
    ----------
    planarfit(indirpf, rawfile, outfile, pfmat='pfitmatrix.csv',
              pf0file='pfitdata0.csv', pf1file='pfitdata1.csv',
              pf2file='pfitdata2.csv', histsteps=50):
    
    Input
    ----- 
    indirpf     str, path of the folder where results will be saved
    rawfile     str, path of the file with raw wind speeds from EddyFlux
    outfile     str, name of the output file 
        
    
    Optional Input
    --------------
    pfmat       str, name of the pfitmatix file, default: 'pfitmatrix.csv'
    pf0file     str, name of the original wind speed file of EDDYPFit, default: 'pfitdata0.csv'
    pf1file     str, name of the one plane fit wind speed file of EDDYPFit, default: 'pfitdata1.csv'
    pf2file     str, name of the sectorial fit wind speed file of EDDYPFit, default: 'pfitdata2.csv'
    histstep    int, histogram steps for plotting (default=50)
                
    
    Output
    ------
    X_pfit.pdf  plot with planar fit
    X_uvw.csv   file with raw wind speeds
    X_wd.pdf    plot with wind rose
    X_wdis.pdf  plot with wind speed distributions 
    
    
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
    ############################################################################
    # reading raw file
    uvw   = np.array(fread.fread('%s' %rawfile, skip=1, cskip=13, nc=3))
    wdhor = np.array(fread.fread('%s' %rawfile, skip=1, cskip=20, nc=1))
    header = np.array(sread.sread('%s' %rawfile, cskip=13, nc=3), dtype='|S5')

    ############################################################################
    # coordinate transformation for u+v
    alpha = 120.*np.pi/180. #rotation angle in rad
    uvw_trans = np.copy(uvw)
    uvw_trans[:,0] = np.where((uvw[:,0]!=-9999) & (uvw[:,1]!=-9999),
                              uvw[:,0]*np.cos(alpha)+uvw[:,1]*np.sin(alpha),
                              -9999)# u
    uvw_trans[:,1] = np.where((uvw[:,0]!=-9999) & (uvw[:,1]!=-9999),
                              -uvw[:,0]*np.sin(alpha)+uvw[:,1]*np.cos(alpha),
                              -9999)# v
    
    ############################################################################
    # writing uvw file with wind speed for EDDYPFit
    file1 = open('%s/%s_uvw.csv' %(indirpf,outfile[:-4]), 'wb')
    output = csv.writer(file1)    
    output.writerow(header[0])                             
    for i in xrange(np.shape(uvw_trans)[0]):
        output.writerow(uvw_trans[i])
    file1.close()
    
    ############################################################################
    # user input to continue
    print "Do EddyPFit with the 'uvw.csv' file now!"
    ui1 = raw_input("Ready or quit (y/n)?: ").lower()
    if ui1 != "y":
        sys.exit()

    ############################################################################
    # reading pfit files
    header0 = np.array(sread.sread('%s/%s' %(indirpf,pf0file)) , dtype='|S1')
    uvw0    = np.array(fread.fread('%s/%s' %(indirpf,pf0file), skip=1), dtype=np.float)
    uvw0_trans = np.copy(uvw0)
    uvw0_trans[:,0] = uvw0[:,0]*np.cos(alpha) - uvw0[:,1]*np.sin(alpha)
    uvw0_trans[:,1] = uvw0[:,0]*np.sin(alpha) + uvw0[:,1]*np.cos(alpha)
    
    header1 = np.array(sread.sread('%s/%s' %(indirpf,pf1file)) , dtype='|S1')
    uvw1    = np.array(fread.fread('%s/%s' %(indirpf,pf1file), skip=1), dtype=np.float)
    uvw1_trans = np.copy(uvw1)
    uvw1_trans[:,0] = uvw1[:,0]*np.cos(alpha) - uvw1[:,1]*np.sin(alpha)
    uvw1_trans[:,1] = uvw1[:,0]*np.sin(alpha) + uvw1[:,1]*np.cos(alpha)
            
    header2 = np.array(sread.sread('%s/%s' %(indirpf,pf2file)) , dtype='|S1')
    uvw2    = np.array(fread.fread('%s/%s' %(indirpf,pf2file), skip=1), dtype=np.float)
    uvw2_trans = np.copy(uvw2)
    uvw2_trans[:,0] = uvw2[:,0]*np.cos(alpha) - uvw2[:,1]*np.sin(alpha)
    uvw2_trans[:,1] = uvw2[:,0]*np.sin(alpha) + uvw2[:,1]*np.cos(alpha)
    
    ############################################################################
    # plots
    # define grid
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec                       
        from matplotlib.mlab import griddata
        import matplotlib.cm as cm
        import matplotlib.mlab as mlab
        import matplotlib.backends.backend_pdf as pdf
        x0 = np.linspace(np.min(np.minimum(uvw0_trans[:,0],uvw0_trans[:,1])),
                         np.max(np.maximum(uvw0_trans[:,0],uvw0_trans[:,1])),500)
        y0 = np.linspace(np.min(np.minimum(uvw0_trans[:,0],uvw0_trans[:,1])),
                         np.max(np.maximum(uvw0_trans[:,0],uvw0_trans[:,1])),500)
        
        x1 = np.linspace(np.min(np.minimum(uvw1_trans[:,0],uvw1_trans[:,1])),
                         np.max(np.maximum(uvw1_trans[:,0],uvw1_trans[:,1])),500)
        y1 = np.linspace(np.min(np.minimum(uvw1_trans[:,0],uvw1_trans[:,1])),
                         np.max(np.maximum(uvw1_trans[:,0],uvw1_trans[:,1])),500)
        
        x2 = np.linspace(np.min(np.minimum(uvw2_trans[:,0],uvw2_trans[:,1])),
                         np.max(np.maximum(uvw2_trans[:,0],uvw2_trans[:,1])),500)
        y2 = np.linspace(np.min(np.minimum(uvw2_trans[:,0],uvw2_trans[:,1])),
                         np.max(np.maximum(uvw2_trans[:,0],uvw2_trans[:,1])),500)
        
        # grid the data.
        z0 = griddata(uvw0_trans[:,0],uvw0_trans[:,1],uvw0_trans[:,2],x0,y0, interp='nn')
        z1 = griddata(uvw1_trans[:,0],uvw1_trans[:,1],uvw1_trans[:,2],x1,y1, interp='nn')
        z2 = griddata(uvw2_trans[:,0],uvw2_trans[:,1],uvw2_trans[:,2],x2,y2, interp='nn')
        
        # plotting contours    
        fig1 = plt.figure(1, figsize=(6,13))
        sub1 = fig1.add_subplot(311, aspect=1)
        fillings = sub1.contourf(x0,y0,z0,20,cmap=plt.cm.jet)
        scat = sub1.scatter(uvw0_trans[:,0],uvw0_trans[:,1],marker='o',c='b',s=0.2,zorder=10)
        cbar = fig1.colorbar(fillings, orientation='vertical')
        xlimits = sub1.get_xlim()
        sub1.plot(np.array([xlimits[0],xlimits[1]]),np.array([0,0]), c='k')
        ylimits = sub1.get_ylim()
        sub1.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='k')
        sub1.set_title('Original wind components\nwith point data')
        cbar.set_label('w [m/s]')
        plt.ylabel('v [m/s]')
    
        sub2 = fig1.add_subplot(312, aspect=1)
        fillings = sub2.contourf(x1,y1,z1,20,cmap=plt.cm.jet)
        cbar = fig1.colorbar(fillings, orientation='vertical')
        xlimits = sub2.get_xlim()
        sub2.plot(np.array([xlimits[0],xlimits[1]]),np.array([0,0]), c='k')
        ylimits = sub2.get_ylim()
        sub2.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='k')
        sub2.set_title('One plane')
        cbar.set_label('w [m/s]')
        plt.ylabel('v [m/s]')
    
        sub3 = fig1.add_subplot(313, aspect=1)
        fillings = sub3.contourf(x2,y2,z2,20,cmap=plt.cm.jet)
        cbar = fig1.colorbar(fillings, orientation='vertical')
        xlimits = sub3.get_xlim()
        sub3.plot(np.array([xlimits[0],xlimits[1]]),np.array([0,0]), c='k')
        ylimits = sub3.get_ylim()
        sub3.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='k')
        sub3.set_title('Sectorial')
        cbar.set_label('w [m/s]')
        plt.xlabel('u [m/s]')
        plt.ylabel('v [m/s]')
    
        # plotting histograms
        mi = np.min(np.array([np.min(uvw0_trans[:,2]),np.min(uvw1_trans[:,2]),
                              np.min(uvw2_trans[:,2])]))
        ma = np.max(np.array([np.max(uvw0_trans[:,2]),np.max(uvw1_trans[:,2]),
                              np.max(uvw2_trans[:,2])]))        
        steps = np.abs((ma-mi)/histsteps)
        bins = np.arange(mi,ma+steps,steps)
        
        fig2 = plt.figure(2, figsize=(6,13))
        fig2.subplots_adjust(hspace=0.3)
        sub4 = fig2.add_subplot(311)
        n0, bins0, patches0 = sub4.hist(uvw0_trans[:,2], bins, color= 'b', histtype='bar')
        ylimits = sub4.get_ylim()
        sub4.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='y', lw=3)
        plt.ylabel('count')
        plt.title('Original w-component:\navg(w)= %.2f, var(w)= %.4f, skew(w)= %.4f'
                  %(np.mean(uvw0_trans[:,2]), np.var(uvw0_trans[:,2]), skew(uvw0_trans[:,2])))
    
        sub5 = fig2.add_subplot(312)
        n1, bins1, patches1 = sub5.hist(uvw1_trans[:,2], bins, color= 'g', histtype='bar')
        ylimits = sub5.get_ylim()
        sub5.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='y', lw=3)
        plt.ylabel('count')
        plt.title('One plane w-component:\navg(w)= %.2f, var(w)= %.4f, skew(w)= %.4f'
                  %(np.mean(uvw1_trans[:,2]), np.var(uvw1_trans[:,2]), skew(uvw1_trans[:,2])))
        
        sub6 = fig2.add_subplot(313)
        n2, bins2, patches2 = sub6.hist(uvw2_trans[:,2], bins, color= 'r', histtype='bar')
        ylimits = sub6.get_ylim()
        sub6.plot(np.array([0,0]),np.array([ylimits[0],ylimits[1]]), c='y', lw=3)
        plt.xlabel('Classes [m/s]')
        plt.ylabel('count')
        plt.title('Sectorial w-component:\navg(w)= %.2f, var(w)= %.4f, skew(w)= %.4f'
                  %(np.mean(uvw2_trans[:,2]), np.var(uvw2_trans[:,2]), skew(uvw2_trans[:,2])))
    
        # wind rose
        fig3 = plt.figure(3, figsize=(6,6))
        pol = fig3.add_subplot(111, polar=True)
        hist, bin_edges= np.histogram(wdhor, bins=36, range=(0,360))
        x = 90-np.arange(5,365,10)
        x = [i*pi/180. for i in x]  # convert to radians
        pol.bar(x, hist, width=10*pi/180)
        pol.set_xticklabels([r'$\sf{90\degree}$',r'$\sf{45\degree}$',r'$\sf{0\degree}$',
                             r'$\sf{315\degree}$',r'$\sf{270\degree}$',r'$\sf{225\degree}$',
                             r'$\sf{180\degree}$',r'$\sf{135\degree}$'], fontsize=15)
        plt.title('Horizontal wind direction frequency')        
        
        plt.show()
    
    ############################################################################
    # user input for saving results
    print "Satisfied with the fit?\ny will save the figures, n will exit without saving!"
    ui2 = raw_input("(y/n)?: ").lower()
    if ui2 != "y":
        sys.exit()
    
    ############################################################################
    # save results
    if plot:
        pp1 = pdf.PdfPages('%s/%s_pfit.pdf'%(indirpf,outfile[:-4]))
        pp2 = pdf.PdfPages('%s/%s_wdis.pdf'%(indirpf,outfile[:-4]))
        pp3 = pdf.PdfPages('%s/%s_wd.pdf'%(indirpf,outfile[:-4]))
        fig1.savefig(pp1, format='pdf')
        fig2.savefig(pp2, format='pdf')
        fig3.savefig(pp3, format='pdf')
        pp1.close()
        pp2.close()
        pp3.close()
    
    print "Rename EddyPFit files?"
    ui3 = raw_input("(y/n)?: ").lower()
    if ui3 != "y":
        sys.exit()
    
    os.rename('%s/%s' %(indirpf,pf0file), '%s/%s_%s.csv' %(indirpf, outfile[:-4], pf0file[:-4]))
    os.rename('%s/%s' %(indirpf,pf1file), '%s/%s_%s.csv' %(indirpf, outfile[:-4], pf1file[:-4]))
    os.rename('%s/%s' %(indirpf,pf2file), '%s/%s_%s.csv' %(indirpf, outfile[:-4], pf2file[:-4]))
    os.rename('%s/%s' %(indirpf,pfmat), '%s/%s_%s.csv' %(indirpf, outfile[:-4], pfmat[:-4]))

if __name__ == '__main__':
    import doctest
    doctest.testmod()