#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.dec2date import dec2date
 
def ustarflag(date, data, flags, isday, outdir, min_thresh=0.1, nboot=100,
              ustar_v=2, plot=False):
    '''
    Friction velocity flagging for Eddy Covariance data. Assesses the threshold
    of friction velocity (ustar) below which a reduction in CO2 flux correlates
    with ustar and returns a flag here. Originally coded by Tino Rau.
    
    
    Definition
    ----------
    ustarflag(date, isday, data, flags, outdir, min_thresh=0.1, nboot=100,
              ustar_v=2, plot=False):
    
    
    Input
    ----- 
    date        np.array(N), julian date
    data        np.array(N,3), data array with CO2 flux (Fco2 [mumol/m2s),
                friction velocity (ustar [m/s]) and air temperature (T [degC])
    inflag      np.array(N,3), dtype=int, quality flag of data, 0 where data is
                good
    isday       np.array(N), dtype=bool, True where it is day and False where
                it is night
    outdir      str, path of the output folder
    
                        
    Optional Input
    --------------
    min_thresh  float, minimum ustar threshold, recommendation: 0.1 for
                overstorey towers, 0.01 for understorey towers (default: 0.1)
    nboot       int, number of boot straps (default: 100)
    ustar_v     int, value which shall be returned when a value is below the 
                ustar threshold (default: 2)
    udef        int/float, missing value of data (default: -9999) NaN values are
                excluded from computations anyhow.
    plot        bool, if True data and spikes are plotted (default: False)
    
    
    Output
    ------
    flag_out    np.array(N), flag array where everything is 0 except where
                values fall below ustar threshold, there it is ustar_v
    
    
    Restrictions
    ------------
    - works ONLY for a data set of ONE FULL year
    - works ONLY for half hourly time steps
    
    
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


    History
    -------
    Written,  AP, Aug 2014
    '''       

    ############################################################################
    # fixed parameters
    nperiods   = 4    # numer of seasons per year
    ntclass    = 6    # number of temperature classes
    corrthresh = 0.3  # correlation coefficient 
    nuclass    = 20   # number of u* classes
    upercent   = 95.  # percentage for average comparison 
    udef       =-9999 # undef of output
     
    ############################################################################
    # data and flags
    Fco2       = data[:,0]
    ustar      = data[:,1]
    T          = data[:,2]
    flag_Fco2  = flags[:,0]
    flag_ustar = flags[:,1]
    flag_T     = flags[:,2]
    
    years, months, days, hours, mins, sc = dec2date(date, fulldate=True)
    nmons = 12/nperiods
    yrmin = np.min(years)
    yrmax = np.max(years)
    nyears = yrmax - yrmin ######## + 1 works only for one full year of data
    if (nyears>1) | (np.unique(months).size<12):
        raise ValueError('ustarflag: only one full year of data can be processed')
    flag = (flag_Fco2==0)&(~isday)&(flag_T==0)&(flag_ustar==0)
    
    ############################################################################
    # prepare bootstrapping
    threshs = np.zeros((nboot,nyears))
    seasonal_threshs = np.zeros((nperiods*nboot,nyears))
    
    ############################################################################
    # calculate thresholds
    for y in xrange(yrmin, yrmax): ######## + 1 works only for one full year of data
        #print 'Year ', y
        periods= np.unique(months[np.where(years==y)])[0:-1:nmons]
        nperiods = len(np.unique(months[np.where(years==y)]))/nmons
        out_threshs = np.array([])
        nstepsyear = len(np.where(years==y)[0])
        #print nstepsyear
        
        ########################################################################
        # start bootstrapping
        for r in xrange(nboot):
            range_min = np.min(np.where(years==y)[0])
            range_max = np.max(np.where(years==y)[0])
            rand      = np.random.randint(0,nstepsyear,size=nstepsyear)
            T_b       = T[range_min+rand]
            flag_b    = flag[range_min+rand]
            years_b   = years[range_min+rand]
            months_b  = months[range_min+rand]
            ustar_b   = ustar[range_min+rand]
            Fco2_b    = Fco2[range_min+rand]
        
            # flag the created bootstrap-arrays
            T_f      = np.extract(flag_b, T_b)
            months_f = np.extract(flag_b, months_b)
            years_f  = np.extract(flag_b, years_b)
            ustar_f  = np.extract(flag_b, ustar_b)
            Fco2_f   = np.extract(flag_b, Fco2_b)
            
            # container for thresholds of the seasons of a year
            period_threshs = np.array([])
            
            ####################################################################
            # loop over seasons    
            for m in periods:            
                period_flag = (months_f>=m) & (months_f<m+nmons)
                
                # flag data --> just night_time and unflagged data is left
                T_m     = np.ma.array(T_f,     mask=~period_flag, fill_value=udef)
                ustar_m = np.ma.array(ustar_f, mask=~period_flag, fill_value=udef)
                Fco2_m  = np.ma.array(Fco2_f,  mask=~period_flag, fill_value=udef)
                
                # sort the array according to the temperature
                ii      = np.ma.argsort(T_m)
                T_m     = T_m[ii]
                ustar_m = ustar_m[ii]
                Fco2_m  = Fco2_m[ii]
                
                # separate into 6 temperature classes with
                # equal size according to quantiles
                class_width = np.ma.count(T_m)/ntclass
                ntlen       = T_m.size
                
                # container for the thresholds of this period
                class_threshs = np.array([])
                
                ################################################################
                # loop for every temperature class
                for i in xrange(ntclass):  
                    T_temp    = T_m[i*class_width:
                                    np.min([(i+1)*class_width, ntlen-1])]
                    u_temp    = ustar_m[i*class_width:
                                        np.min([(i+1)*class_width, ntlen-1])]
                    Fco2_temp = Fco2_m[i*class_width:
                                       np.min([(i+1)*class_width, ntlen-1])]
                    
                    # calculate r 
                    #try: #AP
                    r_abs = np.abs(np.corrcoef(T_temp,u_temp)[0,1])
                    #except IndexError: #AP
                    #    print 'IndexError'#AP
                    #    r_abs = 1 #AP
                        
                    # neglect threshold, if correlation is strong 
                    # online gap tool uses 0.3 as r-threshold, 
                    # original papale paper: 0.4
                    if r_abs >= corrthresh:
                        continue
        
                    # sort the data of the temperature class according
                    # to ustar 
                    ii        = np.ma.argsort(u_temp)
                    T_temp    = T_temp[ii]
                    u_temp    = u_temp[ii]
                    Fco2_temp = Fco2_temp[ii]
                    
                    # set the class width of the ustar classes; class
                    # with equal size??
                    u_class_width = np.ma.count(T_temp)/nuclass
                    nulen         = T_temp.size
                    
                    ############################################################
                    # loop over all ustar classses-1
                    # threshold = if class average > 99% of average of
                    # all classes above
                    # loop over all ustar classes
                    for u in xrange(nuclass-1):
                        T_u    = T_temp[u*u_class_width:
                                        np.min([(u+1)*u_class_width,nulen-1])]
                        u_u    = u_temp[u*u_class_width:
                                        np.min([(u+1)*u_class_width,nulen-1])]
                        Fco2_u = Fco2_temp[u*u_class_width:
                                           np.min([(u+1)*u_class_width,nulen-1])]
                        
                        rest_Fco2 = Fco2_temp[np.min([(u+1)*u_class_width,nulen-1]):]
                        rest_avg  = np.ma.mean(rest_Fco2)
                        avg       = np.ma.mean(Fco2_u)
                        
                        if np.abs(avg) >= np.abs(upercent/100.*rest_avg):
                            class_threshs = np.append(class_threshs, np.ma.mean(u_u))
                            break
                        
                # median of thresholds of all temperature classes is
                # set as the threshold of the period
                if ~np.isnan(np.ma.median(class_threshs)):
                    #print r
                    period_threshs = np.append(period_threshs,
                                               np.ma.median(class_threshs))
                else:
                    #print 'All are nans',class_threshs
                    period_threshs = np.append(period_threshs, udef)
                    
            # threshholds of all periods
            out_threshs = np.append(out_threshs,period_threshs) 
            # if there are values in period_thresholds take the maximum, 
            # else take 90th percentil
            if period_threshs == []:
                threshs[r,y-yrmin] = np.sort(ustar_f)[np.len(ustar_f)*9./10.]
            else:
                threshs[r,y-yrmin] = np.max(period_threshs)
        
        seasonal_threshs[0:nperiods*nboot,y-yrmin] = out_threshs    
    
    ############################################################################
    # set the flags
    
    # assign calculated threshold:
    # for each year separately,
    # for night and day-time data
    # take the median of the nboot bootstrapped thresholds for each year
    # "out_threshs" used for export and later plotting
    
    flag_out = np.zeros_like(flag_Fco2, dtype=np.int)
    
    for y in range(yrmin, yrmax): ######## + 1 works only for one full year data
        if np.any(threshs[:,y-yrmin] > 10):
            print('Hu', y, threshs[:,y-yrmin])
        med = np.median(threshs[:,y-yrmin])
        if (med < min_thresh):
            med = min_thresh
        elif np.isnan(med):
            print('NaN-WARNING!')
            ii = np.argsort(ustar[years==y])    
            oo =np.where(flag_Fco2[ii]==0)[0]
            flag_out[years==y][ii[oo][::np.int(len(oo)*0.9)]]\
                    += ustar_v
        else:
            ii = (years==y).flatten()
            yy = np.where((ustar[ii]<med) & (flag_Fco2[ii]==0))
            uu = np.clip(np.unique(np.concatenate((yy[0], yy[0]+1))), 0,
                         len(flag_Fco2[(years==y).flatten()])-1)
            flag_out[np.where(years==y)[0][uu]]+=ustar_v
            
            #print 'Year: %i, u_star thresh: %5.3f'%(y, med)
            
            ####################################################################
            # plot
            if plot:
                import matplotlib.pyplot as plt
                import matplotlib.backends.backend_pdf as pdf
                fig1 = plt.figure(1)
                sub1 = fig1.add_subplot(111)
                sub1.plot(ustar[np.where(flag & (years==y).flatten())],\
                          Fco2[np.where(flag & (years==y).flatten())], 'bo')
                sub1.axvline(x=med, linewidth=0.75, color='r')
                plt.ylabel('F CO2')
                plt.xlabel('u_star')
                plt.title('u_star thresh: %5.3f'%med)
                plt.show()

                pp1 = pdf.PdfPages(outdir+'/ustar.pdf')
                fig1.savefig(pp1, format='pdf')
                pp1.close()

    ############################################################################
    # save output    
    np.savetxt('%s/u_star_thresh.csv' %outdir, threshs,
               fmt='%4.4f', delimiter=',')
    np.savetxt('%s/u_star_thresh_seasonal.csv' %outdir, seasonal_threshs,
               fmt='%4.4f', delimiter=',')
    
    return flag_out

if __name__ == '__main__':
    import doctest
    doctest.testmod()
