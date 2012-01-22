#!/usr/bin/env python

import sys
import numpy as np
#import pdb
from ufz import dec2date

###############################################################
###############################################################
###############################################################
def gap_filling(tofill,Rg,Tair,vpd,dates,flag_tofill=[],\
        Rg_flag=[],Tair_flag=[],vpd_flag=[],Rg_dev=50.,\
        Tair_dev=2.5,vpd_dev=5,t_int = 48,TrueFilling=True):

    '''
        PURPOSE
        -------

        Fills gaps of flux data derived with Eddy-Covariance (EC) 
        technique according to:
        Reichstein et al. (2005):    
        On the separation of net ecosystem exchange into 
        assimilation and ecosystem respiration: review and 
        improved algorithm. Global Change Biology,11,9,
        p. 1424-1439.
        
        INPUT
        -----
        
        tofill     -> input array fluxes 
        Rg         -> input array gloabal radiation in Wm-2
        Tair       -> input array Air Temperature degree Celsius
        vpd        -> input array Vapour Pressure Deficit in hPa
        flag_tofill-> input array flags of fluxes (0: good
                        quality; >0 bad data)
        Rg_flag    -> input array: flags of global radiation 
        Tair_flag  -> input array: flags of air temperature 
        vpd_flag   -> input array: flags of vapour pressure deficit
        dates      -> date in decimal julian days
        
        PARAMETERS
        ----------
       
        TrueFilling -> True: real filling is performed; just 
                        gaps are filled; gaps are flags_tofill > 0
                       False: every data point is "filled";
        
        OPTIONAL ARGUMENTS
        ------------------
 
        Rg_dev  -> deviation threshold of global radiation (see "*")
        Tair_dev-> deviation threshold of Air Temperature (see "*")
        vpd_dev -> deviation threshold of vpd (see "*")
        t_ind   -> number of data points per day; default 48 

        * Deviations are used as follows: Meteorological condition at 
        time of gap is taken +- the deviation; in a time window 
        fluxes with these values (condition +- deviation) are averaged.
        
        OUTPUT
        ------
       
        filled -> Output numpy array filled data
        qc_gf -> quality classes of fills 
        
        Requires the 'numpy' package available at:
        
        http://numpy.scipy.org/ 

        History
        -------
        Written  TR, 2011
        Modified TR, Sep 2011 - loop over single years removed
                              - loops now over range of input data
        Modified MG, Sep 2011 - include standard deviation output
    '''    
# --------------------------------------------------------------------
# sets the quality class of the gap-filling (Reichstein et al. 2005)

    def set_qual(iteration,case):
        if case ==1:        
            if iteration <= 2:
                return 1
            elif iteration > 4:
                return 3
            else:
                return 2
        elif case == 2:
            if iteration ==0:
                return 1
            elif iteration <= 2:
                return 2
            else:
                return 3

    print 'GAP FILLING'
    # -------------------------------------------------------------
    # PARAMETERS

    # defaults
    # Algorithm checks for similar meteorological conditions 
    # (Global Radiation, Air Temperature and Vapour Pressure
    # Deficit) at which the gap occured e. g. it takes those
    # conditions and checks for flux values with the same
    # conditions in a time window
    # selection: all conditions +- deviation 
    
    # a measure to avoid extraploation into a long gap at 
    # the beginning and end of a year
    longestmarginalgap = 60

    udef = -9999.

    # epsilon: Machine variable for floating point calculations
    eps = sys.float_info.epsilon

    # number of data points per week; basic factor of the time 
    # window
    week = t_int * 7

    (years,months,days,hours,minutes,sc) = dec2date(dates,\
            fulldate=True)
    dec_time = hours + minutes/60
    
    # container for gap quality classes
    qc_gf = np.zeros((len(tofill),1))
    
    # array for filled values
    tofill_gf  = np.array(len(tofill)*[udef],dtype=np.float)
    # MGMG
    tofill_st  = np.array(len(tofill)*[udef],dtype=np.float)
    
    #------------------------------------------------------------------
    # GAP FILLING PROCEDURE
    
    # flag for all meteorological conditions
    flag_met = (Tair_flag == 0) & (vpd_flag == 0) & (Rg_flag == 0) 
    
    # flag for all meteorological conditions and NEE; needed to flag data\
    # to associate useful data (with flag_CO2 ==0)
    flag_total = (flag_met == True) & (flag_tofill == 0) 

    # flag for the case, that just global radiation is available 
    flag_sec = (Rg_flag == 0) & (flag_tofill == 0) 
    
    
    # separte calculation for each year
    largemargin = flag_tofill*0

    firstvalid = np.min(np.where(flag_tofill==0)[0])
    lastvalid  = np.max(np.where(flag_tofill==0)[0])
    
    if firstvalid > (t_int*longestmarginalgap):
        largemargin[ 0 :(firstvalid-t_int*longestmarginalgap)]=1
        
    if lastvalid < (len(largemargin)-t_int*longestmarginalgap):
        largemargin[(lastvalid+t_int*longestmarginalgap)::]=1
    
    # in order to separate years, get start and end indices
    # end index is used for computing the maximal possible index
    # in window size slicing
    #start_y = np.min(np.where(years==year))
    #end_y   = np.max(np.where(years==year))
            
    for j in range(len(tofill_gf)):      
    
        #if j == 12943 :
            #pdb.set_trace()
            # CHECK: NEE AVAILABLE? 
        if TrueFilling:
            if flag_tofill[j] == 0 or largemargin[j] == 1 :
                # no reason to go further, no gap--> continue
                continue
       ##################################################################### 
        # NO NEE AVAILABLE --> 3 METHODS OF GAP FILLING
       
        # 1. Tair, vpd and global radiation; 
        # 2. just global radiation; 
        # 3. none meteorolgical conditions: take the mean of\
        #                                   +- hour 
        # 4. take diurnal mean of adjacent days 
        #
        
        # METHOD ONE: ALL MET CONDITIONS #######################
        if (flag_met[j] == 1):
            #print 'all_met available' 
            
            # get the meteorological conditions 
            rg_cond     = Rg[j]
            tair_cond   = Tair[j]
            vpd_cond    = vpd[j]
            
            # search for values around the met-conditions in a window of time
            # (one week in the first iteration and odd weeks in the next)
            # at the edges of the years (first and last week) a negative index 
            # would be produced --> max function to set min to 0 and max to end
            j1 = [j-np.arange(1,week+1)+1] 
            j2 = [j+np.arange(1,week)]
            appended_jays = np.append(j1,j2)
            win_ind = np.sort(np.clip(appended_jays,0,len(flag_met)-1))
            
            # for better overview: dynamic calculation of radiation threshold
            # minimum 20; maximum 50 [Wm-2] according to private correspondence
            # with Markus Reichstein
            rg_dev = max(20,min(rg_cond,Rg_dev))

            # get boolean array where meteo-conditions are in a given width
            conditions = ((np.abs(rg_cond-Rg[win_ind])-rg_dev)<(-rg_dev*eps))&\
                         ((np.abs(tair_cond-Tair[win_ind])-Tair_dev)\
                         <(-Tair_dev*eps*10))&\
                         ((np.abs(vpd_cond-vpd[win_ind])-vpd_dev)<(-vpd_dev*eps))
            number_4avg = len(conditions[conditions&flag_total[win_ind]])
            # we need at least two samples with similar conditions
            
            if number_4avg>=2:
                # strip NEE values with window size
                period = tofill[win_ind]
                #print 'window',len(period)/48
                # calculate the average for gap filling of the point
                avg = np.mean(period[(conditions==True)&\
                        (flag_total[win_ind]==True)])
                                    #print avg
                tofill_gf[j] = avg
                #### MGMG
                std = np.std(period[(conditions==True)&\
                        (flag_total[win_ind]==True)])
                tofill_st[j] = std
                # MGMG
                # assign also quality category of gap filling
                qc_gf[j] = 1 
                #print 'did filling with met_all and one iteration'                
                #print 'my number',number_4avg 
               
                continue               
                                  
            # If less than three samples with similar meteo conditions are found  
            # --> extend time window to two weeks

            else:
                j1 = [j-np.arange(1,2*week+1)+1] 
                j2 = [j+np.arange(1,2*week)]
                appended_jays = np.append(j1,j2)
                win_ind = np.sort(np.clip(appended_jays,0,len(flag_met)-1))
                
                conditions = ((np.abs(rg_cond-Rg[win_ind])-rg_dev)<(-rg_dev*eps))&\
                         ((np.abs(tair_cond-Tair[win_ind])-Tair_dev)\
                         <(-Tair_dev*eps*10))&\
                         ((np.abs(vpd_cond-vpd[win_ind])-vpd_dev)<(-vpd_dev*eps))
    
                number_4avg = len(conditions[conditions&flag_total[win_ind]])
                                    
                if number_4avg>=2:
                    period = tofill[win_ind]
                    avg = np.mean(period[(conditions==True)&\
                            (flag_total[win_ind]==True)])    
                    tofill_gf[j] = avg
                    #### MGMG
                    std = np.std(period[(conditions==True)&\
                        (flag_total[win_ind]==True)])
                    tofill_st[j] = std
                    ## MGMG
                    #print 'window',len(period)/48
                    #print 'filled with two weeks'
                    qc_gf[j] =  2 
                    continue
            
            # if nothing is found under similar meteo within two weeks, just
            # look for global radiation within one week
            
        # METHOD TWO: just global radiation available  ############ 
         
        # check: wether global R is available
        if (Rg_flag[j] == 0) & TrueFilling:
            #print j,'radiation'
            rg_cond = Rg[j]
            j1 = [j-np.arange(1,week+1)+1] 
            j2 = [j+np.arange(1,week)]
            app_jays = np.append(j1,j2)
            win_ind = np.clip(app_jays,0,len(flag_met)-1)
            
            # for better overview: dynamic calculation of radiation threshold
            # minimum 20; maximum 50 [Wm-2] according to private correspondence
            # with Markus Reichstein
            rg_dev = max(20,min(rg_cond,Rg_dev))
            
            # check for the same radation conditions within the same window
            conditions = ((np.abs(rg_cond-Rg[win_ind])-rg_dev)<(-rg_dev*eps))
            number_4avg = len(conditions[conditions&flag_sec[win_ind]])
            
            if number_4avg>= 2:
                period = tofill[win_ind]
                avg = np.mean(period[(conditions==True) \
                        &(flag_sec[win_ind]==True)])    
                tofill_gf[j] = avg
                ## MGMG
                std = np.std(period[(conditions==True) \
                        &(flag_sec[win_ind]==True)])    
                tofill_st[j] = std
                ## MGMG
                qc_gf[j] = 3
                #print 'window',len(period)/48
                #print number_4avg 
                #print 'filled with radiation one week'
                continue
       
        for size in range(0,3):
            t_wind = t_int*(2*size+1)/2
            j1 = [j-np.arange(1,t_wind+1)+1] 
            j2 = [j+np.arange(1,t_wind)]
            appended_jays = np.append(j1,j2)
            win_ind = np.sort(np.clip(appended_jays,0,len(flag_met)-1))
            dec_time_window = dec_time[win_ind]
            period    = tofill[win_ind]
            conditions = np.abs(dec_time_window-dec_time[j])<1.1
            #print len(period[(conditions==True)&(period.flatten()!=udef)])
            ff = flag_tofill.flatten()[win_ind]
            if len(period[(conditions==True) & (ff==0)])>=2: 
                avg = np.mean(period[(ff == 0) & (conditions==True)])
                tofill_gf[j] = avg
                ## MGMG
                tofill_st[j] = np.std(period[(ff == 0) \
                        & (conditions==True)])
                ## MGMG
                qc_gf[j] = 4
                #print 'filled with 2.5 days'
                #print 'window',len(period)/48
                #print len(period[(conditions==True)&(period.flatten()!=udef)])
                break
    
        if tofill_gf[j]!= udef:
            continue
        
        # start a new cycle with increased window 
        
        if (flag_met[j] == 1):
            #print j,'all_met' 
            for multi in range(3,12):
                #print 'multi in second main iteration',multi
                # search for values around the met-conditions in a window of time
                # (one week in the first iteration and odd weeks in the next)
                j1 = [j-np.arange(1,multi*week+1)+1]
                j2 = [j+np.arange(1,multi*week)]
                comb_jays = np.append(j1,j2)
                win_ind = np.clip(comb_jays,0,len(flag_met)-1)
                           
                # for better overview: dynamic calculation of radiation threshold
                # minimum 20; maximum 50 [Wm-2] according to private correspondence
                # with Markus Reichstein
                rg_dev = max(20,min(rg_cond,Rg_dev))
    
                # get boolean array where meteo-conditions are in a given width
                conditions = ((np.abs(rg_cond-Rg[win_ind])-rg_dev)<(-rg_dev*eps))&\
                         ((np.abs(tair_cond-Tair[win_ind])-Tair_dev)\
                         <(-Tair_dev*eps))&\
                         ((np.abs(vpd_cond-vpd[win_ind])-vpd_dev)<(-vpd_dev*eps))
                
                if len(conditions[(conditions==True)&\
                        (flag_sec[win_ind]==True)])>= 2:
                    period = tofill[win_ind]
                    avg = np.mean(period[(conditions==True)&\
                            (flag_total[win_ind]==True)])
                    tofill_gf[j] = avg
                    ##MGMG
                    std = np.std(period[(conditions==True)&\
                            (flag_total[win_ind]==True)])
                    tofill_st[j] = std
                    ##MGMG
                    
                    qc_gf[j] = set_qual(multi,1)
                    #print 'multi filling with met_all and one iteration'                
                    break 
            
            if tofill_gf[j]!=udef: 
                continue
    
        for multi in range(2,12):
            # 2nd case: just global radiation available
            j1 = [j-np.arange(1,multi*week+1)+1]
            j2 = [j+np.arange(1,multi*week)]
            comb_jays = np.append(j1,j2)
            win_ind = np.clip(comb_jays,0,len(flag_met)-1)
            # check: wether global R is available
            if Rg_flag[j] == 0:
                #print j,'radiation'
                rg_cond = Rg[j]
                
                # get boolean array where meteo-conditions are in a given width
                conditions = ((np.abs(rg_cond-Rg[win_ind])-rg_dev)<(-rg_dev*eps))
                
                if len(conditions[(conditions==True)&\
                        (flag_sec[win_ind]==True)]) >= 2:
                    period = tofill[win_ind]
                    avg = np.mean(period[(conditions==True) & (flag_sec[win_ind]==True)])    
                    tofill_gf[j] = avg
                    ## MGMG
                    std = np.std(period[(conditions==True) & (flag_sec[win_ind]==True)])    
                    tofill_st[j] = std
                    ## MGMG
                    qc_gf[j] = set_qual(multi,2)
                    #print 'filled with radiation: '+str(multi)+'week'
                    break
         
        if tofill_gf[j]!=udef:
            continue
        
        #### 3. METHOD  #######################################################  
        
        # extend the window to two days
        
        for size in range(3,120):
            t_wind = t_int*(2*size+1)/2
            j1 = [j-np.arange(1,t_wind+1)+1] 
            j2 = [j+np.arange(1,t_wind)]
            appended_jays = np.append(j1,j2)
            win_ind = np.sort(np.clip(appended_jays,0,len(flag_met)-1))
            dec_time_window = dec_time[win_ind]
            period    = tofill[win_ind]
            conditions = np.abs(dec_time_window-dec_time[j])<1.1
            #print len(period[(conditions==True)&(period.flatten()!=udef)])
            # Flags are in a matrix array -> flat it so it can be compared
            ff = flag_tofill.flatten()[win_ind]
            if len(period[(conditions==True) & (ff==0)])>=2: 
                ind_4_avg = (period.flatten()!=udef)&(conditions==True) 
                tofill_gf[j] = np.mean(period[ind_4_avg])
                ## MGMG
                tofill_st[j] = np.std(period[ind_4_avg])
                ##MGMG
                qc_gf[j] = 3
                #print 'filled with '+str(len(period)/48)+' days'
                #print 'window',len(period)/48
                #print len(period[(conditions==True)&(period.flatten()!=udef)])
                break
    
        #print 'end' 
    
    if TrueFilling: 
        # Want to have a gap filled array with no udefs in it. So combine
        tofill_gf = np.where(flag_tofill.flatten()==0,tofill.flatten(),\
            tofill_gf.flatten())
    tofill_gf = np.reshape(tofill_gf,(len(tofill_gf),1))
    tofill_st = np.reshape(tofill_st,(len(tofill_st),1))
    
    return tofill_gf,qc_gf,tofill_st    
