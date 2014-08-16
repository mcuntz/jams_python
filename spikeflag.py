import numpy as np
from mad import mad

def spikeflag(date, data, inflag, isday, window=13, iter=1,
              fill_days=1, t_int=48, z=7, deriv=0, udef=-9999, spike_v=2,
              plot=False):
    '''
    Spike detection for Eddy Covariance data (and basically all other data)
    using a moving median absolute difference filter. Multiple iterations
    possible. Originally coded by Tino Rau.
    
    
    Definition
    ----------
    spikeflag(date, data, inflag, isday, window=13, iter=1,
              fill_days=1, t_int=48, z=5.5, deriv=0, udef=-9999, spike_v=2,
              plot=False):
    
    
    Input
    ----- 
    date        np.array(N), julian date (used only for plotting)
    data        np.array(N,M), data array where spike detection is applied on
                each column (M)
    inflag      np.array(N,M), dtype=int, quality flag of data, spike detection
                is only applied where inflag=0, all other data is ignored
    isday       np.array(N), dtype=bool, True where it is day and False where
                it is night
                  
                        
    Optional Input
    --------------
    window      int, size of the moving window where mad is calculated in days
                (default: 13)
    iter        int, how often the running window mad shall be applied
                (default: 1)
    fill_days   int, number of days where mad is applied within moving window
                (default: 1)
    t_int       int, number of data points within one day (default: 48)
    z           int/float, data is allowed to deviate maximum z standard
                deviations from the median (default: 7)
    deriv       int, 0: Act on raw data; 1: use first derivatives;
                2: use 2nd derivatives (default: 0)
    udef        int/float, missing value of data (default: -9999) NaN values are
                excluded from computations anyhow.
    spike_v     int, spike value which shall be returned when a spike is
                detected (default: 2)
    plot        bool, if True data and spikes are plotted (default: False)
    
    
    Output
    ------
    flag        np.array(N), flag array where everything is 0 except where
                spikes were detected, there it is spike_v.
    
    
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
    rows, cols = np.shape(data)
    flag       = np.zeros_like(inflag).astype(np.int)
    
    # mad window length and flag window length 
    period   = np.int(window*t_int)/2
    fill_win = np.int(fill_days*t_int)/2
    
    # calculate dusk and dawn times and separate in day and night
    isdawn      = np.zeros(rows,dtype=np.bool)
    isdusk      = np.zeros(rows,dtype=np.bool)
    dis         = isday - np.roll(isday,-1)
    isdawn[:-1] = np.where(dis[:-1] == -1, True, False)
    isdusk[:-1] = np.where(dis[:-1] == 1, True, False)
    isddday     = isdawn
    tmp         = np.roll(isdusk,1)
    isddday[1:] += tmp[1:]
    isddnight   = isdusk
    tmp         = np.roll(isdawn,1)
    isddnight[1:] += tmp[1:]
    
    # iterate over each column of data
    for col in xrange(cols):
        # iterate as much as iter
        for i in xrange(iter):
            # get day and night data
            day_data   = np.where((isday | isddday) & (inflag[:,col]==0) &
                                  (data[:,col]!=udef | ~np.isnan(data[:,col])),
                                  data[:,col], np.nan)
            night_data = np.where(((~isday) | isddnight) & (inflag[:,col]==0) &
                                  (data[:,col]!=udef | ~np.isnan(data[:,col])),
                                  data[:,col], np.nan)       

            # iterate over flag window
            fill_points = xrange(fill_win, isday.size-1, 2*fill_win)
            for j in fill_points:
                j1 = np.max([ j - period - 1,0])
                j2 = np.min([ j + period + 1,isday.size])
                fill_start = np.max([ j - fill_win,1])
                fill_end   = np.min([ j + fill_win,isday.size-1])
                
                day_flag = mad(np.ma.masked_array(data=day_data[j1:j2],\
                               mask=(np.isnan(day_data[j1:j2]))), z=z,\
                               deriv=deriv)
                
                flag[fill_start:fill_end,col] += np.where\
                        (day_flag[fill_start-j1-1:fill_end-j1-1] == 1,\
                         spike_v, 0)
                
                night_flag = mad(np.ma.masked_array(data=night_data[j1:j2],\
                                 mask=(np.isnan(night_data[j1:j2]))), z=z,\
                                 deriv=deriv)
                
                flag[fill_start:fill_end,col] += np.where\
                        (night_flag[fill_start-j1-1:fill_end-j1-1] == 1,\
                         spike_v, 0)
                
            if plot:
                import matplotlib.pyplot as plt
                fig1 = plt.figure(1)
                sub1 = fig1.add_subplot(111)
                valid = (inflag[:,col]==0) &\
                        (data[:,col]!=udef | ~np.isnan(data[:,col]))
                l1 =sub1.plot(date, data[valid,col], '-b')
                l2 =sub1.plot(date, data[flag!=0,col], 'or')
                plt.show()
    
    return flag

if __name__ == '__main__':
    import doctest
    doctest.testmod()