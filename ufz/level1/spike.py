
#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def spike(datin, thresh=7, toler=0, length=15):
    """
        spike filter which looks for spikes which have a first step > thresh and then a 'plateau' for <= length time steps and come back
        to original value within a tolerance of toler. (something like a rectangular shape, but the 'plateau' does not have to be flat,
        it just needs to exceed the threshold without crossing the value before the spike).
        Returns list with indices of detected spikes.


        Definition
        ----------
        def spike(datin, thresh=7, toler=0, length=15):


        Input
        -----
        datin      1d array or list


        Optional Input
        --------------
        thresh     spike threshold
        toler      tolerance value for the abs difference between first value before spike and first value after spike
        length     maximum length of spike plateau


        Output
        ------
        list with indices of all spikes detected


        Restrictions
        ------------
        --
        

        Examples
        --------
        >>> y = [1, 1, 20,1, 2, 1,15,1, 10, 11, 13, 10, 1.5, 1, 2, 1, 1, -7, -8, 1.4, 3, 1, 1, 1, -50, 1, 1, 1]
        >>> spike(y,thresh=2,toler=0.6,length=10)
        no. of spikes: 5
        no. of data points: 9
        [2, 6, 8, 9, 10, 11, 17, 18, 24]


        >>> y = [0.1, 0.1, 0.2, 0.1, 0.2, 0.1, 1.5, 0.1, 1.0, 1.1, 1.3, 1.0, 0.15, 0.1, 0.2, 0.1, 0.1, -0.7, -0.8, 0.14, 0.3, 0.1, 0.1, 0.1, -5.0, 0.1, 0.1, 0.1]
        >>> spike(y,thresh=0.2,toler=0.06,length=3)
        no. of spikes: 3
        no. of data points: 4
        [6, 17, 18, 24]


        >>> y = [0.1, 0.1, 0.2, 0.1, 0.2, 0.1, 1.5, 0.2, 1.0, 1.1, 1.3, 1.0, 0.15, 0.1, 0.2, 0.1, 0.1, -0.7, 0.8, 0.14, 0.3, 0.1, 0.1, 0.1, -5.0, 0.1, 0.1, 0.1]
        >>> spike(y,thresh=0.2,toler=0.06,length=3)
        no. of spikes: 1
        no. of data points: 1
        [24]

        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2016 Benjamin Dechant


        History
        -------
        Written,  BD, Sep 2016
    """
    debug = False

    diff = [datin[k+1]-datin[k] for k in range(len(datin)-1)]   # differences
    ipos0 = [k for k,a in enumerate(diff) if abs(a) > thresh]      # get index position of spikes  
    if len(ipos0)>0: 
        ipos = ipos0[0]                                          # select first spike

    maxlen=length                                                  # set max. length of spike plateau
    spike_pos_all = []
    while ipos < len(datin)-1:
        if debug == True:
            print('ipos=',ipos)
        for i in range(1,maxlen+1):
            if debug == True:
                print('i=',i)
            spike_pos = []
            tm = [abs(datin[ipos]-datin[ipos+v]) for v in  list(range(2,i+2))]  #diff between first val before spike and all candidates for first val after spike
            tms = [datin[ipos]-datin[ipos+v] for v in  list(range(2,i+1))]  #diff between first val before spike and all candidates for first val after spike (with sign)
            if len(tm)==1 and tm[0] < toler:
                spike_pos = [ipos+1]     
                spike_pos_all.append(spike_pos) 
                if debug == True:
                    print('found peak 1',spike_pos)
                break
            elif len(tm)>1 and tm[0:-2] > thresh and [cmp(val,0) for val in tms] == [cmp(datin[ipos] - datin[ipos+1],0) for n in range(len(tms))] and tm[-1] < toler:   # check thresh, tolerance and no switching of sign
                spike_pos = list(range(ipos+1,ipos+i+1))
                spike_pos_all.append(spike_pos)
                if debug == True:
                    print('found peak 2',spike_pos)
                break
        ipos1 = [k for k,a in enumerate(diff) if abs(a) > thresh  and k >= ipos + i +1] # get index position of next spikes
        if len(ipos1)> 0:
            ipos = ipos1[0]
        else:
            break 

    print('no. of spikes: '+str(len(spike_pos_all)))
    spike_pos_all = [item for sublist in spike_pos_all for item in sublist]    # create list without sublists
    print('no. of data points: '+str(len(spike_pos_all)))
    return spike_pos_all



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
