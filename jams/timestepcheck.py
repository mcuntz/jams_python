#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import os as os
from jams.date2dec import date2dec
from jams.dec2date import dec2date
import re

def timestepcheck(indir, pat, outfile, begin, end, numhead=1, timeint=30,
                  format='ascii', empty='-9999', delimiter=',', skiprows=0):
    '''
    Checks ascii data files with first column being a time stamp in ascii
    of eng style for correctness. If for a given time interval time steps
    are missing, they will will be filled with the correct time stamp and
    empty values. A beginning and end for the output files can be set to which
    data should be kept or filled. Pat defines the file name pattern to consider
    in the check. Multiple files, which match the pat are concatenated. 
    
    
    Definition
    ----------
    timestepcheck(indir, pat, outfile, begin, end, numhead=1, timeint=30,
                  format='ascii', empty='-9999', delimiter=',', skiprows=0):

    
    Input
    ----- 
    indir       str, path of the folder where the input files are
    pat         str, name or regular expression of the input files
    outfile     str, path and name of the output file 
    begin       str, start time of the output file, must be in the same format
                as time stamps in the input files 
    end         str, end time of the output file, must be in the same format
                as time stamps in the input files  
        
    
    Optional Input
    --------------
    numhead     int, number of header lines in the input files (default: 1)
    timeint     int, time interval of the input file in minutes (default: 30)
    format      str, format of time stamps in input files. 'ascii' or 'eng' is
                possible (default: 'ascii')
    empty       str, value for missing values (default: '-9999')
    delimiter   str, delimiter of the input files (default: ',')
    skiprows    int, rows to skip in input files, e.g. logger fuzzle before
                actual data header starts (default: 0)
                
    
    Output
    ------
    outfile     file with missing time steps filled with empty values cut from
                begin to end 
    
    
    Restrictions
    ------------
    TODO: tested thoroughly only for timeint=30
    TODO: more bad value checks can be included, see sternchen and tuedelchen
    
    
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
    Written,  AP, Aug 2014
    '''            
    ###########################################################################
    # time interval list
    interval = range(0,60,timeint)
    #jdint    = date2dec(yr=-4712,mo=1,dy=1,hr=12,mi=timeint)
    #jdmin    = date2dec(yr=-4712,mo=1,dy=1,hr=12,mi=0,sc=30)# (=precision)
    jdint    = date2dec(yr=1,mo=1,dy=1,hr=12,mi=timeint) % 1
    jdmin    = date2dec(yr=1,mo=1,dy=1,hr=12,mi=0,sc=30) % 1# (=precision)
    if format == 'ascii':
        jdbegin  = date2dec(ascii = np.array([begin]))
        jdend    = date2dec(ascii = np.array([end]))
    elif format == 'eng':
        jdbegin  = date2dec(eng = np.array([begin]))
        jdend    = date2dec(eng = np.array([end]))
    
    ###########################################################################
    # reading input directory
    pat = re.compile(pat)
    new = True
    
    filelist = os.listdir(indir)
    for file in filelist:
        if re.search(pat, file):
            if new:
                data = np.loadtxt('./%s/%s'%(indir, file), dtype='|S21',\
                                  delimiter=delimiter, skiprows=skiprows)
                if data.shape[0] == 0:
                    print('Warning: File %s is empty!' %(file))
                else:
                    if np.shape(data.shape)[0] == 1:
                        data = data.reshape((1,-1))
                    new = False
            else:
                add_data = np.loadtxt('./%s/%s'%(indir, file), dtype='|S21',
                                      delimiter=delimiter,
                                      skiprows=numhead+skiprows)
                if add_data.shape[0] == 0:
                    print('Warning: File %s is empty!' %(file))
                elif np.shape(add_data.shape)[0] == 1:
                    add_data = add_data.reshape((1,-1))
                    data = np.append(data, add_data, 0)
                else:
                    data = np.append(data, add_data, 0)

    ###########################################################################
    # sternchen check :-D
    # replace with regular expression check
    data[data=='***********'] = empty #!!! uberprufen auf lange und re
    data[data=='********'] = empty #!!! uberprufen auf lange und re
    
    ###########################################################################
    # tuedelchen check :-D
    # replace with regular expression check
    if data[numhead,0][0] == '"':
        data[numhead:,0] = np.array([x[1:-1] for x in data[numhead:,0]])
    
    ###########################################################################
    # "NAN" check :-D
    # replace with regular expression check
    data[data=='"NAN"'] = empty #!!! uberprufen auf lange und re
    
    ###########################################################################
    # leerzeilen check
    blankline = np.where(data[0:2,0]=='')[0]
    data = np.delete(data, blankline, 0)

    ###########################################################################
    # missing values check
    data[data==''] = empty
    data[data=='""'] = empty

    columns = np.shape(data)[1]-1

    ###########################################################################
    # calculate julian date
    if format == 'ascii':
        import time
        jd = date2dec(ascii = data[numhead:,0])
    elif format == 'eng':
        jd = date2dec(eng = data[numhead:,0])

    ###########################################################################
    # wrong time stamp check
    diff = jd[1:] - jd[:-1]
    minute = np.array([x.split()[1][3:5] for x in data[numhead:,0]]).astype(int)
    nii    = np.nonzero(~np.in1d(minute, interval))[0]
    ts     = np.nonzero(np.less(diff, jdint-jdmin))[0]
    wrong  = np.unique(np.append(nii, ts))
    if data.shape[0]-numhead-2 in wrong:
        wrong = np.append(wrong, [data.shape[0]-numhead-1], 0)
    delete = []
    
    for i in wrong:
        print('\nHERE IS SOMETHING WRONG:\n')
        print('BOF' if numhead+i-2<0 else data[numhead+i-2,:4])
        print('BOF' if numhead+i-1<0 else data[numhead+i-1,:4])
        print('-----------------------------------------------')
        print(data[numhead+i,:4])
        print('-----------------------------------------------')
        print('EOF' if numhead+i+1>=np.shape(data)[0] else data[numhead+i+1,:4])
        print('EOF' if numhead+i+2>=np.shape(data)[0] else data[numhead+i+2,:4])
        
        do = raw_input("\n(d)elete entry, (s)et to empty, (t)ype in date, (i)gnore: ")
        
        if   do == 'd':
            delete += [numhead+i]
        elif do == 's':
            data[numhead+i,1:] = empty
        elif do == 't':
            newdate = str(raw_input("\nreplace with: "))
            data[numhead+i,0] = newdate
#            newmin = str(raw_input("\n%s"%(data[numhead+i,0][:-2])))
#            data[numhead+i,0] = data[numhead+i,0][:-2] + newmin
        elif do == 'i':
            pass
    
    data = np.delete(data, delete, 0)
        
    ###########################################################################
    # calculate julian date again
    if format == 'ascii':
        jd = date2dec(ascii = data[numhead:,0])
    elif format == 'eng':
        jd = date2dec(eng = data[numhead:,0])
    
    ###########################################################################
    # check time step
    diff = jd[1:] - jd[:-1]
    ingap = np.where(np.greater(diff, jdint+jdmin))[0]
    nugap = np.rint((diff[ingap]/jdint)-1)
    
    ###########################################################################
    # insert missing time steps
    for i in range(np.size(ingap))[::-1]:
        where = np.ones(nugap[i], dtype=int)*(ingap[i]+1+numhead)
        if format == 'ascii':
            span  = np.arange(1,nugap[i]+1)*jdint + jd[ingap[i]]
            what  = dec2date(span.astype('|S16').astype(float), ascii=True)
        elif format == 'eng':
            span  = np.arange(1,nugap[i]+1)*jdint + jd[ingap[i]]
            what  = dec2date(span.astype('|S16').astype(float), eng=True)
        what    = np.array([x[:-3] for x in what])
        
        miss    = np.empty((nugap[i],columns), dtype='|S11')
        miss[:] = empty
        what    = np.append(np.reshape(what, (-1,1)), miss, 1)
        data    = np.insert(data, where, what, 0)
    
    ###########################################################################
    # fill/cut up/off beginning and end
    start = np.where(data[:,0]==begin)[0]
    if start == numhead:
        pass
    elif start > numhead:
        data = np.delete(data, np.arange(numhead, start), 0)
    else:
        if format == 'ascii':
            tofill = int((date2dec(ascii = data[numhead,0]) - jdbegin)/jdint)
            span   = np.arange(0,tofill)*jdint + jdbegin # tofill+1 weggenommen!!!
            what   = dec2date(span.astype('|S16').astype(float), ascii=True)
        elif format == 'eng':
            tofill = int((date2dec(eng = data[numhead,0]) - jdbegin)/jdint)-1
            span   = np.arange(0,tofill)*jdint + jdbegin # tofill+1 weggenommen!!!
            what   = dec2date(span.astype('|S16').astype(float), eng=True)

        what  = np.array([x[:-3] for x in what])
        where = np.ones(tofill, dtype=int)*numhead # tofill+1 weggenommen!!!

        #if doidate:
        #    miss    = np.empty((tofill,columns-2), dtype='|S11') # tofill+1 weggenommen!!!
        #else:
        miss    = np.empty((tofill,columns), dtype='|S11') # tofill+1 weggenommen!!!
        miss[:] = empty
        what    = np.append(np.reshape(what, (-1,1)), miss, 1)
        data    = np.insert(data, where, what, 0)
        
    stop = np.where(data[:,0]==end)[0]
    maxind = np.shape(data)[0]-1
    if stop == maxind:
        pass
    elif stop < maxind:
        data = data[:stop+1,:]
    else:
        if format == 'ascii':
            tofill = int((jdend - date2dec(ascii = data[-1,0]))/jdint)
            span   = np.arange(1,tofill+1)*jdint + date2dec(ascii = data[-1,0])
            what   = dec2date(span.astype('|S16').astype(float), ascii=True)
        elif format == 'eng':
            tofill = int((jdend - date2dec(eng = data[-1,0]))/jdint)
            span   = np.arange(1,tofill+1)*jdint + date2dec(eng = data[-1,0])
            what   = dec2date(span.astype('|S16').astype(float), eng=True)

        what  = np.array([x[:-3] for x in what])
        #if doidate:
        #    miss    = np.empty((tofill,columns-2), dtype='|S11')
        #else:
        miss    = np.empty((tofill,columns), dtype='|S11')
        miss[:] = empty
        what    = np.append(np.reshape(what, (-1,1)), miss, 1)
        data    = np.append(data, what, 0)
        
    ###########################################################################
    # save data to txt file
    np.savetxt('%s'%outfile, data, fmt='%s', delimiter=delimiter)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
