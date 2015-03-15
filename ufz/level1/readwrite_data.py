#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.ascii2ascii import ascii2eng
from ufz.fsread import fsread

__all__ = ['read_data', 'write_flag']

# --------------------------------------------------------------------

def read_data(files):
    """
        Read and concatenate data from CHS level1 data files.


        Definition
        ----------
        def read_data(files):


        Input
        -----
        files     (nfile,)-list with CHS data level1 file names

        
        Output
        ------
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags
            where
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile+1,)-list with start indices of the input files in the output arrays
        hdate     date/time header
        hrecord   record header
        hdat      data headers
        hflags    flags headers


        Examples
        --------
        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 3, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags)


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
    """
    for cff, ff in enumerate(files):
        # date, data
        isdat, ssdat = fsread(ff, skip=1, snc=[0], nc=-1) # array
        isdate  = ssdat[:,0]
        irecord = isdat[:,0]
        idat    = isdat[:,1::2]
        iflags  = isdat[:,2::2].astype(np.int)

        # date header, data header
        ihead, shead = fsread(ff, skip=1, snc=[0], nc=-1, header=True) # list
        ihdate   = shead[0]
        ihrecord = ihead[0]
        ihdat    = ihead[1::2]
        ihflags  = ihead[2::2]

        # Concatenate arrays, check that headers are the same
        if cff == 0:
            check_hdat = ihdat # save 1st header for check of following headers 
            if isdate.size > 1:
                sdate   = isdate
                record  = irecord
                dat     = idat
                flags   = iflags
                hdate   = ihdate
                hrecord = ihrecord
                hdat    = ihdat
                hflags  = ihflags
            else: # assure array and header in case of only one input line
                sdate   = np.array([isdate])
                record  = np.array([irecord])
                dat     = np.array([idat])
                flags   = np.array([iflags])
                hdate   = list(ihdate)
                hrecord = list(ihrecord)
                hdat    = list(ihdat)
                hflags  = list(ihflags)
            iidate = [0, sdate.size] # list with start and end indices in total arrays
        else:
            # Check that the headers are the same
            if ihdat != check_hdat: raise ValueError('read_data: names in headers are not the same.')
            # append date and data
            sdate  = np.append(sdate,  isdate,  axis=0)
            record = np.append(record, irecord, axis=0)
            dat    = np.append(dat,    idat,    axis=0)
            flags  = np.append(flags,  iflags,  axis=0)
            iidate.append(sdate.size) # append start/end index list

    # assure YYYY-MM-DD hh:mm:ss format even if files had DD.MM.YYYY hh:m:ss format
    sdate = ascii2eng(sdate, full=True)
    
    return sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags

# --------------------------------------------------------------------

def write_data(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags):
    """
        Write concatenated data back to individual CHS level1 data files.


        Definition
        ----------
        def write_data(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile+1,)-list with start indices of the input files in the output arrays
        hdate     date/time header
        hrecord   record header
        hdat      data headers
        hflags    flags headers

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 3, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags)


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
    """
    # Assure iterable infiles
    if not isinstance(infiles, (list, tuple, np.ndarray)): infiles = [infiles]

    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(infiles)   == len(iidate)-1, 'File list and index list do not conform.'
    assert sdate.size     == ntime,         'Not enough dates.'
    assert record.size    == ntime,         'Not enough record numbers.'
    assert flags.shape[0] == ntime,         'Not enough flag time steps.'
    assert flags.shape[1] == ncol,          'Not enough flag columns.'
    assert len(hdat)      == ncol,          'Not enough data headers.'
    assert len(hflags)    == ncol,          'Not enough flag headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ascii2eng(sdate, full=True)

    # Write individual files
    for ff in range(len(infiles)):
        f = open(infiles[ff], 'wb')
        # header
        hstr = hdate+','+hrecord
        for i in range(len(hdat)):
            hstr += ','+hdat[i]+','+hflags[i]
        print(hstr, file=f)
        # data
        for j in range(iidate[ff], iidate[ff+1]):
            dstr = isdate[j]+','+str(record[j])
            for i in range(len(hdat)):
                #For test: restate -9999 for all flags: dstr += ','+str(dat[j,i])+',-9999'
                dstr += ','+str(dat[j,i])+','+str(flags[j,i])
            print(dstr, file=f)
        f.close()

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
