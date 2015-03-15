#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.ascii2ascii import ascii2eng
from ufz.fsread import fsread

__all__ = ['get_flag', 'read_data', 'set_flag']

# --------------------------------------------------------------------

def get_flag(flags, n):
    """
        Get the flags at position n from CHS data flag vector.


        Definition
        ----------
        def get_flag(flags, n):


        Input
        -----
        flags   ND-array of integers
        n       position to extract

        
        Output
        ------
        ND-array of integers with flag at position n of flags
        returns -1 if input integer is less than 10**(n+1)


        Examples
        --------
        >>> flags = np.array([3, 30, 31, 300, 301, 3001, 3201, 312121212])
        >>> print(get_flag(flags, 0))
        [3 3 3 3 3 3 3 3]

        >>> print(get_flag(flags, 1))
        [-1  0  1  0  0  0  2  1]

        >>> print(get_flag(flags, 2))
        [-1 -1 -1  0  1  0  0  2]

        >>> print(get_flag(flags, 3))
        [-1 -1 -1 -1 -1  1  1  1]

        >>> print(get_flag(flags, 4))
        [-1 -1 -1 -1 -1 -1 -1  2]


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
    ilog10 = np.log10(flags).astype(np.int)         # how many digits
    out = np.ones(flags.shape, dtype=np.int)*(-1)   # if not enough digits -> -1
    ii = np.where((ilog10-n) >= 0)[0]               # otherwise extract n-th digit
    if ii.size > 0:
        out[ii] = (flags[ii] // 10**(ilog10[ii]-n)) % 10

    return out

# --------------------------------------------------------------------

def read_data(files):
    """
        Read and concatenate data from CHS level1 data files.


        Definition
        ----------
        def read_data(files):


        Input
        -----
        files   list with CHS data level1 file names

        
        Output
        ------
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags
            where
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    list of nfile+1 entries with start and end indeces in output arrays of the input files
        hdate     date/time header
        hrecord   record header
        hdat      data headers
        hflags    flags headers


        Examples
        --------
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags = ufz.level1.read_data(files)


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
            iidate = [0, sdate.size] # list with start and end indeces in total arrays
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

def set_flag(flags, n, iflag, ii=None):
    """
        Set the flags at position n to iflag at indeces ii of CHS data flag vector.


        Definition
        ----------
        def set_flag(flags, n, iflag, ii=None):


        Input
        -----
        flags   CHS data flag vector, ND-array of integers
        n       position in flag to set, missing positions will be cretaed and filled with 0
        iflag   set flag to this value


        Optional Input
        --------------
        ii      indices at which to set flag (default: None, i.e. each entry)

        
        Output
        ------
        ND-array of integers with flag at position n of flags set to iflag


        Examples
        --------
        >>> flags = np.array([3, 30, 301, 3101, 312121212])
        >>> print(set_flag(flags, 1, 2, [0,1,2]))
        [       32        32       321      3101 312121212]

        >>> print(set_flag(flags, 3, 2, [0,1,2,3]))
        [     3002      3002      3012      3102 312121212]

        >>> print(set_flag(flags, 1, 2))
        [       32        32       321      3201 322121212]


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
    fflags = flags.copy()
    # extend flag vector if needed
    ilog10 = np.log10(fflags).astype(np.int)  # are there enough flag positions (i)
    jj = np.where(ilog10 < n)[0]
    if jj.size > 0:                          # increase number (filled with 0)
        fflags[jj] *= 10**(n-ilog10[jj])
        ilog10 = np.log10(fflags).astype(np.int)
        
    # get entries in flag vector
    isflags = get_flag(fflags, n)

    # set entries in flag vector
    if ii is None:                                  # set all to iflag
        fflags += 10**(ilog10-n) * (iflag-isflags)
    else:
        if np.size(ii) > 0:                         # set ii to iflag
            fflags[ii] += 10**(ilog10[ii]-n) * (iflag-isflags[ii])

    return fflags

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
