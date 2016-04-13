#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import ufz
from collections import OrderedDict

__all__ = ['read_data', 'write_data', 'write_data_dmp', 'write_data_one_file']

# --------------------------------------------------------------------

def read_data(files, undef=-9999., strip=None, norecord=False, nofill=False):
    """
        Read and concatenate data from CHS level1 data files.


        Definition
        ----------
        def read_data(files, undef=-9999., strip=None, norecord=False, nofill=False):


        Input
        -----
        files     (nfile,)-list with CHS data level1 file names


        Optional Input
        --------------
        undef     fill value for data and flags if non-existant (default: -9999.)
        strip     Strip strings with str.strip(strip) during read.
                  If None then strip quotes " and ' (default).
                  If False then no strip (30% faster).
                  Otherwise strip character given by strip.
        norecord  True: Do not assume that second column is record number.
        nofill    True: do not fill-up equal distant time steps

        
        Output
        ------
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead
            where
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        if norecord==False: record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        if norecord==False: hrecord   record header
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead)


        # Read data without record
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True, nofill=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead)

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

        Copyright 2015-2016 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
        Modified, MC, May 2015 - different variable in different input files
                               - strip
                               - norecord
                               - nofill
    """

    debug = False
    iundef = np.int(undef)

    # Get unique header and time stamps of all input file
    for cff, ff in enumerate(files):
        ihead      = ufz.fread(ff, skip=1, strip=strip, header=True)        # head with TIMESTAMP and RECORD
        if not norecord:
            if (ihead[1].split()[0]).lower() != 'record':
                raise ValueError('Assume the following structure: Date, Record, data, flag, data, flag, ...')
        idate      = ufz.sread(ff, skip=1, nc=1, squeeze=True, strip=strip) # TIMESTAMP
        idate      = ufz.ascii2eng(idate, full=True)
        if cff == 0:
            if norecord:
                hdat   = ihead[1::2]
                hflags = ihead[2::2]
            else:
                hdat   = ihead[2::2]
                hflags = ihead[3::2]
            adate  = idate
        else:
            if norecord:
                hdat   += ihead[1::2]
                hflags += ihead[2::2]
            else:
                hdat   += ihead[2::2]
                hflags += ihead[3::2]
            adate  += idate
    hdat   = list(OrderedDict.fromkeys(hdat))   # unique data head
    hflags = list(OrderedDict.fromkeys(hflags)) # unique flags head
    adate  = list(OrderedDict.fromkeys(adate))  # unique dates
    ii     = ufz.argsort(hdat)
    hdat   = [ hdat[i] for i in ii ]
    hflags = [ hflags[i] for i in ii ]
    adate.sort()

    # Fill missing time steps in all time steps
    if not nofill:
        date  = ufz.date2dec(eng=adate)
        dd    = np.round(np.diff(date)*24.*60.).astype(np.int) # minutes between time steps
        dmin  = np.amin(dd)                                    # time step in minutes
        dt    = np.float(dmin)/(24.*60.)                       # time step in fractional day
        igaps = np.where(dd != dmin)[0]                        # indexes of gaps
        for i in igaps[::-1]:
            nt      = np.round((date[i+1]-date[i])/dt).astype(np.int) # # of missing dates
            newdate = (date[i]+np.arange(1,nt)*dt)[::-1]              # the missing dates in reverse order
            for j in range(nt-1):
                date.insert(i+1, newdate[j]) # fill in missing dates, last one first
        adate = ufz.dec2date(date, eng=True)

    # Read files and fill in output array
    nrow   = len(adate)
    ncol   = len(hdat)
    dat    = np.ones((nrow,ncol))*undef # output array without
    flags  = np.ones((nrow,ncol), dtype=np.int)*iundef           # output array
    if not norecord: record = np.ones(nrow, dtype=np.int)*iundef # output array
    iidate = list()                     # list with indices of dates in dat/flag/record arrays
    iihead = list()                     # list with indices of header in hdat/hflags lists
    for cff, ff in enumerate(files):
        if debug: print('File name: ', ff)
        # date, data
        isdat, ssdat = ufz.fsread(ff, skip=1, snc=[0], nc=-1, strip=strip, strarr=True) # array
        isdate = ssdat[:,0]
        if norecord:
            idat    = isdat[:,0::2]
            iflags  = isdat[:,1::2].astype(np.int)
        else:
            irecord = isdat[:,0].astype(np.int)
            idat    = isdat[:,1::2]
            iflags  = isdat[:,2::2].astype(np.int)
        isdate  = ufz.ascii2eng(isdate, full=True)
        # date header, data header
        ihead, shead = ufz.fsread(ff, skip=1, snc=[0], nc=-1, strip=strip, header=True) # list
        ihdate   = shead[0]
        if norecord:
            ihdat    = ihead[0::2]
            ihflags  = ihead[1::2]
        else:
            ihrecord = ihead[0]
            ihdat    = ihead[1::2]
            ihflags  = ihead[2::2]
        # date and record header
        if cff == 0:
            hdate   = ihdate
            if not norecord: hrecord = ihrecord
        else:
            if norecord:
                if (ihdate != hdate):
                    raise ValueError('Assume the same date headers.')
            else:
                if (ihdate != hdate) or (ihrecord != hrecord):
                    raise ValueError('Assume the same date and record headers.')

        # fill output arrays and index lists
        iiidate = np.where(np.in1d(adate,isdate))[0]  # indexes in dat/flags of time steps in current file
        if debug:
            print('Numbers should match: ', iiidate.size, idat.shape[0])
            print('Current time steps:  ', len(isdate), [ aa for aa in isdate ])
            print('Selected time steps: ', len(iiidate), [ adate[aa] for aa in iiidate ])
        iidate.append(iiidate)                       
        if not norecord: record[iiidate] = irecord[:] # write at appropriate places in record
        iiihead = list()
        for i, h in enumerate(ihdat):
            hh = hdat.index(h)
            iiihead.extend([hh])
            dat[iiidate, hh]   = idat[:,i]            # write at appropriate places in dat
            flags[iiidate, hh] = iflags[:,i]          # write at appropriate places in flags
        iihead.append(iiihead)
    
    if norecord:
        return adate, dat, flags, iidate, hdate, hdat, hflags, iihead
    else:
        return adate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead

# --------------------------------------------------------------------

def write_data(*args):
    """
        Write concatenated data back to individual CHS level1 data files.

        Wrapper to write_data_norecord and write_data_record.


        Definition
        ----------
        def write_data(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead):
        or
        def write_data(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        [record   (n,)-array of record number]
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        [hrecord  record header]
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead)

        
        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data(files, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead)

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
        Written,  MC, May 2015
    """
    if len(args) == 9:
        write_data_norecord(*args)
    elif len(args) == 11:
        write_data_record(*args)
    else:
        raise ValueError('Must have 9 or 11 arguments.')

# --------------------------------------------------------------------

def write_data_dmp(*args):
    """
        Write data to individual Tereno Level2b files.

        Wrapper to write_data_norecord_dmp.


        Definition
        ----------
        def write_data_dmp(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead, hdmp):
        or
        def write_data_dmp(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead, hdmp):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        [record   (n,)-array of record number]
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        [hrecord  record header]
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files
        hdmp      data headers in Data Management Portal (DMP)

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofiles = [ f.replace('level2','level2b') for f in files ]
        hdmp = ufz.level1.get_value_excel(chsxlsfile, sheet, hdat, 'headerout (DB)')
        ufz.level1.write_data_dmp(ofiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead, hdmp)

        
        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofiles = [ f.replace('level2','level2b') for f in files ]
        hdmp = ufz.level1.get_value_excel(chsxlsfile, sheet, hdat, 'headerout (DB)')
        ufz.level1.write_data_dmp(ofiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead, hdmp)

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
        Written,  MC, May 2015
    """
    if len(args) == 10:
        write_data_norecord_dmp(*args)
    elif len(args) == 12:
        write_data_norecord_dmp(args[0],args[1],args[3],args[4],args[5],args[6],args[8],args[9],args[10],args[11])
    else:
        raise ValueError('Must have 10 or 12 arguments.')

# --------------------------------------------------------------------

def write_data_one_file(*args):
    """
        Write concatenated data back to one CHS data file.

        Wrapper to write_data_norecord_one_file and write_data_record_one_file.


        Definition
        ----------
        def write_data_one_file(infiles, sdate, record, dat, flags, hdate, hrecord, hdat, hflags):
        or
        def write_data_one_file(infiles, sdate, dat, flags, hdate, hdat, hflags):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        [record   (n,)-array of record number]
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        hdate     date/time header
        [hrecord  record header]
        hdat      data headers
        hflags    flags headers

        
        Output
        ------
        file will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofile = files[0].replace('level1', 'level2')        
        ufz.level1.write_data_one_file(ofile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags)

        
        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofile = files[0].replace('level1', 'level2')        
        ufz.level1.write_data_one_file(ofile, sdate, dat, flags, hdate, hdat, hflags)

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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, May 2016
    """
    if len(args) == 7:
        write_data_norecord_one_file(*args)
    elif len(args) == 9:
        write_data_record_one_file(*args)
    else:
        raise ValueError('Must have 7 or 9 arguments.')

# --------------------------------------------------------------------

def write_data_norecord(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead):
    """
        Write concatenated data back to individual CHS level1 data files without record column.


        Definition
        ----------
        def write_data_norecord(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data_norecord(files, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead)


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
        Written,  MC, Mar 2015 - from write_data
    """
    # Assure iterable infiles
    if not isinstance(infiles, (list, tuple, np.ndarray)): infiles = [infiles]

    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(infiles)   == len(iidate),  'File list and date index list do not conform.'
    assert len(infiles)   == len(iihead),  'File list and header index list do not conform.'
    assert len(sdate)     == ntime,        'Not enough dates.'
    assert flags.shape[0] == ntime,        'Not enough flag time steps.'
    assert flags.shape[1] == ncol,         'Not enough flag columns.'
    assert len(hdat)      == ncol,         'Not enough data headers.'
    assert len(hflags)    == ncol,         'Not enough flag headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ufz.ascii2eng(sdate, full=True)

    # Write individual files
    for cff, ff in enumerate(infiles):
        f = open(ff, 'wb')
        # header
        hstr = hdate
        ihead = iihead[cff]
        for i in ihead:
            hstr += ','+hdat[i]+','+hflags[i]
        print(hstr, file=f)
        # data
        idate = iidate[cff]
        for j in idate:
            dstr = isdate[j]
            for i in ihead:
                # dstr += ','+str(dat[j,i])+',-9999'   # for test: restate -9999 for all flags:
                dstr += ','+str(dat[j,i])+','+str(flags[j,i])
            print(dstr, file=f)
        f.close()

# --------------------------------------------------------------------

def write_data_norecord_dmp(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead, hdmp):
    """
        Write data to individual Tereno Level2b files.


        Definition
        ----------
        def write_data_norecord_dmp(infiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead, hdmp):


        Input
        -----
        infiles   (nfile,)-list with CHS data level2b file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files
        hdmp      data headers in Data Management Portal (DMP)

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofiles = [ f.replace('level2','level2b') for f in files ]
        hdmp = ufz.level1.get_value_excel(chsxlsfile, sheet, hdat, 'headerout(DB)')
        ufz.level1.write_data_norecord_dmp(ofiles, sdate, dat, flags, iidate, hdate, hdat, hflags, iihead, hdmp)


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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2016 - from write_data_norecord
    """
    # Assure iterable infiles
    if not isinstance(infiles, (list, tuple, np.ndarray)): infiles = [infiles]

    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(infiles)   == len(iidate),  'File list and date index list do not conform.'
    assert len(infiles)   == len(iihead),  'File list and header index list do not conform.'
    assert len(sdate)     == ntime,        'Not enough dates.'
    assert flags.shape[0] == ntime,        'Not enough flag time steps.'
    assert flags.shape[1] == ncol,         'Not enough flag columns.'
    assert len(hdat)      == ncol,         'Not enough data headers.'
    assert len(hflags)    == ncol,         'Not enough flag headers.'
    assert len(hdmp)      == ncol,         'Not enough database headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ufz.ascii2eng(sdate, full=True)
        
    # Tereno flags: OK, DOUBTFUL or BAD
    if np.any(flags > 2): # level1 flags
        oflags = np.maximum(ufz.level1.get_maxflag(flags), 0)
    else:                 # level2 flags
        oflags = flags
    strflags = np.zeros(oflags.shape, dtype='S'+str(len('DOUBTFUL,Other,NIL')))
    ii = np.where(oflags == 0)
    if ii[0].size>0: strflags[ii] = 'OK,NIL,NIL'
    ii = np.where(oflags == 1)
    if ii[0].size>0: strflags[ii] = 'DOUBTFUL,Other,Other'
    ii = np.where(oflags == 2)
    if ii[0].size>0: strflags[ii] = 'BAD,Other,Other'

    # Write individual files
    for cff, ff in enumerate(infiles):
        f = open(ff, 'wb')
        # header
        hstr = 'timestamp,sensorname,value,quality_flag,quality_cause,quality_comment'
        print(hstr, file=f)
        # data
        ihead = iihead[cff]
        idate = iidate[cff]
        for j in idate:
            for i in ihead:
                dstr = isdate[j]+','+hdmp[i]+','+str(dat[j,i])+','+strflags[j,i]
                print(dstr, file=f)
        f.close()

# --------------------------------------------------------------------

def write_data_norecord_one_file(infile, sdate, dat, flags, hdate, hdat, hflags):
    """
        Write concatenated data back to one CHS data file without record column.


        Definition
        ----------
        def write_data_norecord_one_file(infile, sdate, dat, flags, hdate, hdat, hflags):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        hdate     date/time header
        hdat      data headers
        hflags    flags headers

        
        Output
        ------
        file will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, dat, flags, iidate, hdate, hdat, hflags, iihead = ufz.level1.read_data(files, norecord=True)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ofile = files[0].replace('level1', 'level2')
        ufz.level1.write_data_norecord_one_file(ofile, sdate, dat, flags, hdate, hdat, hflags)


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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2016 - from write_data_norecord
    """
    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(sdate)     == ntime,        'Not enough dates.'
    assert flags.shape[0] == ntime,        'Not enough flag time steps.'
    assert flags.shape[1] == ncol,         'Not enough flag columns.'
    assert len(hdat)      == ncol,         'Not enough data headers.'
    assert len(hflags)    == ncol,         'Not enough flag headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ufz.ascii2eng(sdate, full=True)

    # Write individual files
    f = open(infile, 'wb')
    # header
    hstr = hdate
    for i in range(ncol):
        hstr += ','+hdat[i]+','+hflags[i]
    print(hstr, file=f)
    # data
    for j in range(nrow):
        dstr = isdate[j]
        for i in range(ncol):
            dstr += ','+str(dat[j,i])+','+str(flags[j,i])
        print(dstr, file=f)
    f.close()

# --------------------------------------------------------------------

def write_data_record(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead):
    """
        Write concatenated data back to individual CHS level1 data files with record number.


        Definition
        ----------
        def write_data_record(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        iidate    (nfile,)-list with indices in the output arrays of the input files
        hdate     date/time header
        hrecord   record header
        hdat      data headers
        hflags    flags headers
        iihead    (nfile,)-list with indices in the output array of headers in the input files

        
        Output
        ------
        files will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data
        ufz.level1.write_data_record(files, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead)


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
        Modified, MC, May 2015 - different variable in different input files
                               - renamed to write_data_record
    """
    # Assure iterable infiles
    if not isinstance(infiles, (list, tuple, np.ndarray)): infiles = [infiles]

    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(infiles)   == len(iidate),  'File list and date index list do not conform.'
    assert len(infiles)   == len(iihead),  'File list and header index list do not conform.'
    assert len(sdate)     == ntime,        'Not enough dates.'
    assert len(record)    == ntime,        'Not enough record numbers.'
    assert flags.shape[0] == ntime,        'Not enough flag time steps.'
    assert flags.shape[1] == ncol,         'Not enough flag columns.'
    assert len(hdat)      == ncol,         'Not enough data headers.'
    assert len(hflags)    == ncol,         'Not enough flag headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ufz.ascii2eng(sdate, full=True)

    # Write individual files
    for cff, ff in enumerate(infiles):
        f = open(ff, 'wb')
        # header
        hstr = hdate+','+hrecord
        ihead = iihead[cff]
        for i in ihead:
            hstr += ','+hdat[i]+','+hflags[i]
        print(hstr, file=f)
        # data
        idate = iidate[cff]
        for j in idate:
            dstr = isdate[j]+','+str(record[j])
            for i in ihead:
                # dstr += ','+str(dat[j,i])+',-9999'   # for test: restate -9999 for all flags:
                dstr += ','+str(dat[j,i])+','+str(flags[j,i])
            print(dstr, file=f)
        f.close()

# --------------------------------------------------------------------

def write_data_record_one_file(infile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags):
    """
        Write concatenated data back to one CHS data file with record number.


        Definition
        ----------
        def write_data_record_one_file(infile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags):


        Input
        -----
        infiles   (nfile,)-list with CHS data level1 file names
        sdate     (n,)-array of ascii dates in format YYYY-MM-DD hh:mm:ss
        record    (n,)-array of record number
        dat       (n,m)-array of data
        flags     (n,m)-array of flags
        hdate     date/time header
        hrecord   record header
        hdat      data headers
        hflags    flags headers

        
        Output
        ------
        file will be overwritten


        Examples
        --------
        --> see __init__.py for full example of workflow

        # Read data
        files = ufz.files_from_gui(title='Choose Level 1 file(s)')
        sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(files)

        # Set flags if variables were not treated yet
        flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

        # Write back data to one file
        ofile = files[0].replace('level1', 'level2')
        ufz.level1.write_data_record_one_file(ofile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags)


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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2016 - from write_data_record
    """
    # Few checks of sizes
    ntime = dat.shape[0]
    ncol  = dat.shape[1]
    assert len(sdate)     == ntime,        'Not enough dates.'
    assert len(record)    == ntime,        'Not enough record numbers.'
    assert flags.shape[0] == ntime,        'Not enough flag time steps.'
    assert flags.shape[1] == ncol,         'Not enough flag columns.'
    assert len(hdat)      == ncol,         'Not enough data headers.'
    assert len(hflags)    == ncol,         'Not enough flag headers.'

    # assure YYYY-MM-DD hh:mm:ss format even if sdate has DD.MM.YYYY hh:m:ss format
    isdate = ufz.ascii2eng(sdate, full=True)

    # Write individual files
    f = open(infile, 'wb')
    # header
    hstr = hdate+','+hrecord
    for i in range(ncol):
        hstr += ','+hdat[i]+','+hflags[i]
    print(hstr, file=f)
    # data
    for j in range(ntime):
        dstr = isdate[j]+','+str(record[j])
        for i in range(ncol):
            dstr += ','+str(dat[j,i])+','+str(flags[j,i])
        print(dstr, file=f)
    f.close()

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
