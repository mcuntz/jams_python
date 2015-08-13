#!/usr/bin/env python
"""
    Helper functions to deal with CHS level1 data and flags.


    Provided functions
    ------------------
    get_flag         Get the flags at position n from CHS data flag vector.
    get_manual_flags Get start and end dates as well as flag values for a specific variable from a manual flag file.
    get_maxflag      Get the maximal flag of the string with the individual flags.
    get_value_excel  Get value in column of sheet in excelfile given variable name.
    read_data        Read and concatenate data from CHS level1 data files.
    set_flag         Set the flags at position n to iflag at indices ii of CHS data flag vector.
    write_data       Write concatenated data back to individual CHS level1 data files.


    Example
    -------
    # get files
    infiles = ufz.files_from_gui(title='Choose Level 1 file(s)')

    # read files
    sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = ufz.level1.read_data(infiles)

    # Get variables and indexes of variables in dat and flags
    myvars = [ v for v in hdat if v.startswith('stemT') ]
    idx    = [ hdat.index(v) for v in myvars ]

    # Set flags if variables were not treated yet
    flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

    # 1st test - check manual flags
    itest = 1
    jisdate = np.array(ufz.date2dec(eng=isdate))
    for v in myvars:
        msdate, medate, mflags = get_manual_flags(manual_flag_file, v)
        if len(msdate) > 0:
            for dd in range(len(msdate)):
                ii = np.where((jsdate>=msdate[dd]) & (jsdate<=medate[dd]))[0]
                if ii.size>0: flags[:,i] = ufz.level1.set_flag(flags[:,i], itest, mflags[dd], ii)

    # 2nd test - set first flag after the initial 9 to 2 if dat is undef
    itest  = 2
    isflag = 2
    for i in idx:
        ii = np.where(dat[:,i]==undef)[0]
        if ii.size>0: flags[:,i] = ufz.level1.set_flag(flags[:,i], itest, isflag, ii)


    # 3rd test - set second flag to 2 if dat is not in [min,max]
    #            treat only data that had no flag==2 before
    itest = 3
    isflag = 2
    for i, v in zip(idx, myvars):
        mini = ufz.level1.get_value_excel('CHS-measurements.xlsx', 'Forest Hohes Holz', v, 'Min')
        maxi = ufz.level1.get_value_excel('CHS-measurements.xlsx', 'Forest Hohes Holz', v, 'Max')
        maxflags = ufz.level1.get_maxflag(flags[:,i])
        ii = np.where((maxflags < 2) & ((dat[:,i] < swdrmin) | (dat[:,i] > swdrmax)))[0]
        if ii.size>0: flags[:,i] = ufz.level1.set_flag(flags[:,i], itest, isflag, ii)

    # write back data and flags
    ufz.level1.write_data(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead)

    # plot data without any flag
    i = 0
    xx = ufz.date2dec(eng=sdate)
    yy = np.ma.array(dat[:,i], mask=(ufz.level1.get_maxflag(flags[:,i])==2))
    mark1 = sub.plot(xx, yy)


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
    Modified  JM, May 2015 - get_maxflag
              MC, May 2015 - excel - get_value_excel
              MC, Aug 2015 - get_manual_flags
"""
from .excel          import get_value_excel
from .getset_flag    import get_flag, set_flag, get_maxflag
from .manual_flags   import get_manual_flags
from .readwrite_data import read_data, write_data

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.2'
__revision__ = "Revision: 2227"
__date__     = 'Date: 12.08.2015'
