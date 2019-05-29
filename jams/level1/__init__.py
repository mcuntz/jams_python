#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Helper functions to deal with CHS level1 data and flags.


    Provided functions
    ------------------
    constant_values  Checks if a given series of data contains consecutive values which are constant over a certain time period.
    get_flag              Get the flags at position n from CHS data flag vector.
    get_manual_flags      Get start and end dates as well as flag values for a specific variable from a manual flag file.
    get_maxflag           Get the maximal flag of the string with the individual flags.
    get_header_excel      Get the header row of an Excel sheet.
    get_value_excel       Get value in column of sheet in excelfile given variable name.
    read_data             Read and concatenate data from CHS level1 data files.
    set_flag              Set the flags at position n to iflag at indices ii of CHS data flag vector.
    write_data            Write concatenated data back to individual CHS level1 data files.
    write_data_one_file   Write concatenated data back to one CHS data file.
    spike                 Spike filter for rectangularly shaped spikes.


    Example
    -------
    # get files
    infiles = jams.files_from_gui(title='Choose Level 1 file(s)')

    # read files
    sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead = jams.level1.read_data(infiles)

    # Get variables and indexes of variables in dat and flags
    myvars = [ v for v in hdat if v.startswith('stemT') ]
    idx    = [ hdat.index(v) for v in myvars ]

    # Set flags if variables were not treated yet
    flags[:,idx] = np.where(flags[:,idx]==np.int(undef), 9, flags[:,idx])

    # 1st test - check manual flags
    itest = 1
    jisdate = np.array(jams.date2dec(eng=isdate))
    for v in myvars:
        msdate, medate, mflags = get_manual_flags(manual_flag_file, v)
        if len(msdate) > 0:
            for dd in range(len(msdate)):
                ii = np.where((jsdate>=msdate[dd]) & (jsdate<=medate[dd]))[0]
                if ii.size>0: flags[:,i] = jams.level1.set_flag(flags[:,i], itest, mflags[dd], ii)

    # 2nd test - set first flag after the initial 9 to 2 if dat is undef
    itest  = 2
    isflag = 2
    for i in idx:
        ii = np.where(dat[:,i]==undef)[0]
        if ii.size>0: flags[:,i] = jams.level1.set_flag(flags[:,i], itest, isflag, ii)


    # 3rd test - set second flag to 2 if dat is not in [min,max]
    #            treat only data that had no flag==2 before
    itest = 3
    isflag = 2
    for i, v in zip(idx, myvars):
        mini = jams.level1.get_value_excel('CHS-measurements.xlsx', 'Forest Hohes Holz', v, 'Min')
        maxi = jams.level1.get_value_excel('CHS-measurements.xlsx', 'Forest Hohes Holz', v, 'Max')
        maxflags = jams.level1.get_maxflag(flags[:,i])
        ii = np.where((maxflags < 2) & ((dat[:,i] < swdrmin) | (dat[:,i] > swdrmax)))[0]
        if ii.size>0: flags[:,i] = jams.level1.set_flag(flags[:,i], itest, isflag, ii)

    # write back data and flags
    jams.level1.write_data(infiles, sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags, iihead)

    # plot data without any flag
    i = 0
    xx = jams.date2dec(eng=sdate)
    yy = np.ma.array(dat[:,i], mask=(jams.level1.get_maxflag(flags[:,i])==2))
    mark1 = sub.plot(xx, yy)

    # Make level2
    flags = jams.maximum(jams.level1.get_maxflag(flags), 0)

    # write level2 data
    ofile = infiles[0].replace('level1', 'level2')
    jams.level1.write_data_onefile(ofile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags)

    # write level2b data for DMP
    ofile = infiles[0].replace('level1', 'level2b')
    hdmp = jams.level1.get_value_excel(chsxlsfile, sheet, hdat, 'headerout (DB)')
    jams.level1.write_data_dmp(ofile, sdate, record, dat, flags, hdate, hrecord, hdat, hflags, hdmp)


    License
    -------
    This file is part of the JAMS Python package.

    Copyright (c) 2015-2016 Matthias Cuntz - mc (at) macu (dot) de

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

    Copyright 2015 Matthias Cuntz


    History
    -------
    Written,  MC, Mar 2015
    Modified  JM, May 2015 - get_maxflag
              MC, May 2015 - excel - get_value_excel
              MC, Aug 2015 - get_manual_flags
              MC, Jan 2016 - get_header_excel
              MC, Mar 2016 - write_data_dmp, write_data_one_file
              BD, Sep 2016 - spike
"""
from .constant_values import constant_values
from .excel           import get_header_excel, get_value_excel
from .getset_flag     import get_flag, set_flag, get_maxflag
from .manual_flags    import get_manual_flags
from .spike    import spike
from .readwrite_data  import read_data, write_data, write_data_dmp, write_data_dmp_size, write_data_one_file

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.3'
__revision__ = "Revision: 2477"
__date__     = 'Date: 14.01.2016'
