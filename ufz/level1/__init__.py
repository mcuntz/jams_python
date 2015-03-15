#!/usr/bin/env python
"""
    Helper functions to deal with CHS level1 data and flags.


    Provided functions
    ------------------
    get_flag         Get the flags at position n from CHS data flag vector.
    read_data        Read and concatenate data from CHS level1 data files.
    set_flag         Set the flags at position n to iflag at indeces ii of CHS data flag vector.


    Example
    -------
    # get files
    files = ufz.files_from_gui(title='Choose Level 1 file(s)')
    # read files
    sdate, record, dat, flags, iidate, hdate, hrecord, hdat, hflags = read_data(files)
    # julian dates
    date  = ufz.date2dec(eng=sdate)
    # 1st test - set first flag after the initial 3 to 2 if dat is undef
    itest = 1
    for i in idx:
        ii = np.where(dat[:,i]==undef)[0]
        flags[:,i] = set_flag(flags[:,i], itest, 2, ii)
    # plot
    i = 0
    xx = date
    # original data without undef
    yy = np.ma.array(dat[:,i], mask=(get_flag(flags[:,i],1)==2))
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
"""
from .level1 import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 2051"
__date__     = 'Date: 15.03.2015'
