#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Functions for interacting with an open FTP connection.


    Provided functions
    ------------------
    get_binary                  - Get a binary file from an ftp connection
    get_check_binary            - Get a binary file from an ftp connection
                                  True if all bytes transfered
                                  False if local bytes differ after transfer
    get_unix_ascii              - Get a Unix/Linux ascii file from an ftp connection
    get_check_unix_ascii        - Get a Unix/Linux ascii file from an ftp connection
                                  True if all bytes transfered
                                  False if local bytes differ after transfer
    get_windows_ascii           - Get a Windows ascii file from an ftp connection
    get_check_windows_ascii     - Get a Windows ascii file from an ftp connection
                                  True if all bytes transfered
                                  False if local bytes differ after transfer
    get_names                   - Get file names in current directory on ftp connection
    get_names_dates             - Get file names and dates in current directory on ftp connection
    get_names_sizes             - Get file names and sizes in current directory on ftp connection
    get_names_times             - Wrapper for get_names_dates
    get_names_dates_sizes       - Get file names, dates and sizes in current directory on ftp connection
    get_names_times_sizes       - Wrapper for get_names_dates_sizes
    get_size                    - Get the size of a file on an ftp connection
    get_sizes                   - Get file sizes in current directory on ftp connection
    set_mtime(ftp, fname)       - Set access and modification time of local file from date information on ftp connection
    put_binary                  - Put a binary file to an ftp connection
    put_check_binary            - Put a binary file to an ftp connection
                                  True if all bytes transfered
                                  False if remote bytes differ after transfer


    Example
    -------
    # Transfer files from FTP connection if they are larger than existing local files
    import ftplib
    import os
    # open ftp connection
    ftp = ftplib.FTP("ftp.server.de")
    ftp.login("user", "password")
    fisdir = ftp.pwd()
    # go to directory
    ftp.cwd(choose/directory)
    # get all dat files with filesizes
    flsdat, flsdatsize = get_names_sizes(ftp, '*.dat')
    # if you want to get more than one file ending, e.g. .zip and .tar.gz
    fls, flssize = get_names_sizes(ftp)
    flszip     = [ i for i in fls if i.endswith('.zip') ]
    flszipsize = [ flssize[ii] for ii, i in enumerate(fls) if i.endswith('.zip') ]
    flstar     = [ i for i in fls if i.endswith('.tar.gz') ]
    flstarsize = [ flssize[ii] for ii, i in enumerate(fls) if i.endswith('.tar.gz') ]
    # local dat file names and sizes
    llsdat     = os.listdir('*.dat')
    llsdatsize = [ int(os.stat(n).st_size) for n in llsdat ]
    ldict = dict(zip(llsdat, llsdatsize))
    fdict = dict(zip(flsdat, fldatsize))
    def local_lt_ftp(name, local_dict, ftp_dict):
        if name not in local_dict:
            return True
        else:
            if local_dict[name] < ftp_dict[name]:
                return True
            else:
                return False
    for i in flsdat:
        if local_lt_ftp(i, ldict, fdict):
            if not get_check_binary(ftp, i):
                if not get_check_binary(ftp, i):
                    print('Transfer of ', i, ' failed twice.')
                    if os.path.exists(i):
                        os.remove(i)
    ftp.quit()


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014-2016 Matthias Cuntz - mc (at) macu (dot) de

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


    History
    -------
    Written,  MC, Jun-Dec 2014
    Modified, MC, Jan 2016
"""
from .ftp import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.1'
__revision__ = "Revision: 2468"
__date__     = 'Date: 06.01.2016'
