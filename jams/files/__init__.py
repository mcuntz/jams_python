#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    File list functions.


    Provided functions
    ------------------
    fullnames                   Filenames with absolute paths in local directories.
    fullnames_dates             Filenames with absolute paths and modification times in local directories.
    fullnames_dates_sizes       Filenames with absolute paths, modification times, and file sizes in local directories.
    fullnames_sizes             Filenames with absolute paths and file sizes in local directories.
    fullnames_times             Wrapper for fullnames_dates.
    fullnames_times_sizes       Wrapper for fullnames_dates_sizes.
    last_fullname               Filename with absolute paths of last file in local directories.
    last_fullname_date          Filename with absolute paths and modification time of last file in local directories.
    last_fullname_date_size     Filename with absolute paths, modification time, and
                                file size of last file in local directories.
    last_fullname_size          Filename with absolute paths and file size of last file in local directories.
    last_fullname_time          Wrapper for last_fullname_date.
    last_fullname_time_size     Wrapper for last_fullname_date_size.
    last_name                   Filename of last file in local directories.
    last_name_date              Filename and modification time of last file in local directories.
    last_name_date_size         Filename, modification time, and file size of last file in local directories.
    last_name_size              Filename and file size of last file in local directories.
    last_name_time              Wrapper for last_name_date.
    last_name_date_size         Wrapper for last_name_date_size.
    names                       Filenames in local directories.
    names_dates                 Filenames and modification times in local directories.
    names_dates_sizes           Filenames, modification times, and file sizes in local directories.
    names_sizes                 Filenames and file sizes in local directories.
    names_times                 Wrapper for names_dates.
    names_times_sizes           Wrapper for names_dates_sizes.
    newest_fullname             Wrapper for last_fullname.
    newest_fullname_date        Wrapper for last_fullname_date.
    newest_fullname_date_size   Wrapper for last_fullname_date_size.
    newest_fullname_size        Wrapper for last_fullname_size.
    newest_fullname_time        Wrapper for last_fullname_date.
    newest_fullname_time_size   Wrapper for last_fullname_date_size.
    newest_name                 Wrapper for last_name.
    newest_name_date            Wrapper for last_name_date.
    newest_name_date_size       Wrapper for last_name_date_size.
    newest_name_size            Wrapper for last_name_size.
    newest_name_time            Wrapper for newest_name_date.
    newest_name_date_size       Wrapper for newest_name_date_size.


    Example
    -------
    import os
    # get .dat filenames with asolute paths in directories 2013, 2014, ...
    fls = fullnames('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

    # get .dat filenames and sizes in current directory
    os.chdir('change/directory')
    fls, flssize = names_sizes('*.dat')

    # get all filenames and modification times in current directory
    fls, flstimes = names_times()


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
"""
from .files import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1923"
__date__     = 'Date: 05.12.2014'
