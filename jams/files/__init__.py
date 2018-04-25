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
    This file is part of the JAMS Python package.

    The JAMS Python package is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The JAMS Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2014 Matthias Cuntz


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
