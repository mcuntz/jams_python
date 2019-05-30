#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import os
import datetime
import glob
from jams.argsort import argsort
import numpy as np

__all__ = ['fullnames', 'fullnames_dates', 'fullnames_dates_sizes', 'fullnames_sizes',
           'fullnames_times', 'fullnames_times_sizes',
           'last_fullname', 'last_fullname_date', 'last_fullname_date_size',
           'last_fullname_size', 'last_fullname_time', 'last_fullname_time_size',
           'last_name', 'last_name_date', 'last_name_date_size', 'last_name_size',
           'last_name_time', 'last_name_time_size', 
           'names', 'names_dates', 'names_dates_sizes', 'names_sizes', 'names_times',
           'names_times_sizes',
           'newest_fullname', 'newest_fullname_date', 'newest_fullname_date_size',
           'newest_fullname_size', 'newest_fullname_time', 'newest_fullname_time_size',
           'newest_name', 'newest_name_date', 'newest_name_date_size',
           'newest_name_size', 'newest_name_time', 'newest_name_date_size']

# --------------------------------------------------------------------

def fullnames(fname=None, dirs=None):
    """
        Filenames with absolute paths in local directories.


        Definition
        ----------
        def fullnames(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames incl. absolute paths.


        Examples
        --------
        import os
        # get .dat filenames with absolute paths in directories 2013, 2014, ...
        fls = fullnames('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames in current directory
        os.chdir('change/directory')
        fls = fullnames('*.dat')

        # get all filenames in current directory
        fls = fullnames()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        fnames = [ os.path.abspath(f) for f in fnames ]
        lls  += fnames
        if i != '.':
            os.chdir(curdir)
    lls.sort()
    return lls

# --------------------------------------------------------------------

def fullnames_dates(fname=None, dirs=None):
    """
        Filenames with absolute paths and modification times in local directories.


        Definition
        ----------
        def fullnames_dates(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames incl. absolute paths, List of modification times


        Examples
        --------
        import os
        # get .dat filenames with absolute paths and modification times in directories 2013, 2014, ...
        fls, flsdate = fullnames_dates('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and modification times in current directory
        os.chdir('change/directory')
        fls, flsdate = fullnames_dates('*.dat')

        # get all filenames and modification times in current directory
        fls, flsdate = fullnames_dates()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llsd = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        fnames = [ os.path.abspath(f) for f in fnames ]
        lls  += fnames
        llsd += [ datetime.datetime.fromtimestamp(os.stat(n).st_mtime) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llsd = [ llsd[i] for i in ii ]
    return lls, llsd

# --------------------------------------------------------------------

def fullnames_dates_sizes(fname=None, dirs=None):
    """
        Filenames with absolute paths, modification times, and file sizes in local directories.


        Definition
        ----------
        def fullnames_dates_sizes(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames incl. absolute paths, List of modification times, List of file sizes


        Examples
        --------
        import os
        # get .dat filenames with absolute paths, modification times, and file sizes in directories 2013, 2014, ...
        fls, flsdate, flssize = fullnames_dates_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames, modification times, and file sizes in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = fullnames_dates_sizes('*.dat')

        # get all filenames, modification times, and file sizes in current directory
        fls, flsdate, flssize = fullnames_dates_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llsd = []
    llss = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        fnames = [ os.path.abspath(f) for f in fnames ]
        lls  += fnames
        llsd += [ datetime.datetime.fromtimestamp(os.stat(n).st_mtime) for n in fnames ]
        llss += [ int(os.stat(n).st_size) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llsd = [ llsd[i] for i in ii ]
    llss = [ llss[i] for i in ii ]
    return lls, llsd, llss


# --------------------------------------------------------------------

def fullnames_sizes(fname=None, dirs=None):
    """
        Filenames with absolute paths and file sizes in local directories.


        Definition
        ----------
        def fullnames_sizes(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames incl. absolute paths, List of file sizes


        Examples
        --------
        import os
        # get .dat filenames with absolute paths and file sizes in directories 2013, 2014, ...
        fls, flssize = fullnames_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and file sizes in current directory
        os.chdir('change/directory')
        fls, flssize = fullnames_sizes('*.dat')

        # get all filenames and file sizes in current directory
        fls, flssize = fullnames_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llss = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        fnames = [ os.path.abspath(f) for f in fnames ]
        lls  += fnames
        llss += [ int(os.stat(n).st_size) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llss = [ llss[i] for i in ii ]
    return lls, llss

# --------------------------------------------------------------------

def fullnames_times(*args, **kwargs):
    """
        Wrapper for fullnames_dates:
            def fullnames_dates(fname=None, dirs=None):


        Examples
        --------
        import os
        # get .dat filenames with absolute paths and modification times in directories 2013, 2014, ...
        fls, flsdate = fullnames_times('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and modification times in current directory
        os.chdir('change/directory')
        fls, flsdate = fullnames_times('*.dat')

        # get all filenames and modification times in current directory
        fls, flsdate = fullnames_times()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return fullnames_dates(*args, **kwargs)

# --------------------------------------------------------------------

def fullnames_times_sizes(*args, **kwargs):
    """
        Wrapper for fullnames_dates_sizes:
            def fullnames_dates_sizes(fname=None, dirs=None):


        Examples
        --------
        import os
        # get .dat filenames with absolute paths, modification times, and file sizes in directories 2013, 2014, ...
        fls, flsdate, flssize = fullnames_times_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames, modification times, and file sizes in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = fullnames_times_sizes('*.dat')

        # get all filenames, modification times, and file sizes in current directory
        fls, flsdate, flssize = fullnames_times_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return fullnames_dates_sizes(*args, **kwargs)

# --------------------------------------------------------------------

def last_fullname(fname=None, dirs=None):
    """
        Filename with absolute paths of last file in local directories.


        Definition
        ----------
        def last_fullname(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename incl. absolute paths of last modified file


        Examples
        --------
        import os
        # get last .dat filename with absolute paths in directories 2013, 2014, ...
        fls = last_fullname('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename in current directory
        os.chdir('change/directory')
        fls = last_fullname('*.dat')

        # get last filename in current directory
        fls = last_fullname()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld = fullnames_dates(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]]

# --------------------------------------------------------------------

def last_fullname_date(fname=None, dirs=None):
    """
        Filename with absolute paths and modification time of last file in local directories.


        Definition
        ----------
        def last_fullname_date(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename incl. absolute paths of last modified file, Modification time


        Examples
        --------
        import os
        # get last .dat filename with absolute paths and modification time in directories 2013, 2014, ...
        fls, flsdate = last_fullname_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flsdate = last_fullname_date('*.dat')

        # get last filename and modification time in current directory
        fls, flsdate = last_fullname_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld = fullnames_dates(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], fld[ii[-1]]

# --------------------------------------------------------------------

def last_fullname_date_size(fname=None, dirs=None):
    """
        Filename with absolute paths, modification time, and file size of last file in local directories.


        Definition
        ----------
        def last_fullname_date_size(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename incl. absolute paths of last modified file, Modification time, File size


        Examples
        --------
        import os
        # get last .dat filename with absolute paths, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = last_fullname_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = last_fullname_date_size('*.dat')

        # get last filename, modification time, and file size in current directory
        fls, flsdate, flssize = last_fullname_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld, flss = fullnames_dates_sizes(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], fld[ii[-1]], flss[ii[-1]]

# --------------------------------------------------------------------

def last_fullname_size(fname=None, dirs=None):
    """
        Filename with absolute paths and file size of last file in local directories.


        Definition
        ----------
        def last_fullname_size(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename incl. absolute paths of last modified file, File size


        Examples
        --------
        import os
        # get last .dat filename with absolute paths and file size in directories 2013, 2014, ...
        fls, flsdate = last_fullname_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and file size in current directory
        os.chdir('change/directory')
        fls, flsdate = last_fullname_date('*.dat')

        # get last filename and file size in current directory
        fls, flsdate = last_fullname_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld, flss = fullnames_dates_sizes(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], flss[ii[-1]]

# --------------------------------------------------------------------

def last_fullname_time(*args, **kwargs):
    """
        Wrapper for last_fullname_date:
            def last_fullname_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get last .dat filename with absolute paths and modification time in directories 2013, 2014, ...
        fls, flstime = last_fullname_time('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flstime = last_fullname_time('*.dat')

        # get last filename and modification time in current directory
        fls, flstime = last_fullname_time()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_date(*args, **kwargs)

# --------------------------------------------------------------------

def last_fullname_time_size(*args, **kwargs):
    """
        Wrapper for last_fullname_date_size:
            def last_fullname_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get last .dat filename with absolute paths, modification time, and file size in directories 2013, 2014, ...
        fls, flstime, flssize = last_fullname_time_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flstime, flssize = last_fullname_time_size('*.dat')

        # get last filename, modification time, and file size in current directory
        fls, flstime, flssize = last_fullname_time_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_date_size(*args, **kwargs)

# --------------------------------------------------------------------

def last_name(fname=None, dirs=None):
    """
        Filename of last file in local directories.


        Definition
        ----------
        def last_name(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename of last modified file


        Examples
        --------
        import os
        # get last .dat filename in directories 2013, 2014, ...
        fls = last_name('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename in current directory
        os.chdir('change/directory')
        fls = last_name('*.dat')

        # get last filename in current directory
        fls = last_name()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld = names_dates(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]]

# --------------------------------------------------------------------

def last_name_date(fname=None, dirs=None):
    """
        Filename and modification time of last file in local directories.


        Definition
        ----------
        def last_name_date(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename of last modified file, Modification time


        Examples
        --------
        import os
        # get last .dat filename and modification time in directories 2013, 2014, ...
        fls, flsdate = last_name_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flsdate = last_name_date('*.dat')

        # get last filename and modification time in current directory
        fls, flsdate = last_name_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld = names_dates(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], fld[ii[-1]]

# --------------------------------------------------------------------

def last_name_date_size(fname=None, dirs=None):
    """
        Filename, modification time, and file size of last file in local directories.


        Definition
        ----------
        def last_name_date_size(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename of last modified file, Modification time, File size


        Examples
        --------
        import os
        # get last .dat filename, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = last_name_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = last_name_date_size('*.dat')

        # get last filename, modification time, and file size in current directory
        fls, flsdate, flssize = last_name_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld, flss = names_dates_sizes(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], fld[ii[-1]], flss[ii[-1]]

# --------------------------------------------------------------------

def last_name_size(fname=None, dirs=None):
    """
        Filename and file size of last file in local directories.


        Definition
        ----------
        def last_name_size(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename of last modified file, File size


        Examples
        --------
        import os
        # get last .dat filename and file size in directories 2013, 2014, ...
        fls, flsdate = last_name_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and file size in current directory
        os.chdir('change/directory')
        fls, flsdate = last_name_date('*.dat')

        # get last filename and file size in current directory
        fls, flsdate = last_name_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld, flss = names_dates_sizes(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return fls[ii[-1]], flss[ii[-1]]

# --------------------------------------------------------------------

def last_name_time(*args, **kwargs):
    """
        Wrapper for last_name_date:
            last_name_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get last .dat filename and modification time in directories 2013, 2014, ...
        fls, flstime = last_name_time('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flstime = last_name_time('*.dat')

        # get last filename and modification time in current directory
        fls, flstime = last_name_time()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_name_date(*args, **kwargs)

# --------------------------------------------------------------------

def last_name_time_size(*args, **kwargs):
    """
        Wrapper for last_name_date_size:
            last_name_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get last .dat filename, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = last_name_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get last .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = last_name_date_size('*.dat')

        # get last filename, modification time, and file size in current directory
        fls, flsdate, flssize = last_name_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_name_date_size(*args, **kwargs)

# --------------------------------------------------------------------

def names(fname=None, dirs=None):
    """
        Filenames in local directories.


        Definition
        ----------
        def names(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames


        Examples
        --------
        import os
        # get .dat filenames in directories 2013, 2014, ...
        fls = names('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames in current directory
        os.chdir('change/directory')
        fls = names('*.dat')

        # get all filenames in current directory
        fls = names()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        lls  += fnames
        if i != '.':
            os.chdir(curdir)
    lls.sort()
    return lls

# --------------------------------------------------------------------

def names_dates(fname=None, dirs=None):
    """
        Filenames and modification times in local directories.


        Definition
        ----------
        def names_dates(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames, List of modification times


        Examples
        --------
        import os
        # get .dat filenames and modification times in directories 2013, 2014, ...
        fls, flsdate = names_dates('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and modification times in current directory
        os.chdir('change/directory')
        fls, flsdate = names_dates('*.dat')

        # get all filenames and modification times in current directory
        fls, flsdate = names_dates()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llsd = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        lls  += fnames
        llsd += [ datetime.datetime.fromtimestamp(os.stat(n).st_mtime) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llsd = [ llsd[i] for i in ii ]
    return lls, llsd

# --------------------------------------------------------------------

def names_dates_sizes(fname=None, dirs=None):
    """
        Filenames, modification times, and file sizes in local directories.


        Definition
        ----------
        def names_dates_sizes(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames, List of modification times, List of file sizes


        Examples
        --------
        import os
        # get .dat filenames, modification times, and file sizes in directories 2013, 2014, ...
        fls, flsdate = names_dates_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames, modification times, and file sizes in current directory
        os.chdir('change/directory')
        fls, flsdate = names_dates_sizes('*.dat')

        # get all filenames, modification times, and file sizes in current directory
        fls, flsdate = names_dates_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llsd = []
    llss = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        lls  += fnames
        llsd += [ datetime.datetime.fromtimestamp(os.stat(n).st_mtime) for n in fnames ]
        llss += [ int(os.stat(n).st_size) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llsd = [ llsd[i] for i in ii ]
    llss = [ llss[i] for i in ii ]
    return lls, llsd, llss


# --------------------------------------------------------------------

def names_sizes(fname=None, dirs=None):
    """
        Filenames and file sizes in local directories.


        Definition
        ----------
        def names_sizes(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        List of filenames, List of file sizes


        Examples
        --------
        import os
        # get .dat filenames and file sizes in directories 2013, 2014, ...
        fls, flssize = names_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and file sizes in current directory
        os.chdir('change/directory')
        fls, flssize = names_sizes('*.dat')

        # get all filenames and file sizes in current directory
        fls, flssize = names_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
        Modified, MC, Jun 2015 - dirs can be single directory
                               - dirs can be empty
    """
    if dirs is None:
        idirs = ['.']
    else:
        if isinstance(dirs, (list, tuple, np.ndarray, set)):
            idirs = dirs
        else:
            idirs = [dirs]
    lls  = []
    llss = []
    curdir = os.path.realpath(os.path.curdir)
    for i in idirs:
        if i != '.':
            os.chdir(str(i))
        if fname is None:
            fnames = os.listdir('.')
        else:
            fnames = glob.glob(fname)
        lls  += fnames
        llss += [ int(os.stat(n).st_size) for n in fnames ]
        if i != '.':
            os.chdir(curdir)
    ii = argsort(lls)
    lls  = [ lls[i] for i in ii ]
    llss = [ llss[i] for i in ii ]
    return lls, llss

# --------------------------------------------------------------------

def names_times(*args, **kwargs):
    """
        Wrapper for names_dates:
            def names_dates(fname=None, dirs=None):


        Examples
        --------
        import os
        # get .dat filenames and modification times in directories 2013, 2014, ...
        fls, flsdate = names_times('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames and modification times in current directory
        os.chdir('change/directory')
        fls, flsdate = names_times('*.dat')

        # get all filenames and modification times in current directory
        fls, flsdate = names_times()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return names_dates(*args, **kwargs)

# --------------------------------------------------------------------

def names_times_sizes(*args, **kwargs):
    """
        Wrapper for names_dates_sizes:
            def names_dates_sizes(fname=None, dirs=None):


        Examples
        --------
        import os
        # get .dat filenames, modification times, and file sizes in directories 2013, 2014, ...
        fls, flsdate, flssize = names_times_sizes('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get .dat filenames, modification times, and file sizes in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = names_times_sizes('*.dat')

        # get all filenames, modification times, and file sizes in current directory
        fls, flsdate, flssize = names_times_sizes()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return names_dates_sizes(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname(*args, **kwargs):
    """
        Wrapper for last_fullname:
            last_fullname(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths in directories 2013, 2014, ...
        fls = newest_fullname('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename in current directory
        os.chdir('change/directory')
        fls = newest_fullname('*.dat')

        # get newest filename in current directory
        fls = newest_fullname()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname_date(*args, **kwargs):
    """
        Wrapper for last_fullname_date:
            last_fullname_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths and modification time in directories 2013, 2014, ...
        fls, flsdate = newest_fullname_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flsdate = newest_fullname_date('*.dat')

        # get newest filename and modification time in current directory
        fls, flsdate = newest_fullname_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_date(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname_date_size(*args, **kwargs):
    """
        Wrapper for last_fullname_date_size:
            last_fullname_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = newest_fullname_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = newest_fullname_date_size('*.dat')

        # get newest filename, modification time, and file size in current directory
        fls, flsdate, flssize = newest_fullname_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_date_size(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname_size(*args, **kwargs):
    """
        Wrapper for last_fullname_size:
            last_fullname_size(fname=None, dirs=None):


        Definition
        ----------
        def newest_fullname_size(fname=None, dirs=None):


        Optional Input
        --------------
        fname        filename, filename globbing is possible such as '*.dat' (all files)
        dirs         list of or single directory names (default: '.')

        
        Output
        ------
        Filename incl. absolute paths of newest modified file, File size


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths and file size in directories 2013, 2014, ...
        fls, flsdate = newest_fullname_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and file size in current directory
        os.chdir('change/directory')
        fls, flsdate = newest_fullname_date('*.dat')

        # get newest filename and file size in current directory
        fls, flsdate = newest_fullname_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_size(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname_time(*args, **kwargs):
    """
        Wrapper for last_fullname_date:
            def last_fullname_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths and modification time in directories 2013, 2014, ...
        fls, flstime = newest_fullname_time('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flstime = newest_fullname_time('*.dat')

        # get newest filename and modification time in current directory
        fls, flstime = newest_fullname_time()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_time(*args, **kwargs)

# --------------------------------------------------------------------

def newest_fullname_time_size(*args, **kwargs):
    """
        Wrapper for last_fullname_date_size:
            def last_fullname_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename with absolute paths, modification time, and file size in directories 2013, 2014, ...
        fls, flstime, flssize = newest_fullname_time_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flstime, flssize = newest_fullname_time_size('*.dat')

        # get newest filename, modification time, and file size in current directory
        fls, flstime, flssize = newest_fullname_time_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_fullname_date_size(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name(*args, **kwargs):
    """
        Wrapper for last_name:
            last_name(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename in directories 2013, 2014, ...
        fls = newest_name('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename in current directory
        os.chdir('change/directory')
        fls = newest_name('*.dat')

        # get newest filename in current directory
        fls = newest_name()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_name(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name_date(*args, **kwargs):
    """
        Wrapper for last_name_date:
            last_name_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename and modification time in directories 2013, 2014, ...
        fls, flsdate = newest_name_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flsdate = newest_name_date('*.dat')

        # get newest filename and modification time in current directory
        fls, flsdate = newest_name_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    # Files
    fls, fld = names_dates(fname, dirs=dirs)
    if len(fls) == 0: # nothing in there
        return None
    # Dates
    fldec = []
    for idat in fld:
        dec = (idat.year*366. + idat.month*31. + idat.day*1. +
               idat.hour/24. + idat.minute/1440. + idat.second/86400.)
        fldec.append(dec)
    ii = argsort(fldec)
    return last_name_date(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name_date_size(*args, **kwargs):
    """
        Wrapper for last_name_date_size:
            last_name_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = newest_name_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = newest_name_date_size('*.dat')

        # get newest filename, modification time, and file size in current directory
        fls, flsdate, flssize = newest_name_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_name_date_size(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name_size(*args, **kwargs):
    """
        Wrapper for last_name_size:
            last_name_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename and file size in directories 2013, 2014, ...
        fls, flsdate = newest_name_date('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and file size in current directory
        os.chdir('change/directory')
        fls, flsdate = newest_name_date('*.dat')

        # get newest filename and file size in current directory
        fls, flsdate = newest_name_date()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return last_name_size(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name_time(*args, **kwargs):
    """
        Wrapper for newest_name_date:
            newest_name_date(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename and modification time in directories 2013, 2014, ...
        fls, flstime = newest_name_time('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename and modification time in current directory
        os.chdir('change/directory')
        fls, flstime = newest_name_time('*.dat')

        # get newest filename and modification time in current directory
        fls, flstime = newest_name_time()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return newest_name_date(*args, **kwargs)

# --------------------------------------------------------------------

def newest_name_date_size(*args, **kwargs):
    """
        Wrapper for newest_name_date_size:
            newest_name_date_size(fname=None, dirs=None):


        Examples
        --------
        import os
        # get newest .dat filename, modification time, and file size in directories 2013, 2014, ...
        fls, flsdate, flssize = newest_name_date_size('*.dat', dirs=glob.glob('[0-9][0-9][0-9][0-9]'))

        # get newest .dat filename, modification time, and file size in current directory
        os.chdir('change/directory')
        fls, flsdate, flssize = newest_name_date_size('*.dat')

        # get newest filename, modification time, and file size in current directory
        fls, flsdate, flssize = newest_name_date_size()


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Dec 2014
    """
    return newest_name_date_size(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
