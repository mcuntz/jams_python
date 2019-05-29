#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import os
import datetime
import glob
import time
from jams.argsort import argsort

__all__ = ['get_binary', 'get_check_binary',
           'get_mac_ascii', 'get_check_mac_ascii',
           'get_unix_ascii', 'get_check_unix_ascii',
           'get_windows_ascii', 'get_check_windows_ascii',
           'get_names', 'get_names_dates', 'get_names_sizes', 'get_names_times',
           'get_names_dates_sizes', 'get_names_times_sizes',
           'get_size', 'get_sizes',
           'set_mtime',
           'put_binary', 'put_check_binary']

# ------------------------------------------------------------------------------------------

def get_binary(ftp, fname):
    """
        Get a binary file from an open FTP connection.


        Definition
        ----------
        def get_binary(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        get_binary(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014 - modified from
                                 http://www.java2s.com/Tutorial/Python/0420__Network/Binaryfiledownload.htm
    """
    ftp.voidcmd("TYPE I")
    datasock, estsize = ftp.ntransfercmd("RETR "+fname)
    transbytes = 0
    fd = open(fname, 'wb')
    while True:
        buf = datasock.recv(2048)
        if not len(buf):
            break
        fd.write(buf)
        transbytes += len(buf)
    fd.close()
    datasock.close()
    ftp.voidresp()
    set_mtime(ftp, fname)

# ------------------------------------------------------------------------------------------

def get_check_binary(ftp, fname):
    """
        Get a binary file from an open FTP connection and
        check that local filesize is the same as the remote file size.


        Definition
        ----------
        def get_check_binary(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        if not get_check_binary(ftp, 'test.dat'):
            print('Error in transfer')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    rsize = get_size(ftp, fname)
    get_binary(ftp, fname)
    lsize = os.stat(fname).st_size
    if rsize == lsize:
        return True
    else:
        return False

# ------------------------------------------------------------------------------------------

def get_mac_ascii(ftp, fname):
    """
        Get a Mac ascii file from an open FTP connection.


        Definition
        ----------
        def get_mac_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        get_mac_ascii(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014 - modified from
                                 http://www.java2s.com/Tutorial/Python/0420__Network/ASCIIfiledownload.htm
    """
    def writeline(data):
        fd.write(data + '\r') # Mac-Format
    fd = open(fname, 'wt')
    ftp.retrlines('RETR '+fname, writeline)
    fd.close()
    set_mtime(ftp, fname)

# ------------------------------------------------------------------------------------------

def get_check_mac_ascii(ftp, fname):
    """
        Get a Mac ascii file from an open FTP connection and
        check that local filesize is the same as the remote file size.


        Definition
        ----------
        def get_check_mac_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        if not get_check_mac_ascii(ftp, 'test.dat'):
            print('Error in transfer')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    rsize = get_size(ftp, fname)
    get_mac_ascii(ftp, fname)
    lsize = os.stat(fname).st_size
    if rsize == lsize:
        return True
    else:
        return False

# ------------------------------------------------------------------------------------------

def get_unix_ascii(ftp, fname):
    """
        Get a Unix/Linux ascii file from an open FTP connection.


        Definition
        ----------
        def get_unix_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        get_unix_ascii(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014 - modified from
                                 http://www.java2s.com/Tutorial/Python/0420__Network/ASCIIfiledownload.htm
    """
    def writeline(data):
        fd.write(data + '\n') # Unix-Format
    fd = open(fname, 'wt')
    ftp.retrlines('RETR '+fname, writeline)
    fd.close()
    set_mtime(ftp, fname)

# ------------------------------------------------------------------------------------------

def get_check_unix_ascii(ftp, fname):
    """
        Get a Unix/Linux ascii file from an open FTP connection and
        check that local filesize is the same as the remote file size.


        Definition
        ----------
        def get_check_unix_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        if not get_check_unix_ascii(ftp, 'test.dat'):
            print('Error in transfer')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    rsize = get_size(ftp, fname)
    get_unix_ascii(ftp, fname)
    lsize = os.stat(fname).st_size
    if rsize == lsize:
        return True
    else:
        return False

# ------------------------------------------------------------------------------------------

def get_windows_ascii(ftp, fname):
    """
        Get a Windows ascii file from an open FTP connection.


        Definition
        ----------
        def get_windows_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        get_windows_ascii(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014 - modified from
                                 http://www.java2s.com/Tutorial/Python/0420__Network/ASCIIfiledownload.htm
    """
    def writeline(data):
        fd.write(data + '\r\n') # Windows-Format
    fd = open(fname, 'wt')
    ftp.retrlines('RETR '+fname, writeline)
    fd.close()
    set_mtime(ftp, fname)

# ------------------------------------------------------------------------------------------

def get_check_windows_ascii(ftp, fname):
    """
        Get a Windows ascii file from an open FTP connection and
        check that local filesize is the same as the remote file size.


        Definition
        ----------
        def get_check_windows_ascii(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        if not get_check_windows_ascii(ftp, 'test.dat'):
            print('Error in transfer')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    rsize = get_size(ftp, fname)
    get_windows_ascii(ftp, fname)
    lsize = os.stat(fname).st_size
    if rsize == lsize:
        return True
    else:
        return False

# ------------------------------------------------------------------------------------------

def get_names(ftp, fname=None):
    """
        Get filename(s) in current directory of open FTP connection.


        Definition
        ----------
        def get_names(ftp, fname=None):


        Input
        -----
        ftp          instance of the ftplib.FTP class


        Optional input
        --------------
        fname        filename, filename globbing is possible such as '*.dat'


        Output
        ------
        Filenames in current ftp directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls = get_names(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    fdir = []
    if fname is not None:
        ftp.dir(fname, fdir.append)
    else:
        ftp.dir(fdir.append)
    names = []
    for f in fdir:
        i = f.split()
        if i[0][0] == '-': # assure file
            names.append(' '.join(i[8:]))      # filename
    names.sort()
    return names

# ------------------------------------------------------------------------------------------

def get_names_dates(ftp, fname=None):
    """
        Get filename(s) and file modification times in current directory of open FTP connection.


        Definition
        ----------
        def get_names_dates(ftp, fname=None):


        Input
        -----
        ftp          instance of the ftplib.FTP class


        Optional input
        --------------
        fname        filename, filename globbing is possible such as '*.dat'


        Output
        ------
        Filenames, Modification Times


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls, flsdate = get_names_dates(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
        Modified, MC, Jan 2016 - return datetime.datetime
    """
    today    = datetime.date.today()
    thisyear = str(today.year)
    fdir = []
    if fname is not None:
        ftp.dir(fname, fdir.append)
    else:
        ftp.dir(fdir.append)
    names = []
    dates = []
    for f in fdir:
        i = f.split()
        if i[0][0] == '-': # assure file
            names.append(' '.join(i[8:]))      # filename
            yr = i[7]               # year/time column
            if ':' in yr:
                hhmm = yr
                yr = thisyear
            else:
                hhmm = '00:00'
            hhmi = hhmm.split(':')
            hh, mi = [ int(hi) for hi in hhmi ]
            yyyy  = int(yr)
            dd    = int(i[6]) # day
            mm    = i[5]      # month as Jan, Feb, ...
            mm    = datetime.datetime.strptime(mm,"%b").month
            idate = datetime.date(yyyy,mm,dd)
            if idate > today: # time instead for files younger than 6 month, even in the last year
                yyyy -= 1
            idate = datetime.datetime(yyyy,mm,dd,hh,mi)
            dates.append(idate) # file date
    ii = argsort(names)
    names = [ names[i] for i in ii ]
    dates = [ dates[i] for i in ii ]
    return names, dates

# ------------------------------------------------------------------------------------------

def get_names_dates_sizes(ftp, fname=None):
    """
        Get filename(s),  file modification times, and filesizes in current directory of open FTP connection.


        Definition
        ----------
        def get_names_dates_sizes(ftp, fname=None):


        Input
        -----
        ftp          instance of the ftplib.FTP class


        Optional input
        --------------
        fname        filename, filename globbing is possible such as '*.dat'


        Output
        ------
        Filenames, Modification Times, filesizes


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls, flsdate, flssize = get_names_dates_sizes(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
        Modified, MC, Jan 2016 - return datetime.datetime
    """
    today    = datetime.date.today()
    thisyear = str(today.year)
    fdir = []
    if fname is not None:
        ftp.dir(fname, fdir.append)
    else:
        ftp.dir(fdir.append)
    names = []
    dates = []
    sizes = []
    for f in fdir:
        i = f.split()
        if i[0][0] == '-': # assure file
            names.append(' '.join(i[8:])) # filename
            sizes.append(int(i[4]))       # file size
            yr = i[7]         # year/time column
            if ':' in yr:
                hhmm = yr
                yr = thisyear
            else:
                hhmm = '00:00'
            hhmi = hhmm.split(':')
            hh, mi = [ int(hi) for hi in hhmi ]
            yyyy  = int(yr)
            dd    = int(i[6]) # day
            mm    = i[5]      # month as Jan, Feb, ...
            mm    = datetime.datetime.strptime(mm,"%b").month
            idate = datetime.date(yyyy,mm,dd)
            if idate > today: # time instead for files younger than 6 month, even in the last year
                yyyy -= 1
            idate = datetime.datetime(yyyy,mm,dd,hh,mi)
            dates.append(idate)           # file date
    ii = argsort(names)
    names = [ names[i] for i in ii ]
    dates = [ dates[i] for i in ii ]
    sizes = [ sizes[i] for i in ii ]
    return names, dates, sizes

# ------------------------------------------------------------------------------------------

def get_names_sizes(ftp, fname=None):
    """
        Get filename(s) and filesizes in current directory of open FTP connection.


        Definition
        ----------
        def get_names_sizes(ftp, fname=None):


        Input
        -----
        ftp          instance of the ftplib.FTP class


        Optional input
        --------------
        fname        filename, filename globbing is possible such as '*.dat'


        Output
        ------
        Filenames, filesizes


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls, flssize = get_names_sizes(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    fdir = []
    if fname is not None:
        ftp.dir(fname, fdir.append)
    else:
        ftp.dir(fdir.append)
    names = []
    sizes = []
    for f in fdir:
        i = f.split()
        if i[0][0] == '-': # assure file
            names.append(' '.join(i[8:]))      # filename
            sizes.append(int(i[4])) # file size
    ii = argsort(names)
    names = [ names[i] for i in ii ]
    sizes = [ sizes[i] for i in ii ]
    return names, sizes

# ------------------------------------------------------------------------------------------

def get_names_times(*args, **kwargs):
    """
        Wrapper for get_names_dates:
            def get_names_dates(ftp, fname=None):


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls, flsdate = get_names_times(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    return get_names_dates(*args, **kwargs)

# ------------------------------------------------------------------------------------------

def get_names_times_sizes(*args, **kwargs):
    """
        Wrapper for get_names_dates_sizes:
            def get_names_dates_sizes(ftp, fname=None):


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fls, flsdate, flssize = get_names_times_sizes(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    return get_names_dates_sizes(*args, **kwargs)

# ------------------------------------------------------------------------------------------

def get_size(ftp, fname):
    """
        Get file size of file on open FTP connection.


        Definition
        ----------
        def get_size(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        Filesize


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        fs = get_size(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    fdir = []
    ftp.dir(fname, fdir.append)
    return int(fdir[0].split()[4])

# ------------------------------------------------------------------------------------------

def get_sizes(ftp, fname=None):
    """
        Get filesizes of files in current directory of open FTP connection.


        Definition
        ----------
        def get_sizes(ftp, fname=None):


        Input
        -----
        ftp          instance of the ftplib.FTP class


        Optional input
        --------------
        fname        filename, filename globbing is possible such as '*.dat'


        Output
        ------
        Filesizes


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        flssize = get_sizes(ftp, '*.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    fdir = []
    if fname is not None:
        ftp.dir(fname, fdir.append)
    else:
        ftp.dir(fdir.append)
    sizes = []
    for f in fdir:
        i = f.split()
        if i[0][0] == '-': # assure file
            sizes.append(int(i[-5])) # file size
    sizes.sort()
    return sizes

# ------------------------------------------------------------------------------------------

def set_mtime(ftp, fname):
    """
        Set the access and modification times of a local file
        to time of the file with the same name in the current directory of an open FTP connection.


        Definition
        ----------
        def set_mtime(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        fname in current local directory has the same access and modification time as remote file


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        set_mtime(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    today     = datetime.date.today()
    thisyear  = str(today.year)
    fdir = []
    ftp.dir(fname, fdir.append)
    yr = fdir[0].split()[-2]         # year/time column
    if ':' in yr:
        hhmm = yr
        yr = thisyear
    else:
        hhmm = '00:00'
    yyyy  = int(yr)
    dd    = int(fdir[0].split()[-3]) # day
    mm    = fdir[0].split()[-4]      # month as Jan, Feb, ...
    mm    = datetime.datetime.strptime(mm,"%b").month
    idate = datetime.date(yyyy,mm,dd)
    if idate > today: # time instead for files younger than 6 month, even in the last year
        yyyy -= 1
    mtime = int(time.mktime(time.strptime(str(dd)+'.'+str(mm)+'.'+str(yyyy)+' '+hhmm, "%d.%m.%Y %H:%M")))
    os.utime(fname, (mtime, mtime))

# ------------------------------------------------------------------------------------------

def put_binary(ftp, fname):
    """
        Put a binary file to an open FTP connection.


        Definition
        ----------
        def put_binary(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        put_binary(ftp, 'test.dat')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2016
    """
    ftp.voidcmd("TYPE I")
    myfile = open(fname, 'rb')
    ftp.storbinary('STOR ' + fname, myfile)
    myfile.close()

# ------------------------------------------------------------------------------------------

def put_check_binary(ftp, fname):
    """
        Put a binary file to an open FTP connection and
        check that local filesize is the same as the remote file size.


        Definition
        ----------
        def put_check_binary(ftp, fname):


        Input
        -----
        ftp          instance of the ftplib.FTP class
        fname        filename


        Output
        ------
        File fname in current local directory


        Examples
        --------
        import ftplib
        import os
        ftp = ftplib.FTP("ftp.server.de")
        ftp.login("user", "password")
        ftp.cwd('ftp/directory')
        os.chdir('local/directory')
        if not put_check_binary(ftp, 'test.dat'):
            print('Error in transfer')
        ftp.close()


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2016
    """
    lsize = os.stat(fname).st_size
    put_binary(ftp, fname)
    rsize = get_size(ftp, fname)
    if rsize == lsize:
        return True
    else:
        return False


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
