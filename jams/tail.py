#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import os

__all__ = ['tail']

def tail(ifile, n=10, noblank=False, comment='', keepnewline=False):
    """
        Return list with last n lines of file.
        Equivalent to Unix tail utility.

        If n is omitted, it defaults to 10.

        Blank (only whitespace) and comment lines can be excluded.


        Definition
        ----------
        def tail(ifile, n=10, noblank=False, comment='', keepnewline=False):


        Input
        -----
        ifile        file handle or filename


        Optional input parameters
        -------------------------
        n            Get the last n lines from file (default: 10)
        comment      Line gets excluded if first character of line is in comment sequence.
                     Sequence can be string, list or tuple.


        Options
        -------
        noblank      if True, excludes all lines that consists only of
                     whitespace characters.
        keepnewline  if True, newline character at end of string will not be removed.


        Output
        ------
        list with strings of last n lines in file


        Restrictions
        ------------
        For efficiency, the utility first seeks the file pointer to file.seek(-n*8192, os.SEEK_END)
        and then starts reading.
        If the lines are too long, 8192 might not be enough.
        The utility then attempts to determine the file seek from the length of the lines and
        redoes the reading. This is still experimental.


        Examples
        --------
        # Create some data
        >>> filename = 'test.dat'
        >>> ifile = open(filename,'w')
        >>> ifile.writelines('First line\\n')
        >>> ifile.writelines(' \\n')
        >>> ifile.writelines('# First comment\\n')
        >>> ifile.writelines('! Second comment\\n')
        >>> ifile.writelines('Last line\\n')
        >>> ifile.close()

        # Read last 10 lines
        >>> print(tail(filename))
        ['First line', ' ', '# First comment', '! Second comment', 'Last line']

        # Read last 4 lines
        >>> print(tail(filename, 4))
        [' ', '# First comment', '! Second comment', 'Last line']

        # Read last 4 lines, keep newlines in strings
        >>> print(tail(filename, n=4, keepnewline=True))
        [' \\n', '# First comment\\n', '! Second comment\\n', 'Last line\\n']

        # Read last 4 lines excluding blank lines
        >>> print(tail(filename, 4, noblank=True))
        ['First line', '# First comment', '! Second comment', 'Last line']

        # Read last 4 lines excluding lines starting with #
        >>> print(tail(filename, 4, comment='#'))
        ['First line', ' ', '! Second comment', 'Last line']

        # Read last 4 lines excluding lines starting with # and !
        >>> print(tail(filename, 4, comment='#!'))
        ['First line', ' ', 'Last line']

        # Read last 4 lines excluding lines starting with # and S
        >>> print(tail(filename, 4, comment='#S'))
        ['First line', ' ', '! Second comment', 'Last line']

        # Read last 4 lines excluding lines starting with # and L
        >>> print(tail(filename, 4, comment=('#','L')))
        ['First line', ' ', '! Second comment']

        # Read last 4 lines excluding lines starting with # and ! but comment given as list not string
        >>> print(tail(filename, 4, comment=['#','!']))
        ['First line', ' ', 'Last line']

        # Read last 4 lines excluding lines starting with # and ! and blank lines
        >>> print(tail(filename, 4, comment='#!', noblank=True))
        ['First line', 'Last line']

        # Read last 2 lines
        >>> print(tail(filename, n=2))
        ['! Second comment', 'Last line']

        # Read last line excluding lines starting with # and ! and blank lines
        >>> print(tail(filename, 1, comment='#!', noblank=True))
        ['Last line']

        # Read last lines from an open file handle rather than a file name
        >>> f = open(filename, 'r')
        >>> print(tail(f, 1, comment='#!', noblank=True))
        ['Last line']
        >>> f.close()

        # Check intermediate line lengths
        >>> n=1000
        >>> f = open(filename, 'w')
        >>> print('a'*n, file=f)
        >>> print('b'*n, file=f)
        >>> print('c'*n, file=f)
        >>> print('d'*n, file=f)
        >>> f.close()
        >>> ll = tail(filename, n=2)
        >>> if (ll[0] != 'c'*n) or (ll[1] != 'd'*n):
        ...     print('failed')

        # Check long line lengths
        >>> n=10000
        >>> f = open(filename, 'w')
        >>> print('a'*n, file=f)
        >>> print('b'*n, file=f)
        >>> print('c'*n, file=f)
        >>> print('d'*n, file=f)
        >>> f.close()
        >>> ll = tail(filename, n=2)
        >>> if (ll[0] != 'c'*n) or (ll[1] != 'd'*n):
        ...     print('failed')

        # Clean up
        >>> import os
        >>> os.remove(filename)


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
        Written,  MC, Dec 2014
                  MC, Nov 2016 - adapted file handling to Python 2 and 3
    """
    import sys
    if sys.version_info > (3,0):
        import io
        f = ifile if isinstance(ifile,io.TextIOWrapper) else open(ifile,"r")
    else:
        f = ifile if isinstance(ifile,file) else open(ifile,"r")

    # Seek
    bufsize = 8192*n
    try:    # might be out of file
        f.seek(-bufsize, os.SEEK_END)
    except: # stay were you are
        pass

    # Read lines
    doit = True
    while doit:
        count = 0
        liste = [[]]*n
        imax  = []
        for line in f:
            if not line:
                break
            if noblank and not line.strip():
                continue
            if line[0] in comment:
                continue
            if not keepnewline:
                line = line.strip("\n")
            liste.insert(n, line)
            tmp = liste.pop(0)
            imax.append(len(line.strip("\n"))) # lines include \0 at the end.
            count += 1
        imax = max(imax)
        if (imax < bufsize//n) or (count > n):
            doit = False
        else:
            bufsize = imax*8*n # 8 is for safety
            try:    # might be out of file
                f.seek(-bufsize, os.SEEK_END)
            except:
                f.seek(0, os.SEEK_SET)

    if f is not ifile:
        f.close()

    if count < n:
        return liste[-count:]
    else:
        return liste


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # filename = 'test.dat'
        
    # n=1000
    # f = open(filename, 'w')
    # print('a'*n, file=f)
    # print('b'*n, file=f)
    # print('c'*n, file=f)
    # print('d'*n, file=f)
    # f.close()
    # ll = tail(filename, n=2)
    # if (ll[0] != 'c'*n) or (ll[1] != 'd'*n):
    #     print('failed')

    # n=10000
    # f = open(filename, 'w')
    # print('a'*n, file=f)
    # print('b'*n, file=f)
    # print('c'*n, file=f)
    # print('d'*n, file=f)
    # f.close()
    # ll = tail(filename, n=2)
    # if (ll[0] != 'c'*n) or (ll[1] != 'd'*n):
    #     print('failed')

    # import os
    # os.remove(filename)
