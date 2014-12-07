#!/usr/bin/env python
from __future__ import print_function
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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
    """
    # Open file
    f = ifile if isinstance(ifile, file) else open(ifile,"r")

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
