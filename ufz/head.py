#!/usr/bin/env python
from __future__ import print_function
import string

def head(ifile, n=10, noblank=False, comment='', skip=0, keepnewline=False):
    """
        Return list with first n lines of file.
        Equivalent to Unix head utility.

        If n is omitted, it defaults to 10.

        Blank (only whitespace) and comment lines can be excluded.


        Definition
        ----------
        def head(ifile, n=10, noblank=False, comment='', skip=0, keepnewline=True):


        Input
        -----
        ifile        file handle or filename


        Optional input parameters
        -------------------------
        n            Get the first n lines from file (default: 10)
        comment      Line gets excluded if first character of line is in comment sequence.
                     Sequence can be string, list or tuple.
        skip         Number of lines to skip at the beginning of
                     the file (default 0).


        Options
        -------
        noblank      if True, excludes all lines that consists only of
                     whitespace characters.
        keepnewline  if True, newline character at end of string will not be removed.


        Output
        ------
        list with strings of first n lines in file


        Restrictions
        ------------
        Only ascii files.


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

        # Read first 10 lines
        >>> print(head(filename))
        ['First line', ' ', '# First comment', '! Second comment', 'Last line']

        # Read first 10 lines, keep newlines in strings
        >>> print(head(filename, keepnewline=True))
        ['First line\\n', ' \\n', '# First comment\\n', '! Second comment\\n', 'Last line\\n']

        # Read first 10 lines excluding blank lines
        >>> print(head(filename, noblank=True))
        ['First line', '# First comment', '! Second comment', 'Last line']

        # Read first 10 lines excluding lines starting with #
        >>> print(head(filename, comment='#'))
        ['First line', ' ', '! Second comment', 'Last line']

        # Read first 10 lines excluding lines starting with # and !
        >>> print(head(filename, comment='#!'))
        ['First line', ' ', 'Last line']

        # Read first 10 lines excluding lines starting with # and S
        >>> print(head(filename, comment='#S'))
        ['First line', ' ', '! Second comment', 'Last line']

        # Read first 10 lines excluding lines starting with # and L
        >>> print(head(filename, comment=('#','L')))
        ['First line', ' ', '! Second comment']

        # Read first 10 lines excluding lines starting with # and ! but comment given as list not string
        >>> print(head(filename, comment=['#','!']))
        ['First line', ' ', 'Last line']

        # Read first 10 lines excluding lines starting with # and ! and blank lines
        >>> print(head(filename, comment='#!', noblank=True))
        ['First line', 'Last line']

        # Read first 10 lines excluding the first 2 lines
        >>> print(head(filename, skip=2))
        ['# First comment', '! Second comment', 'Last line']

        # Read first 3 lines
        >>> print(head(filename, n=3))
        ['First line', ' ', '# First comment']

        # Read first 1 line excluding lines starting with # and ! and blank lines while skipping 1 line
        >>> print(head(filename, 1, skip=1, comment='#!', noblank=True))
        ['Last line']

        # Read first 10 lines from an open file handle rather than a file name
        >>> f = open(filename, 'r')
        >>> print(head(f, 1, skip=1, comment='#!', noblank=True))
        ['Last line']
        >>> f.close()

        # Clean up
        >>> import os
        >>> os.remove(filename)


        Notes
        -----
        if file handle is given, the file pointer will point to the next line after the n-th line read.


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
        Written,  MC, Nov 2014
    """
    # Open file
    handle = False
    if isinstance(ifile, file):
        handle = True
        f = ifile
    else:
        try:
            f = open(ifile, 'r')
        except IOError:
            raise ValueError('Cannot open file '+ifile)
    # Read lines
    count = 0
    header = []
    if skip > 0:
        iskip = 0
        while iskip < skip:
            l = f.readline()
            iskip += 1
    if noblank and (comment != ''):        # exclude blank, exclude comment
        while count <= n:
            ll = f.readline()
            if len(ll) == 0: break
            l  = ll.strip(string.whitespace)
            if (l != ''):
                if (l[0] not in comment):
                    if keepnewline:
                        header.append(ll)
                    else:
                        if ll.endswith('\n'):
                            header.append(ll[:-1])
                        else:
                            header.append(ll)
                    count += 1
    elif noblank and (comment == ''):      # exclude blank, include comment
        while count <= n:
            ll = f.readline()
            if len(ll) == 0: break
            l     = ll.strip(string.whitespace)
            if (l != ''):
                if keepnewline:
                    header.append(ll)
                else:
                    if ll.endswith('\n'):
                        header.append(ll[:-1])
                    else:
                        header.append(ll)
                count += 1
    elif (not noblank) and (comment != ''):# include blank, exclude comment
        while count <= n:
            ll = f.readline()
            if len(ll) == 0: break
            l  = ll.strip(string.whitespace)
            if (l == ''):
                if keepnewline:
                    header.append(ll)
                else:
                    if ll.endswith('\n'):
                        header.append(ll[:-1])
                    else:
                        header.append(ll)
                count += 1
            else:
                if (l[0] not in comment):
                    if keepnewline:
                        header.append(ll)
                    else:
                        if ll.endswith('\n'):
                            header.append(ll[:-1])
                        else:
                            header.append(ll)
                    count += 1
    else:                                  # include blank, include comment
        for i in range(n):
            ll = f.readline()
            if len(ll) == 0: break
            if keepnewline:
                header.append(ll)
            else:
                if ll.endswith('\n'):
                    header.append(ll[:-1])
                else:
                    header.append(ll)
    if not handle:
        f.close()

    return header


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)