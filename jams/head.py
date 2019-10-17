#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

__all__ = ['head']

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
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014-2016 Matthias Cuntz, David Schaefer - mc (at) macu (dot) de

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
        Written,  MC, Nov 2014
        Modified, DS, Nov 2014
                  MC, Nov 2016 - adapted file handling to Python 2 and 3
    """
    # Open file
    import sys
    if sys.version_info > (3,0):
        import io
        f = ifile if isinstance(ifile,io.TextIOWrapper) else open(ifile,"r")
    else:
        f = ifile if isinstance(ifile,file) else open(ifile,"r")

    # consume lines
    for _ in range(skip):
        f.readline()

    # Read lines
    count = 0
    header = []
    while count < n:
        line = f.readline()
        if not line:
            break
        if noblank and not line.strip():
            continue
        if line[0] in comment:
            continue
        if not keepnewline:
            line = line.strip("\n")
        header.append(line)
        count += 1

    if f is not ifile:
        f.close()

    return header


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
