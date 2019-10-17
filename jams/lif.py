#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

def lif(ifile, noblank=False, comment='', skip=0, maxcol=False):
    """
        Count number of lines in a file.

        Blank (only whitespace) and comment lines can be excluded.


        Definition
        ----------
        def lif(ifile, noblank=False, comment='', skip=0):


        Input
        -----
        ifile         source file name


        Optional input parameters
        -------------------------
        comment      line gets excluded if first character of line is
                     in comment sequence.
                     Sequence can be string, list or tuple
        skip         number of lines to skip at the beginning of
                     the file (default 0)


        Options
        -------
        noblank      if True, exclude all lines that consists only of
                     whitespace characters
        maxcol       if True, return also maximum amount of characters in one line


        Output
        ------
        if maxcol:
            number of lines in file, maximum characters in a line
        else:
            number of lines in file


        Restrictions
        ------------
        Only ascii files.


        Examples
        --------
        # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('First line\\n')
        >>> file.writelines(' \\n')
        >>> file.writelines('# First comment\\n')
        >>> file.writelines('! Second comment\\n')
        >>> file.writelines('Last line\\n')
        >>> file.close()

        # Count lines
        >>> print(lif(filename))
        5
        >>> print(lif(filename,noblank=True))
        4
        >>> print(lif(filename,comment='#'))
        4
        >>> print(lif(filename,comment='#!'))
        3
        >>> print(lif(filename,comment='#S'))
        4
        >>> print(lif(filename,comment=('#','L')))
        3
        >>> print(lif(filename,comment=['#','!']))
        3
        >>> print(lif(filename,comment='#!',noblank=True))
        2
        >>> print(lif(filename,skip=2))
        3
        >>> print(lif(filename,skip=2,maxcol=True))
        (3, 16)

        # Clean up
        >>> import os
        >>> os.remove(filename)


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2009-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jul 2009
        Modified, MC, Nov 2012 - maxcol
                  MC, Feb 2013 - ported to Python 3
                  MC, Dec 2014 - changed similar to elegant code of David for head.py
                  MC, Nov 2016 - adapted file handling to Python 2 and 3
    """
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
    imax = []
    for line in f:
        if not line:
            break
        if noblank and not line.strip():
            continue
        if line[0] in comment:
            continue
        imax.append(len(line.strip("\n"))) # lines include \0 at the end.
        count += 1

    if f is not ifile:
        f.close()

    if maxcol:
        return count, max(imax)
    else:
        return count


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
