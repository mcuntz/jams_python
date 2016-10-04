#!/usr/bin/env python
from __future__ import print_function

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

        Copyright 2009-2014 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Nov 2012 - maxcol
                  MC, Feb 2013 - ported to Python 3
                  MC, Dec 2014 - changed similar to elegant code of David for head.py
    """
    # Open file
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
