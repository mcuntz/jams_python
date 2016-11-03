#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def fread(file, nc=0, skip=0, cskip=0, hskip=0, separator=None,
          squeeze=False, reform=False, skip_blank=False, comment=None,
          fill=False, fill_value=0, strip=None,
          header=False, full_header=False,
          transpose=False, strarr=False):
    """
        Read numbers into 2D-array from a file.

        This routines is exactly the same as sread but transforms
        everything to floats, handling NaN and Inf.


        Definition
        ----------
        def fread(file, nc=0, skip=0, cskip=0, hskip=0, separator=None,
                  squeeze=False, reform=False, skip_blank=False, comment=None,
                  fill=False, fill_value=0, strip=None,
                  header=False, full_header=False,
                  transpose=False, strarr=False):


        Input
        -----
        file         source file name


        Optional Input Parameters
        -------------------------
        nc           number of columns to be read (default: all (nc<=0))
                     nc can be a vector of column indeces,
                     starting with 0; cskip will be ignored then.
        skip         number of lines to skip at the beginning of file (default: 0)
        cskip        number of columns to skip at the beginning of each line (default: 0)
        hskip        number of lines in skip that do not belong to header (default: 0)
        separator    column separator
                     If not given, columns separator are (in order):
                     comma, semi-colon, whitespace
        comment      line gets excluded if first character of line is in comment sequence
                     sequence can be e.g. string, list or tuple
        fill_value   value to fill in if not enough columns in line
                     and fill=True (default: 0, and '' for header)
        strip        Strip strings with str.strip(strip).
                     If None then strip quotes " and ' (default).
                     If False then no strip (30% faster).
                     Otherwise strip character given by strip.


        Options
        -------
        squeeze      True:  2-dim array will be cleaned of degenerated
                            dimension, i.e. results in vector
                     False: array will be two-dimensional as read (default)
        reform       Same as squeeze.
        skip_blank   True:  continues reading after blank line
                     False: stops reading at first blank line (default)
        fill         True:  fills in fill_value if not enough columns in line
                     False: stops execution and returns None if not enough
                            columns in line (default)
        header       True:  header strings will be returned
                     False  numbers in file will be returned (default)
        full_header  True:  header is a string vector of the skipped rows
                     False: header will be split in columns,
                            exactly as the data, and will hold only the
                            selected columns (default)
        transpose    True:  column-major format output(0:ncolumns,0:nlines)
                     False: row-major format output(0:nlines,0:ncolumns) (default)
        strarr       True:  return header as numpy array of strings
                     False: return header as list


        Output
        ------
        Depending on options:
            2D-array of floats if header=False
            String array of file header if header=True
            String vector of file header if header=True and full_header=True


        Restrictions
        ------------
        If header=True then skip is counterintuitive because it is
          actally the number of header rows to be read. This is to
          be able to have the exact same call of the function, once
          with header=False and once with header=True.
        If fill=True, blank lines are not filled but are expected end of file.
        transpose=True has no effect on 1D output such as 1 header line


        Examples
        --------
        >>> # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 2.2 2.3 2.4\\n')
        >>> file.close()

        >>> # Read sample file in different ways
        >>> # data
        >>> from autostring import astr
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=2), 1, pp=True))
        [['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, cskip=1), 1, pp=True))
        [['1.2' '1.3' '1.4']
         ['2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, nc=2, skip=1, cskip=1), 1, pp=True))
        [['1.2' '1.3']
         ['2.2' '2.3']]
        >>> print(astr(fread(filename, nc=[1,3], skip=1), 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> print(astr(fread(filename, nc=1, skip=1), 1, pp=True))
        [['1.1']
         ['2.1']]
        >>> print(astr(fread(filename, nc=1, skip=1, reform=True), 1, pp=True))
        ['1.1' '2.1']

        >>> # header
        >>> print(fread(filename, nc=2, skip=1, header=True))
        ['head1', 'head2']
        >>> print(fread(filename, nc=2, skip=1, header=True, full_header=True))
        ['head1 head2 head3 head4']
        >>> print(fread(filename, nc=1, skip=2, header=True))
        [['head1'], ['1.1']]
        >>> print(fread(filename, nc=1, skip=2, header=True, squeeze=True))
        ['head1', '1.1']
        >>> print(fread(filename, nc=1, skip=2, header=True, strarr=True))
        [['head1']
         ['1.1']]

        >>> # skip blank lines
        >>> file = open(filename, 'a')
        >>> file.writelines('\\n')
        >>> file.writelines('3.1 3.2 3.3 3.4\\n')
        >>> file.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']]

        >>> # skip comment lines
        >>> file = open(filename, 'a')
        >>> file.writelines('# First comment\\n')
        >>> file.writelines('! Second 2 comment\\n')
        >>> file.writelines('4.1 4.2 4.3 4.4\\n')
        >>> file.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, nc=[2], skip_blank=True, comment='#'), 1, pp=True))
        [['1.3']
         ['2.3']
         ['3.3']
         ['2.0']
         ['4.3']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment='#!'), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment=('#','!')), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment=['#','!']), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]

        >>> # fill missing columns
        >>> file = open(filename, 'a')
        >>> file.writelines('5.1 5.2\\n')
        >>> file.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment='#!', fill=True, fill_value=-1), 1, pp=True))
        [[' 1.1' ' 1.2' ' 1.3' ' 1.4']
         [' 2.1' ' 2.2' ' 2.3' ' 2.4']
         [' 3.1' ' 3.2' ' 3.3' ' 3.4']
         [' 4.1' ' 4.2' ' 4.3' ' 4.4']
         [' 5.1' ' 5.2' '-1.0' '-1.0']]

        >>> # transpose
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, transpose=True), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' '2.2']
         ['1.3' '2.3']
         ['1.4' '2.4']]

        >>> # Create some more data with Nan and Inf
        >>> filename1 = 'test1.dat'
        >>> file = open(filename1, 'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 nan Inf "NaN"\\n')
        >>> file.close()

        >>> # Treat Nan and Inf with automatic strip of " and '
        >>> print(astr(fread(filename1, skip=1, transpose=True), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' 'nan']
         ['1.3' 'inf']
         ['1.4' 'nan']]

        >>> # Create some more data with escaped numbers
        >>> filename2 = 'test2.dat'
        >>> file = open(filename2, 'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('"1.1" "1.2" "1.3" "1.4"\\n')
        >>> file.writelines('2.1 nan Inf "NaN"\\n')
        >>> file.close()

        >>> # Strip
        >>> print(astr(fread(filename2,  skip=1,  transpose=True,  strip='"'), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' 'nan']
         ['1.3' 'inf']
         ['1.4' 'nan']]

        >>> # Create some more data with an extra (shorter) header line
        >>> filename3 = 'test3.dat'
        >>> file = open(filename3, 'w')
        >>> file.writelines('Extra header\\n')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 2.2 2.3 2.4\\n')
        >>> file.close()

        >>> print(astr(fread(filename3, skip=2, hskip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(fread(filename3, nc=2, skip=2, hskip=1, header=True))
        ['head1', 'head2']

        >>> # Clean up doctest
        >>> import os
        >>> os.remove(filename)
        >>> os.remove(filename1)
        >>> os.remove(filename2)
        >>> os.remove(filename3)


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

        Copyright 2009-2015 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Feb 2012 - transpose
                  MC, Feb 2013 - ported to Python 3
                  MC, Nov 2014 - bug when nc is list and contains 0
                  MC, Nov 2014 - hskip
                  MC, Feb 2015 - speed: everything list until very end
    """
    #
    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        raise IOError('Cannot open file '+file)
    #
    # Read header and Skip lines
    if hskip > 0:
        ihskip = 0
        while ihskip < hskip:
            tmp = f.readline().rstrip()
            ihskip += 1
    if skip > 0:
        head = ['']*(skip-hskip)
        iskip = 0
        while iskip < (skip-hskip):
            head[iskip] = f.readline().rstrip()
            iskip += 1
    #
    # read first line to determine nc and separator (if not set)
    split = -1
    while True:
        s = f.readline().rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        break
    if separator is None:
        sep = ','
        res = s.split(sep)
        nres = len(res)
        if nres == 1:
            sep = ';'
            res = s.split(sep)
            nres = len(res)
            if nres == 1:
                sep = None
                res = s.split(sep)
                nres = len(res)
    else:
        sep = separator
        res = s.split(sep)
        nres = len(res)
    #
    # Determine indices
    if isinstance(nc, (list, tuple, np.ndarray)):
        nnc = len(nc)
        iinc = nc
    else:
        if nc <= 0:
            nnc = nres-cskip
            iinc = np.arange(nnc, dtype='int') + cskip
        else:
            nnc = nc
            iinc = np.arange(nnc, dtype='int') + cskip
    miinc = max(iinc)
    #
    # Header
    if header:
        # Split header
        var = None
        if (skip-hskip) > 0:
            if full_header:
                var = head
            else:
                var = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miinc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    k += 1
                if (skip-hskip) == 1: var = var[0]
        f.close()
        if strarr:
            var = np.array(var, dtype=np.str)
            if transpose: var = var.T
            if squeeze or reform: var = var.squeeze()
            if fill: var = np.where(var=='', fill_value, var)
        else:
            if fill:
                var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
            if transpose and isinstance(var[0], list):
                var = [list(i) for i in zip(*var)] # transpose
        return var
    #
    # Values - first line
    if (miinc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index: '+s)
    var = list()
    null = line2var(res, var, iinc, strip)
    #
    # Values - rest of file
    for line in f:
        s = line.rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        res = s.split(sep)
        nres = len(res)
        if (miinc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index: '+s)
        null = line2var(res, var, iinc, strip)

    f.close()
    # list -> array
    if fill:
        var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
    var = np.array(var, dtype=np.float)
    if squeeze or reform: var = var.squeeze()
    if transpose: var = var.T

    return var


# Helper for append var with current line already splitted into list
def line2var(res, var, iinc, strip):
    nres = len(res)
    if strip is None:
        tmp = [res[i].strip('"').strip("'") for i in iinc if i < nres]
    elif not strip:
        tmp = [res[i] for i in iinc if i < nres]
    else:
        tmp = [res[i].strip(strip) for i in iinc if i < nres]
    rest = len([ i for i in iinc if i >= nres ])
    if rest > 0:
        tmp.extend(['']*rest)
    var.append(tmp)
    return


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # filename = 'test.dat'
    # file = open(filename,'w')
    # file.writelines('head1 head2 head3 head4\n')
    # file.writelines('1.1 1.2 1.3 1.4\n')
    # file.writelines('2.1 2.2 2.3 2.4\n')
    # file.writelines('\n')
    # file.writelines('3.1 3.2 3.3 3.4\n')
    # file.writelines('# First comment\n')
    # file.writelines('! Second 2 comment\n')
    # file.writelines('4.1 4.2 4.3 4.4\n')
    # file.close()
    # # print(astr(fread(filename, skip=1), 1, pp=True))
    # # print(astr(fread(filename, skip=1, nc=[2], skip_blank=True, comment='#'), 1, pp=T
    #            rue))
    # print(astr(fread(filename, skip=1, skip_blank=True, comment='#!'), 1, pp=True))
