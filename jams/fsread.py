#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from jams.fread import fread
from jams.sread import sread

def fsread(file, nc=0, snc=0, skip=0, cskip=0, hskip=0, separator=None,
           squeeze=False, reform=False, skip_blank=False, comment=None,
           fill=False, fill_value=0, sfill_value='', strip=None,
           header=False, full_header=False,
           transpose=False, strarr=False):
    """
        Read from a file numbers into 2D float array as well as characters into 2D string array.


        Definition
        ----------
        def fsread(file, nc=0, snc=0, skip=0, cskip=0, hskip=0, separator=None,
                   squeeze=False, reform=False, skip_blank=False, comment=None,
                   fill=False, fill_value=0, sfill_value='', strip=None,
                   header=False, full_header=False,
                   transpose=False, strarr=False):


        Input
        -----
        file         source file name


        Optional Input Parameters
        -------------------------
        nc           number of columns to be read into float array (default: all (nc<=0))
                     nc can be a vector of column indeces,
                     starting with 0; cskip will be ignored then.
                     if snc!=0: nc must be iterable or -1 for all other columns.
        snc          number of columns to be read into string array (default: none (snc=0))
                     snc can be a vector of column indeces,
                     starting with 0; cskip will be ignored then.
                     if nc!=0: snc must be iterable or -1 for all other columns.
        skip         number of lines to skip at the beginning of file (default: 0)
        cskip        number of columns to skip at the beginning of each line (default: 0)
        hskip        number of lines in skip that do not belong to header (default: 0)
        separator    column separator
                     If not given, columns separator are (in order):
                     comma, semi-colon, whitespace
        comment      line gets excluded if first character of line is in comment sequence
                     sequence can be e.g. string, list or tuple
        fill_value   value to fill float arry in empty cells or if not enough columns per line
                     and fill=True (default 0)
        sfill_value  value to fill in header and in string array in empty cells
                     or if not enough columns per line and fill=True (default '')
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
                            if nc==0 and snc!=0: return strarr
                     False: return header as list
                            if nc==0 and snc!=0: return list


        Output
        ------
        Depending on options:
            2D-array of floats (snc=0)
            list/2D-array of strings (nc=0 and snc!=0)
            2D-array of floats, 2D-array of strings (nc!=0 and snc!=0)
            list/string array of header ((nc=0 or snc=0) and header=True)
            list/string array of header for float array, list/string array of header for strarr ((nc!=0 and snc!=0) and header=True)
            String vector of full file header (header=True and full_header=True)


        Restrictions
        ------------
        If header=True then skip is counterintuitive because it is
          actally the number of header rows to be read. This is to
          be able to have the exact same call of the function, once
          with header=False and once with header=True.
        If fill=True, blank lines are not filled but are taken as end of file.


        Examples
        --------
        >>> # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 2.2 2.3 2.4\\n')
        >>> file.close()

        >>> # Read sample with fread - see fread for more examples
        >>> from autostring import astr
        >>> print(astr(fsread(filename, nc=[1,3], skip=1), 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> print(fsread(filename, nc=2, skip=1, header=True))
        ['head1', 'head2']

        >>> # Read sample with sread - see sread for more examples
        >>> print(fsread(filename, snc=[1,3], skip=1))
        [['1.2', '1.4'], ['2.2', '2.4']]

        >>> # Some mixed data
        >>> file = open(filename,'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('01.12.2012 1.2 name1 1.4\\n')
        >>> file.writelines('01.01.2013 2.2 name2 2.4\\n')
        >>> file.close()

        >>> # Read columns
        >>> print(astr(fsread(filename, nc=[1,3], skip=1), 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> a, sa = fsread(filename, nc=[1,3], snc=[0,2], skip=1)
        >>> print(astr(a, 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> print(sa[0][0])
        01.12.2012
        >>> print(sa[0][1])
        name1
        >>> print(sa[1][0])
        01.01.2013
        >>> print(sa[1][1])
        name2
        >>> a, sa = fsread(filename, nc=[1,3], snc=-1, skip=1)
        >>> print(astr(a, 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> print(sa[0][0])
        01.12.2012
        >>> print(sa[0][1])
        name1
        >>> print(sa[1][0])
        01.01.2013
        >>> print(sa[1][1])
        name2
        >>> a, sa = fsread(filename, nc=-1, snc=[0,2], skip=1)
        >>> print(astr(a, 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]

        >>> # Read header
        >>> a, sa = fsread(filename, nc=[1,3], snc=[0,2], skip=1, header=True)
        >>> print(a)
        ['head2', 'head4']
        >>> print(sa)
        ['head1', 'head3']

        >>> # Some mixed data with missing values
        >>> file = open(filename,'w')
        >>> file.writelines('head1,head2,head3,head4\\n')
        >>> file.writelines('01.12.2012,1.2,name1,1.4\\n')
        >>> file.writelines('01.01.2013,,name2,2.4\\n')
        >>> file.close()

        >>> print(astr(fsread(filename, nc=[1,3], skip=1, fill=True, fill_value=-1), 1, pp=True))
        [[' 1.2' ' 1.4']
         ['-1.0' ' 2.4']]

        >>> # Clean up doctest
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

        Copyright 2009-2015 Matthias Cuntz


        History
        -------
        Written,  MC, Feb 2015 - modified fread
    """

    # Input error
    if (nc <= -1) and (snc <= -1):
        raise ValueError('nc and snc must numbers or list of indices; -1 means read the rest. nc and snc cannot both be -1.')

    # wrap to fread/sread
    if not isinstance(snc, (list, tuple, np.ndarray)):
        if snc==0:
            return fread(file, nc=nc, skip=skip, cskip=cskip, hskip=hskip, separator=separator,
                         squeeze=squeeze, reform=reform, skip_blank=skip_blank, comment=comment,
                         fill=fill, fill_value=fill_value, strip=strip,
                         header=header, full_header=full_header,
                         transpose=transpose, strarr=strarr)
    # snc!=0
    if not isinstance(nc, (list, tuple, np.ndarray)):
        if nc==0:
            return sread(file, nc=snc, skip=skip, cskip=cskip, hskip=hskip, separator=separator,
                         squeeze=squeeze, reform=reform, skip_blank=skip_blank, comment=comment,
                         fill=fill, fill_value=sfill_value, strip=strip,
                         header=header, full_header=full_header,
                         transpose=transpose, strarr=strarr)

    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        raise IOError('Cannot open file '+file)

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

    # Determine indices
    if isinstance(nc, (list, tuple, np.ndarray)) and isinstance(snc, (list, tuple, np.ndarray)):
        if np.in1d(nc, snc, assume_unique=True).any():
            raise ValueError('float and string indices overlap.')
        iinc  = nc
        iisnc = snc
    elif isinstance(nc, (list, tuple, np.ndarray)):
        iinc   = nc
        iirest = list(np.delete(np.arange(nres), iinc))
        if snc <= -1:
            iisnc = iirest
        else:
            iisnc = iirest[:snc]
    elif isinstance(snc, (list, tuple, np.ndarray)):
        iisnc  = snc
        iirest = list(np.delete(np.arange(nres), iisnc))
        if nc <= -1:
            iinc = iirest
        else:
            iinc = iirest[:nc]
    else:
        # cannot be nc=-1 and snc=-1
        if nc <= -1:
            iisnc = range(snc)
            iinc  = range(snc,nres)
        else:
            if snc <= -1:
                iinc = range(nc)
                iisnc  = range(nc,nres)
            else:
                # red snc first then nc
                iisnc = range(snc)
                iinc  = range(snc,snc+nc)
    nnc    = len(iinc)
    nsnc   = len(iisnc)
    miinc  = max(iinc)
    miisnc = max(iisnc)
    aiinc  = list()
    aiinc.extend(iinc)
    aiinc.extend(iisnc)
    miianc = max(aiinc)
    #
    # Header
    if header:
        # Split header
        var = None
        if (skip-hskip) > 0:
            if full_header:
                var = head
            else:
                var  = list()
                svar = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miianc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    null = line2var(hres, svar, iisnc, False if strip is None else strip)
                    k += 1
                if (skip-hskip) == 1:
                    var  = var[0]
                    svar = svar[0]
        f.close()
        if strarr:
            var  = np.array(var, dtype=np.str)
            svar = np.array(svar, dtype=np.str)
            if transpose:
                var  = var.T
                svar = svar.T
            if squeeze or reform:
                var  = var.squeeze()
                svar = svar.squeeze()
            if fill:
                var  = np.where(var=='', fill_value, var)
                svar = np.where(svar=='', sfill_value, svar)
        else:
            if fill:
                var  = [ [ fill_value  if i=='' else i for i in row ] for row in var ]
                svar = [ [ sfill_value if i=='' else i for i in row ] for row in svar ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
                maxi = max([ len(i) for i in svar])
                if maxi==1: svar = [ i[0] for i in svar ]
            if transpose and isinstance(var[0], list):
                var = [list(i) for i in zip(*var)] # transpose
            if transpose and isinstance(svar[0], list):
                svar = [list(i) for i in zip(*svar)] # transpose
        return var, svar
    #
    # Values - first line
    if (miianc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index: '+s)
    var  = list()
    svar = list()
    null = line2var(res, var, iinc, strip)
    null = line2var(res, svar, iisnc, False if strip is None else strip)
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
        if (miianc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index: '+s)
        null = line2var(res, var, iinc, strip)
        null = line2var(res, svar, iisnc, False if strip is None else strip)
    f.close()
    # list -> array
    if fill:
        var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
    # var  = np.array(var, dtype=np.float)
    # svar = np.array(svar, dtype=np.str)
    # if fill:
    #     svar = np.where(svar=='', sfill_value, svar)
    # if squeeze or reform:
    #     var  = var.squeeze()
    #     svar = svar.squeeze()
    # if transpose:
    #     var  = var.T
    #     svar = svar.T
    var  = np.array(var, dtype=np.float)
    if squeeze or reform: var  = var.squeeze()
    if transpose: var  = var.T
    if strarr:
        svar = np.array(svar, dtype=np.str)
        if transpose: svar = svar.T
        if squeeze or reform: svar = svar.squeeze()
        if fill: svar = np.where(svar=='', sfill_value, svar)
    else:
        if fill:
            svar = [ [ fill_value if i=='' else i for i in row ] for row in svar ]
        if squeeze or reform:
            maxi = max([ len(i) for i in svar])
            if maxi==1: svar = [ i[0] for i in svar ]
        if transpose and isinstance(svar[0], list):
            svar = [list(i) for i in zip(*svar)] # transpose

    return var, svar


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
    if rest > 0: tmp.extend(['']*rest)
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

    # # Some mixed data with missing values
    # file = open(filename,'w')
    # file.writelines('head1,head2,head3,head4\n')
    # file.writelines('01.12.2012,1.2,name1,1.4\n')
    # file.writelines('01.01.2013,,name2,2.4\n')
    # file.close()

    # print(astr(fsread(filename, nc=[1,3], skip=1, fill=True, fill_value=-1), 1, pp=True))
