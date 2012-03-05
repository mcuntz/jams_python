#!/usr/bin/env python
import numpy as np
from lif import * # from ufz

def fread(file, nc=0, skip=0, cskip=0, separator='',
          squeeze=False, reform=False, skip_blank=False, comment='',
          fill=False, fill_value=0,
          header=False, full_header=False,
          quiet=False, transpose=False, strarr=False):
    """
        Read numbers into float array from a file.
        Lines or columns can be skipped.
        Columns can also be picked specifically.
        Blank (only whitespace) and comment lines can be excluded.
        The header of the file can be read separately.

        Definition
        ----------
        def fread(file, nc=0, skip=0, cskip=0, separator='',
                  squeeze=False, reform=False, skip_blank=False, comment='',
                  fill=False, fill_value=0,
                  header=False, full_header=False,
                  quiet=False, transpose=False, strarr=False):


        Input
        -----
        file         source file name


        Optional Input Parameters
        -------------------------
        nc           number of columns to be read (default: all)
                     nc can be a vector of column indeces,
                     starting with 0; cskip will be ignored then.
        separator    columns separator
                     If not given, columns separator are (in order):
                     comma, semi-colon, whitespace
        skip         number of lines to skip at the beginning of
                     the file (default 0)
        cskip        number of columns to skip at the beginning of
                     each line (default 0)
        comment      line gets excluded if first character of line is
                     in comment sequence
                     sequence can be e.g. string, list or tuple
        fill_value   value to fill in if not enough columns line
                     and fill=True (default 0 and '' for header)


        Options
        -------
        squeeze      True:  2-dim array will be cleaned of degenerated
                            dimension, i.e. results in vector
                     False: array will be two-dimensional as read (default)
        reform       Same as squeeze (looks for 'squeeze or reform')
        skip_blank   True:  continues reading after blank line
                     False: stops reading at first blank line (default)
        fill         True:  fills in fill_value if not enough columns in line
                     False: stops execution and returns None if not enough
                            columns in line (default)
        header       True:  header strings will be returned
                     False  numbers in file will be returned (default)
        full_header  True:  header is a string vector of the skipped rows
                     False: header will be split in columns, exactly as the
                            data, and will hold only the selected columns (default)
        quiet        True:  do not show reason if read fails and returns None
                     False: show error for failed read (default)
        transpose    True:  column-major format output(0:ncolumns,0:nlines)
                     False: row-major format output(0:nlines,0:ncolumns) (default)
        strarr       True:  return header as numpy array of strings
                     False: return header as list


        Output
        ------
        Either float array containing numbers from file if header=False
        or string array of file header if header=True
        or string vector of file header if header=True, full_header=True


        Restrictions
        ------------
        If header=True then skip is counterintuitive because it is
          actally the number of header rows to be read. This is to
          be able to have the exact same call of the function, once
          with header=False and once with header=True.
        If fill=True, blank lines are not filled but are expected
          end of file.
        Tested with python >= 2.5


        Examples
        --------
        # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 2.2 2.3 2.4\\n')
        >>> file.close()

        # Read sample file in different ways
        # data
        >>> fread(filename,skip=1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=2)
        array([[ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=1,cskip=1)
        array([[ 1.2,  1.3,  1.4],
               [ 2.2,  2.3,  2.4]])
        >>> fread(filename,nc=2,skip=1,cskip=1)
        array([[ 1.2,  1.3],
               [ 2.2,  2.3]])
        >>> fread(filename,nc=[1,3],skip=1)
        array([[ 1.2,  1.4],
               [ 2.2,  2.4]])
        >>> fread(filename,nc=1,skip=1)
        array([[ 1.1],
               [ 2.1]])
        >>> fread(filename,nc=1,skip=1,reform=True)
        array([ 1.1,  2.1])

        # header
        >>> fread(filename,nc=2,skip=1,header=True)
        ['head1', 'head2']
        >>> fread(filename,nc=2,skip=1,header=True,full_header=True)
        ['head1 head2 head3 head4']
        >>> fread(filename,nc=1,skip=2,header=True)
        [['head1'], ['1.1']]
        >>> fread(filename,nc=1,skip=2,header=True,squeeze=True)
        ['head1', '1.1']
        >>> fread(filename,nc=1,skip=2,header=True,strarr=True)
        array([['head1'],
               ['1.1']], 
              dtype='|S5')

        # skip blank lines
        >>> file = open(filename,'a')
        >>> file.writelines('\\n')
        >>> file.writelines('3.1 3.2 3.3 3.4\\n')
        >>> file.close()
        >>> fread(filename,skip=1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=1,skip_blank=True)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4],
               [ 3.1,  3.2,  3.3,  3.4]])

        # skip comment lines
        >>> file = open(filename,'a')
        >>> file.writelines('# First comment\\n')
        >>> file.writelines('! Second 2 comment\\n')
        >>> file.writelines('4.1 4.2 4.3 4.4\\n')
        >>> file.close()
        >>> fread(filename,skip=1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=1,skip_blank=True)
        FREAD: Line has not enough columns to be indexed: # First comment
        >>> fread(filename,skip=1,skip_blank=True,quiet=True)
        >>> fread(filename,skip=1,skip_blank=True,comment='#')
        FREAD: Requested elements not all numbers: ['!', 'Second', '2', 'comment']
        >>> fread(filename,skip=1,nc=[2],skip_blank=True,comment='#')
        array([[ 1.3],
               [ 2.3],
               [ 3.3],
               [ 2. ],
               [ 4.3]])
        >>> fread(filename,skip=1,skip_blank=True,comment='#!')
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4],
               [ 3.1,  3.2,  3.3,  3.4],
               [ 4.1,  4.2,  4.3,  4.4]])
        >>> fread(filename,skip=1,skip_blank=True,comment=('#','!'))
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4],
               [ 3.1,  3.2,  3.3,  3.4],
               [ 4.1,  4.2,  4.3,  4.4]])
        >>> fread(filename,skip=1,skip_blank=True,comment=['#','!'])
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4],
               [ 3.1,  3.2,  3.3,  3.4],
               [ 4.1,  4.2,  4.3,  4.4]])

        # fill missing columns
        >>> file = open(filename,'a')
        >>> file.writelines('5.1 5.2\\n')
        >>> file.close()
        >>> fread(filename,skip=1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=1,skip_blank=True,comment='#!')
        FREAD: Line has not enough columns to be indexed: 5.1 5.2
        >>> fread(filename,skip=1,skip_blank=True,comment='#!',fill=True,fill_value=-1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4],
               [ 3.1,  3.2,  3.3,  3.4],
               [ 4.1,  4.2,  4.3,  4.4],
               [ 5.1,  5.2, -1. , -1. ]])

        # transpose
        >>> fread(filename,skip=1)
        array([[ 1.1,  1.2,  1.3,  1.4],
               [ 2.1,  2.2,  2.3,  2.4]])
        >>> fread(filename,skip=1,transpose=True)
        array([[ 1.1,  2.1],
               [ 1.2,  2.2],
               [ 1.3,  2.3],
               [ 1.4,  2.4]])

        # Clean up doctest
        >>> import os
        >>> os.remove(filename)


        History
        -------
        Written, MC, Jul 2009
                 MC, Apr 2011 - transpose
    """
    #
    # Determine number of lines in file.
    nr = lif(file, skip=skip, noblank=skip_blank, comment=comment)
    if nr <= 0:
        if not quiet:
            print "FREAD: Empty file %s." % file
        return None
    #
    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        if not quiet:
            print "FREAD: Cannot open file %s for reading." % file
        return None
    #
    # Read header and Skip lines
    count = 0
    if skip > 0:
        head = ['']*skip
        iskip = 0
        while iskip < skip:
            head[iskip] = f.readline().rstrip()
            iskip += 1
    #
    # read first line to determine nc and separator (if not set)
    split = -1
    while True:
        s = f.readline().rstrip()
        if comment != '':
            if s != '':
                if s[0] not in comment:
                    break
                else:
                    continue
        if skip_blank:
            if s != '':
                break
            else:
                continue
        break
    if separator == '':
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
    count += 1
    #
    # Determine indeces
    if nc == 0:
        nnc = nres-cskip
        iinc = np.arange(nnc, dtype='int') + cskip
    else:
        if type(nc) == type(0):
            nnc = nc
            iinc = np.arange(nnc, dtype='int') + cskip
        else:
            nnc = len(nc)
            iinc = nc
    miinc = max(iinc)
    #
    # Header
    if header:
        # Split header
        if fill_value == 0:
            fillit = ''
        else:
            fillit = fill_value
        var = None
        if skip > 0:
            if full_header:
                var = head
            else:
                if skip == 1:
                    hres = head[0].split(sep)
                    nhres = len(hres)
                    if miinc >= nhres:
                        if fill:
                            var = list()
                            m = 0
                            for i in iinc:
                                if iinc[i] < nhres:
                                    var.append(hres[i])
                                    m += 1
                            for i in range(nnc-m):
                                var.append(fillit)
                        else:
                            if not quiet:
                                print ('FREAD: First header line has not enough '
                                       'columns to be indexed: %s' % head[0])
                            f.close()
                            return None
                    else:
                        var = [hres[i] for i in iinc]
                else:
                    var = ['']
                    k = 0
                    while k < skip:
                        hres = head[k].split(sep)
                        nhres = len(hres)
                        if miinc >= nhres:
                            if fill:
                                htmp = list()
                                m = 0
                                for i in iinc:
                                    if iinc[i] < nhres:
                                        htmp.append(hres[i])
                                        m += 1
                                for i in range(nnc-m):
                                    htmp.append(fillit)
                            else:
                                if not quiet:
                                    print ('FREAD: Header line has not enough '
                                           'columns to be indexed: %s' % head[k])
                                f.close()
                                return None
                        else:
                            htmp = [hres[i] for i in iinc]
                        if (squeeze or reform) and (len(htmp)==1):
                            var.extend(htmp)
                        else:
                            var.append(htmp)
                        k += 1
                    del var[0]
        f.close()
        if transpose:
            var = np.transpose(var, tuple(reversed(range(np.ndim(var)))))
        if strarr:
            var = np.array(var,dtype=np.str)
        return var
    #
    # Values
    iline = 0
    var = np.empty([nr,nnc], dtype='float')
    if miinc >= nres:
        if fill:
            m = 0
            for i in iinc:
                if iinc[i] < nres:
                    var[iline,m] = np.float(res[i])
                    m += 1
            var[iline,m:miinc+1] = fill_value
        else:
            if not quiet:
                print 'FREAD: First line has not enough columns to be indexed: %s' % s
            f.close()
            return None
    else:
        var[iline,0:nnc] = np.array([float(res[i]) for i in iinc],
                                    dtype='float')
    # Read rest of file
    while True:
        s = f.readline().rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment != '':
            if (s[0] in comment): continue
        res = s.split(sep)
        nres = len(res)
        if miinc >= nres:
            if fill:
                iline += 1
                m = 0
                for i in iinc:
                    if iinc[i] < nres:
                        try:
                            var[iline,m] = np.float(res[i])
                        except ValueError:
                            if not quiet:
                                print 'FREAD: Tried to convert "%s"  from Line: %s' % (res[i], s)
                            f.close()
                            return None
                        m += 1
                var[iline,m:miinc+1] = fill_value
            else:
                if not quiet:
                    print 'FREAD: Line has not enough columns to be indexed: %s' % s
                f.close()
                return None
        else:
            iline += 1
            try:
                var[iline,0:nnc] = np.array([float(res[i]) for i in iinc],
                                            dtype='float')
            except ValueError:
                if not quiet:
                    print 'FREAD: Requested elements not all numbers: %s' % ([res[i] for i in iinc])
                f.close()
                return None
        count += 1
        if count == nr: break
    # remove allocated elements if skip_blank=True or comment!=''
    var = var[0:iline+1,0:nnc]
    if transpose:
        var = np.transpose(var, tuple(reversed(range(np.ndim(var)))))
    if (squeeze or reform):
        var = np.squeeze(var)
    f.close()
    return var

if __name__ == '__main__':
    import doctest
    doctest.testmod()
