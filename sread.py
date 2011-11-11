#!/usr/bin/env python
from lif import *
import numpy as np

def sread(file, nc=0, skip=0, cskip=0, separator='',
          squeeze=False, reform=False, skip_blank=False, comment='',
          fill=False, fill_value='',
          header=False, full_header=False,
          quiet=False):
    """
        Read strings into string array from a file.
        Lines or columns can be skipped.
        Columns can also be picked specifically.
        Blank (only whitespace) and comment lines can be excluded.
        The header of the file can be read separately.

        This routines is exactly the same as fread but reads
        everything as strings except of floats.
        
        Definition
        ----------
        def sread(file, nc=0, skip=0, cskip=0, separator='',
                  squeeze=False, reform=False, skip_blank=False, comment='',
                  fill=False, fill_value='',
                  header=False, full_header=False,
                  quiet=False):
                     
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
                      and fill=True (default '')
                       
        Options
        -------
        squeeze      True:  2-dim array will be cleaned of degenerated
                            dimension, i.e. results in vector
                     False: array will be two-dimensional as read (default)
        reform       Same as squeeze (looks for squeeze or reform)
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
        quiet        True:  do not show reason if read fails and returns None
                     False: show error for failed read (default)
                            
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
        
        Example
        ----------
        # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('head1 head2 head3 head4\\n')
        >>> file.writelines('1.1 1.2 1.3 1.4\\n')
        >>> file.writelines('2.1 2.2 2.3 2.4\\n')
        >>> file.close()
        
        # Read sample file in different ways
        # data
        >>> sread(filename,skip=1)
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> sread(filename,skip=2)
        ['2.1', '2.2', '2.3', '2.4']
        >>> sread(filename,skip=1,cskip=1)
        [['1.2', '1.3', '1.4'], ['2.2', '2.3', '2.4']]
        >>> sread(filename,nc=2,skip=1,cskip=1)
        [['1.2', '1.3'], ['2.2', '2.3']]
        >>> sread(filename,nc=[1,3],skip=1)
        [['1.2', '1.4'], ['2.2', '2.4']]
        >>> sread(filename,nc=1,skip=1)
        [['1.1'], ['2.1']]
        >>> sread(filename,nc=1,skip=1,reform=True)
        ['1.1', '2.1']

        # header
        >>> sread(filename,nc=2,skip=1,header=True)
        ['head1', 'head2']
        >>> sread(filename,nc=2,skip=1,header=True,full_header=True)
        ['head1 head2 head3 head4']
        >>> sread(filename,nc=1,skip=2,header=True)
        [['head1'], ['1.1']]
        >>> sread(filename,nc=1,skip=2,header=True,squeeze=True)
        ['head1', '1.1']

        # skip blank lines
        >>> file = open(filename,'a')
        >>> file.writelines('\\n')
        >>> file.writelines('3.1 3.2 3.3 3.4\\n')
        >>> file.close()
        >>> sread(filename,skip=1)
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> sread(filename,skip=1,skip_blank=True)
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4']]

        # skip comment lines
        >>> file = open(filename,'a')
        >>> file.writelines('# First\\n')
        >>> file.writelines('! Second second comment\\n')
        >>> file.writelines('4.1 4.2 4.3 4.4\\n')
        >>> file.close()
        >>> sread(filename,skip=1)
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> sread(filename,skip=1,skip_blank=True)
        SREAD: Line has not enough columns to be indexed: # First
        >>> sread(filename,skip=1,skip_blank=True,quiet=True)
        >>> sread(filename,skip=1,skip_blank=True,fill=True,fill_value='M')
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['#', 'First', 'M', 'M'], ['!', 'Second', 'second', 'comment'], ['4.1', '4.2', '4.3', '4.4']]
        >>> sread(filename,skip=1,skip_blank=True,comment='#')
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['!', 'Second', 'second', 'comment'], ['4.1', '4.2', '4.3', '4.4']]
        >>> sread(filename,skip=1,skip_blank=True,comment='#!')
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
        >>> sread(filename,skip=1,skip_blank=True,comment=('#','!'))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
        >>> sread(filename,skip=1,skip_blank=True,comment=['#','!'])
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]

        # Clean up doctest
        >>> import os
        >>> os.remove(filename)

        History
        -------
        Written, MC, Jul. 2009
    """
    #
    # Determine number of lines in file.
    nr = lif(file, skip=skip, noblank=skip_blank, comment=comment)
    if nr <= 0:
        if not quiet:
            print "SREAD: Empty file %s." % file
        return None
    #
    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        if not quiet:
            print "SREAD: Cannot open file %s for reading." % file
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
                                var.append(fill_value)
                        else:
                            if not quiet:
                                print ('SREAD: First header line has not enough '
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
                                    htmp.append(fill_value)
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
        return var
    #
    # Values
    if nr == 1:
        if miinc >= nres:
            if fill:
                var = list()
                m = 0
                for i in iinc:
                    if iinc[i] < nres:
                        var.append(res[i])
                        m += 1
                for i in range(nnc-m):
                    var.append(fill_value)
            else:
                if not quiet:
                    print ('SREAD: First line has not enough '
                           'columns to be indexed: %s' % s)
                f.close()
                return None
        else:
            var = [res[i] for i in iinc]
    else:                    
        var = ['']
        if miinc >= nres:
            if fill:
                tmp = list()
                m = 0
                for i in iinc:
                    if iinc[i] < nres:
                        tmp.append(res[i])
                        m += 1
                for i in range(nnc-m):
                    tmp.append(fill_value)
            else:
                if not quiet:
                    print ('SREAD: Line has not enough '
                           'columns to be indexed: %s' % s)
                f.close()
                return None
        else:
            tmp = [res[i] for i in iinc]
        if (squeeze or reform) and (len(tmp)==1):
            var.extend(tmp)
        else:
            var.append(tmp)
        k = 1
        while k < nr:
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
                    tmp = list()
                    m = 0
                    for i in iinc:
                        if iinc[i] < nres:
                            tmp.append(res[i])
                            m += 1
                    for i in range(nnc-m):
                        tmp.append(fill_value)
                else:
                    if not quiet:
                        print ('SREAD: Line has not enough '
                               'columns to be indexed: %s' % s)
                    f.close()
                    return None
            else:
                tmp = [res[i] for i in iinc]
            if (squeeze or reform) and (len(tmp)==1):
                var.extend(tmp)
            else:
                var.append(tmp)
            k += 1
        del var[0]
    f.close()
    return var

if __name__ == '__main__':
    import doctest
    doctest.testmod()
