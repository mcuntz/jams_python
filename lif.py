#!/usr/bin/env python
import string

def lif(file, noblank=False, comment='', skip=0):
    """
        Counts the numer of lines in a file.
        Blank (only whitespace) and comment lines can be excluded.

        Definition
        ----------
        def lif(file, noblank=False, comment='', skip=0):


        Input
        -----
        file         source file name


        Optional input parameters
        -------------------------
        comment      line gets excluded if first character of line is
                      in comment sequence
                     sequence can be e.g. string, list or tuple
        skip         number of lines to skip at the beginning of
                       the file (default 0)


        Options
        -------
        noblank      excludes all lines that consists only of
                       whitespace characters


        Output
        ------
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
        >>> lif(filename)
        5
        >>> lif(filename,noblank=True)
        4
        >>> lif(filename,comment='#')
        4
        >>> lif(filename,comment='#!')
        3
        >>> lif(filename,comment='#S')
        4
        >>> lif(filename,comment=('#','L'))
        3
        >>> lif(filename,comment=['#','!'])
        3
        >>> lif(filename,comment='#!',noblank=True)
        2
        >>> lif(filename,skip=2)
        3

        # Clean up
        >>> import os
        >>> os.remove(filename)


        History
        -------
        Written, MC, Jul. 2009
    """
    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        print 'LIF: Cannot open file %s for reading.' % file
        return 0
    # Count lines
    count = 0
    if skip > 0:
        iskip = 0
        while iskip < skip:
            l = f.readline()
            iskip += 1
    if noblank and (comment != ''):        # exclude blank, exclude comment
        for line in f:
            l = line.strip(string.whitespace)
            if (l != ''):
                if (l[0] not in comment): count += 1
    elif noblank and (comment == ''):      # exclude blank, include comment
        for line in f:
            l = line.strip(string.whitespace)
            if (l != ''): count += 1
    elif (not noblank) and (comment != ''):# include blank, exclude comment
        for line in f:
            l = line.strip(string.whitespace)
            if (l == ''):
                count += 1
            else:
                if (l[0] not in comment): count += 1
    else:                                  # include blank, include comment
        count = sum(1 for line in f)
    f.close()
    return count

if __name__ == '__main__':
    import doctest
    doctest.testmod()
