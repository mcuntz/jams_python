#!/usr/bin/env python
from __future__ import print_function

__all__ = ['tee']

def tee(*args, **kwargs):
    """
        Prints arguments on screen and in file.
        
        Calls print function twice, once with the keyword file and once without, i.e. prints on sys.stdout.


        Definition
        ----------
        def tee(*args, **kwargs):


        Input
        -----
        All arguments of the print function.


        Optional Input
        --------------
        All keyword arguments of the print function.

        The file argument must be an object with a write(string) method.
        If the file argument is None, then there will be no second print.


        Output
        ------
        Print on sys.stdout (screen) and in file if given.


        Restrictions
        ------------
        If the keyword file is not given, arguments and keyword arguments are passed
        to the print function only once, i.e. tee is then just a wrapper for print.


        Examples
        --------
        >>> st = 'Output on screen and in log file'
        >>> ff = 'tee_log.txt'

        # write
        >>> f = open(ff, 'w')
        >>> tee(st, file=f)
        Output on screen and in log file
        >>> f.close()

        # test
        >>> f = open(ff, 'r')
        >>> test = f.readline()
        >>> f.close()
        >>> test = test[:-1] # rm trailing newline character
        >>> if test == st:
        ...     print('Yes')
        ... else:
        ...     print('No')
        Yes

        >>> import os
        >>> os.remove(ff)

        >>> f=None
        >>> st = 'Output only on screen'
        >>> tee(st, file=f)
        Output only on screen


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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Oct 2014
        Modified, Mc, Nov 2016 - file=None
    """
    if 'file' in kwargs:
        if kwargs['file'] is not None:
            print(*args, **kwargs) # file
        del kwargs['file']
    print(*args, **kwargs)         # screen


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # st = 'Output on screen and in log file'
    # ff = 'tee_log.txt'
    # # write
    # f = open(ff, 'w')
    # tee(st, file=f)
    # f.close()

    # # test
    # f = open(ff, 'r')
    # test = f.readline()
    # f.close()
    # # rm trailing newline character
    # test = test[:-1]
    # if test == st:
    #     print('Yes')
    # else:
    #     print('No')

    # import os
    # os.remove(ff)
