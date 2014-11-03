#!/usr/bin/env python
from __future__ import print_function

def tee(*args, **kwargs):
    """
        Prints arguments on screen and in file.
        
        Calls print function twice, once with the keyword file and once without, i.e. prints on sys.stdout.


        Definition
        ----------
        def tee(*args, **kwargs):


        Input
        -----
        all arguments of the print function.


        Optional Input
        --------------
        all keyword arguments of the print function.

        The file argument must be an object with a write(string) method.


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
        Written,  MC, Oct 2014
    """
    print(*args, **kwargs)
    if 'file' in kwargs:
        del kwargs['file']
        print(*args, **kwargs)


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
