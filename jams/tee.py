#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

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
