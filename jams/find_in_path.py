#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import os
import sys

def find_in_path(filename):
    """
        Find the file named filename in the sys.path.
        Returns the full path name if found, None if not found.


        Definition
        ----------
        def find_in_path(filename):


        Input
        -----
        filename    str; name of searched file

        Output
        ------
        Full path of file. None if not found.


        Examples
        --------
        >>> datei = 'find_in_path.py'
        >>> isdatei = find_in_path(datei)
        >>> if isdatei is None:
        ...     print('No')
        ... else:
        ...     print('Yes')
        Yes

        >>> datei = 'gapfill.py'
        >>> isdatei = find_in_path(datei)
        >>> if isdatei is None:
        ...     print('No')
        ... else:
        ...     print('Yes')
        No


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2013 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Mar 2013
    """
    for dirname in sys.path:
        possible = os.path.join(dirname, filename)
        if os.path.isfile(possible): return possible
    return None


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # datei = 'wiki.pdf'
    # isdatei = find_in_path(datei)
    # if isdatei is None:
    #     print('No')
    # else:
    #     print('Yes')
    # # Yes

    # datei = 'Humor-Sans.ttf'
    # isdatei = find_in_path(datei)
    # if isdatei is None:
    #     print('No')
    # else:
    #     print('Yes')
    # # No
