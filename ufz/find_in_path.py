#!/usr/bin/env python
from __future__ import print_function
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
        >>> datei = 'Humor-Sans.ttf'
        >>> isdatei = find_in_path(datei)
        >>> if isdatei is None:
        ...     print('No')
        ... else:
        ...     print('Yes')
        Yes

        >>> datei = 'Humor_Sans.ttf'
        >>> isdatei = find_in_path(datei)
        >>> if isdatei is None:
        ...     print('No')
        ... else:
        ...     print('Yes')
        No


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

        Copyright 2013 Matthias Cuntz


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

    # datei = 'Humor-Sans.ttf'
    # isdatei = find_in_path(datei)
    # if isdatei is None:
    #     print('No')
    # else:
    #     print('Yes')
    # # Yes

    # datei = 'Humor_Sans.ttf'
    # isdatei = find_in_path(datei)
    # if isdatei is None:
    #     print('No')
    # else:
    #     print('Yes')
    # # No
