#!/usr/bin/env python
from __future__ import print_function
def tsym(name):
    """
        Returns raw unicodes for common symbols.

        Definition
        ----------
        def tsym(name):


        Input
        -----
        name        Symbol name, case insensitive


        Output
        ------
        Raw unicode string


        Restrictions
        ------------
        Only knows a very limited number of names
        Use empty function to return names


        Known symbols
        ------------
        degree      Degree symbol
        degreec     Degree Celcius
        degree c    Degree Celcius
        mu          Lowercase greek mu
        peclet      Peclet symbol
        permil      Permille sign
        permille    Permille sign
        per mil     Permille sign
        per mille   Permille sign


        Examples
        --------
        # \00 can become \\x in docstring
        >>> import sys
        >>> if sys.hexversion > int('0x3000000',base=16):
        ...     if tsym('degree') == r'\\u00B0': print('Good')
        ... else:
        ...     if tsym('degree') == r'\u00B0': print('Good')
        Good
        >>> print(tsym('degreec'))
        \u2103
        >>> print(tsym('degree C'))
        \u2103
        >>> if sys.hexversion > int('0x3000000',base=16):
        ...     if tsym('mu') == r'\\u00B5': print('Good')
        ... else:
        ...     if tsym('mu') == r'\u00B5': print('Good')
        Good
        >>> print(tsym('peclet'))
        \u2118
        >>> print(tsym('permil'))
        \u2030
        >>> print(tsym('Permil'))
        \u2030
        >>> print(tsym('permille'))
        \u2030
        >>> print(tsym('per mil'))
        \u2030
        >>> print(tsym('per mille'))
        \u2030


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2011
        Modified, MC, Feb 2013 - ported to Python 3
    """
    #
    # Define symbol dictionary
    symdict = ({
               'degree'    : r'\u00B0',
               'degreec'   : r'\u2103',
               'degree c'  : r'\u2103',
               'mu'        : r'\u00B5',
               'peclet'    : r'\u2118',
               'permil'    : r'\u2030',
               'permille'  : r'\u2030',
               'per mil'   : r'\u2030',
               'per mille' : r'\u2030' 
              })

    #
    # lookup symbol
    try:
        out = symdict[name.lower()]
    except KeyError:
        print("TSYM: Symbol not known: %s" % name)
        return None

    return out

if __name__ == '__main__':
    import doctest
    doctest.testmod()
