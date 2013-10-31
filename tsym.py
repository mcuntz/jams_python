#!/usr/bin/env python
from __future__ import print_function
def tsym(name):
    """
        Returns unicodes for common symbols.


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
        deg         Degree symbol
        degree      Degree symbol
        degc        Degree Celcius
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
        >>> print(tsym('deg'))
        \u00B0
        >>> print(tsym('degree'))
        \u00B0
        >>> print(tsym('degC'))
        \u2103
        >>> print(tsym('degreec'))
        \u2103
        >>> print(tsym('degree C'))
        \u2103
        >>> print(tsym('mu'))
        \u00B5
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
                  MC, Mar 2013 - removed raw
    """
    #
    # Define symbol dictionary
    symdict = ({
               'deg'       : '\u00B0',
               'degree'    : '\u00B0',
               'degc'      : '\u2103',
               'degreec'   : '\u2103',
               'degree c'  : '\u2103',
               'mu'        : '\u00B5',
               'peclet'    : '\u2118',
               'permil'    : '\u2030',
               'permille'  : '\u2030',
               'per mil'   : '\u2030',
               'per mille' : '\u2030'
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
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
