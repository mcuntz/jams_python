#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

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
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2011-2013 Matthias Cuntz - mc (at) macu (dot) de

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
