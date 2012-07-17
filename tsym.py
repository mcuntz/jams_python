#!/usr/bin/env python
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
        # \00 becomes \\x in docstring
        >>> tsym('degree')
        u'\\xb0'
        >>> tsym('degreec')
        u'\u2103'
        >>> tsym('degree C')
        u'\u2103'
        >>> tsym('mu')
        u'\\xb5'
        >>> tsym('peclet')
        u'\u2118'
        >>> tsym('permil')
        u'\u2030'
        >>> tsym('Permil')
        u'\u2030'
        >>> tsym('permille')
        u'\u2030'
        >>> tsym('per mil')
        u'\u2030'
        >>> tsym('per mille')
        u'\u2030'


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

        Copyright 2011 Matthias Cuntz


        History
        -------
        Written, MC, Jun 2011
    """
    #
    # Define symbol dictionary
    symdict = ({
               'degree'    : ur'\u00B0',
               'degreec'   : ur'\u2103',
               'degree c'  : ur'\u2103',
               'mu'        : ur'\u00B5',
               'peclet'    : ur'\u2118',
               'permil'    : ur'\u2030',
               'permille'  : ur'\u2030',
               'per mil'   : ur'\u2030',
               'per mille' : ur'\u2030' 
              })

    #
    # lookup symbol
    try:
        out = symdict[name.lower()]
    except KeyError:
        print "TSYM: Symbol not known: %s" % name
        return None

    return out

if __name__ == '__main__':
    import doctest
    doctest.testmod()
