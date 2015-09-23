#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['ascii2ascii', 'eng2ascii', 'ascii2eng']

def ascii2ascii(edate, full=False, eng=False):
    """
        Convert date notationas between English YYYY-MM-DD hh:mm:ss and ascii DD.MM.YYYY hh:mm:ss.


        Definition
        ----------
        def ascii2ascii(edate, full=False, eng=False, ascii=None):


        Input
        -----
        list/ND-array of date strings


        Optional Input
        --------------
        full    True:  output dates arr all in full format DD.MM.YYYY hh:mm:ss; missing time inputs are 00 on output
                False: output dates are as long as input dates,
                e.g. [YYYY-MM-DD, YYYY-MM-DD hh:mm] gives [DD.MM.YYYY, DD.MM.YYYY hh:mm]
        eng     True:  output format is English YYYY-MM-DD hh:mm:ss
                False: output format is ascii DD.MM.YYYY hh:mm:ss (default)


        Output
        ------
        list/ND-array of date strings in chosen date format (default: ascii)


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(ascii2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']

        >>> print(ascii2ascii(edate, eng=True))
        ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

        >>> print(ascii2ascii(edate, eng=True, full=True))
        ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']

        >>> print(ascii2ascii(list(edate)))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(ascii2ascii(tuple(edate)))
        ('12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786')

        >>> print(ascii2ascii(np.array(edate)))
        ['12.11.2014 12:00' '01.03.2015 17:56:00' '01.12.1990' '04.05.1786']

        >>> print(ascii2ascii(edate[0]))
        12.11.2014 12:00


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
        Written,  MC, Feb 2015
        Modified, MC, Sep 2015 - removed date2dec and dec2date
    """

    # Input type and shape
    if isinstance(edate, list):
        idate  = np.array(edate)
    elif isinstance(edate, tuple):
        idate  = np.array(edate)
    elif isinstance(edate, np.ndarray):
        idate   = edate.flatten()
    else:
        idate  = np.array([edate])
    ndate = idate.size

    # Convert to given output type
    odate = np.zeros((ndate,), dtype='|S19') # 'DD.MM.YYYY hh:mm:ss' are 19 characters
    for i, d in enumerate(idate):
        if eng:
            if '-' in d:
                dd = d
            else:
                dd = d[6:10]+'-'+d[3:5]+'-'+d[0:2]+d[10:] # ascii -> eng
        else:
            if '-' not in d:
                dd = d
            else:
                dd = d[8:10]+'.'+d[5:7]+'.'+d[0:4]+d[10:] # eng -> ascii
        if not full:
            odate[i] = dd
        else:
            if len(d) < 11:
                dd += ' 00:00:00'
            else:
                dd += ':00:00'
            odate[i] = dd[:19]

    # Return right type
    if isinstance(edate, list):
        return list(odate)
    elif isinstance(edate, tuple):
        return tuple(odate)
    elif isinstance(edate, np.ndarray):
        return odate.reshape(edate.shape)
    else:
        return odate[0]


def eng2ascii(edate, full=False):
    """
        Wrapper function for ascii2ascii with ascii date format output, i.e. eng=False.
        def ascii2ascii(edate, full=False, eng=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(eng2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(eng2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']
    """
    return ascii2ascii(edate, full=full, eng=False)


def ascii2eng(edate, full=False):
    """
        Wrapper function for ascii2ascii with english date format output, i.e. eng=True.
        def ascii2ascii(edate, full=False, eng=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2eng(edate))
        ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

        >>> print(ascii2eng(edate, full=True))
        ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']
    """
    return ascii2ascii(edate, full=full, eng=True)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
    # print(edate)
    # print(ascii2ascii(edate))
    # print(ascii2ascii(edate, full=True))
    # print(ascii2ascii(edate, eng=True))
    # print(ascii2ascii(edate, eng=True, full=True))
    # print(type(ascii2ascii(edate)))
    # print(type(ascii2ascii(list(edate))))
    # print(type(ascii2ascii(tuple(edate))))
    # print(type(ascii2ascii(np.array(edate))))
    # print(type(ascii2ascii(edate[0])))
    # print(ascii2ascii(edate))
    # print(ascii2ascii(list(edate)))
    # print(ascii2ascii(tuple(edate)))
    # print(ascii2ascii(np.array(edate)))
    # print(ascii2ascii(edate[0]))
