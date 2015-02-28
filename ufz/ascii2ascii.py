#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.date2dec import date2dec
from ufz.dec2date import dec2date

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
    """
    
    # Input size and shape
    islist   = False
    istuple  = False
    isarray  = False
    isscalar = False
    if isinstance(edate, list):
        islist = True
        idate  = np.array(edate)
    elif isinstance(edate, tuple):
        istuple = True
        idate  = np.array(edate)
    elif isinstance(edate, np.ndarray):
        isarray = True
        ishape  = edate.shape
        idate   = edate.flatten()
    else:
        isscalar = True
        idate  = np.array([edate])
    ndate = idate.size
        
    # Search eng and ascii dates
    iieng = [ i for i, d in enumerate(idate) if '-' in d ]
    ieng  = [ d for d in idate if '-' in d ]
    iiasc = [ i for i, d in enumerate(idate) if '-' not in d ]
    iasc  = [ d for d in idate if '-' not in d ]

    # Copy ascii to ascii and transform eng to ascii
    odate = np.zeros((ndate,), dtype='|S19')
    if eng:
        if len(iieng) > 0:
            odate[iieng] = dec2date(date2dec(eng=ieng),   eng=True)
        if len(iiasc) > 0:
            odate[iiasc] = dec2date(date2dec(ascii=iasc), eng=True)
    else:
        if len(iieng) > 0:
            odate[iieng] = dec2date(date2dec(eng=ieng),   ascii=True)
        if len(iiasc) > 0:
            odate[iiasc] = dec2date(date2dec(ascii=iasc), ascii=True)
    
    # Cut to input lengths
    if not full:
        for i, d in enumerate(odate):
            ni = len(idate[i])
            no = len(d)
            if ni != no:
                odate[i] = d[:ni]
                
    # Return
    if islist:
        return list(odate)
    elif istuple:
        return tuple(odate)
    elif isarray:
        return odate.reshape(ishape)
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
