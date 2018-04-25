#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['ascii2ascii',
           'ascii2en', 'ascii2fr', 'ascii2us', 'ascii2eng',
           'en2ascii', 'fr2ascii', 'us2ascii', 'eng2ascii']

def ascii2ascii(edate, full=False, en=False, fr=False, us=False, eng=False):
    """
        Convert date notationas between ascii DD.MM.YYYY hh:mm:ss, English YYYY-MM-DD hh:mm:ss,
        and American MM/DD/YYYY hh:mm:ss.
        Input can only be ascii, English and American. Output can also be French DD/MM/YYYY hh:mm:ss.
        Use fr2ascii first for French input formats.


        Definition
        ----------
        def ascii2ascii(edate, full=False, en=False, fr=False, us=False):


        Input
        -----
        list/ND-array of date strings in ascii, English or American format.


        Optional Input
        --------------
        full    True:  output dates arr all in full format DD.MM.YYYY hh:mm:ss; missing time inputs are 00 on output
                False: output dates are as long as input dates,
                e.g. [YYYY-MM-DD, YYYY-MM-DD hh:mm] gives [DD.MM.YYYY, DD.MM.YYYY hh:mm]
        en      True:  output format is English YYYY-MM-DD hh:mm:ss
                False: output format is ascii DD.MM.YYYY hh:mm:ss (default)
        fr      True:  output format is French DD/MM/YYYY hh:mm:ss
                False: output format is ascii DD.MM.YYYY hh:mm:ss (default)
        us      True:  output format is American MM/DD/YYYY hh:mm:ss
                False: output format is ascii DD.MM.YYYY hh:mm:ss (default)
        eng     Same as en: obsolete.


        Output
        ------
        list/ND-array of date strings in chosen date format (default: ascii)


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(ascii2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']

        >>> print(ascii2ascii(edate, en=True))
        ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

        >>> print(ascii2ascii(edate, en=True, full=True))
        ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']

        >>> print(ascii2ascii(list(edate)))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(ascii2ascii(tuple(edate)))
        ('12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786')

        >>> print(ascii2ascii(np.array(edate)))
        ['12.11.2014 12:00' '01.03.2015 17:56:00' '01.12.1990' '04.05.1786']

        >>> print(ascii2ascii(edate[0]))
        12.11.2014 12:00

        >>> print(ascii2ascii(edate, us=True))
        ['11/12/2014 12:00', '03/01/2015 17:56:00', '12/01/1990', '05/04/1786']

        >>> print(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True))
        ['11/12/2014 12:00:00', '03/01/2015 17:56:00', '12/01/1990 00:00:00', '05/04/1786 00:00:00']

        >>> print(ascii2ascii(edate, fr=True))
        ['12/11/2014 12:00', '01/03/2015 17:56:00', '01/12/1990', '04/05/1786']

        >>> print(ascii2ascii(edate, fr=True, full=True))
        ['12/11/2014 12:00:00', '01/03/2015 17:56:00', '01/12/1990 00:00:00', '04/05/1786 00:00:00']


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014-2018 Matthias Cuntz


        History
        -------
        Written,  MC, Feb 2015
        Modified, MC, Sep 2015 - removed date2dec and dec2date
                  MC, Nov 2016 - adapted docstring to Python 2 and 3
                  MC, Mar 2018 - us, eng->en, fr
    """
    assert (en+fr+us <= 1), 'en, fr and us keywords mutually exclusive.'

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
    odate = np.zeros((ndate,), dtype='|U19') # 'DD.MM.YYYY hh:mm:ss' are 19 characters
    for i, d in enumerate(idate):
        if en:
            if '-' in d:
                dd = d                                    # en -> en
            elif '/' in d:
                dd = d[6:10]+'-'+d[0:2]+'-'+d[3:5]+d[10:] # us -> en
            else:
                dd = d[6:10]+'-'+d[3:5]+'-'+d[0:2]+d[10:] # ascii -> en
        elif fr:
            if '-' in d:
                dd = d[8:10]+'/'+d[5:7]+'/'+d[0:4]+d[10:] # en -> fr
            elif '/' in d:
                dd = d[3:5]+'/'+d[0:2]+'/'+d[6:10]+d[10:] # us -> fr
            else:
                dd = d[0:2]+'/'+d[3:5]+'/'+d[6:10]+d[10:] # ascii -> fr
        elif us:
            if '-' in d:
                dd = d[5:7]+'/'+d[8:10]+'/'+d[0:4]+d[10:] # en -> us
            elif '/' in d:
                dd = d                                    # us -> us
            else:
                dd = d[3:5]+'/'+d[0:2]+'/'+d[6:10]+d[10:] # ascii -> us
        else:
            if '-' in d:
                dd = d[8:10]+'.'+d[5:7]+'.'+d[0:4]+d[10:] # en -> ascii
            elif '/' in d:
                dd = d[3:5]+'.'+d[0:2]+'.'+d[6:10]+d[10:] # us -> ascii
            else:
                dd = d                                    # ascii -> ascii
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


def ascii2en(edate, full=False):
    """
        Wrapper function for ascii2ascii with English date format output, i.e. en=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2en(edate))
        ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

        >>> print(ascii2en(edate, full=True))
        ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']
    """
    return ascii2ascii(edate, full=full, en=True)


def ascii2fr(edate, full=False):
    """
        Wrapper function for ascii2ascii with French date format output, i.e. fr=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2fr(edate))
        ['12/11/2014 12:00', '01/03/2015 17:56:00', '01/12/1990', '04/05/1786']

        >>> print(ascii2fr(edate, full=True))
        ['12/11/2014 12:00:00', '01/03/2015 17:56:00', '01/12/1990 00:00:00', '04/05/1786 00:00:00']
    """
    return ascii2ascii(edate, full=full, fr=True)


def ascii2us(edate, full=False):
    """
        Wrapper function for ascii2ascii with American date format output, i.e. us=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2ascii(edate, us=True))
        ['11/12/2014 12:00', '03/01/2015 17:56:00', '12/01/1990', '05/04/1786']

        >>> print(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True))
        ['11/12/2014 12:00:00', '03/01/2015 17:56:00', '12/01/1990 00:00:00', '05/04/1786 00:00:00']
    """
    return ascii2ascii(edate, full=full, us=True)


def ascii2eng(edate, full=False):
    """
        Wrapper function for ascii2ascii with English date format output, i.e. en=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(ascii2eng(edate))
        ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

        >>> print(ascii2eng(edate, full=True))
        ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']
    """
    return ascii2ascii(edate, full=full, en=True)


def en2ascii(edate, full=False):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, en=True)
        >>> print(en2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(en2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']
    """
    return ascii2ascii(edate, full=full)


def fr2ascii(edate, full=False):
    """
        Convert French date notation DD/MM/YYYY hh:mm:ss to ascii notation DD.MM.YYYY hh:mm:ss.
        Simply replaces / with ., assuring iterable type


        Definition
        ----------
        def fr2ascii(edate, full=False):


        Input
        -----
        list/ND-array of date strings in French format DD/MM/YYYY hh:mm:ss


        Optional Input
        --------------
        full    True:  output dates arr all in full format DD.MM.YYYY hh:mm:ss; missing time inputs are 00 on output
                False: output dates are as long as input dates,
                e.g. [DD/MM/YYYY, DD/MM/YYYY hh:mm] gives [DD.MM.YYYY, DD.MM.YYYY hh:mm]


        Output
        ------
        list/ND-array of date strings in ascii format DD.MM.YYYY hh:mm:ss


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['12/11/2014 12:00', '01/03/2015 17:56:00', '01/12/1990', '04/05/1786']
        >>> print(fr2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(fr2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']

        >>> print(fr2ascii(list(edate)))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(fr2ascii(tuple(edate)))
        ('12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786')

        >>> print(fr2ascii(np.array(edate)))
        ['12.11.2014 12:00' '01.03.2015 17:56:00' '01.12.1990' '04.05.1786']

        >>> print(fr2ascii(edate[0]))
        12.11.2014 12:00


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2018 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2018
    """

    # Input type and shape
    if isinstance(edate, list):
        idate  = edate
    elif isinstance(edate, tuple):
        idate  = list(edate)
    elif isinstance(edate, np.ndarray):
        idate   = list(edate.flatten())
    else:
        idate  = [edate]
    odate = [ d.replace('/','.') for d in idate ]

    if full:
        for i, d in enumerate(odate):
            if len(d) < 11:
                d += ' 00:00:00'
            else:
                d += ':00:00'
            odate[i] = d[:19]

    # Return right type
    if isinstance(edate, list):
        return odate
    elif isinstance(edate, tuple):
        return tuple(odate)
    elif isinstance(edate, np.ndarray):
        return np.array(odate).reshape(edate.shape)
    else:
        return odate[0]


def us2ascii(edate, full=False):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, us=True)
        >>> print(us2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(us2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']
    """
    return ascii2ascii(edate, full=full)


def eng2ascii(edate, full=False):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> import sys
        >>> pyver = sys.version_info
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, en=True)
        >>> print(eng2ascii(edate))
        ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

        >>> print(eng2ascii(edate, full=True))
        ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']
    """
    return ascii2ascii(edate, full=full)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import sys
    # pyver = sys.version_info
    # edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
    # print(ascii2ascii(edate))
    # # ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

    # print(ascii2ascii(edate, full=True))
    # # ['12.11.2014 12:00:00', '01.03.2015 17:56:00', '01.12.1990 00:00:00', '04.05.1786 00:00:00']

    # print(ascii2ascii(edate, en=True))
    # # ['2014-11-12 12:00', '2015-03-01 17:56:00', '1990-12-01', '1786-05-04']

    # print(ascii2ascii(edate, en=True, full=True))
    # # ['2014-11-12 12:00:00', '2015-03-01 17:56:00', '1990-12-01 00:00:00', '1786-05-04 00:00:00']

    # print(ascii2ascii(list(edate)))
    # # ['12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786']

    # print(ascii2ascii(tuple(edate)))
    # # ('12.11.2014 12:00', '01.03.2015 17:56:00', '01.12.1990', '04.05.1786')

    # print(ascii2ascii(np.array(edate)))
    # # ['12.11.2014 12:00' '01.03.2015 17:56:00' '01.12.1990' '04.05.1786']

    # print(ascii2ascii(edate[0]))
    # # 12.11.2014 12:00

    # print(ascii2ascii(edate, us=True))
    # # ['11/12/2014 12:00', '03/01/2015 17:56:00', '12/01/1990', '05/04/1786']

    # print(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True))
    # # ['11/12/2014 12:00:00', '03/01/2015 17:56:00', '12/01/1990 00:00:00', '05/04/1786 00:00:00']

    # print(ascii2ascii(edate, fr=True))
    # # ['12/11/2014 12:00', '01/03/2015 17:56:00', '01/12/1990', '04/05/1786']

    # print(ascii2ascii(edate, fr=True, full=True))
    # # ['12/11/2014 12:00:00', '01/03/2015 17:56:00', '01/12/1990 00:00:00', '04/05/1786 00:00:00']
