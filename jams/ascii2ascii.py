#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['ascii2ascii',
           'ascii2en', 'ascii2fr', 'ascii2us', 'ascii2eng',
           'en2ascii', 'fr2ascii', 'us2ascii', 'eng2ascii']

def ascii2ascii(edate, full=False, en=False, fr=False, us=False, eng=False, YY=False):
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
        YY      boolean (default: False)
                Year in input file is 2-digit year. Every year above current year will be taken in 1900.


        Output
        ------
        list/ND-array of date strings in chosen date format (default: ascii)


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(" ".join(ascii2ascii(edate)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(ascii2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

        >>> print(" ".join(ascii2ascii(edate, en=True)))
        2014-11-12 12:00 2015-03-01 17:56:00 1990-12-01 1786-05-04

        >>> print(" ".join(ascii2ascii(edate, en=True, full=True)))
        2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00 1786-05-04 00:00:00

        >>> print(" ".join(ascii2ascii(list(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(ascii2ascii(tuple(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(ascii2ascii(np.array(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(ascii2ascii(edate[0]))
        12.11.2014 12:00

        >>> print(" ".join(ascii2ascii(edate, us=True)))
        11/12/2014 12:00 03/01/2015 17:56:00 12/01/1990 05/04/1786

        >>> print(" ".join(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True)))
        11/12/2014 12:00:00 03/01/2015 17:56:00 12/01/1990 00:00:00 05/04/1786 00:00:00

        >>> print(" ".join(ascii2ascii(edate, fr=True)))
        12/11/2014 12:00 01/03/2015 17:56:00 01/12/1990 04/05/1786

        >>> print(" ".join(ascii2ascii(edate, fr=True, full=True)))
        12/11/2014 12:00:00 01/03/2015 17:56:00 01/12/1990 00:00:00 04/05/1786 00:00:00

        # YY=True
        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> print(" ".join(ascii2ascii(edate, YY=True)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990
        >>> print(" ".join(ascii2ascii(edate, en=True, YY=True)))
        2014-11-12 12:00 2015-03-01 17:56:00 1990-12-01
        >>> print(" ".join(ascii2ascii(edate, us=True, YY=True)))
        11/12/2014 12:00 03/01/2015 17:56:00 12/01/1990
        >>> print(" ".join(ascii2ascii(ascii2ascii(edate, en=True, YY=True), us=True, full=True)))
        11/12/2014 12:00:00 03/01/2015 17:56:00 12/01/1990 00:00:00
        >>> print(" ".join(ascii2ascii(edate, fr=True, full=True, YY=True)))
        12/11/2014 12:00:00 01/03/2015 17:56:00 01/12/1990 00:00:00

        
        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015-2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Feb 2015
        Modified, MC, Sep 2015 - removed date2dec and dec2date
                  MC, Nov 2016 - adapted docstring to Python 2 and 3
                  MC, Mar 2018 - us, eng->en, fr
                  MC, Nov 2018 - YY
                  MC, Jun 2019 - eng->en working again
    """
    if eng: en = True
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
    odate = list()
    if YY:
        import time as ptime
        iyr2 = int(ptime.asctime()[-2:])
        for i, d in enumerate(idate):
            if en:
                if '-' in d:
                    if int(d[0:2]) > iyr2:
                        dd = '19'+d                                  # en -> en
                    else:
                        dd = '20'+d
                elif '/' in d:
                    if int(d[6:8]) > iyr2:
                        dd = '19'+d[6:8]+'-'+d[0:2]+'-'+d[3:5]+d[8:] # us -> en
                    else:
                        dd = '20'+d[6:8]+'-'+d[0:2]+'-'+d[3:5]+d[8:]
                else:
                    if int(d[6:8]) > iyr2:
                        dd = '19'+d[6:8]+'-'+d[3:5]+'-'+d[0:2]+d[8:] # ascii -> en
                    else:
                        dd = '20'+d[6:8]+'-'+d[3:5]+'-'+d[0:2]+d[8:]
            elif fr:
                if '-' in d:
                    if int(d[0:2]) > iyr2:
                        dd = d[6:8]+'/'+d[3:5]+'/19'+d[0:2]+d[8:] # en -> fr
                    else:
                        dd = d[6:8]+'/'+d[3:5]+'/20'+d[0:2]+d[8:]
                elif '/' in d:
                    if int(d[6:8]) > iyr2:
                        dd = d[3:5]+'/'+d[0:2]+'/19'+d[6:8]+d[8:] # us -> fr
                    else:
                        dd = d[3:5]+'/'+d[0:2]+'/20'+d[6:8]+d[8:]
                else:
                    if int(d[6:8]) > iyr2:
                        dd = d[0:2]+'/'+d[3:5]+'/19'+d[6:8]+d[8:] # ascii -> fr
                    else:
                        dd = d[0:2]+'/'+d[3:5]+'/20'+d[6:8]+d[8:]
            elif us:
                if '-' in d:
                    if int(d[0:2]) > iyr2:
                        dd = d[3:5]+'/'+d[6:8]+'/19'+d[0:2]+d[8:] # en -> us
                    else:
                        dd = d[3:5]+'/'+d[6:8]+'/20'+d[0:2]+d[8:]
                elif '/' in d:
                    if int(d[6:8]) > iyr2:
                        dd = d[0:6]+'19'+d[6:]                    # us -> us
                    else:
                        dd = d[0:6]+'20'+d[6:]
                else:
                    if int(d[6:8]) > iyr2:
                        dd = d[3:5]+'/'+d[0:2]+'/19'+d[6:8]+d[8:] # ascii -> us
                    else:
                        dd = d[3:5]+'/'+d[0:2]+'/20'+d[6:8]+d[8:]
            else:
                if '-' in d:
                    if int(d[0:2]) > iyr2:
                        dd = d[6:8]+'.'+d[3:5]+'.19'+d[0:2]+d[8:] # en -> ascii
                    else:
                        dd = d[6:8]+'.'+d[3:5]+'.20'+d[0:2]+d[8:]
                elif '/' in d:
                    if int(d[6:8]) > iyr2:
                        dd = d[3:5]+'.'+d[0:2]+'.19'+d[6:8]+d[8:] # us -> ascii
                    else:
                        dd = d[3:5]+'.'+d[0:2]+'.20'+d[6:8]+d[8:] # us -> ascii
                else:
                    if int(d[6:8]) > iyr2:
                        dd = d[0:6]+'19'+d[6:]                    # ascii -> ascii
                    else:
                        dd = d[0:6]+'20'+d[6:]
            odate.append(dd)
    else:
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
            odate.append(dd)

    if full:
        odate = [ (d+' 00:00:00')[:19] if len(d) < 11 else (d+':00:00')[:19] for d in odate ]


    # Return right type
    if isinstance(edate, list):
        return odate
    elif isinstance(edate, tuple):
        return tuple(odate)
    elif isinstance(edate, np.ndarray):
        return np.array(odate).reshape(edate.shape)
    else:
        return odate[0]


def ascii2en(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with English date format output, i.e. en=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(" ".join(ascii2en(edate)))
        2014-11-12 12:00 2015-03-01 17:56:00 1990-12-01 1786-05-04

        >>> print(" ".join(ascii2en(edate, full=True)))
        2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00 1786-05-04 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> print(" ".join(ascii2en(edate, full=True, YY=True)))
        2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00
    """
    return ascii2ascii(edate, en=True, **kwarg)


def ascii2fr(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with French date format output, i.e. fr=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(" ".join(ascii2fr(edate)))
        12/11/2014 12:00 01/03/2015 17:56:00 01/12/1990 04/05/1786

        >>> print(" ".join(ascii2fr(edate, full=True)))
        12/11/2014 12:00:00 01/03/2015 17:56:00 01/12/1990 00:00:00 04/05/1786 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> print(" ".join(ascii2fr(edate, full=True, YY=True)))
        12/11/2014 12:00:00 01/03/2015 17:56:00 01/12/1990 00:00:00
    """
    return ascii2ascii(edate, fr=True, **kwarg)


def ascii2us(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with American date format output, i.e. us=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(" ".join(ascii2ascii(edate, us=True)))
        11/12/2014 12:00 03/01/2015 17:56:00 12/01/1990 05/04/1786

        >>> print(" ".join(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True)))
        11/12/2014 12:00:00 03/01/2015 17:56:00 12/01/1990 00:00:00 05/04/1786 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> print(" ".join(ascii2ascii(ascii2ascii(edate, en=True, YY=True), us=True, full=True)))
        11/12/2014 12:00:00 03/01/2015 17:56:00 12/01/1990 00:00:00

    """
    return ascii2ascii(edate, us=True, **kwarg)


def ascii2eng(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with English date format output, i.e. en=True.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> print(" ".join(ascii2eng(edate)))
        2014-11-12 12:00 2015-03-01 17:56:00 1990-12-01 1786-05-04

        >>> print(" ".join(ascii2eng(edate, full=True)))
        2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00 1786-05-04 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> print(" ".join(ascii2eng(edate, full=True, YY=True)))
        2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00
    """
    return ascii2ascii(edate, en=True, **kwarg)


def en2ascii(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, en=True)
        >>> print(" ".join(en2ascii(edate)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(en2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> edate = ascii2ascii(edate, en=True, YY=True)
        >>> print(" ".join(en2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00
    """
    return ascii2ascii(edate, **kwarg)


def fr2ascii(edate, full=False, YY=False):
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
        YY      boolean (default: False)
                Year in input file is 2-digit year. Every year above current year will be taken in 1900.


        Output
        ------
        list/ND-array of date strings in ascii format DD.MM.YYYY hh:mm:ss


        Examples
        --------
        >>> edate = ['12/11/2014 12:00', '01/03/2015 17:56:00', '01/12/1990', '04/05/1786']
        >>> print(" ".join(fr2ascii(edate)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(fr2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

        >>> print(" ".join(fr2ascii(list(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(fr2ascii(tuple(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(fr2ascii(np.array(edate))))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(fr2ascii(edate[0]))
        12.11.2014 12:00

        # YY=True
        >>> edate = ['12/11/14 12:00', '01/03/15 17:56:00', '01/12/90']
        >>> print(" ".join(fr2ascii(edate, YY=True)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990

        >>> print(" ".join(fr2ascii(edate, full=True, YY=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Mar 2018
        Modified, MC, Nov 2018 - YY
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

    if YY:
        import time as ptime
        iyr2 = int(ptime.asctime()[-2:])
        odate = [ d[0:6]+'19'+d[6:] if int(d[6:8]) > iyr2 else d[0:6]+'20'+d[6:] for d in odate ] # ascii -> ascii  

    if full:
        odate = [ (d+' 00:00:00')[:19] if len(d) < 11 else (d+':00:00')[:19] for d in odate ]

    # Return right type
    if isinstance(edate, list):
        return odate
    elif isinstance(edate, tuple):
        return tuple(odate)
    elif isinstance(edate, np.ndarray):
        return np.array(odate).reshape(edate.shape)
    else:
        return odate[0]


def us2ascii(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, us=True)
        >>> print(" ".join(us2ascii(edate)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(us2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> edate = ascii2ascii(edate, us=True, YY=True)
        >>> print(" ".join(us2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00
    """
    return ascii2ascii(edate, **kwarg)


def eng2ascii(edate, **kwarg):
    """
        Wrapper function for ascii2ascii with ascii date format output.
        def ascii2ascii(edate, full=False, en=False, us=False):


        Examples
        --------
        >>> edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
        >>> edate = ascii2ascii(edate, en=True)
        >>> print(" ".join(eng2ascii(edate)))
        12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

        >>> print(" ".join(eng2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

        >>> edate = ['14-11-12 12:00', '01.03.15 17:56:00', '90-12-01']
        >>> edate = ascii2ascii(edate, en=True, YY=True)
        >>> print(" ".join(eng2ascii(edate, full=True)))
        12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00
    """
    return ascii2ascii(edate, **kwarg)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # edate = ['2014-11-12 12:00', '01.03.2015 17:56:00', '1990-12-01', '04.05.1786']
    # print(" ".join(ascii2ascii(edate)))
    # # 12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

    # print(" ".join(ascii2ascii(edate, full=True)))
    # # 12.11.2014 12:00:00 01.03.2015 17:56:00 01.12.1990 00:00:00 04.05.1786 00:00:00

    # print(" ".join(ascii2ascii(edate, en=True)))
    # # 2014-11-12 12:00 2015-03-01 17:56:00 1990-12-01 1786-05-04

    # print(" ".join(ascii2ascii(edate, en=True, full=True)))
    # # 2014-11-12 12:00:00 2015-03-01 17:56:00 1990-12-01 00:00:00 1786-05-04 00:00:00

    # print(" ".join(ascii2ascii(list(edate))))
    # # 12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

    # print(" ".join(ascii2ascii(tuple(edate))))
    # # 12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

    # print(" ".join(ascii2ascii(np.array(edate))))
    # # 12.11.2014 12:00 01.03.2015 17:56:00 01.12.1990 04.05.1786

    # print(ascii2ascii(edate[0]))
    # # 12.11.2014 12:00

    # print(" ".join(ascii2ascii(edate, us=True)))
    # # 11/12/2014 12:00 03/01/2015 17:56:00 12/01/1990 05/04/1786

    # print(" ".join(ascii2ascii(ascii2ascii(edate, en=True), us=True, full=True)))
    # # 11/12/2014 12:00:00 03/01/2015 17:56:00 12/01/1990 00:00:00 05/04/1786 00:00:00

    # print(" ".join(ascii2ascii(edate, fr=True)))
    # # 12/11/2014 12:00 01/03/2015 17:56:00 01/12/1990 04/05/1786

    # print(" ".join(ascii2ascii(edate, fr=True, full=True)))
    # # 12/11/2014 12:00:00 01/03/2015 17:56:00 01/12/1990 00:00:00 04/05/1786 00:00:00
