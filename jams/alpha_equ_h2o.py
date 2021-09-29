#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['alpha_equ_h2o']


def alpha_equ_h2o(temp, isotope=None, undef=-9999., eps=False, greater1=True):
    """
    Calculates fractionation factors during liquid-water vapour
    equilibration at temperature temp [K].
    It does not use the atmospheric convention, i.e. factor<1, but defaults to
    >1 (greater1=True).


    Definition
    ----------
    def alpha_equ_h2o(temp, isotope=None, undef=-9999., eps=False,
                      greater1=True):


    Input
    -----
    temp       temperature (array) [K]


    Optional Input
    --------------
    isotope
        which water isotopologue: 1: HDO; 2: H218O; else: no fractionation,
        return 1 (default)
    undef
        exclude from calculations (default: -9999)
    eps
        reports epsilon=alpha-1 instead of alpha
    greater1
        not atmospheric convention, i.e. alpha > 1 (default)


    Output
    ------
    Equilibrium fractionation factors


    Restrictions
    ------------
    None


    Literature
    ----------
    Majoube, M. (1971)
        Fractionnement en oxygene-18 entre la glace et la vapeur d'eau
        Journal De Chimie Physique Et De Physico-Chimie Biologique,
        68(4), 625-636.


    Examples
    --------
    >>> from autostring import astr
    >>> T0 = 273.15
    >>> T  = np.array([0, 10., 15., 25.])
    >>> print(astr(alpha_equ_h2o(T+T0, isotope=0), 4))
    ['1.0000' '1.0000' '1.0000' '1.0000']

    >>> print(astr(alpha_equ_h2o(T+T0, isotope=2), 4))
    ['1.0117' '1.0107' '1.0102' '1.0094']

    >>> print(astr(alpha_equ_h2o(np.ma.array(T+T0, mask=(T==0.)), isotope=2, greater1=False), 4))
    ['-9999.0000' '    0.9894' '    0.9899' '    0.9907']

    >>> print(astr(alpha_equ_h2o(T+T0, isotope=1, eps=True)*1000., 4))
    ['112.3194' ' 97.6829' ' 91.1296' ' 79.3443']

    >>> print(astr(alpha_equ_h2o(0.+T0, isotope=2, eps=True)*1000., 4))
    11.7187


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python
    library, Department of Computational Hydrosystems, Helmholtz Centre for
    Environmental Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Matthias Cuntz, Sep 2014
    Modified, Matthias Cuntz, Sep 2021 - code refactoring
    """
    # Check scalar, list or array
    islist = False
    if isinstance(temp, list):
        islist = True
        temp = np.array(temp)
    isarr = np.ndim(temp)
    if (isarr == 0):  # scalar
        mtemp = temp
    else:
        # mask undefined
        mtemp = np.ma.array(temp, mask=(temp == undef))

    # Coefficients of exponential function
    if (isotope == 1):    # HDO
        a = +2.4844e+4
        b = -7.6248e+1
        c = +5.261e-2
    elif (isotope == 2):  # H218O
        a = +1.137e+3
        b = -4.156e-1
        c = -2.067e-3
    else:
        a = 0.
        b = 0.
        c = 0.

    # alpha+
    out = np.ma.exp( (a / mtemp + b) / mtemp + c)

    # alpha-
    if not greater1:
        out = 1./out

    # epsilon
    if eps:
        out -= 1.

    # return same as input type
    if islist:
        # fill undefined
        return list(out.filled(undef))
    elif isarr:
        return out.filled(undef)
    else:
        if temp == undef:
            return undef
        else:
            return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # T0 = 273.15
    # T  = np.array([0, 10., 15., 25.])
    # print(astr(alpha_equ_h2o(T+T0, isotope=0), 4))
    # # ['1.0000' '1.0000' '1.0000' '1.0000']

    # print(astr(alpha_equ_h2o(T+T0, isotope=1), 4))
    # # ['1.0117' '1.0107' '1.0102' '1.0094']

    # print(astr(alpha_equ_h2o(np.ma.array(T+T0, mask=(T==0.)), isotope=1, greater1=False), 4))
    # # ['-9999.0000' '    0.9894' '    0.9899' '    0.9907']

    # print(astr(alpha_equ_h2o(T+T0, isotope=2, eps=True)*1000., 4))
    # # ['112.3194' ' 97.6829' ' 91.1296' ' 79.3443']

    # print(astr(alpha_equ_h2o(0.+T0, isotope=2, eps=True)*1000., 4))
    # # 11.7187
