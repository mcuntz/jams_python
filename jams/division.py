#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

def division(a, b, otherwise=np.nan, prec=0.):
    """
        Divide two arrays, return "otherwise" if division by 0.


        Definition
        ----------
        def division(a, b, otherwise=np.nan, prec=0.):
          There is a wrapper function for convenience with the short name 'div'.


        Input
        -----
        a            enumerator
	b            denominator


        Optional Input
        --------------
	otherwise    value to return if b=0 (default: np.nan)
	prec         if |b|<|prec| then otherwise


        Output
        ------
        a/b          if |b|>|prec|
        otherwise    if |b|<=|prec|


        Restrictions
        ------------
        None.


        Examples
        --------
        >>> from autostring import astr
        >>> print(astr(division([1., 2., 3.], 2.),1,pp=True))
        ['0.5' '1.0' '1.5']

        >>> print(astr(division([1., 1., 1.], [2., 1., 0.]),1,pp=True))
        ['0.5' '1.0' 'nan']

        >>> print(astr(division([1., 1., 1.], [2., 1., 0.], 0.),1,pp=True))
        ['0.5' '1.0' '0.0']

        >>> print(astr(division([1., 1., 1.], [2., 1., 0.], otherwise=0.),1,pp=True))
        ['0.5' '1.0' '0.0']

        >>> print(astr(division([1., 1., 1.], [2., 1., 0.], prec=1.),1,pp=True))
        ['0.5' 'nan' 'nan']


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

        Copyright 2012-2014 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, May 2012 - div
                  MC, Feb 2013 - ported to Python 3
                  MC, Oct 2014 - do not return masked array if no masked array given
                  MC, Sep 2015 - bug: returned non-masked array in case of masked array input
    """
    oldsettings = np.geterr()
    np.seterr(divide='ignore')

    if isinstance(a, np.ma.masked_array) or isinstance(b, np.ma.masked_array):
        out = np.ma.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)
    else:
        out = np.where(np.abs(np.array(b)) > np.abs(prec), np.array(a)/np.array(b), otherwise)

    np.seterr(**oldsettings)

    return out



def div(*args, **kwargs):
    """
        Wrapper function for division
        def division(a, b, otherwise=np.nan, prec=0.):


        Examples
        --------
        >>> from autostring import astr
        >>> print(astr(div([1., 2., 3.], 2.),1,pp=True))
        ['0.5' '1.0' '1.5']

        >>> print(astr(div([1., 1., 1.], [2., 1., 0.]),1,pp=True))
        ['0.5' '1.0' 'nan']

        >>> print(astr(div([1., 1., 1.], [2., 1., 0.], 0.),1,pp=True))
        ['0.5' '1.0' '0.0']

        >>> print(astr(div([1., 1., 1.], [2., 1., 0.], otherwise=0.),1,pp=True))
        ['0.5' '1.0' '0.0']

        >>> print(astr(div([1., 1., 1.], [2., 1., 0.], prec=1.),1,pp=True))
        ['0.5' 'nan' 'nan']
    """
    return division(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

