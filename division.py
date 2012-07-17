#!/usr/bin/env python
import numpy as np

def division(a, b, otherwise=np.nan, prec=0.):
    """
        Divide two arrays, return "otherwise" if division by 0.

        Definition
        ----------
        def division(a, b, otherwise=np.nan, prec=0.):
          There is a wrapper function for convenience with the short name 'div'
        def div(a, b, otherwise=np.nan, prec=0.):


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
        >>> print division([1., 2., 3.], 2.)
        [ 0.5  1.   1.5]

        >>> print division([1., 1., 1.], [2., 1., 0.])
        [ 0.5  1.   nan]

        >>> print division([1., 1., 1.], [2., 1., 0.], 0.)
        [ 0.5  1.   0. ]

        >>> print division([1., 1., 1.], [2., 1., 0.], otherwise=0.)
        [ 0.5  1.   0. ]

        >>> print division([1., 1., 1.], [2., 1., 0.], prec=1.)
        [ 0.5  nan  nan]

        >>> print div([1., 1., 1.], [2., 1., 0.], prec=1.)
        [ 0.5  nan  nan]


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

        Copyright 2012 Matthias Cuntz

        
        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, May 2012 - div
    """
    return np.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)

def div(a, b, otherwise=np.nan, prec=0.):
    return np.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)
 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
