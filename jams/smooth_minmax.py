#!/usr/bin/env python
from __future__ import print_function


def smax(x, y, eta=1. - 1.e-5):
    """
        Calculate smooth maximum of two numbers.


        Definition
        ----------
        def around(x, y, eta=1.):


        Input
        -----
        x        ND array of floats
        y        ND array of floats


        Optional Input
        --------------
        eta        Shape of smoothing function (default: 1._dp)


        Output
        ------
        Smoothed maximum


        Restrictions
        ------------
        If eta is set to 1, then this function becomes the standard maximum function.


        Examples
        --------
        >>> from jams.autostring import astr
        >>> print(astr(smax(3.5967, 3.6),3,pp=True))
            3.610


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

        Copyright 2016-2016 Stephan Thober


        History
        -------
        Written,  ST, Mar 2016
        Modified, 
    """
    #
    if (eta > 1.):
        raise ValueError('eta is {:3.4f}, but cannot be larger than 1.'.format(eta))
    #
    out = (x + y + ((x + y)**2. - 4 * eta * x * y)**0.5) / (2. * eta)
    return out


def smin(x, y, eta=1. - 1.e-5):
    """
        Calculate smooth minimum of two numbers.


        Definition
        ----------
        def around(x, y, eta=1. - 1.e-5):


        Input
        -----
        x        ND array of floats
        y        ND array of floats


        Optional Input
        --------------
        eta        Shape of smoothing function (default: 1._dp)


        Output
        ------
        Smoothed minimum


        Restrictions
        ------------
        If eta is set to 1, then this function becomes the standard minimum function.


        Examples
        --------
        >>> from jams.autostring import astr
        >>> print(astr(smin(3.5967, 3.6),3,pp=True))
            3.587


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

        Copyright 2016-2016 Stephan Thober


        History
        -------
        Written,  ST, Mar 2016
        Modified, 
    """
    #
    if (eta > 1.):
        raise ValueError('eta is {:3.4f}, but cannot be larger than 1.'.format(eta))
    #
    out = (x + y - ((x + y)**2. - 4 * eta * x * y)**0.5) / (2. * eta)
    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
