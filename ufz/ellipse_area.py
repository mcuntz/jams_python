#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def ellipse_area(a, b=None):
    """
        Area of ellipse with two axes a and b.


        Definition
        ----------
        def ellipse_area(a, b):


        Input
        -----
        a          main axis


        Optional Input
        --------------
        b          secondary axis


        Output
        ------
        area of ellipse or circle if secondary axis not given


        Restrictions
        ------------
        None


        Examples
        --------
        >>> from autostring import astr
        >>> area = ellipse_area(1, 1)
        >>> print(astr(area, 3))
        3.142

        >>> area = ellipse_area(2)
        >>> print(astr(area, 3))
        12.566
        

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

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Oct 2014
    """
    if b is None:
        return np.pi*a*a
    else:
        return np.pi*a*b


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # area = ellipse_area(1, 1)
    # print(astr(area, 3))

    # area = ellipse_area(2)
    # print(astr(area, 3))
