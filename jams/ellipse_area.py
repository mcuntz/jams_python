#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

__all__ = ['ellipse_area']

def ellipse_area(a, b=None):
    """
        Area of ellipse with major and minor axes a and b.


        Definition
        ----------
        def ellipse_area(a, b=None):


        Input
        -----
        a          semi-major axis


        Optional Input
        --------------
        b          semi-minor axis


        Output
        ------
        area of ellipse or circle if semi-minor axis not given


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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Oct 2014 - in Python seminar
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
