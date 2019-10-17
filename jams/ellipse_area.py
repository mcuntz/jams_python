#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
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
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
