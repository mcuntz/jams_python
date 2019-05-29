#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.const import T25

__all__ = ['dielectric_water']

def dielectric_water(T):
    """
        Dielectric constant of liquid water in F m^-1.


        Definition
        ----------
        def dielectric_water(T):


        Input
        -----
        list/ND-array of temperatures [K]

        
        Output
        ------
        list/ND-array of dielectric constants of liquid water [F m^-1]


        Restrictions
        ------------
        Valid between 0 degC and 60 degC


        References
        ----------
        Kaatze, U (2007) Reference liquids for the calibration of dielectric
                         sensors and measurement instruments,
                         Measurement Science and Technology, 18


        Examples
        --------
        >>> from autostring import astr
        >>> print(astr(dielectric_water(300.), 3))
        77.693


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  AP, Mar 2014 - in const
                  MC, Mar 2015 - own subroutine
    """
    if isinstance(T, np.ma.masked_array):
        return 78.35*np.ma.exp(-4.55e-3*(T-T25))
    else:
        return 78.35*np.exp(-4.55e-3*(T-T25))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
