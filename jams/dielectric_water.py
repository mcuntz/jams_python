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
