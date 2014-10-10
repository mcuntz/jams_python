#!/usr/bin/env python
"""
    Provides physical, mathematical, computational, and isotope constants.


    Definition
    ----------
    Pi = 3.141592653589793238462643383279502884197
    ...
    Define the following constants:
        Mathematical
            Pi, Pi2, TwoPi, Sqrt2
        Physical
            Gravity, T0, P0, T25, sigma, R, Na, REarth
        Isotope
            RPDB
        Computational
            tiny


    Examples
    --------
    area = ufz.const.Pi * radius**2


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

    Copyright 2014 Matthias Cuntz


    History
    -------
    Written,  MC, Oct 2014
"""
from .const import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1796"
__date__     = 'Date: 05.10.2014'
