#!/usr/bin/env python
import numpy as np
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
    >>> print Pi
    3.14159265359

    >>> print Sqrt2
    1.41421356237

    >>> print Gravity
    9.81

    >>> print T0
    273.15

    >>> print sigma
    5.67e-08

    >>> print RPDB
    0.0112372

    >>> print tiny
    1e-06

    >>> print REarth
    6371000.


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
    Written, MC, Jan 2012
"""
# Mathematical
Pi    = 3.141592653589793238462643383279502884197    # Pi
Pi2   = 1.57079632679489661923132169163975144209858  # Pi/2
TwoPi = 6.283185307179586476925286766559005768394    # 2*Pi
Sqrt2 = 1.41421356237309504880168872420969807856967  # Sqrt(2)

# Physical
Gravity = 9.81          # Standard average Earth's gravity [m^2/s]
T0      = 273.15        # Celcius <-> Kelvin [K]
P0      = 101325.       # Standard pressure [Pa]
T25     = 298.15        # Standard ambient temperature [K]
sigma   = 5.67e-08      # Stefan-Boltzmann constant [W/m^2/K^4]
R       = 8.3144621     # Ideal gas constant [J/K/mol]
Na      = 6.02214129e23 # Avogrado number [mol^-1]
REarth  = 6371009.      # Radius of Earth [m]

# Isotope
RPDB = 0.0112372 # Isotope ratio of VPDB-CO2

# Computational
tiny = 1e-6
eps  = np.finfo(np.float).eps


if __name__ == '__main__':
    import doctest
    doctest.testmod()
