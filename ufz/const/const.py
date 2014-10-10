#!/usr/bin/env python
from __future__ import print_function
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
    >>> from autostring import astr
    >>> print(astr(Pi,3,pp=True))
    3.142

    >>> print(astr(Sqrt2,3,pp=True))
    1.414

    >>> print(astr(Gravity,3,pp=True))
    9.810

    >>> print(astr(T0,3,pp=True))
    273.150

    >>> print(astr(sigma,3,pp=True))
    5.670e-08

    >>> print(astr(RPDB,3,pp=True))
    0.011

    >>> print(astr(tiny,3,pp=True))
    1.000e-06

    >>> print(astr(REarth,3,pp=True))
    6371000.000

    >>> print(astr(M_WV,3,pp=True))
    18.015

    >>> print(astr(M_DAIR,3,pp=True))
    28.964


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

    Copyright 2012-2014 Matthias Cuntz


    History
    -------
    Written,  MC, Jan 2012
    Modified, MC, Feb 2013 - ported to Python 3
    Modified, AP, Mar 2014 - add: dielectric constant H2O
    Modified, AP, Sep 2014 - add: heat capacity of quartz, air and water, density of quartz
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
M_CO2   = 44.01         # Molar mass CO2 [g*mol^-1]
M_WV    = 18.01528      # Molar mass water vapour [g*mol^-1]
M_DAIR  = 28.9644       # Molar mass of dry air [g*mol^-1]
# from Cambell, G., Soil Physics with BASIC, Elsevier Science, 1985:
rhoq    = 2.65          # density of quartz [g cm-3]
cqua    = 0.80e3        # heat capacity of quartz [J kg-1 K-1]
cwat    = 4.18e3        # heat capacity of water [J kg-1 K-1]
cair    = 1.01e3        # heat capacity of air [J kg-1 K-1]
lam     = 2.45e6        # specific heat of vaporization of water [J/kg]

def dielH2O(T):         # dielectric constant of water [F/m]
    '''
    in:     temperature T [K]
    out:    dielectric constant of H2O [F/m]
    source: Kaatze, U.: Reference liquids for the calibration of dielectric
            sensors and measurement instruments,
            Measurement Science and Technology, 2007, 18
            valid for 0 - 60 degC
    '''
    return 78.35*np.ma.exp(-4.55e-3*(T-298.15))

# Isotope
RPDB = 0.0112372 # Isotope ratio of VPDB-CO2

# Computational
tiny = 1e-6
eps  = np.finfo(np.float).eps


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
