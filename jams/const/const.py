#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
"""
    Provides physical, mathematical, computational, and isotope constants.


    Definition
    ----------
    Pi = 3.141592653589793238462643383279502884197
    ...
    Define the following constants:
        Mathematical
            Pi, Pi2, Pi3, TwoPi, Sqrt2
        Physical
            Gravity, T0, P0, T25, sigma, R, R_air, R_H2O, Na, REarth
        Isotope
            R13VPDB, R18VSMOW, R2VSMOW
        Computational
            tiny, huge, eps
        Material
            mmol_co2, mmol_h2o, mmol_air
            density_quartz, cheat_quartz, cheat_water, cheat_air, latentheat_vaporization


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

    >>> print(astr(R13VPDB,3,pp=True))
    0.011

    >>> print(astr(tiny,3,pp=True))
    1.000e-06

    >>> print(astr(REarth,3,pp=True))
    6371000.000

    >>> print(astr(mmol_h2o,3,pp=True))
    18.015

    >>> print(astr(mmol_air,3,pp=True))
    28.964

    >>> print(astr(density_quartz,3,pp=True))
    2.650

    >>> print(astr(cheat_quartz,3,pp=True))
    800.000

    >>> print(astr(cheat_water,3,pp=True))
    4180.000

    >>> print(astr(cheat_air,3,pp=True))
    1010.000

    >>> print(astr(latentheat_vaporization,3,pp=True))
    2.450e+06


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2012-2017 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Jan 2012
    Modified, MC, Feb 2013 - ported to Python 3
              AP, Mar 2014 - dielectric constant H2O
              AP, Sep 2014 - heat capacity of quartz, air and water, density of quartz
              MC, Mar 2015 - Pi3=Pi/3
                           - rename heat capacities, molar masses, density of quartz
                           - moved dielH2O to own routine/file
                           - R13VPDB, R18VSMOW, R2VSMOW
              MC, Nov 2016 - tiny->np.finfo(np.float).tiny, huge
              RL, May 2017 - gas constants for dry air and water vapour
"""

__all__ = ['Pi', 'Pi2', 'Pi3', 'TwoPi', 'Sqrt2',
           'Gravity', 'T0', 'P0', 'T25', 'sigma', 'R','R_air', 'R_H2O', 'Na', 'REarth',
           'mmol_co2', 'mmol_h2o', 'mmol_air',
           'density_quartz', 'cheat_quartz', 'cheat_water', 'cheat_air', 'latentheat_vaporization',
           'R13VPDB', 'R18VSMOW', 'R2VSMOW',
           'tiny', 'huge', 'eps']

# Mathematical
Pi    = 3.141592653589793238462643383279502884197    # Pi
Pi2   = 1.57079632679489661923132169163975144209858  # Pi/2
Pi3   = 1.0471975511965977461542144610931676280656   # Pi/3
TwoPi = 6.283185307179586476925286766559005768394    # 2*Pi
Sqrt2 = 1.41421356237309504880168872420969807856967  # Sqrt(2)

# Physical
Gravity  = 9.81          # Standard average Earth's gravity [m^2 s^-1]
T0       = 273.15        # Celcius <-> Kelvin [K]
P0       = 101325.       # Standard pressure [Pa]
T25      = 298.15        # Standard ambient temperature [K]
sigma    = 5.67e-08      # Stefan-Boltzmann constant [W m^-2 K^-4]
R        = 8.3144621     # Ideal gas constant [J K^-1 mol^-1]
R_air    = 287.06        # gas constant of dry air [J K^-1 kg^-1]
R_H2O    = 461.4         # gas constant of water vapour [J K^-1 kg^-1]
Na       = 6.02214129e23 # Avogrado number [mol^-1]
REarth   = 6371009.      # Radius of Earth [m]

# Material
mmol_co2 = 44.01         # Molar mass CO2 [g mol^-1]
mmol_h2o = 18.01528      # Molar mass water [g mol^-1]
mmol_air = 28.9644       # Molar mass of dry air [g mol^-1]
# from Cambell G (1985) Soil Physics with BASIC, Elsevier Science
density_quartz = 2.65    # density of quartz [g cm^-3]
cheat_quartz   = 800.    # heat capacity of quartz [J kg^-1 K^-1]
cheat_water    = 4180.   # heat capacity of water [J kg^-1 K^-1]
cheat_air      = 1010.   # heat capacity of air [J kg^-1 K^-1]
latentheat_vaporization = 2.45e6 # latent heat of vaporization of water [J kg^-1]

# Isotope
R13VPDB  = 0.0112372     # 13C isotope ratio of VPDB
R18VSMOW = 2005.2e-6     # 18O isotope ratio of VSMOW
R2VSMOW  = 155.76e-6     # 2H  isotope ratio of VSMOW

# Computational
eps  = np.finfo(np.float).eps
huge = np.finfo(np.float).max
tiny = np.finfo(np.float).tiny


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
