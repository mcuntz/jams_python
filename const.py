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
            Gravity, T0, P0, T25, sigma, R, Na
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

# Isotope
RPDB = 0.0112372 # Isotope ratio of VPDB-CO2

# Computational
tiny = 1e-6

if __name__ == '__main__':
    import doctest
    doctest.testmod()
