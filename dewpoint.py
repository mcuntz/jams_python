#!/usr/bin/env python
import numpy as np

def dewpoint(pres, Celsius=False):
    """
        Calculates the dew point [K] from ambient humidity [Pa].

        Definition
        ----------
        def dewpoint(pres):


        Input
        -----
        pres       actual vapour pressure [Pa]


        Optional Input
        --------------
        Celsius    If True, return degree Celsius instead of Kelvin


        Output
        ------
        Dew point in Kelvin [K]


        Restrictions
        ------------
        None.


        References
        ---------
        Uses Bucks original vapour pressure formulation based on Tetens formula
        Buck, A. L., New equations for computing vapour pressure and enhancement factor,
                     J. Appl. Meteorol., 20, 1527-1532, 1981.

        Examples
        --------
        >>> Ta = 20. + 273.15
        >>> from esat import * # from ufz
        >>> es = esat(Ta, formula='Buck_original')
        >>> print dewpoint(es)
        293.15
        
        >>> print dewpoint(es, Celsius=True)
        20.0


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
    pw = 611.21
    c1 = 240.97
    c2 = 17.502
    w   = np.log(pres/pw)
    out = w*c1/(c2-w)
    if Celsius:
      return out
    else:
      import const
      return out+const.T0

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
