#!/usr/bin/env python
import numpy as np

def closest(vec, num, value=False):
    """
        Determine the array index in a vector of the element that is closest in value to a given number.

        Definition
        ----------
        def closest(vec, num, value=False):


        Input
        -----
        vec        vector
        num        number to compare


        Optional Input
        --------------
        value      If true, give closest number instead of index


        Output
        ------
        Index of element closest to given number num.


        Restrictions
        ------------
        None.


        Examples
        --------
        >>> vec = np.arange(100)/99.*5.
        >>> print closest(vec, 3.125)
        62
        
        >>> print closest(vec, 3.125, value=True)
        3.13131313131


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
    out = np.ma.argmin(np.ma.abs(np.ma.array(vec)-num))
    if value:
      return vec[out]
    else:
      return out

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
