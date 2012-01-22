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
        >>> from ufz import esat
        >>> es = esat(Ta, formula='Buck_original')
        >>> print dewpoint(es)
        293.15
        
        >>> print dewpoint(es, Celsius=True)
        20.0
        

        History
        -------
        Written, MC, Jan 2012
    """
    T0 = 273.15
    pw = 611.21
    c1 = 240.97
    c2 = 17.502
    w   = np.log(pres/pw)
    out = w*c1/(c2-w)
    if Celsius:
      return out
    else:
      return out+T0

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
