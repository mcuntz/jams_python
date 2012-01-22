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
