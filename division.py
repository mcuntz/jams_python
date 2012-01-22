#!/usr/bin/env python
import numpy as np

def division(a, b, otherwise=np.nan, prec=0.):
    """
        Divide two arrays, return "otherwise" if division by 0.

        Definition
        ----------
        def division(a, b, otherwise=0., prec=0.):


        Input
        -----
        a            enumerator
	b            denominator


        Optional Input
        --------------
	otherwise    value to return if b=0 (default: np.nan)
	prec         if |b|<|prec| then otherwise


        Output
        ------
        a/b          if |b|>|prec|
        otherwise    if |b|<=|prec|


        Restrictions
        ------------
        None.


        Examples
        --------
        >>> print division([1., 2., 3.], 2.)
        [ 0.5  1.   1.5]

        >>> print division([1., 1., 1.], [2., 1., 0.])
        [ 0.5  1.   nan]

        >>> print division([1., 1., 1.], [2., 1., 0.], 0.)
        [ 0.5  1.   0. ]

        >>> print division([1., 1., 1.], [2., 1., 0.], otherwise=0.)
        [ 0.5  1.   0. ]

        >>> print division([1., 1., 1.], [2., 1., 0.], prec=1.)
        [ 0.5  nan  nan]


        History
        -------
        Written, MC, Jan 2012
    """
    return np.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
