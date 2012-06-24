#!/usr/bin/env python
import numpy as np

def heaviside(x, value=1., unitstep=False, zero=False):
    """
        Heaviside (or unit step) operator
            H =  0  for x <  0
            H = 1/2 for x =  0
            H =  1  for x >  0
          if unitstep
            H =  0  for x <  0
            H =  1  for x >= 0
          if zero
            H =  0  for x =< 0
            H =  1  for x >  0

        Definition
        ----------
        def heaviside(x, value=1., unitstep=False, zero=False):


        Input
        -----
        x            value or array of values


        Optional Input
        --------------
        value        output is heaviside *= value
	unitstep     If True, H(0)=1 instead of 1/2
	zero         If True, H(0)=0 instead of 1/2


        Output
        ------
        Heaviside function 0, 1/2, and 1


        Restrictions
        ------------
        Returns False if error.


        Examples
        --------
        >>> print heaviside([-1,0.,1.])
        [0.0 0.5 1.0]

        >>> print heaviside([-1,0.,1.], zero=True)
        [ 0.  0.  1.]

        >>> print heaviside([-1,0.,1.], unitstep=True)
        [ 0.  1.  1.]

        >>> print heaviside([-1,0.,1.], value=2)
        [0.0 1.0 2.0]

        >>> print heaviside([-1,0.,1.], zero=True, unitstep=True)
        HEAVISIDE: unitstep and zero mutually exclusive.
        False


        History
        -------
        Written, MC, Jan 2012
    """
    if zero and unitstep:
        raise ValueError('unitstep and zero mutually exclusive.')

    if zero:
        out = np.where(np.ma.array(x) > 0., 1., 0.)
    elif unitstep:
        out = np.where(np.ma.array(x) >= 0., 1., 0.)
    else:
        out = (np.where(np.ma.array(x) > 0., 1., 0.) - np.where(np.ma.array(x) < 0., 1., 0.) + 1.) * 0.5
    out *= value
    return out

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
