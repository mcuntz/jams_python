#!/usr/bin/env python
import numpy as np

def around(num, powten, ceil=False, floor=False):
    """
        Round to the passed power of ten.

        Definition
        ----------
        def around(num, powten=None, ceil=False, floor=False):


        Input
        -----
        num        number array


        Optional Input
        --------------
        powten     Power of ten array
                   If missing, simple round (ceil, floor) is taken.
        ceil       ceil instead of round to the nearest power of ten
        floor      floor instead of round to the nearest power of ten


        Output
        ------
        Rounded values.


        Restrictions
        ------------
        Powten is exactly opposite of decimal keyword of numpy.around.

        From numpy.around documentation:
        'For values exactly halfway between rounded decimal values,
        Numpy rounds to the nearest even value. Thus 1.5 and 2.5 round to 2.0,
        -0.5 and 0.5 round to 0.0, etc. Results may also be surprising due to the
        inexact representation of decimal fractions in the IEEE floating point
        standard and errors introduced when scaling by powers of ten.'


        Examples
        --------
        >>> around(np.array([3.5967,345.5967]), -3)
        array([   3.597,  345.597])

        >>> around(np.array([1994344,345.5967]), [3,-3])
        array([  1.99400000e+06,   3.45597000e+02])

        >>> around(np.array([1994344,345.5967]), [3,-3], ceil=True)
        array([  1.99500000e+06,   3.45597000e+02])

        >>> around(np.array([1994344,345.5967]), [3,-3], floor=True)
        array([  1.99400000e+06,   3.45596000e+02])

        >>> around(np.array([3.5967,345.5967]), 3)
        array([ 0.,  0.])

        >>> around(np.array([3.5967,345.5967]), 3, ceil=True)
        array([ 1000.,  1000.])


        History
        -------
        Written, MC, Jun 2011
    """
    #
    # Check input
    if (ceil and floor):
        raise ValueError('ceil and floor keywords given.')
    if (powten == None):
        ipowten = 0
    else:
        ipowten = np.array(powten)
    nnum = np.size(num)
    npowten = np.size(ipowten)
    if ((npowten != nnum) and (npowten != 1)):
        raise ValueError('powten must be scalar or have array size of input numbers.')
    #
    # Shift decimal point
    #out = num * 10.**(-ipowten)
    out = num * np.exp(-ipowten*np.log(10.))
    # Round/ceil/floor
    #eps = np.MachAr().eps
    eps = np.finfo(np.float).eps
    if (ceil):
        out = np.ceil(out-10.*eps)
    elif (floor):
        out = np.floor(out+10.*eps)
    else:
       out = np.around(out)
    # Shift back decimal point
    #out *= 10.**ipowten
    out *= np.exp(ipowten*np.log(10.))

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod()
