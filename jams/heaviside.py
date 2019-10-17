#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
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
            H =  0  for x <= 0
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
        >>> from autostring import astr
        >>> print(astr(heaviside([-1,0.,1.]),1,pp=True))
        ['0.0' '0.5' '1.0']

        >>> print(astr(heaviside([-1,0.,1.], zero=True),1,pp=True))
        ['0.0' '0.0' '1.0']

        >>> print(astr(heaviside([-1,0.,1.], unitstep=True),1,pp=True))
        ['0.0' '1.0' '1.0']

        >>> print(astr(heaviside([-1,0.,1.], value=2),1,pp=True))
        ['0.0' '1.0' '2.0']


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Matthias Cuntz - mc (at) macu (dot) de

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
                  MC, Apr 2014 - assert
    """
    assert (zero+unitstep) < 2, 'unitstep and zero mutually exclusive.'

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
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

