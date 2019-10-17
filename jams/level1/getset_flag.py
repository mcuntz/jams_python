#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['get_flag', 'get_maxflag', 'set_flag']

# --------------------------------------------------------------------

def get_flag(flags, n):
    """
        Get the flags at position n from CHS data flag vector.


        Definition
        ----------
        def get_flag(flags, n):


        Input
        -----
        flags   ND-array of integers
        n       position to extract


        Output
        ------
        ND-array of integers with flag at position n of flags
        returns -1 if input integer is less than 10**(n+1)


        Examples
        --------
        --> see __init__.py for full example of workflow

        >>> flags = np.array([9, 90, 91, 900, 901, 9001, 9201, 912121212])
        >>> print(get_flag(flags, 0))
        [9 9 9 9 9 9 9 9]

        >>> print(get_flag(flags, 1))
        [-1  0  1  0  0  0  2  1]

        >>> print(get_flag(flags, 2))
        [-1 -1 -1  0  1  0  0  2]

        >>> print(get_flag(flags, 3))
        [-1 -1 -1 -1 -1  1  1  1]

        >>> print(get_flag(flags, 4))
        [-1 -1 -1 -1 -1 -1 -1  2]

        >>> flags = np.array([[9, 90, 91, 900, 901, 9001,  9201, 912121212],
        ...                   [9, 91, 92, 901, 912, 9101, -9999, 912121212]])
        >>> print(get_flag(flags, 1))
        [[-1  0  1  0  0  0  2  1]
         [-1  1  2  0  1  1 -2  1]]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Mar 2015
    """

    out = np.ones(flags.shape, dtype=np.int)*(-1)   # if not enough digits -> -1

    # deal with -9999
    jj = np.where(flags <= 0)
    if jj[0].size>0: out[jj] = -2

    # normal flags starting with 9
    jj = np.where(flags > 0)
    if jj[0].size>0:
        oo = out[jj]
        ilog10 = np.log10(flags[jj]).astype(np.int) # how many digits
        ii = np.where((ilog10-n) >= 0)              # otherwise extract n-th digit
        if ii[0].size > 0:
            oo[ii] = (flags[jj][ii] // 10**(ilog10[ii]-n)) % 10
        out[jj] = oo

    return out


# --------------------------------------------------------------------

def set_flag(flags, n, iflag, ii=None):
    """
        Set the flags at position n to iflag at indeces ii of CHS data flag vector.


        Definition
        ----------
AL        def set_flag(flags, n, iflag, ii=None):


        Input
        -----
        flags   CHS data flag vector, ND-array of integers
        n       position in flag to set, missing positions will be cretaed and filled with 0
        iflag   set flag to this value


        Optional Input
        --------------
        ii      indices at which to set flag (default: None, i.e. each entry)


        Output
        ------
        ND-array of integers with flag at position n of flags set to iflag


        Examples
        --------
        --> see __init__.py for full example of workflow

        >>> flags = np.array([9, 90, 901, 9101, 912121212, -9999])
        >>> print(set_flag(flags, 1, 2, [0,1,2]))
        [       92        92       921      9101 912121212        90]

        >>> print(set_flag(flags, 3, 2, [0,1,2,3]))
        [     9002      9002      9012      9102 912121212      9000]

        >>> print(set_flag(flags, 1, 2))
        [       92        92       921      9201 922121212        92]

        >>> flags = np.array([[9, 90, 901, 9101, 912121212, -9999],
        ...                   [9, 90, 901, 9101, 912121212, -9999]])
        >>> print(set_flag(flags, 1, 1, [[0,0,1],[0,1,2]]))
        [[       91        91       901      9101 912121212        90]
         [       90        90       911      9101 912121212        90]]

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Mar 2015
        Modified, MC, Jan 2016 - ND arrays possible
    """
    fflags = flags.copy()
    uu = np.where(fflags < 0)                        # if undef, set to treated that mean 9
    if uu[0].size > 0: fflags[uu] = 9
    # extend flag vector if needed
    ilog10 = np.log10(fflags).astype(np.int)         # are there enough flag positions (i)
    jj = np.where(ilog10 < n)
    if jj[0].size > 0:                               # increase number of entries (filled with 0)
        fflags[jj] *= 10**(n-ilog10[jj])
        ilog10      = np.log10(fflags).astype(np.int)

    # get entries in flag vector
    isflags = get_flag(fflags, n)

    # set entries in flag vector
    if ii is None:                                  # set all to iflag
        fflags += 10**(ilog10-n) * (iflag-isflags)
    else:
        if np.size(ii) > 0:                         # set ii to iflag
            fflags[ii] += 10**(ilog10[ii]-n) * (iflag-isflags[ii])

    return fflags


# --------------------------------------------------------------------

def get_maxflag(flags, no=0):
    """
        Returns the overall flag for each flag from CHS data flag vector.


        Definition
        ----------
        def get_maxflag(flags, no):


        Input
        -----
        flags   ND-array of integers
        no      column or list of columns that will be excluded from the determination of the maximum flag


        Output
        ------
        ND-array of integers with overall flag
        returns:
            2, if at least one flag is 2
            1, if at least one flag is 1
            0, if all flags are 0
           -1, if flag is '9'
           -2, if flag is still -9999


        Examples
        --------
        --> see __init__.py for full example of workflow

        >>> flags = np.array([9, 90, 91, 900, 901, 9001, 9201, 912121212, -9999])
        >>> print(get_maxflag(flags))
        [-1  0  1  0  1  1  2  2 -2]

        >>> print(get_maxflag(flags, no=1))
        [-1 -1 -1  0  1  1  1  2 -2]

        >>> print(get_maxflag(flags, no=[1,3]))
        [-1 -1 -1  0  1  0  0  2 -2]

        >>> flags = np.array([[9, 90, 91, 900, 901, 9001,  9201, 912121212],
        ...                   [9, 91, 92, 901, 912, 9101, -9999, 912121212]])
        >>> print(get_maxflag(flags))
        [[-1  0  1  0  1  1  2  2]
         [-1  1  2  1  2  1 -2  2]]

        >>> flags = np.array([[9, 90, 91, 900, 901, 9001,  9201, 912121212],
        ...                   [9, 91, 92, 901, 912, 9101, -9999, 912121212]])
        >>> print(get_maxflag(flags, no=[1,3]))
        [[-1 -1 -1  0  1  0  0  2]
         [-1 -1 -1  1  2  0 -2  2]]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  JM, Mar 2015
        Modified, JV, Oct 2015 - no
                  MC, Jan 2016 - no iterable
    """
    # assue no as list
    try:
        nno = len(no)
    except:
        no = [no]

    # Determine overall flag:
    #   maxflag=-2, if flag is '-9999'
    #   maxflag=-1, if flag is '9'
    #   maxflag=0,  if all individual flags are 0, e.g. '9000000'
    #   maxflag=1,  if at least one flag is 1,     e.g. '9000100'
    #   maxflag=2,  if at least one flag is 2,     e.g. '9002010'
    fmax = np.amax(flags)
    if fmax == -9999:
        # return -2 if no flags set yet
        return np.ones(np.shape(flags), dtype=np.int)*(-2)
    elif fmax == 9:
        # return -1 if no flags set yet but flags were already initialised
        return np.ones(np.shape(flags), dtype=np.int)*(-1)
    else:
        maxflag = np.ones(np.shape(flags), dtype=np.int)*(-2)
        n_flags = np.log10(fmax).astype(np.int)
        for ii in range(n_flags):
            iii = ii+1
            if iii not in no:
                maxflag = np.maximum(maxflag, get_flag(flags,iii))
        return maxflag


# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
