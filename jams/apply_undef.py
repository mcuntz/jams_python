#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['apply_undef']

def apply_undef(func, undef, *args, **kwargs):
    """
        Use a function on arguments masked with undef.

        Arguments will be masked with undef (masked arrays) and then
        passed to the function. The return value will be filled with
        undef at the masked entries.

        For example, using np.ma.add two arrays a and b with undef
            c = apply_undef(np.ma.add, undef, a, b)
        is equivalent to
            c = (np.ma.array(a, mask=(a==undef)) + np.ma.array(b, mask=(a==undef))).filled(undef)


        Definition
        ----------
        def apply_undef(func, undef, *args, **kwargs)


        Input
        -----
        func        numpy.ma.MaskedArray aware function
        undef       float
                    value to be masked
        args        list of arguments that will be masked before passing them to func


        Optional Input
        --------------
        kwargs     dictionary with keyword parameters passed to func


        Output
        ------
        func(*args) where args were masked and func(*args) filled with undef


        Restrictions
        ------------
        func must be aware of masked arrays.


        Examples
        --------
        >>> undef = -999

        # test scalars
        >>> print(apply_undef(np.ma.add, undef, 1, 2))
        3
        >>> print(apply_undef(np.ma.add, undef, undef, 2))
        -999.0
        
        # test iterables: can be list or tuple as well
        >>> a = np.array([1, 2, 3])
        >>> b = [2, undef, 4]
        >>> print(apply_undef(np.ma.add, undef, a, b))
        [   3 -999    7]
        >>> print(apply_undef(np.add, undef, a, b)) # numpy functions are mostly masked-aware now
        [   3 -999    7]

        # test function keywords
        >>> print(apply_undef(np.mean, undef, np.vstack([a, b])))
        2.4
        >>> print(apply_undef(np.mean, undef, np.vstack([a, b]), axis=0))
        [1.5 2.  3.5]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jan 2018
    """
    #
    # check input
    assert len(args) > 0, 'function arguments must be given.'
    #
    # prep arguments
    arg = []
    for a in args:
        if np.iterable(a):
            if isinstance(a, (list,tuple)):
                arg.append(np.ma.array(a, mask=(np.array(a)==undef)))
            elif isinstance(a, (np.ndarray,np.ma.MaskedArray)):
                arg.append(np.ma.array(a, mask=(a==undef)))
            else:
                raise IOError('Iterable not known.')
        else:
            arg.append(np.ma.array(a, mask=(a==undef)))
    # apply function
    res = func(*arg, **kwargs)

    # fill result
    if isinstance(res, (np.ma.MaskedArray,np.ma.core.MaskedIterator)):
        res = res.filled(undef)
    elif isinstance(res, np.ma.core.MaskedConstant):
        if res is np.ma.masked:
            res = undef
        else:
            res = res.data.reshape((1,))[0]

    return res


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # undef = -999
    # a = np.array([1, 2, 3])
    # b = [2, undef, 4]
    # print(apply_undef(np.ma.add, undef, a, b))
    # # [   3 -999    7]
    # print(apply_undef(np.add, undef, a, b))
    # # [   3 -999    7]

    # print(apply_undef(np.mean, undef, b))
    # # 3.0

    # print(apply_undef(np.mean, undef, np.vstack([a, b]), axis=0))
    # # [ 1.5  2.   3.5]

    # print(apply_undef(np.ma.add, undef, 1, 2))
    # # 3

    # print(apply_undef(np.ma.add, undef, undef, 2))
    # # -999.0
