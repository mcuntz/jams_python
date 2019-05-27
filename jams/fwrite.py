#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def fwrite(fname, arr, header=None, precision='10.0', delimiter=' '):
    """
        Write numbers of 2D-array to a file.

        A header can be given as optional as well as a precision.

        
        Definition
        ----------
        def fwrite(fname, data, header=None, precision='10.0'):

        
        Input
        -----
        fname        target file name
        arr          2d numpy array to write


        Optional Input Parameters
        -------------------------
        header       list of header elements: a header element is a list of two strings:
                     the first entry is the header argument, the second the value
        precision    floating point precision of array to write
        delimiter    delimiter to separate values, default ' '

        
        Examples
        --------
        >>> from jams import fread
        >>> # Clean up doctest
        >>> filename = 'fwrite.test'
        >>> header = [['Description', 'testing'], ['author', 'ST']]
        >>> data = np.arange(10).reshape(2, 5)
        >>> fwrite(filename, data, header=header)

        >>> fread(filename, nc=2, skip=2, header=True)
        [['Description', 'testing'], ['author', 'ST']]
        >>> fread(filename, skip=2)
        array([[0., 1., 2., 3., 4.],
               [5., 6., 7., 8., 9.]])

        >>> import os
        >>> os.remove(filename)


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2009-2015 Stephan Thober


        History
        -------
        Written,  ST, Feb 2016
        Modified, ST, Aug 2017 - delimiter
    """
    if not type(arr) is np.ndarray:
        raise ValueError('function fwrite: argument arr must be numpy.ndarray')
    fo = open(fname, 'w')
    if not header is None:
        # write header
        for ll in np.arange(len(header)):
            write_str = str(header[ll][0]) + ' ' + str(header[ll][1]) + '\n'
            fo.write(write_str)
    # write arr
    for ll in np.arange(arr.shape[0]):
        # format is mRM compatible
        write_str = delimiter.join(['{:' + precision + 'f}'] * arr.shape[1]).format(*arr[ll, :]) + '\n'
        fo.write(write_str)
    fo.close()    


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
