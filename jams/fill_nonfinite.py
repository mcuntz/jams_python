#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from jams.kernel_regression import kernel_regression

def nan_helper(y, nan=False, inf=False, undef=None):
    """
        Helper to handle indices and logical indices of NaNs, Infs or undefs.


        Definition
        ----------
        def nan_helper(y, nan=False, inf=False, undef=None):


        Input
        -----
        y        1d numpy array with possible missing values


        Optional Input
        --------------
        At least one of the following has to be given
        nan        if True, check only for NaN and not Inf.
        inf        if True, check only for Inf and not NaN.
        undef      if given then check for undef value rather than NaN and Inf.


        Output
        ------
        ind      logical indices of missing values
        find     function, with signature indices = find(ind),
                 to convert logical indices of NaNs to 'equivalent' indices


        Examples
        --------
        >>> # linear interpolation of NaNs
        >>> y       = np.array([1, np.nan, 3])
        >>> nans, z = nan_helper(y, nan=True)
        >>> y[nans] = np.interp(z(nans), z(~nans), y[~nans])


        History
        -------
        Written,  MC, Jul 2013 - modified from
                                 http://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
                  MC, Apr 2014 - assert
    """
    assert not ((not nan) & (not inf) & (undef is None)), 'at least one of nan, inf or undef has to be given.'

    out = np.zeros(y.shape, dtype=np.bool)
    if nan:
        out = out | np.isnan(y)
    if inf:
        out = out | np.isinf(y)
    if undef!=None:
        out = out | (y==undef)

    return out, lambda ind: ind.nonzero()[0]


def fill_nonfinite(xin, yin=None, nan=None, inf=None, undef=None, method='interp', **kwargs):
    """
        Fill missing values by interpolation.

        Checks for NaN and Inf in xin and interpolates along column (1st dim).
        If yin is None, then interpolation distance is determined by index number.
        If undef=number is given, it will be interpolated as well.
        If one of nan of inf is given, then only that will be checked,
        i.e. no keyword or only undef given: nan=True, inf=True.
             only nan and undef given: nan=True, inf=False.


        Definition
        ----------
        def fill_nonfinite(xin, yin=None, nan=None, inf=None, undef=None, method='interp', **kwargs):


        Input
        -----
        xin     ND array


        Optional Input
        --------------
        yin        ND array (default: None)
                   If not given, then y=xin and x=np.arange(len(xin[:,0]))
        nan        if True, interpolate NaN (default: None, see above).
        inf        if True, interpolate Inf (default: None, see above).
        undef      if given, interpolate xin==undef (default: None).
        method     Interpolation method (default: 'interp'):
                   'interp': linear interpolation with np.interp.
                   'linear': same as 'interp'
                   'kernel_regression': non-linear interpolation with jams.kernel_regression.
        **kwargs   All other keyword arguments will be passed to interpolation routine.


        Output
        ------
        xin filled at NaN and Inf values.


        Examples
        --------
        >>> from autostring import astr

        >>> a = np.arange(3)+1.
        >>> print(astr(fill_nonfinite(a),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a[1] = np.nan
        >>> print(astr(fill_nonfinite(a),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a[1] = np.inf
        >>> print(astr(fill_nonfinite(a,inf=True),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a = np.arange(3)+1.
        >>> print(astr(fill_nonfinite(a,undef=2),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> print(astr(fill_nonfinite(a,undef=2,method='Kernel_regression'),1,pp=True))
        ['1.0' '2.0' '3.0']
        >>> print(astr(fill_nonfinite(a,undef=2,method='Kernel_regression',silverman=True),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a = np.arange(10)+1.
        >>> x = a[:]
        >>> x[1::3] -= 0.5
        >>> print(astr(fill_nonfinite(x,a,undef=2),1,pp=True))
        [' 1.0' ' 1.5' ' 3.0' ' 4.0' ' 4.5' ' 6.0' ' 7.0' ' 7.5' ' 9.0' '10.0']
        >>> print(astr(fill_nonfinite(x,a,undef=2,method='Kernel_regression',silverman=True),1,pp=True))
        [' 1.0' ' 1.5' ' 3.0' ' 4.0' ' 4.5' ' 6.0' ' 7.0' ' 7.5' ' 9.0' '10.0']


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2013-2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jul 2013
        Modified, MC, Jan 2018 - yin
                  MC, Jan 2018 - added kernel_regression and **kwargs
    """
    #
    # yin
    if yin is not None:
        xx = xin
        yy = yin
    else:
        xx = None
        yy = xin
    #
    # Assure ND-array
    isone = False
    if np.ndim(yy) == 1:
        isone = True
        yy = yy[:,np.newaxis]
        ny = 1
    else:
        ny = yy.shape[1]
    if xx is None:
        nx = 0
    else:
        if np.ndim(xx) == 1:
            xx = xx[:,np.newaxis]
            nx = 1
        else:
            nx = xx.shape[1]
            assert nx == ny, 'x and y ND-arrays must have same dimensions if x is not 1D-array.'
    #
    # Check options - if none given then nan and inf are True
    inan = False
    iinf = False
    if (nan==None) & (inf==None):
        inan = True
        iinf = True
    else:
        if (nan!=None): inan=nan
        if (inf!=None): iinf=inf

    # interpolate along each column
    yout = yy.astype(np.float) # np.interp always returns float
    for i in range(ny):
        yy1 = yout[:,i]
        if xx is not None:
            xx1 = xx[:,min(i,nx-1)].astype(np.float)
        nans, ind = nan_helper(yy1, nan=inan, inf=iinf, undef=undef)            
        if np.sum(nans) > 0:
            if (method.lower() == 'interp') or (method.lower() == 'linear'):
                if xx is None:
                    yy1[nans] = np.interp(ind(nans), ind(~nans), yy1[~nans], **kwargs)
                else:
                    yy1[nans] = np.interp(xx1[nans], xx1[~nans], yy1[~nans], **kwargs)
            elif method.lower() == 'kernel_regression':
                if xx is None:
                    yy1[nans] = kernel_regression(ind(~nans), yy1[~nans], xout=ind(nans), **kwargs)
                else:
                    yy1[nans] = kernel_regression(xx1[~nans], yy1[~nans], xout=xx1[nans], **kwargs)
            else:
                raise ValueError('Interpolation method unknown:', method)
            yout[:,i] = yy1

    if isone:
        return np.squeeze(yout)
    else:
        return yout


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # a = np.arange(3)+1.
    # print(astr(fill_nonfinite(a),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a[1] = np.nan
    # print(astr(fill_nonfinite(a),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a[1] = np.inf
    # print(astr(fill_nonfinite(a,inf=True),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a = np.arange(3)+1.
    # print(astr(fill_nonfinite(a,undef=2),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # print(astr(fill_nonfinite(a,undef=2,method='Kernel_regression'),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # print(astr(fill_nonfinite(a,undef=2,method='Kernel_regression',silverman=True),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a = np.arange(10)+1.
    # x = a[:]
    # x[1::3] -= 0.5
    # print(astr(fill_nonfinite(x,a,undef=2),1,pp=True))
    # #[' 1.0' ' 1.5' ' 3.0' ' 4.0' ' 4.5' ' 6.0' ' 7.0' ' 7.5' ' 9.0' '10.0']
    # print(astr(fill_nonfinite(x,a,undef=2,method='Kernel_regression',silverman=True),1,pp=True))
    # #[' 1.0' ' 1.5' ' 3.0' ' 4.0' ' 4.5' ' 6.0' ' 7.0' ' 7.5' ' 9.0' '10.0']
