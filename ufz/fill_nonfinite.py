#!/usr/bin/env python
from __future__ import print_function
import numpy as np

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
        find     function, with signature indices= find(ind),
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


def fill_nonfinite(xin, nan=None, inf=None, undef=None):
    """
        Fill missing values by linear interpolation.

        Checks for NaN and Inf in xin and interpolates along column (1st dim).
        If undef=number is given, it will be interpolated as well.
        If one of nan of inf is given, then only that will be checked,
        i.e. no keyword or only undef given: nan=True, inf=True.
             only nan and undef given: nan=True, inf=False.


        Definition
        ----------
        def fill_nonfinite(xin, nan=None, inf=None, undef=None):


        Input
        -----
        xin     ND array


        Optional Input
        --------------
        nan        if True, interpolate NaN.
        inf        if True, interpolate Inf.
        undef      if given, interpolate xin==undef.


        Output
        ------
        xin filled at NaN and Inf values.


        Examples
        --------
        >>> from autostring import astr
        >>> a = np.array([1, 2, 3])
        >>> print(astr(fill_nonfinite(a),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a = np.array([1, np.nan, 3])
        >>> print(astr(fill_nonfinite(a),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a = np.array([1, np.inf, 3])
        >>> print(astr(fill_nonfinite(a,inf=True),1,pp=True))
        ['1.0' '2.0' '3.0']

        >>> a = np.array([1, 2, 3])
        >>> print(astr(fill_nonfinite(a,undef=2),1,pp=True))
        ['1.0' '2.0' '3.0']


        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2013
    """
    #
    # Assure ND-array
    isone = False
    if np.ndim(xin) == 1:
        isone = True
        xin = xin[:,np.newaxis]
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
    xin = xin.astype(np.float) # np.interp always returns float
    for i in range(xin.shape[1]):
        yy = xin[:,i]
        nans, ind = nan_helper(yy, nan=inan, inf=iinf, undef=undef)
        if ind(nans).size > 0:
            yy[nans] = np.interp(ind(nans), ind(~nans), yy[~nans])
            xin[:,i] = yy

    if isone:
        return np.squeeze(xin)
    else:
        return xin


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # a = np.array([1, 2, 3])
    # print(astr(fill_nonfinite(a),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a = np.array([1, np.nan, 3])
    # print(astr(fill_nonfinite(a),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a = np.array([1, np.inf, 3])
    # print(astr(fill_nonfinite(a,inf=True),1,pp=True))
    # #['1.0' '2.0' '3.0']
    # a = np.array([1, 2, 3])
    # print(astr(fill_nonfinite(a,undef=2),1,pp=True))
    # #['1.0' '2.0' '3.0']
