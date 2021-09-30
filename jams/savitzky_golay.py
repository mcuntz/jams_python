#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from math import factorial
import scipy.signal as sps


__all__ = ['savitzky_golay', 'sg', 'savitzky_golay2d', 'sg2d']


def savitzky_golay(y, window, order, deriv=0, rate=1):
    """
    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.


    Definition
    ----------
    def savitzky_golay(y, window, order, deriv=0, rate=1):


    Input
    -----
    y         array_like, shape (N,)
              the values of the time history of the signal.
    window    int
              the length of the window. Must be an odd integer number.
    order     int
              the order of the polynomial used in the filtering.
              Must be less then `window` - 1.


    Optional Input
    --------------
    deriv : int
        the order of the derivative to compute
        (default = 0 means only smoothing)
    rate : int
        ??? output will be multiplied by rate**deriv

    Output
    ------
    ndarray : shape(N,)
        the smoothed signal (or it's n-th derivative).


    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.


    References
    ----------
    [1] A. Savitzky, M. J. E. Golay
        Smoothing and Differentiation of Data by Simplified Least Squares
        Procedures. Analytical Chemistry, 1964, 36 (8), pp 1627-1639.
    [2] W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Numerical Recipes 3rd Edition: The Art of Scientific Computing
        Cambridge University Press ISBN-13: 9780521880688


    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(1)
    >>> t = np.linspace(-4, 4, 10)
    >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    >>> from autostring import astr
    >>> print(astr(savitzky_golay(y, window=3, order=1),3,pp=True))
    [' 0.007' ' 0.010' ' 0.022' ' 0.320' ' 0.562' ' 0.609' ' 0.310' ' 0.080' '-0.009' ' 0.007']

    >>> print(astr(savitzky_golay(y, 3, 1, deriv=1),3,pp=True))
    [' 0.000' '-0.050' ' 0.073' ' 0.442' ' 0.295' '-0.304' '-0.368' '-0.120' ' 0.009' ' 0.000']

    >>> print(astr(savitzky_golay(y, 3, 1, deriv=1, rate=10),3,pp=True))
    [' 0.000' '-0.502' ' 0.729' ' 4.416' ' 2.952' '-3.039' '-3.683' '-1.201' ' 0.092' ' 0.000']


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python
    library, Department of Computational Hydrosystems, Helmholtz Centre for
    Environmental Research - UFZ, Leipzig, Germany.

    Copyright (c) 2012-2021 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Matthias Cuntz, Oct 2012
                  - from SciPy cookbook:
                    http://www.scipy.org/Cookbook/SavitzkyGolay
    Modified, Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Apr 2014 - assert
              Matthias Cuntz, Sep 2021 - code refactoring
    """
    #
    # Check input
    try:
        window = np.abs(int(window))
        order  = np.abs(int(order))
    except ValueError as msg:
        raise ValueError("window and order have to be of type int.")
    assert ((window % 2 == 1) and (window >= 1)), "window size must be a positive odd number."
    assert window >= (order + 2), "window is too small for polynomial order."
    order_range = list(range(order+1))
    half_window = (window-1) // 2
    #
    # precompute coefficients (m)
    b = np.mat([[k**i for i in order_range]
                for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    #
    # pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0]  - np.abs(y[1:half_window+1][::-1]   - y[0])
    lastvals  = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    #
    # Filter with pre-computed coefficients m
    return np.convolve(m[::-1], y, mode='valid')


def sg(y, window, order, deriv=0, rate=1):
    """
    Wrapper function for savitzky_golay
    def sg(y, window, order, deriv=0, rate=1):


   Examples
   --------
   >>> import numpy as np
   >>> np.random.seed(1)
   >>> t = np.linspace(-4, 4, 10)
   >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
   >>> from autostring import astr
   >>> print(astr(sg(y, window=3, order=1),3,pp=True))
   [' 0.007' ' 0.010' ' 0.022' ' 0.320' ' 0.562' ' 0.609' ' 0.310' ' 0.080' '-0.009' ' 0.007']

   >>> print(astr(sg(y, 3, 1, deriv=1),3,pp=True))
   [' 0.000' '-0.050' ' 0.073' ' 0.442' ' 0.295' '-0.304' '-0.368' '-0.120' ' 0.009' ' 0.000']

   >>> print(astr(sg(y, 3, 1, deriv=1, rate=10),3,pp=True))
   [' 0.000' '-0.502' ' 0.729' ' 4.416' ' 2.952' '-3.039' '-3.683' '-1.201' ' 0.092' ' 0.000']
    """
    return savitzky_golay(y, window, order, deriv=deriv, rate=rate)


def savitzky_golay2d(z, window, order, deriv=None):
    """
    Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.


    Definition
    ----------
    def savitzky_golay2d(z, window, order, deriv=None)


    Input
    -----
    z         array_like, shape (N,M)
              the value array.
    window    int
              the length of the window. Must be an odd integer number.
    order     int
              the order of the polynomial used in the filtering.
              Must be less then `window` - 1.


    Optional Input
    --------------
    deriv     None=0, 'both=1, 'col'=2, 'row'=3
              None or 0: normal smoothing
              'both' or 1: returns the gradient (first derivatives)
              'col' or 2, or 'row' or 3: indicates the direction of
                                         the derivative

    Output
    ------
    ndarray : shape(N,M) if deriv=0, 2 or 3;
              tuple(shape(N,M),shape(N,M)) if deriv=1
        the smoothed signal, or it's 1st derivative in one direction or both.


    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    Savitsky-Golay filters can also be used to smooth two dimensional
    data affected by noise. The algorithm is exactly the same as for
    the one dimensional case, only the math is a bit more tricky. The
    basic algorithm is as follow:
    1. for each point of the two dimensional matrix extract a sub-matrix,
       centered at that point and with a size equal to an odd number
       "window_size".
    2. for this sub-matrix compute a least-square fit of a polynomial surface,
       defined as
         p(x,y) = a0 + a1*x + a2*y + a3*x2 + a4*y2 + a5*x*y + ... .
       Note that x and y are equal to zero at the central point.
    3. replace the initial central point with the value computed with the fit.


    References
    ----------
    [1] A. Savitzky, M. J. E. Golay
        Smoothing and Differentiation of Data by Simplified Least Squares
        Procedures. Analytical Chemistry, 1964, 36 (8), pp 1627-1639.
    [2] W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Numerical Recipes 3rd Edition: The Art of Scientific Computing
        Cambridge University Press ISBN-13: 9780521880688


    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(1)
    >>> x = np.linspace(-3,3,5)
    >>> y = np.linspace(-3,3,3)
    >>> X, Y = np.meshgrid(x,y)
    >>> Z = np.exp( -(X**2+Y**2))
    >>> Zn = Z + np.random.normal(0, 0.2, Z.shape)
    >>> from autostring import astr
    >>> print(astr(savitzky_golay2d(Zn, window=3, order=2),3,pp=True))
    [[' 0.080' ' 0.052' '-0.106' '-0.165' ' 0.161']
     ['-0.204' ' 0.170' ' 0.664' ' 0.181' ' 0.048']
     [' 0.111' '-0.230' ' 0.094' '-0.170' ' 0.240']]

    >>> s2 = savitzky_golay2d( Zn, 3, 2, deriv=1)
    >>> print(astr(s2[0],3,pp=True))
    [[' 0.156' ' 0.046' '-0.046' ' 0.065' ' 0.257']
     [' 0.305' ' 0.087' '-0.007' '-0.055' ' 0.230']
     [' 0.305' ' 0.126' '-0.039' '-0.158' ' 0.120']]

    >>> s3 = savitzky_golay2d( Zn, 3, 2, deriv='both')
    >>> if np.any(s2[0]!=s3[0]): print('Error 01')
    >>> if np.any(s2[1]!=s3[1]): print('Error 02')

    >>> s2 = savitzky_golay2d( Zn, 3, 2, deriv=2)
    >>> print(astr(s2,3,pp=True))
    [['-0.069' ' 0.510' ' 0.638' ' 0.446' ' 0.128']
     ['-0.102' '-0.047' '-0.018' ' 0.039' ' 0.027']
     [' 0.556' ' 0.251' ' 0.000' ' 0.092' ' 0.111']]

    >>> s3 = savitzky_golay2d( Zn, 3, 2, deriv='col')
    >>> if np.any(s2!=s3): print('Error 03')


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2012-2021 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Matthias Cuntz, Oct 2012
                  - from SciPy cookbook:
                    http://www.scipy.org/Cookbook/SavitzkyGolay
    Modified, Matthias Cuntz, Nov 2012
                  - replaced fftconvolve by convolve because of crash on Mac
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Apr 2014 - assert
              Matthias Cuntz, Sep 2021 - code refactoring
    """
    #
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0
    assert window % 2 == 1, 'window must be odd.'
    assert window**2 >= n_terms, 'order is too high for the window size.'
    half_size = window // 2
    #
    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]
    #
    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window )
    dy = np.tile( ind, [window, 1]).reshape(window**2, )
    #
    # build matrix of system of equation
    A = np.empty( (window**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:, i] = (dx**exp[0]) * (dy**exp[1])
    #
    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs(
        np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs( np.flipud(
        z[-half_size-1:-1, :] ) - band )
    # left band
    band = np.tile( z[:, 0].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs(
        np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:, -1].reshape(-1, 1), [1, half_size] )
    Z[half_size:-half_size, -half_size:] = band + np.abs(
        np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    #
    # pad input array with appropriate values at the four corners
    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs(
        np.flipud(np.fliplr(z[1:half_size+1, 1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs(
        np.flipud(np.fliplr(z[-half_size-1:-1, -half_size-1:-1]) ) - band )
    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs(
        np.flipud(Z[half_size+1:2*half_size+1, -half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:, half_size].reshape(-1, 1)
    Z[-half_size:, :half_size] = band - np.abs(
        np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )
    #
    # solve system and convolve
    # crashed on example at Mac OSX 10.7.5 with Python 2.7.1
    if (deriv is None) | (deriv == 0):
        m = np.linalg.pinv(A)[0].reshape((window, -1))
        # return sps.fftconvolve(Z, m, mode='valid')
        return sps.convolve(Z, m, mode='valid')
    elif (deriv == 'both') | (deriv == 1):
        c = np.linalg.pinv(A)[1].reshape((window, -1))
        r = np.linalg.pinv(A)[2].reshape((window, -1))
        # return (sps.fftconvolve(Z, -r, mode='valid'),
        #         sps.fftconvolve(Z, -c, mode='valid'))
        return (sps.convolve(Z, -r, mode='valid'),
                sps.convolve(Z, -c, mode='valid'))
    elif (deriv == 'col') | (deriv == 2):
        c = np.linalg.pinv(A)[1].reshape((window, -1))
        # return sps.fftconvolve(Z, -c, mode='valid')
        return sps.convolve(Z, -c, mode='valid')
    elif (deriv == 'row') | (deriv == 3):
        r = np.linalg.pinv(A)[2].reshape((window, -1))
        # return sps.fftconvolve(Z, -r, mode='valid')
        return sps.convolve(Z, -r, mode='valid')


def sg2d(*args, **kwargs):
    """
    Wrapper function for savitzky_golay2d
    def savitzky_golay2d(z, window, order, deriv=None):


    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(1)
    >>> x = np.linspace(-3,3,5)
    >>> y = np.linspace(-3,3,3)
    >>> X, Y = np.meshgrid(x,y)
    >>> Z = np.exp( -(X**2+Y**2))
    >>> Zn = Z + np.random.normal(0, 0.2, Z.shape)
    >>> from autostring import astr
    >>> print(astr(sg2d(Zn, window=3, order=2),3,pp=True))
    [[' 0.080' ' 0.052' '-0.106' '-0.165' ' 0.161']
     ['-0.204' ' 0.170' ' 0.664' ' 0.181' ' 0.048']
     [' 0.111' '-0.230' ' 0.094' '-0.170' ' 0.240']]

    >>> s2 = sg2d( Zn, 3, 2, deriv=1)
    >>> print(astr(s2[0],3,pp=True))
    [[' 0.156' ' 0.046' '-0.046' ' 0.065' ' 0.257']
     [' 0.305' ' 0.087' '-0.007' '-0.055' ' 0.230']
     [' 0.305' ' 0.126' '-0.039' '-0.158' ' 0.120']]

    >>> s3 = sg2d( Zn, 3, 2, deriv='both')
    >>> if np.any(s2[0]!=s3[0]): print('Error 01')
    >>> if np.any(s2[1]!=s3[1]): print('Error 02')

    >>> s2 = sg2d( Zn, 3, 2, deriv=2)
    >>> print(astr(s2,3,pp=True))
    [['-0.069' ' 0.510' ' 0.638' ' 0.446' ' 0.128']
     ['-0.102' '-0.047' '-0.018' ' 0.039' ' 0.027']
     [' 0.556' ' 0.251' ' 0.000' ' 0.092' ' 0.111']]

    >>> s3 = sg2d( Zn, 3, 2, deriv='col')
    >>> if np.any(s2!=s3): print('Error 03')
    """
    return savitzky_golay2d(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import matplotlib.pyplot as plt
    # #
    # # 1D
    # t = np.linspace(-4, 4, 500)
    # y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    # ysg = sg(y, window=31, order=4)

    # plt.figure(1)
    # plt.plot(t, y, label='Noisy signal')
    # plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    # plt.plot(t, ysg, 'r', label='Filtered signal')
    # plt.legend()

    # #
    # # 2D
    # x = np.linspace(-3,3,100)
    # y = np.linspace(-3,3,200)
    # X, Y = np.meshgrid(x,y)
    # Z = np.exp( -(X**2+Y**2))
    # Zn = Z + np.random.normal(0, 0.2, Z.shape)
    # Zf = sg2d(Zn, window=29, order=4)

    # plt.matshow(Z)
    # plt.matshow(Zn)
    # plt.matshow(Zf)

    # plt.show()

    # np.random.seed(1)
    # t = np.linspace(-4, 4, 10)
    # y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    # print savitzky_golay(y, window=3, order=1)
    # #[ 0.0067223   0.01048351  0.02186601  0.32004949  0.56168932  0.60865232
    # #  0.31034614  0.08045836 -0.00911575  0.00652022]

    # print savitzky_golay(y, 3, 1, deriv=1)
    # #[  0.00000000e+00  -5.02294966e-02   7.29450608e-02   4.41633401e-01
    # #   2.95156495e-01  -3.03885643e-01  -3.68285622e-01  -1.20119683e-01
    # #   9.21248188e-03  -1.73472348e-18]

    # print savitzky_golay(y, 3, 1, deriv=1, rate=10)
    # #[ -2.77555756e-17  -5.02294966e-01   7.29450608e-01   4.41633401e+00
    # #  2.95156495e+00  -3.03885643e+00  -3.68285622e+00  -1.20119683e+00
    # #  9.21248188e-02  -1.38777878e-17]

    # np.random.seed(1)
    # x = np.linspace(-3,3,5)
    # y = np.linspace(-3,3,3)
    # X, Y = np.meshgrid(x,y)
    # Z = np.exp( -(X**2+Y**2))
    # Zn = Z + np.random.normal(0, 0.2, Z.shape)
    # print sg2d( Zn, window=3, order=2)
    # #[[ 0.07980947  0.05211804 -0.10551094 -0.16506245  0.16105317]
    # # [-0.20426342  0.17009221  0.66354671  0.18103661  0.04756387]
    # # [ 0.11059771 -0.22965185  0.09383012 -0.17015538  0.23954209]]

    # s2 = sg2d( Zn, 3, 2, deriv=1)
    # print s2[0]
    # #[[ 0.15577951  0.04649446 -0.04612122  0.06501884  0.25715001]
    # # [ 0.30484864  0.08679688 -0.00702995 -0.05463381  0.23040468]
    # # [ 0.30484864  0.12564873 -0.03918197 -0.15848258  0.12036621]]
    # print s2[1]
    # #[[-0.06945119  0.5099898   0.63791906  0.44568577  0.12792925]
    # # [-0.10196687 -0.04682891 -0.01845718  0.03876769  0.02683618]
    # # [ 0.55571728  0.25086864  0.          0.09216819  0.11135048]]

    # s3 = sg2d( Zn, 3, 2, deriv='both')
    # if np.any(s2[0]!=s3[0]): print 'Error 01'
    # if np.any(s2[1]!=s3[1]): print 'Error 02'

    # s2 = sg2d( Zn, 3, 2, deriv=2)
    # print s2
    # #[[-0.06945119  0.5099898   0.63791906  0.44568577  0.12792925]
    # # [-0.10196687 -0.04682891 -0.01845718  0.03876769  0.02683618]
    # # [ 0.55571728  0.25086864  0.          0.09216819  0.11135048]]

    # s3 = sg2d( Zn, 3, 2, deriv='col')
    # if np.any(s2!=s3): print 'Error 03'
