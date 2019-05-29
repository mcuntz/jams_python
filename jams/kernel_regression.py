#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.optimize as opt # fmin_tnc
from jams.division import division

def kernel_regression(x, y, h=None, silverman=False, xout=None):
    """
        Multi-dimensional non-parametric kernel regression.

        Optimal bandwidth can be estimated by cross-validation
        or by using Silverman''s rule-of-thumb.


        Definition
        ----------
        def kernel_regression(x, y, h=None, silverman=False):


        Input
        -----
        x          ndarray(n,k) of independent values
        y          array(n) of dependent values


        Optional Input
        --------------
        h          None:  determine optimal h
                   float > 0: use for calculating regression values
        silverman  False: determine h via cross-validation
                   True:  use Silverman''s rule-of-thumb


        Output
        ------
        Fitted values at x


        References
        ----------
        Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
            In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
            application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> h = kernel_regression_h(x,y)
        >>> print(np.allclose(h, [0.172680, 9.516907], atol=0.0001))
        True

        >>> yk = kernel_regression(x,y,h)
        >>> print(np.allclose(yk[0:6], [0.52241, 0.52570, 0.54180, 0.51781, 0.47644, 0.49230], atol=0.0001))
        True

        >>> yk = kernel_regression(x,y)
        >>> print(np.allclose(yk[0:6], [0.52241, 0.52570, 0.54180, 0.51781, 0.47644, 0.49230], atol=0.0001))
        True

        >>> h = kernel_regression_h(x,y,silverman=True)
        >>> print(np.allclose(h, [0.229190, 1.903381], atol=0.0001))
        True

        >>> yk = kernel_regression(x,y,h)
        >>> print(np.allclose(yk[0:6], [0.691153, 0.422809, 0.545844, 0.534315, 0.521494, 0.555426], atol=0.0001))
        True

        >>> ss = np.shape(x)
        >>> nn = 5
        >>> xx = np.empty((nn,ss[1]))
        >>> xx[:,0] = np.amin(x[:,0]) + (np.amax(x[:,0])-np.amin(x[:,0])) * np.arange(nn,dtype=np.float)/np.float(nn)
        >>> xx[:,1] = np.amin(x[:,1]) + (np.amax(x[:,1])-np.amin(x[:,1])) * np.arange(nn,dtype=np.float)/np.float(nn)
        >>> yk = kernel_regression(x,y,h,xout=xx)
        >>> print(np.allclose(yk, [0.605485, 0.555235, 0.509529, 0.491191, 0.553325], atol=0.0001))
        True


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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2012 - inspired by Matlab routines
                                 of Yingying Dong, Boston College and Yi Cao, Cranfield University
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    #
    # Check input
    ss = np.shape(x)
    n  = ss[0]
    assert n == np.size(y), 'size(x,0) != size(y): '+str(n)+' != '+str(np.size(y))
    if np.size(ss) == 1: # to deal with 1d-arrays
        xx = x[:,np.newaxis]
    else:
        xx = x
    ss = np.shape(xx)
    d  = ss[1]
    #
    # determine h
    if h is None:
        hh = kernel_regression_h(xx,y,silverman=silverman)
    else:
        if np.size(np.shape(h))==0:
            hh = np.repeat(h,d)
        else:
            hh = np.array(h)
        assert np.size(hh) == d, 'size(h) must be 1 or size(x,1): '
    #
    # Calc regression
    if xout is None:
        xxout = xx
    else:
        if np.size(np.shape(xout)) == 1:
            xxout = xout[:,np.newaxis]
        else:
            xxout = xout
    ssout = np.shape(xxout)
    nout  = ssout[0]
    dout  = ssout[1]
    assert d == dout, 'size(x,1) != size(xout,1): '+str(d)+' != '+str(dout)
    # allocate output
    out = np.empty(nout)
    # Loop through each regression point
    for i in range(nout):
        # scaled deference from regression point
        z = (xx - xxout[i,:]) / hh
        # nadaraya-watson estimator of gaussian multivariate kernel
        out[i] = nadaraya_watson(z, y)
    #
    return out


def kernel_regression_h(x, y, silverman=False):
    """
        Optimal bandwidth for multi-dimensional non-parametric kernel regression
        using cross-validation or Silverman''s rule-of-thumb.

        Definition
        ----------
        def kernel_regression_h(x, y):


        Input
        -----
        x          ndarray(n,k) of independent values
        y          array(n) of dependent values


        Optional Input
        --------------
        silverman  False: determine h via cross-validation
                   True:  use Silverman''s rule-of-thumb


        Output
        ------
        Optimal bandwidth. If multidimensional regression then h is vector,
        assuming diagonal bandwith matrix.


        References
        ----------
        Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
            In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
            application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> h = kernel_regression_h(x,y)
        >>> print(np.allclose(h, [0.172680, 9.516907], atol=0.0001))
        True

        >>> h = kernel_regression_h(x,y,silverman=True)
        >>> print(np.allclose(h, [0.229190, 1.903381], atol=0.0001))
        True


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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2012 - inspired by Matlab routines
                                 of Yingying Dong, Boston College and Yi Cao, Cranfield University
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
                  MC, Jan 2018 - bug in boot_h: x.size->x.shape[0]
    """
    #
    # Check input
    ss = np.shape(x)
    n  = ss[0]
    assert x.shape[0] == np.size(y), 'size(x,0) != size(y): '+str(n)+' != '+str(np.size(y))
    if np.ndim(x) == 1: # to deal with 1d-arrays
        xx = x[:,np.newaxis]
    else:
        xx = x
    d = xx.shape[1]
    #
    # Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    # Very similar to stats.gaussian_kde
    h = (4./np.float(d+2)/np.float(n))**(1./np.float(d+4)) * np.std(xx,axis=0,ddof=1)
    #
    if not silverman:
        # Find the optimal h
        bounds = [(0.2*i,5.0*i) for i in h]
        if n<=100:
            h, nfeval, rc = opt.fmin_tnc(cross_valid_h, h, bounds=bounds,
                                         args=(xx, y), approx_grad=True, disp=False,
                                         maxfun=1000, xtol=1e-10, ftol=1e-10)
        else:
            h, nfeval, rc = opt.fmin_tnc(boot_h, h, bounds=bounds,
                                         args=(xx, y), approx_grad=True, disp=False,
                                         maxfun=1000, xtol=1e-10, ftol=1e-10)
    #
    return h


def cross_valid_h(h, x, y):
    """
        Helper function that calculates cross-validation function for the
        Nadaraya-Watson estimator, which is basically the mean square error
        where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).
    """
    n = x.shape[0]
    # allocate output
    out = np.empty(n)
    # Loop through each regression point
    for i in range(n):
        # all-1 points
        xx     = np.delete(x,i,axis=0)
        yy     = np.delete(y,i,axis=0)
        z      = (xx - x[i,:]) / h
        out[i] = nadaraya_watson(z, yy)
    cv = np.sum((y-out)**2) / np.float(n)
    #
    return cv


def boot_h(h, x, y):
    """
        Helper function that calculates bootstrap function for the
        Nadaraya-Watson estimator, which is basically the mean square error
        where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).
        This does basically cross_valid_h for 100 random points.
    """
    n   = 100
    ind = np.random.randint(x.shape[0],size=n)
    # allocate output
    out = np.empty(n)
    # Loop through each bootstrap point
    for i in range(n):
        # all-1 points
        xx     = np.delete(x,i,axis=0)
        yy     = np.delete(y,i,axis=0)
        z      = (xx - x[i,:]) / h
        out[i] = nadaraya_watson(z, yy)
    cv = np.sum((y[ind]-out)**2) / np.float(n)
    #
    return cv


def nadaraya_watson(z, y):
    """
        Helper function that calculates the Nadaraya-Watson estimator for a given kernel.
        Until now there is only the gaussian kernel.
    """
    kerf   = (1./np.sqrt(2.*np.pi)) * np.exp(-0.5*z*z)
    w      = np.prod(kerf,1)
    out    = division(np.dot(w,y), np.sum(w), np.nan)
    #
    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    # nn = 1000
    # x = np.zeros((nn,2))
    # x[:,0] = np.arange(nn,dtype=np.float)/float(nn-1)
    # x[:,1] = 1./(x[:,0]+0.1)
    # y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
    # h = kernel_regression_h(x,y)
    # print(h)
    # yy = kernel_regression(x,y,h,xout=x)
    # print(yy[0], yy[-1])
    # print(kernel_regression(x,y))
    # h = kernel_regression_h(x,y,silverman=True)
    # print(h)
    # print(kernel_regression(x,y,h))
    # ss = np.shape(x)
    # nn = 5
    # xx = np.empty((nn,ss[1]))
    # xx[:,0] = np.amin(x[:,0]) + (np.amax(x[:,0])-np.amin(x[:,0])) * np.arange(nn,dtype=np.float)/np.float(nn)
    # xx[:,1] = np.amin(x[:,1]) + (np.amax(x[:,1])-np.amin(x[:,1])) * np.arange(nn,dtype=np.float)/np.float(nn)
    # print(kernel_regression(x,y,h,xout=xx))

