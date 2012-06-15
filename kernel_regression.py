#!/usr/bin/env python
import numpy as np
import scipy.optimize as opt # fmin_tnc
from division import *

def kernel_regression(x, y, h=None, silverman=False):
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
        Hardle, W., & Muller, M. (2000). Multivariate and semiparametric kernel regression.
            In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
            application (pp. 357–392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> h = kernel_regression_h(x,y)
        >>> print h
        [ 0.17267991  9.51690721]

        >>> print kernel_regression(x,y,h)
        [ 0.52240942  0.52569869  0.54179588  0.51780751  0.47644246  0.49230198
          0.60344807  0.77747176  0.95450466  1.09603888]

        >>> print kernel_regression(x,y)
        [ 0.52240942  0.52569869  0.54179588  0.51780751  0.47644246  0.49230198
          0.60344807  0.77747176  0.95450466  1.09603888]

        >>> h = kernel_regression_h(x,y,silverman=True)
        >>> print h
        [ 0.22919046  1.90338144]

        >>> print kernel_regression(x,y,h)
        [ 0.69115273  0.42280858  0.54584447  0.53431539  0.52149406  0.55542563
          0.64206536  0.76189995  0.88777986  1.00014619]


        History
        -------
        Written in MC, Jun 2011 - inspired by Matlab routine of Yingying Dong, Boston College
    """
    #
    # Check input
    ss = np.shape(x)
    n  = ss[0]
    if n != np.size(y):
        raise ValueError('size(x,0) != size(y): '+str(n)+' != '+str(np.size(y)))
    if np.size(ss) == 1: # to deal with 1d-arrays
        xx = x[:,np.newaxis]
    else:
        xx = x
    ss = np.shape(xx)
    d  = ss[1]
    #
    # determine h
    if h == None:
        hh = kernel_regression_h(xx,y,silverman=silverman)
    else:
        if np.size(np.shape(h))==0:
            hh = np.repeat(h,d)
        else:
            hh = np.array(h)
        if (np.size(hh)!=d):
            raise ValueError('size(h) must be 1 or size(x,1): ')
    #
    # Calc regression
    # Scaling first, make diagonal matrix
    hh1 = 1. / np.prod(hh)
    # allocate output
    out = np.empty(n)
    # Loop through each regression point
    for i in xrange(n):
        # scaled deference from regression point
        z      = (xx - xx[i,:]) / hh
        # gaussian multivariate kernel
        kerf   = (1./np.sqrt(2.*np.pi)) * np.exp(-0.5*z*z)
        # multiplicative kernel
        w      = np.prod(kerf,1)
        out[i] = division(np.dot(w,y), np.sum(w), np.nan)
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
        Hardle, W., & Muller, M. (2000). Multivariate and semiparametric kernel regression.
            In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
            application (pp. 357–392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> print kernel_regression_h(x,y)
        [ 0.17267991  9.51690721]

        >>> print kernel_regression_h(x,y,silverman=True)
        [ 0.22919046  1.90338144]


        History
        -------
        Written in MC, Jun 2011 - inspired by Matlab routine of Yingying Dong, Boston College
    """
    #
    # Check input
    ss = np.shape(x)
    n  = ss[0]
    if n != np.size(y):
        raise ValueError('size(x,0) != size(y): '+str(n)+' != '+str(np.size(y)))
    if np.size(ss) == 1: # to deal with 1d-arrays
        xx = x[:,np.newaxis]
    else:
        xx = x
    ss = np.shape(xx)
    d  = ss[1]
    #
    # Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    # Very similar to stats.gaussian_kde
    h = (4./np.float(d+2)/np.float(n))**(1./np.float(d+4)) * np.std(xx,axis=0,ddof=1)
    #
    if not silverman:
        # Find the optimal h
        bounds = [(0.2*i,5.0*i) for i in h]
        h, nfeval, rc  = opt.fmin_tnc(cross_valid_h, h, bounds=bounds,
                                      args=(xx, y), approx_grad=True, disp=False,
                                      maxfun=1000, xtol=1e-10, ftol=1e-10)
    #
    return h


def cross_valid_h(h, x, y):
    """
        Helper function that calculates cross-validation function for the
        Nadaraya-Watson estimator, which is basically the mean square error
        where model estimate is replaced by the jackknife estimate (Hardle et al. 2000).
    """
    n = np.size(x[:,0])
    # allocate output
    out = np.empty(n)
    # Loop through each regression point
    for i in xrange(n):
        # all-1 points
        xx     = np.delete(x,i,axis=0)
        yy     = np.delete(y,i,axis=0)
        z      = (xx - x[i,:]) / h
        kerf   = (1./np.sqrt(2.*np.pi)) * np.exp(-0.5*z*z)
        w      = np.prod(kerf,1)
        out[i] = division(np.dot(w,yy), np.sum(w), np.nan)
    cv = np.sum((y-out)**2) / np.float(n)
    #
    return cv


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # x = np.zeros((10,2))
    # x[:,0] = np.arange(10,dtype=np.float)/9.
    # x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
    # y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
    # print y
    # h = kernel_regression_h(x,y)
    # print h
    # print kernel_regression(x,y,h)
    # print kernel_regression(x,y)
    # h = kernel_regression_h(x,y,silverman=True)
    # print h
    # print kernel_regression(x,y,h)
