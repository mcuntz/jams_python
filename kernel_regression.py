#!/usr/bin/env python
import numpy as np
import scipy.optimize as opt # fmin_tnc
from division import *

def kernel_regression(x, y, h=None, linear=False):
    """
        Multi-dimensional non-parametric regression with either kernel regression
        or local linear regression using the Quartic kernel.

        Optimal bandwidth can be estimated by cross-validation.

        Definition
        ----------
        def kernel_regression(x, y, linear=False, h=None):


        Input
        -----
        x          ndarray(n,k) of independent values
        y          array(n) of dependent values


        Optional Input
        --------------
        linear     False: local constant kernel regression
                   True:  local linear regression using the Quartic kernel
        h          None:  determine optimal h
                   float > 0: use for calculating regression values

        Output
        ------
        Fitted values at x


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> h = kernel_regression_h(x,y)
        >>> print kernel_regression(x,y,h)
        [ 0.70404103  0.01294352  1.04548174  0.5591734   0.27353146  0.31533275
          0.51589033  0.78084832  1.07240574  1.36797238]

        >>> print kernel_regression(x,y)
        [ 0.70404103  0.01294352  1.04548174  0.5591734   0.27353146  0.31533275
          0.51589033  0.78084832  1.07240574  1.36797238]

        >>> print kernel_regression(x,y,h=0.5)
        [ 0.70404103  0.01294352  0.98649551  0.56013557  0.32292917  0.34133147
          0.52673436  0.78545776  1.07487758  1.30446591]


        History
        -------
        Written in Matlab Yingying Dong, Boston College, July, 2008
        Transferred to Python, MC, Jun 2011
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
    if linear:
        raise ValueError('linear=True not working yet')
    #
    # determine h
    if h == None:
        h = kernel_regression_h(xx,y,linear=linear)
    #
    # Calc regression
    out = np.empty(n)
    for i in xrange(n):
        dis    = xx - xx[i,:]
        u      = division(dis, np.std(dis,axis=0,ddof=1)/h, np.nan)
        Kernel = 15./16. * (1.-u**2)**2 * (np.abs(u)<=1.)
        w      = np.prod(Kernel,1)
        sw     = np.sum(w)
        if linear:
            t1     = np.transpose(dis * w[:,np.newaxis])
            t2     = np.sum(t1,1)
            lhs    = sw*np.dot(t1,dis) - np.outer(t2,t2)
            rhs    = sw*np.dot(t1,y)  - np.dot(np.outer(t2,w),y)
            out[i] = np.nan
            try:
                # it works until here but the linalg is not working yet
                # in Matlab it is: b = lhs\rhs
                b      = np.transpose(np.dot(np.matrix(lhs).I,rhs))
                out[i] = division(np.dot(w,(yo-np.dot(dis,b))), sw, np.nan)
            except:
                pass
        else:
            out[i] = division(np.dot(w,y), sw, np.nan)

    return out



def kernel_regression_h(x, y, h0=None, linear=False):
    """
        Optimal bandwidth for a multi-dimensional non-parametric regression using cross-validation.
        
        Either kernel regression or local linear regression using the Quartic kernel.

        Definition
        ----------
        def kernel_regression_h(x, y, linear=False, h0=None):


        Input
        -----
        x          ndarray(n,k) of independent values
        y          array(n) of dependent values


        Optional Input
        --------------
        linear     False: local constant kernel regression
                   True:  local linear regression using the Quartic kernel
        h0         None:  determine starting value of h in cross-validation
                   float > 0: use as starting h in cross-validation

        Output
        ------
        Optimal bandwidth


        Examples
        --------
        >>> import numpy as np
        >>> x = np.zeros((10,2))
        >>> x[:,0] = np.arange(10,dtype=np.float)/9.
        >>> x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
        >>> y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
        >>> print kernel_regression_h(x,y)
        0.365148365498

        >>> print kernel_regression_h(x,y,h0=0.5)
        0.472301961489

        History
        -------
        Written in Matlab Yingying Dong, Boston College, July, 2008
        Transferred to Python, MC, Jun 2011 - changed initial search of h0
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
    if linear:
        raise ValueError('linear=True not working yet')
    #
    # determine h0
    eps = np.finfo(np.float).eps
    if h0 == None:
        h00    = 10.
        nhLin  = 9
        nrange = np.arange(nhLin,dtype=np.float)/(np.float(nhLin)-1.)
        # search in interval [h00-10^(-k+1),h00+(10^-k+2)]
        niter  = 2
        for k in xrange(niter):
            n5   = 10.**(-(k-1))
            hLin = np.maximum(h00 - n5 + (2.*n5*nrange) * np.where(k==0, 0.4, 1.), eps)
            mseFun = np.zeros(nhLin)
            for i in xrange(nhLin):
                mseFun[i] = kernel_regression_MSE(hLin[i],xx,y,linear)
                if np.sum(mseFun==np.min(mseFun))==1:
                    h00 = hLin[mseFun==np.min(mseFun)]
                else:
                    h00 = np.mean(hLin[mseFun==np.min(mseFun)])
        h0 = np.float(h00)
        h0min = h0 - 10.**(-niter)
        h0max = h0 + 10.**(-niter)
    else:
        h0 = np.float(h0)
        if (h0<0.) | (h0>10.):
            raise ValueError('h0 should be 0<=h0<=10')
        h0min = eps
        h0max = 10.
    #
    # Find the optimal h
    h, nfeval, rc  = opt.fmin_tnc(kernel_regression_MSE, [h0], bounds=[(h0min,h0max)],
                                  args=(xx, y, linear), approx_grad=True, disp=False,
                                  maxfun=1000, xtol=1e-10, ftol=1e-10)
    
    return np.float(h)


def kernel_regression_MSE(h, x, y, linear=False):
    """
        Helper function for kernel_regression_h that returns
        mean square error of obs and the cross-validate kernels,
        i.e. with one obs removed.
    """
    n, k = np.shape(x)
    crossval  = np.ones(n)*np.nan
    # calculate the cross validation criterion function
    for i in xrange(n):
        # remove ith observation
        xo = np.delete(x,i,axis=0)
        yo = np.delete(y,i)
        # calculate kernel function
        dis    = xo - x[i,:]
        u      = division(dis, np.std(dis,axis=0,ddof=1)/h, np.nan)
        Kernel = 15./16. * (1.-u**2)**2 * (np.abs(u)<=1.)
        # calculate weights
        w      = np.prod(Kernel,1)
        sw     = np.sum(w)
        # calculate yhat, using kernel or local linear regression
        if linear:
            t1     = np.transpose(dis * w[:,np.newaxis])
            t2     = np.sum(t1,1)
            lhs    = sw*np.dot(t1,dis) - np.outer(t2,t2)
            rhs    = sw*np.dot(t1,yo)  - np.dot(np.outer(t2,w),yo)
            crossval[i] = np.nan
            try:
                # it works until here but the linalg is not working yet
                # in Matlab it is: b = lhs\rhs
                b           = np.transpose(np.dot(np.matrix(lhs).I,rhs))
                crossval[i] = division(np.dot(w,(yo-np.dot(dis,b))), sw, np.nan)
            except:
                pass
        else:
            crossval[i] = division(np.dot(w,yo), sw, np.nan)
    # MSE
    flag     = np.isfinite(crossval)
    crossval = np.ma.array(crossval, mask=(~flag))
    yy       = np.ma.array(y, mask=(~flag))
    MSE = np.ma.mean((yy-crossval)*(yy-crossval))
    MSE = np.ma.filled(MSE, fill_value=np.finfo(np.float).max)

    return MSE


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # x = np.zeros((10,2))
    # x[:,0] = np.arange(10,dtype=np.float)/9.
    # x[:,1] = 1./(np.arange(10,dtype=np.float)/9.+0.1)
    # y      = 1. + x[:,0]**2 - np.sin(x[:,1])**2
    # h = kernel_regression_h(x,y)
    # yy = kernel_regression(x,y,h,linear=False)
    # print h
    # print y
    # print yy
    # print kernel_regression(x,y)
    # print kernel_regression_h(x,y,h0=0.5)
    # print kernel_regression(x,y,h=0.5)
