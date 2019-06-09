#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
"""
    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2016-2017 Matthias Cuntz - mc (at) macu (dot) de

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

    Modified, MC, Dec 2017 - multinormal
"""

__all__ = ['exponential', 'laplace',
           'gauss', 'norm', 'normal',                           # Normal (Gauss)
           'multigauss', 'multinorm', 'multinormal',            # Multivariate Normal (Gauss)
           'ep', 'sep', 'sep_fs', 'sep_fs_mean', 'sep_fs_std',  # (Skew) Exponential Power
           'st', 'st_fs', 'st_fs_mean', 'st_fs_std', 't']       # (Skew) Student-t


# --------------------------------------------------------------------


def ep(x, loc=0., sca=1., beta=0., sig=None):
    """
        The exponential power distribution with given location, scale, and kurtosis.


        Definition
        ----------
        def ep(x, loc=0., sca=1., beta=0., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        beta       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Exponential power pdf at x


        Examples
        --------
        >>> print(np.allclose(sep(1., 2., 2., 0.5, 2.), sep((1.-2.)/2., xi=0.5, beta=2.)/2.))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 0.), normal(1.3, 0., 1.)))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 1.), laplace(1.3, 0., 1./np.sqrt(2.))))
        True

        >>> print(np.allclose(sep(1.3, 0., np.sqrt(2.), 1., 1.), laplace(1.3, 0., 1.)))
        True


        History
        -------
        Written,  MC, May 2016
    """
    if sig is not None: sca = sig
    return ep01((x-loc)/sca, beta)/sca


def ep01(x, beta=0.):
    """
        The exponential power distribution with given kurtosis parameter, location zero and unit scale.


        Definition
        ----------
        def ep01(x, beta=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        beta       kurtosis parameter


        Output
        ------
        Exponential power pdf with loc=0, sca=1 at x


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    if beta != -1.0:
        b1 = 0.5*(1.0 + beta)
        b3 = 1.5*(1.0 + beta)
        g1 = ss.gamma(b1)
        g3 = ss.gamma(b3)
        # -> 0
        c_beta = (g3/g1)**(1.0/(1.0+beta))
        # -> sqrt(1/12)
        om_beta = np.sqrt(g3)/((1.0+beta)*np.sqrt(g1**3))
    else:
        c_beta  = 0.0
        om_beta = np.sqrt(1.0/12.0)

    # pdf
    if np.abs(beta+1.0) < 0.003: # 2/(1-0.997) ~ 666
        # Uniform between [-x1,x1]
        height = om_beta
        x1 = 0.5/height # int(pdf) = 1 = 2*x1*height
        pdf = np.where((x > (-x1)) & (x < x1), height, 0.0)
        if not np.iterable(x): pdf = pdf[0]
        return pdf
    else:
        return om_beta * np.exp(-c_beta*np.abs(x)**(2.0/(1.0+beta)))


# --------------------------------------------------------------------


def exponential(x, loc=0., sca=1., sig=None):
    """
        Exponential probability density function (pdf).


        Definition
        ----------
        def exponential(x, loc=0., sca=1., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        sig        standard deviation, overwrites scale


        Output
        ------
        Exponential pdf at x


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return exponential01((x-loc)/sca)/sca


def exponential01(x):
    """
        Exponential probability density function (pdf) at location zero and with unit scale.


        Definition
        ----------
        def exponential01(x):


        Input
        -----
        x          array_like quantiles


        Output
        ------
        Exponential pdf with loc=0, sca=1. at x


        History
        -------
        Written,  MC, May 2016
    """

    iexp = np.where(x<0.0, 0., np.exp(-x))
    if np.iterable(x):
        return iexp
    else:
        return float(iexp)


# --------------------------------------------------------------------


def gauss(*args, **kwargs):
    """
        Wrapper for normal

        def normal(x, loc=0., sca=1., sig=None):
    """
    return normal(*args, **kwargs)


# --------------------------------------------------------------------


def laplace(x, loc=0., sca=1., sig=None):
    """
        Laplace probability density function (pdf).


        Definition
        ----------
        def laplace(x, loc=0., sca=1., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        sig        standard deviation, overwrites scale


        Output
        ------
        Laplace pdf at x


        Examples
        --------
        >>> print(str(laplace(0.)))
        0.5

        >>> print(str(laplace(1.) - 0.5/np.e))
        0.0

        >>> print(str(laplace(0., 0., 2.)))
        0.25

        >>> print(str(laplace(0., 2., 2.) - 0.25/np.e))
        0.0

        >>> print(np.allclose(laplace(1., 2., 2.), laplace((1.-2.)/2.)/2.))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig/np.sqrt(2.)
    return laplace01((x-loc)/sca)/sca


def laplace01(x):
    """
        Laplace probability density function (pdf) with at location zero and unit scale.


        Definition
        ----------
        def laplace01(x):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        None


        Output
        ------
        Laplace pdf with loc=0 and sca=1 at x


        Examples
        --------
        >>> print(str(laplace01(0.)))
        0.5

        >>> print(str(laplace01(1.) - 0.5/np.e))
        0.0


        History
        -------
        Written,  MC, May 2016
    """

    return 0.5 * np.exp(-np.abs(x))


# --------------------------------------------------------------------


def multigauss(*args, **kwargs):
    """
        Wrapper for multinormal

        def multinormal(x, loc=0., sca=1., sig=None):
    """
    return multinormal(*args, **kwargs)


# --------------------------------------------------------------------


def multinorm(*args, **kwargs):
    """
        Wrapper for multinormal

        def multinormal(x, loc=0., sca=1., sig=None):
    """
    return multinormal(*args, **kwargs)


# --------------------------------------------------------------------


def multinormal(x, loc=0., sca=1., cov=None):
    """
        Multivariate normal (Gauss) probability density function (pdf).


        Definition
        ----------
        def multinormal(x, loc=0., sca=1., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        array_like location
        sca        array_like scale
                   Scalar: single variance for all dimensions
                   1D: variances for all dimensions
                   2D: covariance matrix
        cov        array_like standard deviation, overwrites scale


        Output
        ------
        Multivariate normal (Gauss) pdf at x


        Examples
        --------
        >>> print(np.allclose(multinormal(0.), 1./np.sqrt(2.*np.pi)))
        True

        >>> print(np.allclose(multinormal(1.), 1./np.sqrt(2.*np.pi*np.e)))
        True

        >>> print(np.allclose(multinormal(0., 0., 2.), 1./np.sqrt(2.*np.pi)/np.sqrt(2.)))
        True

        >>> print(np.allclose(multinormal(0., np.sqrt(2.), 1.)*np.sqrt(2.*np.pi), 1./np.e))
        True

        >>> print(np.allclose(multinormal(1., 2., 3.), multinormal((1.-2.)/np.sqrt(3.))/np.sqrt(3.)))
        True


        History
        -------
        Written,  MC, Dec 2017
    """

    if cov is not None: sca = cov
    k = np.size(x)
    if k==1:
        ix = np.array([x])
    else:
        ix = np.array(x)
    l = np.size(loc)
    if l==1:
        iloc = np.array([loc])
    else:
        iloc = np.array(loc)
    if np.ndim(sca) == 0:
        icov = np.diag(np.full(k, sca)) # np.full((k,k), sca)
    elif np.ndim(sca) == 1:
        icov = np.diag(sca)
    elif np.ndim(sca) == 2:
        icov = sca
    detcov = np.linalg.det(icov)
    assert detcov != 0., 'Degenerate case: determinant of covariance matrix == 0.'
    part1 = np.sqrt((2.*np.pi)**k * detcov)
    # transpose different to matrix formula because of broadcasting rules in Python
    # x (ndim,) or (nsample,ndim); loc scalar or (ndim,); cov (ndim,ndim)
    part2 = -0.5 * np.dot(np.dot((ix-iloc).T, np.linalg.inv(icov)), (ix-iloc))
    return np.exp(part2)/part1


# --------------------------------------------------------------------


def norm(*args, **kwargs):
    """
        Wrapper for normal

        def normal(x, loc=0., sca=1., sig=None):
    """
    return normal(*args, **kwargs)


# --------------------------------------------------------------------


def normal(x, loc=0., sca=1., sig=None):
    """
        Normal (Gauss) probability density function (pdf).


        Definition
        ----------
        def normal(x, loc=0., sca=1., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        sig        standard deviation, overwrites scale


        Output
        ------
        Normal (Gauss) pdf at x


        Examples
        --------
        >>> print(np.allclose(normal(0.), 1./np.sqrt(2.*np.pi)))
        True

        >>> print(np.allclose(normal(1.), 1./np.sqrt(2.*np.pi*np.e)))
        True

        >>> print(np.allclose(normal(0., 0., 2.), 0.5/np.sqrt(2.*np.pi)))
        True

        >>> print(np.allclose(normal(0., np.sqrt(2.), 1.)*np.sqrt(2.*np.pi), 1./np.e))
        True

        >>> print(np.allclose(normal(1., 2., 3.), normal((1.-2.)/3.)/3.))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return normal01((x-loc)/sca)/sca


def normal01(x):
    """
        Normal (Gauss) probability density function (pdf) at location zero and unit scale.


        Definition
        ----------
        def normal01(x):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        None


        Output
        ------
        Normal pdf with loc=0 and sca=1 at x


        Examples
        --------
        >>> print(np.allclose(normal01(0.), 1./np.sqrt(2.*np.pi)))
        True

        >>> print(np.allclose(normal01(1.), 1./np.sqrt(2.*np.pi*np.e)))
        True


        History
        -------
        Written,  MC, May 2016
    """

    return 1./np.sqrt(2.*np.pi) * np.exp(-0.5*x**2)


# --------------------------------------------------------------------


def sep(x, loc=0., sca=1., xi=1., beta=0., sig=None):
    """
        The skew exponential power distribution with given location, scale, skewness, and kurtosis.


        Definition
        ----------
        def sep(x, loc=0., sca=1., xi=1., beta=0., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        beta       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skew exponential power pdf at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.


        Examples
        --------
        >>> print(np.allclose(sep(1., 2., 2., 0.5, 2.), sep((1.-2.)/2., xi=0.5, beta=2.)/2.))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 0.), normal(1.3, 0., 1.)))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 1.), laplace(1.3, 0., 1./np.sqrt(2.))))
        True

        >>> print(np.allclose(sep(1.3, 0., np.sqrt(2.), 1., 1.), laplace(1.3, 0., 1.)))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return sep01((x-loc)/sca, xi, beta)/sca


def sep01(x, xi=1., beta=0.):
    """
        The skew exponential power distribution with given skewness and kurtosis parameters,
        location zero and unit scale.


        Definition
        ----------
        def sep01(x, xi=1., beta=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        xi         skewness parameter
        beta       kurtosis parameter


        Output
        ------
        Skew exponential power pdf with loc=0, sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.


        History
        -------
        Written,  MC, May 2016
    """

    mu  = sep01_fs_mean(xi, beta)
    sig = sep01_fs_std(xi, beta)
    z   = mu + sig*x

    return sig * sep01_fs(z, xi, beta)


# --------------------------------------------------------------------


def sep_fs(x, loc=0., sca=1., xi=1., beta=0., sig=None):
    """
        The skew exponential power distribution Fernandez C & Steel M (1998)
        with given location, scale, skewness, and kurtosis.


        Definition
        ----------
        def sep_fs(x, loc=0., sca=1., xi=1., beta=0., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        beta       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skew exponential power pdf at x after Fernandez and Steel (1998)


        Examples
        --------
        >>> print(np.allclose(sep(1., 2., 2., 0.5, 2.), sep((1.-2.)/2., xi=0.5, beta=2.)/2.))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 0.), normal(1.3, 0., 1.)))
        True

        >>> print(np.allclose(sep(1.3, 0., 1., 1., 1.), laplace(1.3, 0., 1./np.sqrt(2.))))
        True

        >>> print(np.allclose(sep(1.3, 0., np.sqrt(2.), 1., 1.), laplace(1.3, 0., 1.)))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return sep01_fs((x-loc)/sca, xi, beta)/sca


def sep01_fs(x, xi=1., beta=0.):
    """
        The skew exponential power distribution after Fernandez and Steel (1998)
        with location zero and unit scale.

        Definition
        ----------
        def sep01_fs(x, xi=1., beta=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        xi         skewness parameter
        beta       kurtosis parameter


        Output
        ------
        Skew exponential power pdf with loc=0, sca=1 at x after Fernandez and Steel (1998)


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    alpha = np.where(x<0.0, xi, 1./xi)
    if not np.iterable(x): alpha = float(alpha)

    return 2.0/(xi+1./xi) * ep01(alpha*x, beta)


# --------------------------------------------------------------------


def sep_fs_mean(loc=0., sca=1., xi=1., beta=0., sig=None):
    """
        Mean of skew exponential power distribution with given skewness and kurtosis parameters,
        location and scale.


        Definition
        ----------
        def sep_fs_mean(loc=0., sca=1., xi=1., beta=0., sig=None):


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        beta       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Mean of skew exponential power pdf at x


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return sep01_fs_mean(xi, beta) * sca + loc


def sep01_fs_mean(xi=1., beta=0.):
    """
        Mean of skew exponential power distribution with given skewness and kurtosis parameters,
        location zero and unit scale.


        Definition
        ----------
        def sep01_fs_mean(xi=1., beta=0.):


        Optional Input
        --------------
        xi         skewness parameter
        beta       kurtosis parameter


        Output
        ------
        Mean of skew exponential power pdf with loc=0, sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.

        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    if beta != -1.0:
        b1 = 0.5*(1.0 + beta)
        b2 =      1.0 + beta
        b3 = 1.5*(1.0 + beta)
        g1 = ss.gamma(b1)
        g2 = ss.gamma(b2)
        g3 = ss.gamma(b3)
        # -> sqrt(3/4)
        M1 = g2 / np.sqrt(g3*g1)
    else:
        M1 = np.sqrt(0.75)

    if xi != 1:
        return M1*(xi-1./xi)
    else:
        return 0.0


# --------------------------------------------------------------------


def sep_fs_std(sca=1., xi=1., beta=0., sig=None):
    """
        Standard deviation of skew exponential power distribution with given skewness and kurtosis parameters,
        location and scale.


        Definition
        ----------
        def sep_fs_std(sca=1., xi=1., beta=0., sig=None):


        Optional Input
        --------------
        sca        scale
        xi         skewness parameter
        beta       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Standard deviation of skew exponential power pdf at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig
    return sep01_fs_std(xi, beta) * sca


def sep01_fs_std(xi=1., beta=0.):
    """
        Standard deviation of skew exponential power distribution with given skewness and kurtosis parameters,
        location zero and unit scale.


        Definition
        ----------
        def sep01_fs_std(xi=1., beta=0.):


        Optional Input
        --------------
        xi         skewness parameter
        beta       kurtosis parameter


        Output
        ------
        Standard deviation of skew exponential power pdf with loc=0, sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.


        History
        -------
        Written,  MC, May 2016
    """

    mu  = sep01_fs_mean(xi, beta)
    if xi != 1.:
        M1 = mu / (xi-1./xi)
    else:
        M1 = 0.
    M2  = 1.0
    var = (M2-M1**2)*(xi**2+1./xi**2) + 2.0*M1**2 - M2
    if var > 0.0:
        sig  = np.sqrt(var)
    else:
        sig  = 0.0

    return sig


# --------------------------------------------------------------------


def st(x, nu, loc=0., sca=1., xi=1., sig=None):
    """
        The standardised skewed Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def st(x, nu, loc=0., sca=1., xi=1., sig=None):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skewed Student t pdf at x


        Examples
        --------
        >>> print(np.allclose(st(1., 3., 2., 0.5, 2.), st((1.-2.)/0.5, 3., xi=2.)/0.5))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)
    return st01((x-loc)/sca, nu, xi)/sca


def st01(x, nu, xi=1.):
    """
        The standardised skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def st01(x, nu, xi=1.):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        xi       skewness parameter


        Output
        ------
        Skewed Student t pdf with loc=0 and sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.


        History
        -------
        Written,  MC, May 2016
    """

    mu  = st01_fs_mean(nu, xi)
    sca = st01_fs_std(nu, xi) * np.sqrt((nu-2.)/nu)
    z   = mu + sca*x

    return sca * st01_fs(z, nu, xi)


# --------------------------------------------------------------------


def st_fs(x, nu, loc=0., sca=1., xi=1., sig=None):
    """
        The skewed Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def st(x, nu, loc=0., sca=1., xi=1., sig=None):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skewed Student t pdf at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)
    return st01_fs((x-loc)/sca, nu, xi)/sca


def st01_fs(x, nu, xi=1.):
    """
        The skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.

        If xi is not 1 then mean is not zero and standard deviation is not 1.


        Definition
        ----------
        def st01_fs(x, nu, xi=1.):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        xi         skewness parameter


        Output
        ------
        Skewed Student t pdf with loc=0 and sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    alpha = np.where(x<0.0, xi, 1./xi)
    if not np.iterable(x): alpha = float(alpha)

    return 2.0/(xi+1./xi) * t01(alpha*x, nu)


def st_fs_mean(nu, loc=0., sca=1., xi=1., sig=None):
    """
        The mean of the skewed Student t distribution with given degrees of freedom,
        location, scale, and skewness.


        Definition
        ----------
        def st_fs_mean(nu, loc=0., sca=1., xi=1., sig=None):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         skewness parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Mean of skewed Student t pdf


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)
    return st01_fs_mean(nu, xi) * sca + loc


def st01_fs_mean(nu, xi=1.):
    """
        The mean of the skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def st01_fs_mean(nu, xi=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        xi         skewness parameter


        Output
        ------
        Mean of skewed Student t pdf with loc=0 and sca=1 at x


        Literature
        ----------
        Scharnagl B et al. (2015) Inverse modelling of in situ soil water dynamics:
            accounting for heteroscedastic, autocorrelated, and non-Gaussian distributed residuals.
            Hydrology And Earth System Sciences Discussions 12, 2155-2199.


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    if nu <= 1.:
        return np.inf
    else:
        mu  = 2.0 * (xi-1./xi) * ss.gamma(0.5*(nu+1.0)) / ss.gamma(0.5*nu)
        mu *= (nu-2.0)/(nu-1.0) * nu/(nu-2.) / np.sqrt(np.pi*nu)

    return mu


def st_fs_std(nu, sca=1., xi=1., sig=None):
    """
        The standard deviation of the skewed Student t distribution with given degrees of freedom,
        scale, and skewness.


        Definition
        ----------
        def st_fs_std(nu, sca=1., xi=1., sig=None):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        sca        scale
        xi         skewness parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Standard deviation of skewed Student t pdf


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)
    return st01_fs_std(nu, xi) * sca


def st01_fs_std(nu, xi=1.):
    """
        The standard deviation of the skewed Student t distribution with given degrees of freedom and skewness,
        and unit scale.


        Definition
        ----------
        def st01_fs_std(nu, xi=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        xi       skewness parameter


        Output
        ------
        Standard deviation of skewed Student t pdf with sca=1 at x


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    if nu <= 2.:
        return np.inf
    else:
        mu = st01_fs_mean(nu, xi)
        if xi != 1.:
            M1 = mu / (xi-1./xi)
        else:
            M1 = 0.
        M2  = nu/(nu-2.)
        var = (M2-M1**2)*(xi**2+1./xi**2)+2.*M1**2-M2
        if var > 0.0:
            return np.sqrt(var)
        else:
            return 0.0


# --------------------------------------------------------------------


def t(x, nu, loc=0., sca=1., sig=None):
    """
        The Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def t(x, nu, loc=0., sca=1., sig=None):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        sig        standard deviation, overwrites scale


        Output
        ------
        Student t pdf at x


        Examples
        --------
        >>> print(np.allclose(st(1., 3., 2., 0.5, 2.), st((1.-2.)/0.5, 3., xi=2.)/0.5))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)
    return t01((x-loc)/sca, nu)/sca


def t01(x, nu):
    """
        The Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def t01(x, nu):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Output
        ------
        Student t pdf with loc=0 and sca=1 at x


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    c = ss.gamma(0.5*(nu+1.0)) / (ss.gamma(0.5*nu) * np.sqrt(np.pi*nu))
    return c * (1.0 + x**2/nu)**(-0.5*(nu+1.0))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # import scipy.special as spec
    # import scipy.stats as sp
    # import jams

    # nn = 10000000
    # xx = np.arange(nn)/float(nn) * 30.*100. - 15.*100.

    # loc=1.1
    # sca=2.2
    # sig=sca
    # xi=3.3
    # beta=0.5
    # nu=4.4

    # dx = xx[2]-xx[1]

    # __all__ = ['ep', 'exponential', 'laplace', 'normal', 'sep', 'st', 't']

    # print('loc/sca/sig/sig**2/xi/beta/nu:', loc, sca, sig, sig**2, xi, beta, nu)
    # for dd in __all__:
    #     # use sca
    #     if dd == 'ep':
    #         pdf = jams.distributions.ep(xx, loc, sca, beta)
    #     elif dd == 'exponential':
    #         pdf = jams.distributions.exponential(xx, loc, sca)
    #     elif dd == 'laplace':
    #         pdf = jams.distributions.laplace(xx, loc, sca)
    #     elif dd == 'normal':
    #         pdf = jams.distributions.normal(xx, loc, sca)
    #     elif dd == 'sep':
    #         pdf = jams.distributions.sep(xx, loc, sca, xi, beta)
    #     elif dd == 'st':
    #         pdf = jams.distributions.st(xx, nu, loc, sca, xi)
    #     elif dd == 't':
    #         pdf = jams.distributions.t(xx, nu, loc, sca)

    #     dmean = np.sum(xx * pdf * dx)
    #     dvar  = np.sum((xx-dmean)**2 * pdf * dx)
    #     dstd  = np.sqrt(dvar)
    #     dskew = np.sum(((xx-dmean)/dstd)**3 * pdf * dx)
    #     dkurt = np.sum(((xx-dmean)/dstd)**4 * pdf * dx)
    #     print(dd, ':', dmean, dstd, dvar, dskew, dkurt)

    #     # use sig
    #     if dd == 'ep':
    #         pdf = jams.distributions.ep(xx, loc, beta=beta, sig=sig)
    #     elif dd == 'exponential':
    #         pdf = jams.distributions.exponential(xx, loc, sig=sig)
    #     elif dd == 'laplace':
    #         pdf = jams.distributions.laplace(xx, loc, sig=sig)
    #     elif dd == 'normal':
    #         pdf = jams.distributions.normal(xx, loc, sig=sig)
    #     elif dd == 'sep':
    #         pdf = jams.distributions.sep(xx, loc, xi=xi, beta=beta, sig=sig)
    #     elif dd == 'st':
    #         pdf = jams.distributions.st(xx, nu, loc, xi=xi, sig=sig)
    #     elif dd == 't':
    #         pdf = jams.distributions.t(xx, nu, loc, sig=sig)

    #     dmean = np.sum(xx * pdf * dx)
    #     dvar  = np.sum((xx-dmean)**2 * pdf * dx)
    #     dstd  = np.sqrt(dvar)
    #     dskew = np.sum(((xx-dmean)/dstd)**3 * pdf * dx)
    #     dkurt = np.sum(((xx-dmean)/dstd)**4 * pdf * dx)
    #     print(dd, ':', dmean, dstd, dvar, dskew, dkurt)

    #     # # xx   = xx[xx>0.]
    #     # # pdf  = jams.distributions.t(xx, nu, loc, sca)
    #     # # M1   = 2.* np.sum(xx * pdf * dx)
    #     # # M2   = nu/(nu-2.) * sca**2
    #     # # dvarfs = (M2-M1**2)*(xi**2+1./xi**2)+2.*M1**2-M2
    #     # # # print('Fernandez var:', dvarfs)

    #     # mu_shmc  = 2.0 * (xi-1./xi) * spec.gamma(0.5*(nu+1.0)) / spec.gamma(0.5*nu)
    #     # mu_shmc *= (nu-2.0)/(nu-1.0)
    #     # mu_shmc *= nu/(nu-2.) / np.sqrt(np.pi*nu)
    #     # mu_shmc *= sca

    #     # # var_shmc = -mu_shmc**2 + (xi**3+1./xi**3)/(xi+1./xi)
    #     # # # print('mu_shmc:', mu_shmc, var_shmc)

    #     # if xi != 1.:
    #     #     M1 = mu_shmc / (xi-1./xi)
    #     # else:
    #     #     M1 = 0.
    #     # M2   = nu/(nu-2.) * sca**2
    #     # dvar_mc = (M2-M1**2)*(xi**2+1./xi**2)+2.*M1**2-M2
    #     # print('dvar_mc:', mu_shmc, dvar_mc)
