#!/usr/bin/env python
from __future__ import print_function
import numpy as np
"""
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

    Copyright 2016 Matthias Cuntz
"""

__all__ = ['ep', 'exponential', 'gauss', 'laplace', 'normal',
           'sep', 'sep_mean', 'sep_std', 'ssep',
           'ssstudentt', 'sstudentt', 'sstudentt_mean', 'sstudentt_std', 'studentt']

def ep(x, loc=0., sca=1., kurt=0., sig=None):
    """
        The exponential power distribution with given location, scale, and kurtosis.


        Definition
        ----------
        def ep(x, loc=0., sca=1., kurt=0., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        kurt       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Exponential power pdf at x
        

        Examples
        --------
        >>> print(np.allclose(sep(1., 2., 2., 0.5, 2.), sep((1.-2.)/2., skew=0.5, kurt=2.)/2.))
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
    if sig is None:
        return ep01((x-loc)/sca, kurt)/sca
    else:
        sca = sig
        return ep01((x-loc)/sca, kurt)/sca


def ep01(x, kurt=0.):
    """
        The exponential power distribution with given skewness and kurtosis, location zero and unit scale.


        Definition
        ----------
        def sep01(x, kurt=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        kurt       kurtosis parameter
        

        Output
        ------
        Exponential power pdf with loc=0, sca=1 at x


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    beta = kurt
    
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


def exponential(x, loc=0., sca=1., theta=1., sig=None):
    """
        Exponential probability density function (pdf).


        Definition
        ----------
        def exponential(x, loc=0., sca=1., theta=1., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        theta      duration (shape) parameter
        sig        standard deviation, overwrites scale
        

        Output
        ------
        Exponential pdf at x


        History
        -------
        Written,  MC, May 2016
    """

    if sig is None:
        return exponential01((x-loc)/sca)/sca
    else:
        sca = sig
        return exponential01((x-loc)/sca)/sca


def exponential01(x, theta=1.):
    """
        Exponential probability density function (pdf) at location zero and with unit scale.


        Definition
        ----------
        def exponential01(x, theta=1.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        theta      duration (shape) parameter
        

        Output
        ------
        Exponential pdf with loc=0, sca=1. at x


        History
        -------
        Written,  MC, May 2016
    """

    return np.exp(-x/theta)/theta

    
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

    if sig is None:
        return laplace01((x-loc)/sca)/sca
    else:
        sca = sig/np.sqrt(2.)
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

        >>> print(np.allclose(normal(1., 2., 2.), normal((1.-2.)/2.)/2.))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is None:
        return normal01((x-loc)/sca)/sca
    else:
        sca = sig
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

    return 1./np.sqrt(2.*np.pi) * np.exp(-0.5*x*x)

    
# --------------------------------------------------------------------


def sep(x, loc=0., sca=1., skew=1., kurt=0.):
    """
        The skew exponential power distribution with given location, scale, skewness, and kurtosis.


        Definition
        ----------
        def sep(x, loc=0., sca=1., skew=1., kurt=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        kurt       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skew exponential power pdf at x
        

        Examples
        --------
        >>> print(np.allclose(sep(1., 2., 2., 0.5, 2.), sep((1.-2.)/2., skew=0.5, kurt=2.)/2.))
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

    return sep01((x-loc)/sca, skew, kurt)/sca


def sep01(x, skew=1., kurt=0.):
    """
        The skew exponential power distribution with given skewness and kurtosis, location zero and unit scale.


        Definition
        ----------
        def sep01(x, skew=1., kurt=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        skew       skewness parameter
        kurt       kurtosis parameter
        

        Output
        ------
        Skew exponential power pdf with loc=0, sca=1 at x
        

        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  MC, May 2016
    """

    alpha = np.where(x<0.0, skew, 1./skew)
    if not np.iterable(x): alpha = alpha[0]

    return 2.0/(skew+1./skew) * ep01(alpha*x, kurt)


def sep_mean(loc=0., sca=1., skew=1., kurt=0.):
    """
        Mean of skew exponential power distribution with given skewness and kurtosis, location and scale.


        Definition
        ----------
        def sep_mean(loc=0., sca=1., skew=1., kurt=0.):


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        kurt       kurtosis parameter
        

        Output
        ------
        Mean of skew exponential power pdf at x


        History
        -------
        Written,  MC, May 2016
    """

    return sep01_mean(skew, kurt) * sca + loc


def sep01_mean(skew=1., kurt=0.):
    """
        Mean of skew exponential power distribution with given skewness and kurtosis, location zero and unit scale.


        Definition
        ----------
        def sep01_mean(skew=1., kurt=0.):


        Optional Input
        --------------
        skew       skewness parameter
        kurt       kurtosis parameter
        

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

    beta = kurt # notation of Schoups and Vrugt (2010)    
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
        M1      = np.sqrt(0.75)

    if skew != 1:
        return M1*(skew-1./skew)
    else:
        return 0.0


def sep_std(sca=1., skew=1., kurt=0.):
    """
        Standard deviation of skew exponential power distribution with given skewness and kurtosis, location and scale.


        Definition
        ----------
        def sep_std(sca=1., skew=1., kurt=0.):


        Optional Input
        --------------
        sca        scale
        skew       skewness parameter
        kurt       kurtosis parameter
        

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

    return sep01_std(skew, kurt) * sca


def sep01_std(skew=1., kurt=0.):
    """
        Standard deviation of skew exponential power distribution with given skewness and kurtosis,
        location zero and unit scale.


        Definition
        ----------
        def sep01_std(skew=1., kurt=0.):


        Optional Input
        --------------
        skew       skewness parameter
        kurt       kurtosis parameter
        

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

    mu  = sep01_mean(skew, kurt)
    if skew != 1.:
        M1 = mu / (skew-1./skew)
    else:
        M1 = 0.
    M2  = 1.0
    var = (M2-M1**2)*(skew**2+1./skew**2) + 2.0*M1**2 - M2
    if var > 0.0:
        sig  = np.sqrt(var)
    else:
        sig  = 0.0

    return sig

    
# --------------------------------------------------------------------


def ssep(x, loc=0., sca=1., skew=1., kurt=0., sig=None):
    """
        The standardise skew exponential power distribution with given location, scale, skewness, and kurtosis.


        Definition
        ----------
        def ssep(x, loc=0., sca=1., skew=1., kurt=0., sig=None):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        kurt       kurtosis parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Standardised skew exponential power pdf at x
        

        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.
        

        Examples
        --------
        >>> print(np.allclose(ssep(1., 2., 2., 0.5, 2.), ssep((1.-2.)/2., skew=0.5, kurt=2.)/2.))
        True

        >>> print(np.allclose(ssep(1.3, 0., 1., 1., 0.), normal(1.3, 0., 1.)))
        True

        >>> print(np.allclose(ssep(1.3, 0., 1., 1., 1.), laplace(1.3, 0., 1./np.sqrt(2.))))
        True

        >>> print(np.allclose(ssep(1.3, 0., np.sqrt(2.), 1., 1.), laplace(1.3, 0., 1.)))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is None:
        return ssep01((x-loc)/sca, skew, kurt)/sca
    else:
        sca = sig
        return ssep01((x-loc)/sca, skew, kurt)/sca


def ssep01(x, skew=1., kurt=0.):
    """
        The standardised skew exponential power distribution with given skewness and kurtosis, location zero and unit scale.


        Definition
        ----------
        def ssep01(x, skew=1., kurt=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        skew       skewness parameter
        kurt       kurtosis parameter
        

        Output
        ------
        Standardised skew exponential power pdf with loc=0, sca=1 at x
        

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

    mu  = sep01_mean(skew, kurt)
    sig = sep01_std(skew, kurt)
    z   = mu + sig*x
    
    return sig * sep01(z, skew, kurt)


# --------------------------------------------------------------------


def ssstudentt(x, nu, loc=0., sca=1., skew=1., sig=None):
    """
        The standardised skewed Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def ssstudentt(x, nu, loc=0., sca=1., skew=1., sig=None):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        sig        standard deviation, overwrites scale


        Output
        ------
        Skewed Student t pdf at x
        

        Examples
        --------
        >>> print(np.allclose(ssstudentt(1., 2., 2., 0.5, 2.), ssstudentt((1.-2.)/0.5, 2., skew=2.)/0.5))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is None:
        return ssstudentt01((x-loc)/sca, nu, skew)/sca
    else:
        sca = sig
        return ssstudentt01((x-loc)/sca, nu, skew)/sca


def ssstudentt01(x, nu, skew=1.):
    """
        The standardised skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def ssstudentt01(x, nu, skew=1.):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        skew       skewness parameter
        

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

    mu  = sstudentt01_mean(nu, skew)
    sig = sstudentt01_std(nu, skew)
    z   = mu + sig*x

    return sig * sstudentt01(z, nu, skew)

    
# --------------------------------------------------------------------


def sstudentt(x, nu, loc=0., sca=1., skew=1.):
    """
        The skewed Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def ssstudentt(x, nu, loc=0., sca=1., skew=1., sig=None):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter


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

    return sstudentt01((x-loc)/sca, nu, skew)/sca


def sstudentt01(x, nu, skew=1.):
    """
        The skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.
        
        If skew is not 1 then mean is not zero and standard deviation is not 1.


        Definition
        ----------
        def sstudentt01(x, nu, skew=1.):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        skew       skewness parameter
        

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

    alpha = np.where(x<0.0, skew, 1./skew)
    if not np.iterable(x): alpha = alpha[0]        

    return 2.0/(skew+1./skew) * studentt01(alpha*x, nu)


def sstudentt_mean(nu, loc=0., sca=1., skew=1.):
    """
        The mean of the skewed Student t distribution with given degrees of freedom,
        location, scale, and skewness.


        Definition
        ----------
        def sstudentt_mean(nu, loc=0., sca=1., skew=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        skew       skewness parameter
        

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

    return sstudentt01_mean(nu, skew) * sca + loc


def sstudentt01_mean(nu, skew=1.):
    """
        The mean of the skewed Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def sstudentt01_mean(nu, skew=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        skew       skewness parameter
        

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

    mu  = 2.0 * (skew-1./skew) * ss.gamma(0.5*(nu+1.0)) / ss.gamma(0.5*nu)
    mu *= (nu-2.0)/(nu-1.0) * nu/(nu-2.) / np.sqrt(np.pi*nu)

    return mu


def sstudentt_std(nu, sca=1., skew=1.):
    """
        The standard deviation of the skewed Student t distribution with given degrees of freedom,
        scale, and skewness.


        Definition
        ----------
        def sstudentt_std(nu, sca=1., skew=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        sca        scale
        skew       skewness parameter
        

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

    return sstudentt01_std(nu, skew) * sca


def sstudentt01_std(nu, skew=1.):
    """
        The standard deviation of the skewed Student t distribution with given degrees of freedom and skewness,
        and unit scale.


        Definition
        ----------
        def sstudentt01_std(nu, skew=1.):


        Input
        -----
        nu         degrees of freedom


        Optional Input
        --------------
        skew       skewness parameter
        

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

    mu = sstudentt01_mean(nu, skew)
    if skew != 1.:
        M1 = mu / (skew-1./skew)
    else:
        M1 = 0.
    M2  = nu/(nu-2.)
    var = (M2-M1**2)*(skew**2+1./skew**2)+2.*M1**2-M2
    if var > 0.0:
        return np.sqrt(var)
    else:
        return 0.0

    
# --------------------------------------------------------------------


def studentt(x, nu, loc=0., sca=1., sig=None):
    """
        The Student t distribution with given degrees of freedom, location, scale, and skewness.


        Definition
        ----------
        def studentt(x, nu, loc=0., sca=1., sig=None):


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
        >>> print(np.allclose(ssstudentt(1., 2., 2., 0.5, 2.), ssstudentt((1.-2.)/0.5, 2., skew=2.)/0.5))
        True


        History
        -------
        Written,  MC, May 2016
    """

    if sig is None:
        return studentt01((x-loc)/sca, nu)/sca
    else:
        sca = np.sqrt((nu-2.)/nu)*sig
        return studentt01((x-loc)/sca, nu)/sca


def studentt01(x, nu):
    """
        The Student t distribution with given degrees of freedom and skewness,
        location zero and unit scale.


        Definition
        ----------
        def studentt01(x, nu):


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

    c      = ss.gamma(0.5*(nu+1.0)) / (ss.gamma(0.5*nu) * np.sqrt(np.pi*nu))
    return c * (1.0 + x**2/nu)**(-0.5*(nu+1.0))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # import scipy.special as spec
    # import scipy.stats as sp
    # import ufz

    # nn = 10000000
    # xx1 = np.arange(nn)/float(nn) * 30.*100. - 15.*100.

    # loc=1.
    # sca=2.
    # skew=2.
    # kurt=0.5
    # nu=3.1
    # theta=1.5

    # dx = xx1[2]-xx1[1]

    # __all__ = ['ep', 'exponential', 'laplace', 'normal', 'sep', 'ssstudentt', 'studentt']
    # __all__ = ['ssep', 'ssstudentt']

    # for dd in __all__:
    #     xx = xx1
    #     if dd == 'ep':
    #         pdf = ufz.distributions.ep(xx, loc, sca, kurt)
    #     elif dd == 'exponential':
    #         xx = xx1[xx1>0.]
    #         pdf = ufz.distributions.exponential(xx, loc, sca, theta)
    #     elif dd == 'laplace':
    #         pdf = ufz.distributions.laplace(xx, loc, sig=sca)
    #     elif dd == 'normal':
    #         pdf = ufz.distributions.normal(xx, loc, sca)
    #     elif dd == 'ssep':
    #         pdf = ufz.distributions.ssep(xx, loc, sca, skew, kurt)
    #     elif dd == 'ssstudentt':
    #         pdf = ufz.distributions.ssstudentt(xx, nu, loc, sca, skew)
    #     elif dd == 'studentt':
    #         pdf = ufz.distributions.studentt(xx, nu, loc, sca)

    #     dmean = np.sum(xx * pdf * dx)
    #     dvar  = np.sum((xx-dmean)**2 * pdf * dx)
    #     dstd  = np.sqrt(dvar)
    #     dskew = np.sum(((xx-dmean)/dstd)**3 * pdf * dx)
    #     dkurt = np.sum(((xx-dmean)/dstd)**4 * pdf * dx)
    #     # print(dd, ':', dmean, dstd)
    #     print(dd, ':', dmean, dvar, dskew, dkurt)

    #     # xx   = xx1[xx1>0.]
    #     # pdf  = ufz.distributions.studentt(xx, nu, loc, sca)
    #     # M1   = 2.* np.sum(xx * pdf * dx)
    #     # M2   = nu/(nu-2.) * sca**2
    #     # dvarfs = (M2-M1**2)*(skew**2+1./skew**2)+2.*M1**2-M2
    #     # # print('Fernandez var:', dvarfs)

    #     xi = skew
    #     mu_shmc  = 2.0 * (xi-1./xi) * spec.gamma(0.5*(nu+1.0)) / spec.gamma(0.5*nu)
    #     mu_shmc *= (nu-2.0)/(nu-1.0)
    #     mu_shmc *= nu/(nu-2.) / np.sqrt(np.pi*nu)
    #     mu_shmc *= sca
        
    #     # var_shmc = -mu_shmc**2 + (xi**3+1./xi**3)/(xi+1./xi)
    #     # # print('mu_shmc:', mu_shmc, var_shmc)

    #     if xi != 1.:
    #         M1 = mu_shmc / (xi-1./xi)
    #     else:
    #         M1 = 0.
    #     M2   = nu/(nu-2.) * sca**2
    #     dvar_mc = (M2-M1**2)*(skew**2+1./skew**2)+2.*M1**2-M2
    #     print('dvar_mc:', mu_shmc, dvar_mc)
