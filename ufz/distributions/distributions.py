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

__all__ = ['ep', 'exponential', 'gauss', 'laplace', 'normal', 'sep', 'ssstudentt', 'sstudentt', 'studentt']

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
        return om_beta # ToDo - vector
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


def sep(x, loc=0., sca=1., skew=1., kurt=0., sig=None):
    """
        The skew exponential power distribution with given location, scale, skewness, and kurtosis.


        Definition
        ----------
        def sep(x, loc=0., sca=1., skew=1., kurt=0., sig=None):


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

    if sig is None:
        return sep01((x-loc)/sca, skew, kurt)/sca
    else:
        sca = sig
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
        

        Examples
        --------
        None


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    xi   = skew
    beta = kurt
    
    if beta != -1.0:
        b1 = 0.5*(1.0 + beta)
        b2 =      1.0 + beta
        b3 = 1.5*(1.0 + beta)
        g1 = ss.gamma(b1)
        g2 = ss.gamma(b2)
        g3 = ss.gamma(b3)
        # -> sqrt(3/4)
        M1 = g2 / np.sqrt(g3*g1)
        # -> 0
        c_beta = (g3/g1)**(1.0/(1.0+beta))
        # -> sqrt(1/12)
        om_beta = np.sqrt(g3)/((1.0+beta)*np.sqrt(g1**3))
    else:
        M1      = np.sqrt(0.75)
        c_beta  = 0.0
        om_beta = np.sqrt(1.0/12.0)
    M2  = 1.0
    xi1 = 1.0/xi

    # coefficients
    mu_xi   = M1*(xi-xi1)
    sig_xi  = (M2-M1*M1)*(xi*xi+xi1*xi1) + 2.0*M1*M1 - M2
    if sig_xi > 0.0:
        sig_xi  = np.sqrt(sig_xi)
    else:
        sig_xi  = 0.0
        
    a_xi = mu_xi+sig_xi*x
    if np.iterable(x):
        for i in range(len(x)):
            if a_xi[i] < 0.0:
                a_xi[i] = a_xi[i] * xi
            else:
                a_xi[i] = a_xi[i] * xi1
    else:
        if a_xi < 0.0:
            a_xi = a_xi * xi
        else:
            a_xi = a_xi * xi1
            
    # pdf
    if np.abs(beta+1.0) < 0.003: # 2/(1-0.997) ~ 666
        return 2.0*sig_xi/(xi+xi1) * om_beta
    else:
        return 2.0*sig_xi/(xi+xi1) * om_beta * np.exp(-c_beta*np.abs(a_xi)**(2.0/(1.0+beta)))

    
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


        History
        -------
        Written,  MC, May 2016
    """
    # import scipy.special as ss

    # xi = skew

    # mu_xi  = 2.0 * (xi-1./xi) * ss.gamma(0.5*(nu+1.0)) / ss.gamma(0.5*nu)
    # mu_xi *= (nu-2.0)/(nu-1.0)
    # mu_xi *= nu/(nu-2.) / np.sqrt(np.pi*nu)
    # # mu_xi *= sca
    
    # if xi != 1.:
    #     M1 = mu_xi / (xi-1./xi)
    # else:
    #     M1 = 0.
    # M2  = nu/(nu-2.)
    # # M2 *= sca**2
    # var = (M2-M1**2)*(skew**2+1./skew**2)+2.*M1**2-M2
    # sig_xi = np.sqrt(var)

    # z = mu_xi + sig_xi*x
    # pdf = sig_xi * sstudentt01(z, nu, xi)

    # return pdf

    mu  = sstudentt01_mean(nu, skew)
    sig = sstudentt01_std(nu, skew)
    z   = mu + sig*x
    return sig * sstudentt01(z, nu, skew)

    
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
