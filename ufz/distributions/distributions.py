#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['gauss', 'laplace', 'normal', 'sep', 'sstudent']

def gauss(*args, **kwargs):
    """
        Wrapper for normal

        def normal(x, mu=0., sig=1.):
    """
    return normal(*args, **kwargs)


def laplace(x, mu=0., sig=1.):
    """
        Laplace probability density function (pdf).


        Definition
        ----------
        def laplace(x, mu=0., sig=1.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        mu         mean
        sig        standard deviation
        

        Output
        ------
        Laplace pdf at x


        Restrictions
        ------------
        None
        

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


        History
        -------
        Written,  MC, May 2016
    """

    return laplace01((x-mu)/sig)/sig


def laplace01(x):
    """
        Laplace probability density function (pdf) with with zero mean and unit standard deviation.


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
        Laplace pdf with mean=0 and stddev=1 at x


        Restrictions
        ------------
        None
        

        Examples
        --------
        >>> print(str(laplace01(0.)))
        0.5

        >>> print(str(laplace01(1.) - 0.5/np.e))
        0.0


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


        History
        -------
        Written,  MC, May 2016
    """

    return 0.5 * np.exp(-np.abs(x))


def normal(x, mu=0., sig=1.):
    """
        Normal (Gauss) probability density function (pdf).


        Definition
        ----------
        def normal(x, mu=0., sig=1.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        mu         mean
        sig        standard deviation
        

        Output
        ------
        Normal (Gauss) pdf at x


        Restrictions
        ------------
        None
        

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


        History
        -------
        Written,  MC, May 2016
    """

    return normal01((x-mu)/sig)/sig


def normal01(x):
    """
        Normal (Gauss) probability density function (pdf) with zero mean and unit standard deviation.


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
        Normal pdf with mean=0 and stddev=1 at x


        Restrictions
        ------------
        None
        

        Examples
        --------
        >>> print(np.allclose(normal01(0.), 1./np.sqrt(2.*np.pi)))
        True

        >>> print(np.allclose(normal01(1.), 1./np.sqrt(2.*np.pi*np.e)))
        True


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


        History
        -------
        Written,  MC, May 2016
    """

    return 1./np.sqrt(2.*np.pi) * np.exp(-0.5*x*x)


def sep(x, mu=0., sig=1., skew=1., kurt=0.):
    """
        The skew exponential power distribution with given mean, standard deviation, skewness, and kurtosis.


        Definition
        ----------
        def sep(x, mu=0., sig=1., skew=1., kurt=0.):


        Input
        -----
        x          array_like quantiles


        Optional Input
        --------------
        mu         mean
        sig        standard deviation
        skew       skewness parameter
        kurt       kurtosis parameter


        Output
        ------
        Skew exponential power pdf at x


        Restrictions
        ------------
        None
        

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


        History
        -------
        Written,  MC, May 2016
    """

    return sep01((x-mu)/sig, skew, kurt)/sig


def sep01(x, skew=1., kurt=0.):
    """
        The skew exponential power distribution with given skewness and kurtosis, zero mean and unit standard deviation.


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
        Skew exponential power pdf with mean=0, stddev=1 at x


        Restrictions
        ------------
        None
        

        Examples
        --------


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


def sstudentt(x, nu, mu=0., sig=1., skew=1.):
    """
        The skewed Student t distribution with given degrees of freedom, mean, standard deviation, and skewness.


        Definition
        ----------
        def sstudentt(x, nu, mu=0., sig=1., skew=1.):


        Input
        -----
        x          array_like quantiles
        nu         degrees of freedom


        Optional Input
        --------------
        mu         mean
        sig        standard deviation
        skew       skewness parameter


        Output
        ------
        Skewed Student t pdf at x


        Restrictions
        ------------
        None
        

        Examples
        --------
        >>> print(np.allclose(sstudentt(1., 2., 2., 0.5, 2.), sstudentt((1.-2.)/0.5, 2., skew=2.)/0.5))
        True


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


        History
        -------
        Written,  MC, May 2016
    """

    return sstudentt01((x-mu)/sig, nu, skew)/sig


def sstudentt01(x, nu, skew=1.):
    """
        The skewed Student t distribution with given degrees of freedom and skewness,
        zero mean and unit standard deviation.


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
        Skewed Student t pdf with mean=0 and stddev=1 at x


        Restrictions
        ------------
        None
        

        Examples
        --------


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


        History
        -------
        Written,  MC, May 2016
    """
    import scipy.special as ss

    xi = skew

    xi1    = 1.0/xi
    c      = ss.gamma(0.5*(nu+1.0)) / (ss.gamma(0.5*nu) * np.sqrt(np.pi*nu))
    mu_xi  = 2.0 * (xi*xi-xi1*xi)/(xi+xi1) * c * (nu-2.0)/(nu-1.0)
    
    sig_xi = -mu_xi*mu_xi + (xi**3+xi**3)/(xi+xi1)
    if sig_xi > 0.0:
        sig_xi = np.sqrt(sig_xi)
    else:
        sig_xi = 0.0
        
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

    return 2.0 * sig_xi / (xi+xi1) * c * (1.0 + a_xi*a_xi/nu)**(-0.5*(nu+1.0))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
