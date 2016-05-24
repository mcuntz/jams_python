#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.special as sp
import scipy.stats as ss

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

    Copyright 2016 Juliane Mai, Dmitri Kavetski
"""

__all__ = ['sample_ep01',      # sample from exponential power function EP(0,1,beta)
           'sample_ep',        # sample from (general) exponential power function EP(loc,sca,beta)
           'sample_sep_fs',    # sample from skew exponential power distribution obtained by using the approach of Fernandez and Steel
           'sample_ssep',      # sample from standardized skew exponential power distribution (mean=zero, std. dev.=1)
           'sample_sep'        # sample from the (general) skew exponential power distribution with
                               # given location, scale, and parameters controlling skewness and kurtosis.
           ]

def sample_ep01(nn, beta=0.):
    """
        Samples from the exponential power function EP(0,1,beta) are drawn.

        This distribution is symmetric.

        The kurtosis of this sample can be calculated by
             bb   = 2./(beta+1.)
             kurt = (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3.
        

        Definition
        ----------
        def sample_ep01(nn, beta=0.):


        Input
        -----
        nn       -> number of samples

        
        Optional Input
        --------------
        beta       parameter which controls the kurtosis


        Output
        ------
        Samples drawn from exponential power function


        Restrictions
        ------------
        Not known


        Examples
        --------
        >>> np.random.seed(seed=12345)
        >>> 

        
        Literature
        --------
        Schoups, G., & Vrugt, J. A. (2010).
        A formal likelihood function for parameter and predictive inference of hydrologic
        models with correlated, heteroscedastic, and non-Gaussian errors.
        Water Resources Research, 46(10).
        --> Steps (1)-(3) described on page 5

        
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

        Copyright 2010-2016 Juliane Mai


        History
        -------
        Written  JM, May 2016
        Modified 
    """

    bb   = (1.+beta)/2.
        
    # (1) Generate a sample gg from gamma distribution with shape
    #     parameter (1+beta)/2 and scale parameter 1.
    gg = np.random.gamma(bb, 1.0, nn)

    # (2) Generate a random sign ss (+1 or -1) with equal probability.
    ss = np.where(np.random.rand(nn) < 0.5, -1., 1.)

    # (3) Compute EP which is a sample from the exponential power distribution EP(0,1,beta)
    EP01 = ss * np.abs(gg)**bb * np.sqrt(sp.gamma(bb)) / np.sqrt(sp.gamma(3.*bb))

    return EP01

def sample_ep(nn, loc=0., sca=1., beta=0.):
    """
        Samples from the (general) exponential power function EP(loc,sca,beta) are drawn.

        The mean is loc, standard deviation is sca.

        This distribution is symmetric.

        The kurtosis of this sample can be calculated by
             bb   = 2./(beta+1.)
             kurt = (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3.

        Definition
        ----------
        def sample_ep(nn, loc=0., sca=1., beta=0.):


        Input
        -----
        nn       -> number of samples

        
        Optional Input
        --------------
        loc        location
        sca        scale
        beta       parameter which controls kurtosis


        Output
        ------
        Samples drawn from exponential power function


        Restrictions
        ------------
        Not known


        Examples
        --------
        >>> np.random.seed(seed=12345)
        >>> 

        
        Literature
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

        Copyright 2010-2016 Juliane Mai


        History
        -------
        Written  JM, May 2016
        Modified 
    """

    EP01 = sample_ep01(nn, beta=beta)
    EP   = loc + sca * EP01

    return EP

def sample_sep_fs(nn, xi=1., beta=0.):
    """
        Samples from the skew exponential power distribution obtained
        by using the approach of Fernandez and Steel.
        Parameters which control skewness and kurtosis needs to be given.
        If xi = 1. then mean=0. and stddev=1.


        Definition
        ----------
        def sample_sep_fs(nn, xi=1., beta=0.):


        Input
        -----
        nn       -> number of samples


        Optional Input
        --------------
        xi         parameter which controls skewness
        beta       parameter which controls the kurtosis


        Output
        ------
        Samples from the skew exponential power distribution


        Examples
        --------
        None


        Literature
        --------
        Schoups, G., & Vrugt, J. A. (2010).
        A formal likelihood function for parameter and predictive inference of hydrologic
        models with correlated, heteroscedastic, and non-Gaussian errors.
        Water Resources Research, 46(10).
        --> Steps (4)-(5) described on page 5


        History
        -------
        Written,  JM, May 2016
    """

    # (1) - (3) 
    EP = sample_ep01(nn, beta=beta)

    # (4) Generate a random sign w_t (+1 or -1) with probabilities 1-xi/(xi+1/xi) and xi/(xi+1/xi)
    #     wt = np.where(np.random.rand(nn) > 1.-xi/(xi+1./xi), 1, -1)
    # (5) Compute SEP_t which is a sample from the skew exponential power distribution
    #     SEP(mu_xi, sigma_xi, xi, beta) with mu_xi and sigma_xi given by (A5) and (A6)
    #     SEPt = wt * np.abs(EPt) * xi**wt
    #
    # TODO: put this 3 lines in a function
    #       sample_skewed = skew_fernandez_steel(sample_symmetric, xi)
    #       --> can be re-used for the student-t<
    rr     = np.random.rand(nn)
    absEP  = np.abs(EP)
    SEP_fs = np.where(rr > 1.-xi/(xi+1./xi), absEP * xi, -absEP / xi)

    return SEP_fs

def sample_ssep(x, xi=1., beta=0.):
    """
        Samples from the standardized skew exponential power distribution, obtained
        by standardizing the distribution resulting from the approach of Fernandez and Steel.
        Mean is zero and standard deviation is one.


        Definition
        ----------
        def sample_ssep(nn, xi=1., beta=0.):


        Input
        -----
        nn       -> number of samples


        Optional Input
        --------------
        xi         parameter which controls the skewness
        beta       parameter which controls the kurtosis


        Output
        ------
        Samples from the standardized skew exponential power distribution


        Examples
        --------
        None


        Literature
        --------
        Schoups, G., & Vrugt, J. A. (2010).
        A formal likelihood function for parameter and predictive inference of hydrologic
        models with correlated, heteroscedastic, and non-Gaussian errors.
        Water Resources Research, 46(10).
        --> Steps (6) described on page 5


        History
        -------
        Written,  JM, May 2016
    """

    SEP_fs = sample_sep_fs(nn, xi=xi, beta=beta)
    
    # (6) Standardize SEP_fs
    mean_sep_fs  = ufz.distributions.sep_mean(skew=xi, kurt=beta)
    std_sep_fs   = ufz.distributions.sep_std(skew=xi, kurt=beta)
    sSEP         = (SEP_fs - mean_sep_fs) / std_sep_fs   # standardized SEP (=Schoups and Vrugt's a_t)

    return sSEP

def sample_sep(nn, loc=0., sca=1., xi=1., beta=0.):
    """
        Samples from the (general) skew exponential power distribution with
        given location, scale, and parameters controlling skewness and kurtosis.
        This distribution is obtained by shifting and scaling the standardized SEP distribution.

        The location loc and the scale sca are the mean and standard deviation
        of the distribution respectively.

        The xi parameter needs to be positive.
             1 = symmetric
            >1 = right skewed
            <1 = left skewed

        The beta parameter needs to be from [-1,1]
             1 = fat tail
             0 = medium tails
            -1 = thin tail
            
        If xi = 1 and beta = 0, the SEP is the Gaussian distribution.
        


        Definition
        ----------
        def sample_sep(nn, loc=0., sca=1., xi=1., beta=0.):


        Input
        -----
        nn       -> number of samples


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         parameter which controls the skewness
        beta       parameter which controls the kurtosis


        Output
        ------
        Samples from the (general) skew exponential power


        Examples
        --------
        None


        Literature
        --------
        None 


        History
        -------
        Written,  JM, May 2016
    """
    
    sSEP = sample_ssep(nn, xi=xi, beta=beta)
    SEP    = loc + sca * sSEP 

    return SEP

# # ------------------------------------------------------------
# # Do some tests for the sampling of distributions
# # ------------------------------------------------------------

# import ufz
# import matplotlib.pylab as plt
# plt.close("all")

# nn    = 1000000   # number of random samples drawn
# loc   = 5.0       # location parameter = mean
# sca   = 0.8       # scale parameer = standard deviation
# xi    = 1.0       # skew parameter (0,inf) 1=symmetric, >1=right skewed, <1=left skewed
# beta  = -0.5       # ex. kurtosis parameter (-1,1] 1=fat tail, -1=thin tail

# print('loc            = ', loc)
# print('sca            = ', sca)
# print('xi             = ', xi)
# print('beta           = ', beta)
# print('')

# # --------------
# # EP tests
# # --------------
# ep_samples = sample_ep(nn, loc=loc, sca=sca, beta=beta)
# plt.hist(ep_samples,bins=100,normed=True)
# print('EP:  mean_samp      = ', np.mean(ep_samples))
# print('EP:  std_samp       = ', np.std(ep_samples))
# print('EP:  skew_samp      = ', ss.skew(ep_samples))
# print('EP:  kurt_samp      = ', ss.kurtosis(ep_samples, fisher=True))
# print('')

# # Compare with pdf of EP
# dx = 0.01
# ep_xx = np.arange(-10,10,dx)
# ep_yy = ufz.distributions.ep(ep_xx, loc=loc, sca=sca, kurt=beta)

# mean_num = np.sum(ep_xx*ep_yy*dx)
# var_num  = np.sum((ep_xx-mean_num)**2*ep_yy*dx)
# std_num  = np.sqrt(var_num)
# skew_num  = np.sum(((ep_xx-mean_num)/std_num)**3*ep_yy*dx)
# kurt_num  = np.sum(((ep_xx-mean_num)/std_num)**4*ep_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
# print('EP:  mean_num       = ', mean_num )
# print('EP:  std_num        = ', std_num  )
# print('EP:  skew_num       = ', skew_num )
# print('EP:  kurt_num       = ', kurt_num )
# print('')

# bb = 2./(beta+1.)
# print('EP:  calc kurt      = ', (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3. )

# plt.plot(ep_xx, ep_yy, 'r-', lw=5, alpha=0.6)
# plt.show()

# # --------------
# # SEP tests
# # --------------
# sep_samples = sample_sep(nn, loc=loc, sca=sca, xi=xi, beta=beta)
# plt.hist(sep_samples,bins=100,normed=True)
# print('SEP: mean_samp      = ', np.mean(sep_samples))
# print('SEP: std_samp       = ', np.std(sep_samples))
# print('SEP: skew_samp      = ', ss.skew(sep_samples))
# print('SEP: kurt_samp      = ', ss.kurtosis(sep_samples, fisher=True))
# print('')

# # Compare with pdf of SEP
# dx = 0.01
# sep_xx = np.arange(-10,10,dx)
# sep_yy = ufz.distributions.ssep(sep_xx, loc=loc, sca=sca, skew=xi, kurt=beta)

# mean_num = np.sum(sep_xx*sep_yy*dx)
# var_num  = np.sum((sep_xx-mean_num)**2*sep_yy*dx)
# std_num  = np.sqrt(var_num)
# skew_num  = np.sum(((sep_xx-mean_num)/std_num)**3*sep_yy*dx)
# kurt_num  = np.sum(((sep_xx-mean_num)/std_num)**4*sep_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
# print('SEP: mean_num       = ', mean_num )
# print('SEP: std_num        = ', np.sqrt(var_num) )
# print('SEP: skew_num       = ', skew_num )
# print('SEP: kurt_num       = ', kurt_num )
# print('')

# plt.plot(sep_xx, sep_yy, 'r-', lw=5, alpha=0.6)
# plt.show()
