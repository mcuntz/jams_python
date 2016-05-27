#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.special as sp
import scipy.stats as ss
from .distributions import sep_mean, sep_std
from .distributions import sstudentt_mean, sstudentt_std

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
           'sample_sep',       # sample from the (general) skew exponential power distribution with
                               # given location, scale, and parameters controlling skewness and kurtosis.
           #
           'sample_studentt01',      # sample from Student-t distribution studentt(loc=0,sca=1)
           'sample_studentt',        # sample from (general) Student-t distribution studentt(loc,sca)
           'sample_sstudentt_fs',    # sample from skew Student-t distribution obtained by using the approach of Fernandez and Steel
           'sample_ssstudentt',      # sample from standardized skew student-t distribution (mean=zero, std. dev.=1)
           'sample_sstudentt'        # sample from the (general) skew Student-t distribution with
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

def sample_studentt01(nn, nu):
    """
        Samples from the Student-t distribution studentt(df=nu, loc=0, sca=1) are drawn.

        This distribution is symmetric.
        

        Definition
        ----------
        def sample_studentt01(nn):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom

        
        Optional Input
        --------------
        None


        Output
        ------
        Samples drawn from Student-t distribution


        Restrictions
        ------------
        Not known


        Examples
        --------
        None

        
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

    studentt01 = ss.t.rvs(nu, loc=0.0, scale=1.0, size=nn)

    return studentt01

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

def sample_studentt(nn, nu, loc=0., sca=1.):
    """
        Samples from the (general) Student-t distribution studentt(nu,loc,sca) are drawn.

        The mean is loc, standard deviation is sca.

        This distribution is symmetric.

        Definition
        ----------
        def sample_studentt(nn, nu, loc=0., sca=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom

        
        Optional Input
        --------------
        loc        location
        sca        scale


        Output
        ------
        Samples drawn from (general) Student-t distribution


        Restrictions
        ------------
        Not known


        Examples
        --------
        None 

        
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

    studentt01 = sample_studentt01(nn, nu)
    studentt   = studentt01 * sca + loc
    # studentt = ss.t.rvs(nu, loc=loc, scale=sca, size=nn)

    return studentt

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

        Fernandez C & Steel M (1998).
        On Bayesian modeling of fat tails and skewness.
        Journal of the American Statistical Association 93, 359-371.


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
    SEP_fs = skew_fernandez_steel(EP, xi)

    return SEP_fs

def sample_sstudentt_fs(nn, nu, xi=1.):
    """
        Samples from the skew Student-t distribution obtained
        by using the approach of Fernandez and Steel.
        Parameters which control skewness needs to be given.


        Definition
        ----------
        def sample_sstudentt_fs(nn, nu, skew=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        xi       parameter which controls skewness


        Output
        ------
        Samples from the skew Student-t distribution


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
    studentt = sample_studentt01(nn, nu)

    # (4) Generate a random sign w_t (+1 or -1) with probabilities 1-xi/(xi+1/xi) and xi/(xi+1/xi)
    #     wt = np.where(np.random.rand(nn) > 1.-xi/(xi+1./xi), 1, -1)
    # (5) Compute SEP_t which is a sample from the skew exponential power distribution
    #     SEP(mu_xi, sigma_xi, xi, beta) with mu_xi and sigma_xi given by (A5) and (A6)
    #     SEPt = wt * np.abs(EPt) * xi**wt
    sstudentt_fs = skew_fernandez_steel(studentt, xi)

    return sstudentt_fs

def sample_ssep(nn, xi=1., beta=0.):
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
    mean_sep_fs  = sep_mean(skew=xi, kurt=beta)
    std_sep_fs   = sep_std(skew=xi, kurt=beta)
    sSEP         = (SEP_fs - mean_sep_fs) / std_sep_fs   # standardized SEP (=Schoups and Vrugt's a_t)

    return sSEP

def sample_ssstudentt(nn, nu, xi=1.):
    """
        Samples from the standardized skew Student-t distribution, obtained
        by standardizing the distribution resulting from the approach of Fernandez and Steel.
        Mean is zero and standard deviation is one.


        Definition
        ----------
        def sample_ssstudentt(nn, nu, xi=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        xi         parameter which controls the skewness


        Output
        ------
        Samples from the standardized skew Student-t distribution


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

    sstudentt_fs = sample_sstudentt_fs(nn, nu, xi=xi)
    
    # (6) Standardize sstudentt_fs
    mean_sstudentt_fs  = sstudentt_mean(nu, skew=xi)
    std_sstudentt_fs   = sstudentt_std(nu, skew=xi)
    ssstudentt         = (sstudentt_fs - mean_sstudentt_fs) / std_sstudentt_fs   # standardized skewed Student-t

    return ssstudentt

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

def sample_sstudentt(nn, nu, loc=0., sca=1., xi=1., sigma=None):
    """
        Samples from the (general) skew Student-t distribution with
        given location, scale, and parameter controlling skewness.
        This distribution is obtained by shifting and scaling the standardized skewed Student-t distribution.

        The location loc and the scale sca are the mean and standard deviation
        of the distribution respectively.

        The xi parameter needs to be positive.
             1 = symmetric
            >1 = right skewed
            <1 = left skewed        


        Definition
        ----------
        def sample_sstudentt(nn, nu, loc=0., sca=1., xi=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        sigma      standard deviation
        xi         parameter which controls the skewness


        Output
        ------
        Samples from the (general) skew Student-t distribution


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

    ssstudentt = sample_ssstudentt(nn, nu, xi=xi)
    
    if (sigma is None):
        sigma = sca*np.sqrt(nu/(nu-2.))
        
    sstudentt  = loc + sigma * ssstudentt
    

    return sstudentt

# ------------------------------------------------------------
# Private routines
# ------------------------------------------------------------
def skew_fernandez_steel(sample_symmetric, xi):
    """
        Skews any given symmetric sample with a skewness factor xi
        following the approch of Fernandez & Steel.     
        

        Definition
        ----------
        def skew_fernandez_steel(sample_symmetric, xi):


        Input
        -----
        sample_symmetric       -> sample drawn from a symmetric distribution
        xi                     -> skewness factor


        Optional Input
        --------------
        None


        Output
        ------
        Samples from the skewed version of symmetric distribution used to generate sample


        Examples
        --------
        None


        Literature
        --------
        Fernandez C & Steel M (1998).
        On Bayesian modeling of fat tails and skewness.
        Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  JM, May 2016
    """
        
    nn       = np.shape(sample_symmetric)[0]
    rr       = np.random.rand(nn)
    absSymm  = np.abs(sample_symmetric)
    sample_fs   = np.where(rr > 1.-xi/(xi+1./xi), absSymm * xi, -absSymm / xi)

    return sample_fs

