#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.special as sp
import scipy.stats as ss
from ufz.distributions import sep_mean, sep_std
from ufz.distributions import sstudentt_mean, sstudentt_std

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
    std_sstudentt_fs   = sstudentt_std(nu, skew=xi) * np.sqrt((nu-2.)/nu)
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

def sample_sstudentt(nn, nu, loc=0., sca=1., xi=1., sig=None):
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
        sig        standard deviation
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
    
    if sig is not None:
        sca = sig * np.sqrt((nu-2.)/nu)
        
    sstudentt  = loc + sca * ssstudentt
    

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


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # ------------------------------------------------------------
    # # Do some tests for the sampling of distributions
    # # ------------------------------------------------------------

    # import ufz
    # import matplotlib.pylab as plt
    # plt.close("all")

    # nn    = 10000000   # number of random samples drawn
    # loc   = 0.0       # location parameter = mean
    # sca   = 2.0       # scale parameer = standard deviation
    # xi    = 1.5       # skew parameter (0,inf) 1=symmetric, >1=right skewed, <1=left skewed
    # beta  = -0.5       # ex. kurtosis parameter (-1,1] 1=fat tail, -1=thin tail
    # nu    = 6        # degrees of freedom for student-t

    # test_EP         = False
    # test_SEP        = False
    # test_Studentt   = False
    # test_ssStudentt = False
    # test_sStudentt  = True

    # print('loc  = ', loc)
    # print('sca  = ', sca)
    # print('xi   = ', xi)
    # print('beta = ', beta)
    # print('nu   = ', nu)
    # print('')

    # if (test_EP):
    #     # --------------
    #     # EP tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     ep_samples = ufz.distributions.sample_ep(nn, loc=loc, sca=sca, beta=beta)
    #     plt.hist(ep_samples,bins=100,normed=True)
    #     print('EP:  mean_samp      = ', np.mean(ep_samples))
    #     print('EP:  std_samp       = ', np.std(ep_samples))
    #     print('EP:  skew_samp      = ', ss.skew(ep_samples))
    #     print('EP:  kurt_samp      = ', ss.kurtosis(ep_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of EP
    #     dx = 0.01
    #     ep_xx = np.arange(-10,10,dx)
    #     ep_yy = ufz.distributions.ep(ep_xx, loc=loc, sca=sca, kurt=beta)

    #     mean_num = np.sum(ep_xx*ep_yy*dx)
    #     var_num  = np.sum((ep_xx-mean_num)**2*ep_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((ep_xx-mean_num)/std_num)**3*ep_yy*dx)
    #     kurt_num  = np.sum(((ep_xx-mean_num)/std_num)**4*ep_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('EP:  mean_num       = ', mean_num )
    #     print('EP:  std_num        = ', std_num  )
    #     print('EP:  skew_num       = ', skew_num )
    #     print('EP:  kurt_num       = ', kurt_num )
    #     print('')

    #     bb = 2./(beta+1.)
    #     print('EP:  calc kurt      = ', (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3. )

    #     plt.plot(ep_xx, ep_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()

    # if (test_SEP):
    #     # --------------
    #     # SEP tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     sep_samples = ufz.distributions.sample_sep(nn, loc=loc, sca=sca, xi=xi, beta=beta)
    #     plt.hist(sep_samples,bins=100,normed=True)
    #     print('SEP: mean_samp      = ', np.mean(sep_samples))
    #     print('SEP: std_samp       = ', np.std(sep_samples))
    #     print('SEP: skew_samp      = ', ss.skew(sep_samples))
    #     print('SEP: kurt_samp      = ', ss.kurtosis(sep_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of SEP
    #     dx = 0.01
    #     sep_xx = np.arange(-10,10,dx)
    #     sep_yy = ufz.distributions.ssep(sep_xx, loc=loc, sca=sca, skew=xi, kurt=beta)

    #     mean_num = np.sum(sep_xx*sep_yy*dx)
    #     var_num  = np.sum((sep_xx-mean_num)**2*sep_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((sep_xx-mean_num)/std_num)**3*sep_yy*dx)
    #     kurt_num  = np.sum(((sep_xx-mean_num)/std_num)**4*sep_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('SEP: mean_num       = ', mean_num )
    #     print('SEP: std_num        = ', np.sqrt(var_num) )
    #     print('SEP: skew_num       = ', skew_num )
    #     print('SEP: kurt_num       = ', kurt_num )
    #     print('')

    #     plt.plot(sep_xx, sep_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()


    # if (test_Studentt):
    #     # --------------
    #     # Student-t tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     studentt_samples = ufz.distributions.sample_studentt(nn, nu, loc=loc, sca=sca)
    #     plt.hist(studentt_samples,bins=100,normed=True)
    #     print('Student-t:  mean_samp      = ', np.mean(studentt_samples),                 ' ~ ',loc,                    ' : theoretical')
    #     print('Student-t:  std_samp       = ', np.std(studentt_samples),                  ' ~ ',sca*np.sqrt(nu/(nu-2.)),' : theoretical')
    #     print('Student-t:  skew_samp      = ', ss.skew(studentt_samples),                 ' ~ ',0.0,                    ' : theoretical')
    #     print('Student-t:  kurt_samp      = ', ss.kurtosis(studentt_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of Student-t
    #     dx = 0.0001
    #     studentt_xx = np.arange(-50+loc,50+loc,dx)
    #     studentt_yy = ufz.distributions.studentt(studentt_xx, nu, loc=loc, sca=sca)

    #     mean_num = np.sum(studentt_xx*studentt_yy*dx)
    #     var_num  = np.sum((studentt_xx-mean_num)**2*studentt_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((studentt_xx-mean_num)/std_num)**3*studentt_yy*dx)
    #     kurt_num  = np.sum(((studentt_xx-mean_num)/std_num)**4*studentt_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('Student-t:  mean_num       = ', mean_num )
    #     print('Student-t:  std_num        = ', std_num  )
    #     print('Student-t:  skew_num       = ', skew_num )
    #     print('Student-t:  kurt_num       = ', kurt_num )
    #     print('')

    #     plt.plot(studentt_xx, studentt_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()

    # if (test_ssStudentt):
    #     sstudentt_fs_samp = ufz.distributions.sample_sstudentt_fs(nn, nu, xi=xi)
    #     mean_sstudentt_fs  = ufz.distributions.sstudentt_mean(nu, skew=xi)
    #     std_sstudentt_fs   = ufz.distributions.sstudentt_std(nu, skew=xi)

    #     print('ssStudent-t:  mean_func      = ', mean_sstudentt_fs)
    #     print('ssStudent-t:  std_func       = ', std_sstudentt_fs)
    #     print('')
    #     print('ssStudent-t:  mean_samp      = ', np.mean(sstudentt_fs_samp))
    #     print('ssStudent-t:  std_samp       = ', np.std(sstudentt_fs_samp))

    # if (test_sStudentt):
    #     # --------------
    #     # skewed Student-t tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     # sstudentt_samples = ufz.distributions.sample_sstudentt(nn, nu, loc=loc, sig=sca*np.sqrt(nu/(nu-2.)), xi=xi)
    #     sstudentt_samples = ufz.distributions.sample_sstudentt(nn, nu, loc=loc, sca=sca, xi=xi)
    #     plt.hist(sstudentt_samples,bins=100,normed=True)
    #     print('sStudent-t:  mean_samp      = ', np.mean(sstudentt_samples))
    #     print('sStudent-t:  std_samp       = ', np.std(sstudentt_samples))
    #     print('sStudent-t:  skew_samp      = ', ss.skew(sstudentt_samples))
    #     print('sStudent-t:  kurt_samp      = ', ss.kurtosis(sstudentt_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of Student-t
    #     dx = 0.0001
    #     sstudentt_xx = np.arange(-30+loc,30+loc,dx)
    #     # sstudentt_yy = ufz.distributions.ssstudentt(sstudentt_xx, nu, loc=loc, sig=sca*np.sqrt(nu/(nu-2.)), skew=xi)
    #     sstudentt_yy = ufz.distributions.ssstudentt(sstudentt_xx, nu, loc=loc, sca=sca, skew=xi)

    #     mean_num = np.sum(sstudentt_xx*sstudentt_yy*dx)
    #     var_num  = np.sum((sstudentt_xx-mean_num)**2*sstudentt_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((sstudentt_xx-mean_num)/std_num)**3*sstudentt_yy*dx)
    #     kurt_num  = np.sum(((sstudentt_xx-mean_num)/std_num)**4*sstudentt_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('sStudent-t:  mean_num       = ', mean_num )
    #     print('sStudent-t:  std_num        = ', std_num  )
    #     print('sStudent-t:  skew_num       = ', skew_num )
    #     print('sStudent-t:  kurt_num       = ', kurt_num )
    #     print('')

    #     plt.plot(sstudentt_xx, sstudentt_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()
