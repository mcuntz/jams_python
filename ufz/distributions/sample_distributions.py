#!/usr/bin/env python
from __future__ import print_function
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

    Copyright 2016 Juliane Mai, Dmitri Kavetski, Matthias Cuntz
"""
import numpy as np

__all__ = ['sample_ep',     # sample from exponential power distribution of Box and Tiao (1992)
           'sample_sep',    # sample from skew exponential power distribution
           'sample_sep_fs', # sample from skew exponential power distribution after Fernandez and Steel (1998)
           'sample_st',     # sample from the skewed Student-t distribution
           'sample_st_fs',  # sample from skewed Student-t distribution after Fernandez and Steel (1998)
           'sample_t']      # sample from Student-t distribution


# --------------------------------------------------------------------


def sample_ep(nn, loc=0., sca=1., beta=0., sig=None):
    """
        Samples from the exponential power function of Box and Tiao (1992) ep(loc,sca,beta) are drawn.

        The distribution is symmetric.

        The kurtosis of this sample can be calculated by
             bb   = 2./(beta+1.)
             kurt = (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3.


        Definition
        ----------
        def sample_ep(nn, loc=0., sca=1., beta=0., sig=None):


        Input
        -----
        nn         number of samples


        Optional Input
        --------------
        loc        location
        sca        scale
        beta       parameter controlling kurtosis
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples drawn from exponential power function


        Examples
        --------


        History
        -------
        Written  JM, May 2016
    """

    if sig is not None: sca = sig

    EP01 = sample_ep01(nn, beta=beta)
    EP   = loc + sca * EP01

    return EP


def sample_ep01(nn, beta=0.):
    """
        Samples from the exponential power function of Box and Tiao (1992) ep(0,1,beta) are drawn.

        The distribution is symmetric.

        The kurtosis of this sample can be calculated by
             bb   = 2./(beta+1.)
             kurt = (sp.gamma(5./bb)*sp.gamma(1./bb))/((sp.gamma(3./bb))**2) - 3.


        Definition
        ----------
        def sample_ep01(nn, beta=0.):


        Input
        -----
        nn         number of samples


        Optional Input
        --------------
        beta       parameter which controls the kurtosis


        Output
        ------
        Samples drawn from exponential power function


        Examples
        --------


        Literature
        --------
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.
        --> Steps (1)-(3) described on page 5


        History
        -------
        Written  JM, May 2016
    """
    import scipy.special as sp

    bb = (1.+beta)/2.

    # (1) Generate a sample gg from gamma distribution with shape
    #     parameter (1+beta)/2 and scale parameter 1.
    gg = np.random.gamma(bb, 1.0, nn)

    # (2) Generate a random sign ss (+1 or -1) with equal probability.
    ss = np.where(np.random.rand(nn) < 0.5, -1., 1.)

    # (3) Compute EP which is a sample from the exponential power distribution EP(0,1,beta)
    EP01 = ss * np.abs(gg)**bb * np.sqrt(sp.gamma(bb)) / np.sqrt(sp.gamma(3.*bb))

    return EP01


# --------------------------------------------------------------------


def sample_sep(nn, loc=0., sca=1., xi=1., beta=0., sig=None):
    """
        Samples from the skew exponential power distribution with
        given location, scale, and parameters controlling skewness and kurtosis.

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
        def sample_sep(nn, loc=0., sca=1., xi=1., beta=0., sig=None):


        Input
        -----
        nn         number of samples


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         parameter which controls the skewness
        beta       parameter which controls the kurtosis
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples from the (general) skew exponential power


        Examples
        --------
        None


        History
        -------
        Written,  JM, May 2016
    """

    if sig is not None: sca = sig

    sSEP = sample_sep01(nn, xi=xi, beta=beta)
    SEP  = loc + sca * sSEP

    return SEP


def sample_sep01(nn, xi=1., beta=0.):
    """
        Samples from the skew exponential power distribution with location zero and scale one.


        Definition
        ----------
        def sample_sep01(nn, xi=1., beta=0.):


        Input
        -----
        nn         number of samples


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
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.
        --> Steps (6) described on page 5


        History
        -------
        Written,  JM, May 2016
    """
    from ufz.distributions import sep_fs_mean, sep_fs_std

    SEP_fs = sample_sep01_fs(nn, xi=xi, beta=beta)

    # (6) Standardize SEP_fs
    mean_sep_fs = sep_fs_mean(skew=xi, kurt=beta)
    std_sep_fs  = sep_fs_std(skew=xi, kurt=beta)
    sSEP        = (SEP_fs - mean_sep_fs) / std_sep_fs   # standardized SEP (=Schoups and Vrugt's a_t)

    return sSEP


# --------------------------------------------------------------------


def sample_sep_fs(nn, loc=0., sca=1., xi=1., beta=0., sig=None):
    """
        Samples from the skew exponential power distribution with
        given location, scale, and parameters controlling skewness and kurtosis,
        using the approach of Fernandez and Steel (1998).

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
        def sample_sep_fs(nn, loc=0., sca=1., xi=1., beta=0., sig=None):


        Input
        -----
        nn         number of samples


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         parameter which controls the skewness
        beta       parameter which controls the kurtosis
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples from the (general) skew exponential power


        Examples
        --------
        None


        History
        -------
        Written,  JM, May 2016
    """

    if sig is not None: sca = sig

    sSEP = sample_sep01_fs(nn, xi=xi, beta=beta)
    SEP  = loc + sca * sSEP

    return SEP


def sample_sep01_fs(nn, xi=1., beta=0.):
    """
        Samples from the skew exponential power distribution with location zero and scale one and
        parameters controlling skewness and kurtosis,
        using the approach of Fernandez and Steel (1998).

        If xi = 1. then mean=0. and stddev=1.


        Definition
        ----------
        def sample01_sep_fs(nn, xi=1., beta=0.):


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
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.
        Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and predictive
            inference of hydrologic models with correlated, heteroscedastic, and non-Gaussian errors.
            Water Resources Research 46, W10531.
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
    SEP_fs = skew_fernandez_steel(EP, xi)

    return SEP_fs



# --------------------------------------------------------------------



def sample_st(nn, nu, loc=0., sca=1., xi=1., sig=None):
    """
        Samples from the skew Student-t distribution with
        given location, scale, and parameter controlling skewness.

        The parameter controlling skewness xi needs to be positive.
             1 = symmetric
            >1 = right skewed
            <1 = left skewed


        Definition
        ----------
        def sample_st(nn, nu, loc=0., sca=1., xi=1., sig=None):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         parameter controlling skewness
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples from the skew Student-t distribution


        Examples
        --------
        None


        History
        -------
        Written,  JM, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)

    sst = sample_st01(nn, nu, xi=xi)
    st  = loc + sca * sst

    return st


def sample_st01(nn, nu, xi=1.):
    """
        Samples from the skew Student-t distribution with loc=0 and scale=1.


        Definition
        ----------
        def sample_st01(nn, nu, xi=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        xi         parameter controlling skewness


        Output
        ------
        Samples from the skew Student-t distribution iwth loc=0, sca=1


        Examples
        --------
        None


        History
        -------
        Written,  JM, May 2016
    """
    from ufz.distributions import st_fs_mean, st_fs_std

    st_fs = sample_st_fs(nn, nu, xi=xi)

    # (6) Standardize st_fs
    mean_st_fs = st_fs_mean(nu, skew=xi)
    std_st_fs  = st_fs_std(nu, skew=xi) * np.sqrt((nu-2.)/nu)
    sst        = (st_fs - mean_st_fs) / std_st_fs   # standardized skew Student-t

    return sst


# --------------------------------------------------------------------


def sample_st_fs(nn, nu, loc=0., sca=1., xi=1., sig=None):
    """
        Samples from the skew Student-t distribution with
        given location, scale, and parameter controlling skewness,
        using the approach of Fernandez and Steel (1998).

        The parameter controlling skewness xi needs to be positive.
             1 = symmetric
            >1 = right skewed
            <1 = left skewed


        Definition
        ----------
        def sample_st_fs(nn, nu, loc=0., sca=1., xi=1., sig=None):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        xi         parameter controlling skewness
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples from the skew Student-t distribution


        Examples
        --------
        None


        History
        -------
        Written,  JM, May 2016
    """

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)

    sst = sample_st01_fs(nn, nu, xi=xi)
    st  = loc + sca * sst

    return st


def sample_st01_fs(nn, nu, xi=1.):
    """
        Samples from the skew Student-t distribution obtained
        by using the approach of Fernandez and Steel (1998).


        Definition
        ----------
        def sample_st01_fs(nn, nu, xi=1.):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        xi       parameter controlling skewness


        Output
        ------
        Samples from the skew Student-t distribution after Fernandez and Steel (1998).


        Examples
        --------
        None


        Literature
        --------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  JM, May 2016
    """

    # (1) - (3)
    t = sample_t01(nn, nu)

    # (4) Generate a random sign w_t (+1 or -1) with probabilities 1-xi/(xi+1/xi) and xi/(xi+1/xi)
    #     wt = np.where(np.random.rand(nn) > 1.-xi/(xi+1./xi), 1, -1)
    # (5) Compute SEP_t which is a sample from the skew exponential power distribution
    #     SEP(mu_xi, sigma_xi, xi, beta) with mu_xi and sigma_xi given by (A5) and (A6)
    #     SEPt = wt * np.abs(EPt) * xi**wt
    st_fs = skew_fernandez_steel(t, xi)

    return st_fs


# --------------------------------------------------------------------


def sample_t(nn, nu, loc=0., sca=1., sig=None):
    """
        Samples from the Student-t distribution t(nu,loc,sca) are drawn.

        The mean is loc, standard deviation is sca*np.sqrt(nu/(nu-2)).

        The distribution is symmetric.


        Definition
        ----------
        def sample_t(nn, nu, loc=0., sca=1., sig=None):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Optional Input
        --------------
        loc        location
        sca        scale
        sig        standard deviation, overwrites sca


        Output
        ------
        Samples drawn from Student-t distribution


        Examples
        --------
        None


        History
        -------
        Written  JM, May 2016
    """
    # import scipy.stats as ss

    if sig is not None: sca = sig * np.sqrt((nu-2.)/nu)

    t01 = sample_t01(nn, nu)
    t   = t01 * sca + loc
    # t = ss.t.rvs(nu, loc=loc, scale=sca, size=nn)

    return t


def sample_t01(nn, nu):
    """
        Samples from the Student-t distribution with loc=0 and scale=1 t(df=nu, loc=0, sca=1) are drawn.

        The distribution is symmetric.


        Definition
        ----------
        def sample_t01(nn, nu):


        Input
        -----
        nn       -> number of samples
        nu       -> degrees of freedom


        Output
        ------
        Samples drawn from Student-t distribution with loc=0, scale=1.


        Examples
        --------
        None


        History
        -------
        Written  JM, May 2016
    """
    import scipy.stats as ss

    t01 = ss.t.rvs(nu, loc=0.0, scale=1.0, size=nn)

    return t01


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


        Output
        ------
        Samples from the skewed version of symmetric distribution used to generate sample


        Examples
        --------
        None


        Literature
        ----------
        Fernandez C & Steel M (1998) On Bayesian modeling of fat tails and skewness.
            Journal of the American Statistical Association 93, 359-371.


        History
        -------
        Written,  JM, May 2016
    """

    rr        = np.random.rand(np.size(sample_symmetric))
    absSymm   = np.abs(sample_symmetric)
    sample_fs = np.where(rr > 1.-xi/(xi+1./xi), absSymm*xi, -absSymm/xi)

    if np.iterable(sample_symmetric):
        return sample_fs
    else:
        return float(sample_fs)


# ------------------------------------------------------------


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # ------------------------------------------------------------
    # # Do some tests for the sampling of distributions
    # # ------------------------------------------------------------

    # import scipy.special as sp
    # import scipy.stats as ss
    # import ufz
    # import matplotlib.pylab as plt
    # plt.close("all")

    # nn    = 10000000   # number of random samples drawn
    # loc   = 0.0       # location parameter = mean
    # sca   = 2.0       # scale parameer = standard deviation
    # xi    = 1.5       # skew parameter (0,inf) 1=symmetric, >1=right skewed, <1=left skewed
    # beta  = -0.5       # ex. kurtosis parameter (-1,1] 1=fat tail, -1=thin tail
    # nu    = 6        # degrees of freedom for student-t

    # test_ep  = True
    # test_sep = True
    # test_t   = True
    # test_st_fs = True
    # test_st  = True

    # print('loc  = ', loc)
    # print('sca  = ', sca)
    # print('xi   = ', xi)
    # print('beta = ', beta)
    # print('nu   = ', nu)
    # print('')

    # if test_ep:
    #     # --------------
    #     # EP tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     ep_samples = ufz.distributions.sample_ep(nn, loc=loc, sca=sca, beta=beta)
    #     plt.figure()
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

    # if test_sep:
    #     # --------------
    #     # SEP tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     sep_samples = ufz.distributions.sample_sep(nn, loc=loc, sca=sca, xi=xi, beta=beta)
    #     plt.figure()
    #     plt.hist(sep_samples,bins=100,normed=True)
    #     print('SEP: mean_samp      = ', np.mean(sep_samples))
    #     print('SEP: std_samp       = ', np.std(sep_samples))
    #     print('SEP: skew_samp      = ', ss.skew(sep_samples))
    #     print('SEP: kurt_samp      = ', ss.kurtosis(sep_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of SEP
    #     dx = 0.01
    #     sep_xx = np.arange(-10,10,dx)
    #     sep_yy = ufz.distributions.sep(sep_xx, loc=loc, sca=sca, skew=xi, kurt=beta)

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


    # if test_t:
    #     # --------------
    #     # Student-t tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     t_samples = ufz.distributions.sample_t(nn, nu, loc=loc, sca=sca)
    #     plt.figure()
    #     plt.hist(t_samples,bins=100,normed=True)
    #     print('Student-t:  mean_samp      = ', np.mean(t_samples),                 ' ~ ',loc,                    ' : theoretical')
    #     print('Student-t:  std_samp       = ', np.std(t_samples),                  ' ~ ',sca*np.sqrt(nu/(nu-2.)),' : theoretical')
    #     print('Student-t:  skew_samp      = ', ss.skew(t_samples),                 ' ~ ',0.0,                    ' : theoretical')
    #     print('Student-t:  kurt_samp      = ', ss.kurtosis(t_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of Student-t
    #     dx = 0.0001
    #     t_xx = np.arange(-50+loc,50+loc,dx)
    #     t_yy = ufz.distributions.t(t_xx, nu, loc=loc, sca=sca)

    #     mean_num = np.sum(t_xx*t_yy*dx)
    #     var_num  = np.sum((t_xx-mean_num)**2*t_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((t_xx-mean_num)/std_num)**3*t_yy*dx)
    #     kurt_num  = np.sum(((t_xx-mean_num)/std_num)**4*t_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('Student-t:  mean_num       = ', mean_num )
    #     print('Student-t:  std_num        = ', std_num  )
    #     print('Student-t:  skew_num       = ', skew_num )
    #     print('Student-t:  kurt_num       = ', kurt_num )
    #     print('')

    #     plt.plot(t_xx, t_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()

    # if (test_st_fs):
    #     st_fs_samp = ufz.distributions.sample_st_fs(nn, nu, xi=xi)
    #     mean_st_fs  = ufz.distributions.st_fs_mean(nu, skew=xi)
    #     std_st_fs   = ufz.distributions.st_fs_std(nu, skew=xi)

    #     print('ssStudent-t:  mean_func      = ', mean_st_fs)
    #     print('ssStudent-t:  std_func       = ', std_st_fs)
    #     print('')
    #     print('ssStudent-t:  mean_samp      = ', np.mean(st_fs_samp))
    #     print('ssStudent-t:  std_samp       = ', np.std(st_fs_samp))

    # if (test_st):
    #     # --------------
    #     # skewed Student-t tests
    #     # --------------
    #     print('')
    #     print('-----------------------------------------------------')
    #     # st_samples = ufz.distributions.sample_st(nn, nu, loc=loc, sig=sca*np.sqrt(nu/(nu-2.)), xi=xi)
    #     st_samples = ufz.distributions.sample_st(nn, nu, loc=loc, sca=sca, xi=xi)
    #     plt.figure()
    #     plt.hist(st_samples,bins=100,normed=True)
    #     print('sStudent-t:  mean_samp      = ', np.mean(st_samples))
    #     print('sStudent-t:  std_samp       = ', np.std(st_samples))
    #     print('sStudent-t:  skew_samp      = ', ss.skew(st_samples))
    #     print('sStudent-t:  kurt_samp      = ', ss.kurtosis(st_samples, fisher=True))
    #     print('')

    #     # Compare with pdf of Student-t
    #     dx = 0.0001
    #     st_xx = np.arange(-30+loc,30+loc,dx)
    #     # st_yy = ufz.distributions.st(st_xx, nu, loc=loc, sig=sca*np.sqrt(nu/(nu-2.)), skew=xi)
    #     st_yy = ufz.distributions.st(st_xx, nu, loc=loc, sca=sca, skew=xi)

    #     mean_num = np.sum(st_xx*st_yy*dx)
    #     var_num  = np.sum((st_xx-mean_num)**2*st_yy*dx)
    #     std_num  = np.sqrt(var_num)
    #     skew_num  = np.sum(((st_xx-mean_num)/std_num)**3*st_yy*dx)
    #     kurt_num  = np.sum(((st_xx-mean_num)/std_num)**4*st_yy*dx) - 3.0   # Fisher definition (Gaussian kurt=0.0)
    #     print('sStudent-t:  mean_num       = ', mean_num )
    #     print('sStudent-t:  std_num        = ', std_num  )
    #     print('sStudent-t:  skew_num       = ', skew_num )
    #     print('sStudent-t:  kurt_num       = ', kurt_num )
    #     print('')

    #     plt.plot(st_xx, st_yy, 'r-', lw=5, alpha=0.6)
    #     plt.show()
