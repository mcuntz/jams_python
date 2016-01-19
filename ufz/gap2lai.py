#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from scipy import integrate as ig
from scipy.stats import beta
from ufz.kernel_regression import kernel_regression, kernel_regression_h

__all__ = ['gap2lai', 'leafprojection']

def gap2lai(tgap, lgap, G, alpha, boot=False):
    """
        Computes effective leaf area index Le, total leaf area index Lt and
        clumping factor omega out of canopy gap fraction observations. When
        boot=False, for each gap fraction the output is calculated. If you give
        boot a number, the input is bootstrapped and mean +- standard error of
        the mean is returned as output. If you don't have large gap fraction
        lgap, you can give an array of zeros assuming random distribution of the
        leaves, Le equals Lt then.
        
        
        Definition
        ----------
        def gap2lai(tgap, lgap, G, alpha, boot=False):


        Input
        -----
        tgap         array(N), total gap fraction, e.g. fraction of the number
                     of all gap pixels in your cover photo divided by the total
                     number of pixels in your image. [-]
        lgap         array(N), large gap fraction, e.g. fraction of the number
                     of pixels only found in the largest gap of your cover photo
                     divided by the total number of pixels in your image. [-]
        G            float, value of the leaf projection function at view zenith
                     angle alpha [-]
        alpha        float, view zenith angle at which your observations are
                     made. 0 equals zenith, 90 equals horizont. [deg]
                     
        
        Optional Input
        --------------
        boot         bool or int
                     if set False (default), for every gap fraction
                     value a separate output of Le, Lt and omega is computed.
                     If set to an int, the gap fractions are bootstrapped int times
                     and a mean +- standard error of the mean for Le, Lt and
                     omega is computed. WARNING: you are not allowed to average
                     the separate Le or Lt values after computation to get a
                     mean Le or Lt since you would introduce apparent clumping
                     correction.
                     If set True then boot=1.
        

        Output
        ------
        Lt           total leaf area index [m2leaf/m2ground]. if boot=False
                     (default), array(N) with Lt for each gap fraction. if
                     boot=int, array(2) with mean Lt and standard error of the
                     mean Lt 
        Le           effective leaf area index [m2leaf/m2ground]. if boot=False
                     (default), array(N) with Le for each gap fraction. if
                     boot=int, array(2) with mean Le and standard error of the
                     mean Le 
        omega        clumping index [-]. if boot=False (default), array(N) with
                     omega for each gap fraction. if boot=int, array(2) with 
                     mean omega and standard error of the mean omega 
        
        
        Literature
        ----------
        Macfarlane et al. (2007)
            Estimation of leaf area index in eucalypt forest using digital
            photography.
            Agricultural and Forest Meteorology 143, 176 - 188
        Ryu et al. (2010)
            On the correct estimation of effective leaf area index: Does it
            reveal information on clumping effects?
            Agricultural and Forest Meteorology 150, 463 - 472 
        
        
        Examples
        --------
        >>> from autostring import astr

        >>> # Create some data
        >>> tgap  = np.arange(0.,1.1,0.1)
        >>> lgap  = np.arange(0.,1.1,0.1)*0.1
        >>> G     = 0.5
        >>> alpha = 53.7
     
        >>> # compute each observation separately
        >>> Lt, Le, omega = gap2lai(tgap, lgap, G, alpha)
        >>> print(astr(Lt, 2, join=True))
        --   2.81 1.97 1.47 1.11 0.84 0.62 0.43 0.27 0.13 0.00
        >>> print(astr(Le, 2, join=True))
        --   2.73 1.91 1.43 1.08 0.82 0.60 0.42 0.26 0.12 0.00
        >>> print('!', astr(omega, 2, join=True), '!')
        ! --   0.97 0.97 0.97 0.97 0.98 0.98 0.98 0.99 0.99 --   !
     
        >>> # compute mean and standard error of the mean
        >>> Lt, Le, omega = gap2lai(tgap, lgap, G, alpha, boot=10000)
        >>> print(astr(Lt, 1, join=True))
        0.9 0.2
        >>> print(astr(Le, 1, join=True))
        0.8 0.2
        >>> print(astr(omega, 1, join=True))
        1.0 0.0

        >>> # if you don't have lgap, Le and Lt are equal
        >>> Lt, Le, omega = gap2lai(tgap, np.zeros_like(tgap), G, alpha)
        >>> print(np.ma.all(Lt==Le))        
        True
        
        
        License
        -------
        This file is part of the UFZ Python package.
    
        It is NOT released under the GNU Lesser General Public License, yet.
        
        If you use this routine, please contact Arndt Piayda.
        
        Copyright 2014 Arndt Piayda


        History
        -------
        Written,  AP, Sep 2014
    """
    
    ###########################################################################
    # safety checks
    tgaps, lgaps = tgap.size, lgap.size
    assert (tgaps==lgaps) & (tgap.ndim==1) & (lgap.ndim==1), "gap2lai: tgap and lgap must be 1D and of same size"
    assert isinstance(G, float) & isinstance(alpha, float), "gap2lai: G and alpha must be of type float"    
    assert (np.ma.max(tgap)<=1.) & (np.ma.max(lgap)<=1.), "gap2lai: max of lgap and tgap must be <=1."
    assert (np.ma.min(tgap)>=0.) & (np.ma.min(lgap)>=0.), "gap2lai: min of lgap and tgap must be >=0."
    assert np.ma.all((tgap-lgap)>=0.), "gap2lai: tgap must be larger or equal lgap"
    assert (np.signbit(G)==False) & (np.signbit(alpha)==False), "gap2lai: G and alpha must be positive"
    
    ###########################################################################
    # bootstrap input if desired
    if boot:
        tgap, lgap = gap2lai_bootstrap2(tgap, lgap, boot)
    #print(tgap)
    # radians from angle in degrees
    rad = np.deg2rad(alpha)
    # fraction of foliage cover, defined as proportion of ground area covered
    # by vertical projection of foliage and branches (crown cover)    
    fc = 1.- lgap
    # proportion of ground area covered by vertical projection of foliage and
    # branches within the perimeter of the crowns of individual plants
    ff = 1.- tgap
    # crown porosity
    phi = 1. - ff/fc
    # effective leaf area index without clumping correction
    Le = -np.ma.log(1.-ff)*np.ma.cos(rad)/G
    Le = np.ma.where(Le==-0.0, 0.0, Le) # safety for case: tgap=1 
    # clumping index
    omega = ((1.-phi)*np.ma.log(1.-ff))/(np.ma.log(phi)*ff)
    # total leaf area index with clumping correction
    Lt = -fc*np.ma.log(phi)*np.ma.cos(rad)/G
    # if bootstrap is not false, mean and standard error of the mean is returned
    if boot:
        Lt    = np.ma.array([np.ma.mean(Lt), np.ma.std(Lt)]).flatten()
        Le    = np.ma.array([np.ma.mean(Le), np.ma.std(Le)]).flatten()
        omega = np.ma.array([np.ma.mean(omega), np.ma.std(omega)]).flatten()
    
    return Lt, Le, omega


def leafprojection(alpha, theta, t360=False, kernel=False, min=0., max=90.,
                   step=5., h=None, boot=False):
    """
        Computes leaf projection funtion G for a certain view zenith angle alpha
        from measurements of leaf zenith angles theta (angle of the leaf surface
        normal to the zenith). When boot is set to an int, the input thetas are
        bootstrapped and mean +- standard error of the mean is returned as
        output. You can coose between the fast beta distribution function
        (default) or the more representative kernel regression.
        
        
        Definition
        ----------
        def leafprojection(alpha, theta, t360=False, kernel=False, min=0., max=90.,
                           step=5., h=None, boot=False):


        Input
        -----
        alpha        float, view zenith angle at which you want the leaf
                     projection value G. 0 equals zenith, 90 equals horizont.
                     [deg]
        theta        np.ma.array(N), leaf inclination zenith angle measurements
                     (angle of the leaf surface normal to the zenith) [deg] -180
                     to 180 with 0 pointing to zenith. If t360=True: leaf zenith
                     angle measurements [deg] 0 to 360 with 0 pointing to
                     zenith, clockwise.        


        Optional Input
        --------------
        t360         bool, if False (default) theta  -180 to 180, else 0 to 360
        kernel       bool, if False (default) beta function is used, if True,
                     kernel regression is used (takes much more time) but closer
                     to measurements 
        min          float, minimum of theta histogram (only used for kernel
                     regression, default: 0.)
        max          float, maximum of theta histogram (only used for kernel
                     regression, default: 90.)
        step         float, step size of theta histogram (only used for kernel
                     regression, default: 5.)
        h            float, if None (default) correlation lenght for kernel
                     regression is calculated from the input, else you can give
                     a float here.
        boot         bool or int, if set False (default), the original theta
                     distribution is used for the computations. If set to an
                     int, the theta distribution is bootstrapped int times and
                     a mean +- standard error of the mean for the projection
                     function G is computed
        

        Output
        ------
        G            float, value of the leaf projection function at the view
                     zenith angle alpha, given the leaf view zenith angle
                     distribution measurements theta [-]. If boot=int, tuple(2)
                     with mean G and standard error of the mean G at alpha
                     
        
        Literature
        ----------
        Goel & Strebel (1984)
            Simple Beta Distribution Representation of Leaf Orientation in
            Vegetation Canopies.
            Agron. J., 1984, 76, 800-802
        Wang et al. (2007)
            Comparison of leaf angle distribution functions: Effects on
            extinction coefficient and fraction of sunlit foliage.
            Agricultural and Forest Meteorology, 2007, 143, 106 - 122
        
        
        Examples
        --------
        >>> # use beta ditribution for a view zenith angle of 57.3 deg and
        >>> # measured leaf inclination in -180 to 180 deg.
        >>> np.random.seed(1)
        >>> theta = np.ma.array(np.random.random(50)*180.)
        >>> print(np.ma.round(leafprojection(57.3, theta), 2))
        0.5
        
        >>> # use beta ditribution for a view zenith angle of 57.3 deg and
        >>> # measured leaf inclination in 0 to 360 deg.
        >>> np.random.seed(1)
        >>> theta = np.ma.array(np.random.random(50)*360.)
        >>> print(np.ma.round(leafprojection(57.3, theta, t360=True), 2))
        0.5
        
        >>> # use kernel regression for a view zenith angle of 57.3 deg and
        >>> # measured leaf inclination in 0 to 360 deg.
        >>> np.random.seed(1)
        >>> theta = np.ma.array(np.random.random(50)*360.)
        >>> print(np.ma.round(leafprojection(57.3, theta, t360=True, kernel=True), 2))
        0.5
        
        >>> # use kernel regression for a view zenith angle of 57.3 deg and
        >>> # measured leaf inclination in 0 to 360 deg. and set h fix to 1.0
        >>> np.random.seed(1)
        >>> theta = np.ma.array(np.random.random(50)*360.)
        >>> print(np.ma.round(leafprojection(57.3, theta, t360=True, kernel=True, h=1.0), 2))
        0.5
        
        >>> # use kernel regression for a view zenith angle of 57.3 deg and
        >>> # measured leaf inclination in 0 to 360 deg. Estimate the error with
        >>> # 10 bootstraps
        >>> np.random.seed(1)
        >>> theta = np.ma.array(np.random.random(50)*360.)
        >>> from autostring import astr
        >>> print(astr(leafprojection(57.3, theta, t360=True, kernel=True, boot=10), 3, joinall=True))
        0.503 0.003
        
        
        License
        -------
        This file is part of the UFZ Python package.
    
        It is NOT released under the GNU Lesser General Public License, yet.
        
        If you use this routine, please contact Arndt Piayda.
        
        Copyright 2014 Arndt Piayda


        History
        -------
        Written,  AP, Sep 2014
    """
    
    assert not np.all(theta.mask) & np.ma.isMA(theta), 'leafprojection: theta is fully masked.'
    
    ###########################################################################
    # normalizing function for angles
    def t(thetal):
        return 2.*thetal/np.pi
    
    # beta distribution
    def f_beta(thetal, mu, nu):
        return beta.pdf(t(thetal), mu, nu)*2./np.pi
        
    # azimuth condition
    def vartheta(thetal, theta):
        return np.ma.arccos((1./np.ma.tan(theta)) * (1./np.ma.tan(thetal)))
    
    #undefined for 90 degrees due to tangens
    def phi(thetal, theta):
        return np.ma.where(np.ma.absolute((1./(np.ma.tan(theta)) * (1./np.ma.tan(thetal))))>1.,
                        np.ma.cos(theta)*np.ma.cos(thetal),
                        np.ma.cos(theta)*np.ma.cos(thetal)*(1. + (2./np.pi)*(np.ma.tan(
                        vartheta(thetal, theta))-vartheta(thetal, theta))))
    
    # thetal must be first argument, because first argument is taken for
    # integration by .quad
    def integ(thetal, theta, mu, nu):
        return phi(thetal, theta)*f_beta(thetal, mu, nu)
    def integ_k(thetal, theta, mbins, hist, h):
        ker = kernel_regression(mbins, hist, h, xout=np.array([t(thetal)]))
        return phi(thetal, theta)*ker*2./np.pi
    
    # leaf projection function with beta function
    def G(theta, mu, nu):
        if isinstance(theta, (int, float)):
            return ig.quad(integ, 0., np.pi/2., args=(theta, mu, nu))[0]
        else:
            out = np.array([ig.quad(integ, 0, np.pi/2.,args=(x, mu, nu))[0]
                             for x in theta])
            return np.ma.array(out, mask=theta.mask)
            
    # leaf projection function with kernel regression
    def G_ker(theta, mbins, hist, h):
        if isinstance(theta, (int, float)):
            return ig.quad(integ_k, 0, np.pi/2., args=(theta, mbins, hist, h))[0]
        else:
            out = np.array([ig.quad(integ_k, 0, np.pi/2.,args=(x, mbins, hist, h))[0]
                             for x in theta])
            return np.ma.array(out, mask=theta.mask)
    
    ###########################################################################
    # converting angles to 0-90 deg span and to rad
    if t360:
        theta = np.deg2rad(np.ma.where(theta<=90., theta,
                           np.ma.where(theta<=180., 180.-theta,
                           np.ma.where(theta<=270., theta-180., theta-270.))))
    else:
        theta = np.deg2rad(np.ma.where(np.ma.absolute(theta)<=90.,
                                       np.ma.absolute(theta),
                                       180.-np.ma.absolute(theta)))
    rad = np.deg2rad(alpha)
    
    if not kernel:
        #######################################################################
        # parameter determination for beta function
        if boot:
            theta = gap2lai_bootstrap(theta, boot)

        sig_0_2 = np.ma.mean(t(theta), axis=0)*(1.-np.ma.mean(t(theta), axis=0))
        sig_t_2 = np.ma.var(t(theta), axis=0)
        mu      = np.ma.mean(t(theta), axis=0)*(sig_0_2/sig_t_2-1.)
        nu      = (1.-np.ma.mean(t(theta), axis=0))*(sig_0_2/sig_t_2-1.)
        
        if boot:
            G = np.ma.array([G(rad, x[0], x[1]) for x in zip(mu, nu)])
            return np.ma.mean(G), np.ma.std(G)
        else:
            return G(rad, mu, nu)

    else:
        #######################################################################
        # parameter determination for kernel regression
        # histogram
        if boot:
            theta = gap2lai_bootstrap(theta, boot)

        min, max, step  = np.deg2rad(min), np.deg2rad(max), np.deg2rad(step)
        hist, bins, mbins = gap2lai_axishist(theta, min, max, step, t)
        
        # kernel h estimation
        if boot:
            G = np.ma.empty(boot)
            for i in xrange(boot):
                if not h:
                    h = kernel_regression_h(mbins[:,i], hist[:,i])
                G[i] = G_ker(rad, mbins[:,i], hist[:,i], h)
            return np.ma.mean(G), np.ma.std(G)
        else:
            h = kernel_regression_h(mbins, hist)
            return G_ker(rad, mbins, hist, h)

def gap2lai_axishist(x, min, max, step, t):
    """
    gap2lai: histogram along 0 axis of an array
    """
    inbins = np.arange(t(min),t(max+step),t(step))
    
    if x.ndim==1:
        hist, bins = np.histogram(t(x.compressed()), bins=inbins, density=True)
        mbins = bins[:-1]+t(step)/2.
    else:
        bs, xs = inbins.size, x.shape[1]
        hist, bins = np.ma.empty((bs-1, xs)), np.ma.empty((bs, xs))
        for i in xrange(xs):
            hist[:,i], bins[:,i] = np.histogram(t(x[:,i].compressed()), bins=inbins, density=True)
        mbins = bins[:-1,:]+t(step)/2.
                        
    return hist, bins, mbins

def gap2lai_bootstrap(x, boot):
    """
    gap2lai: bootstrap a 1D array
    """
    xs = x.size
    assert x.ndim==1, "bootstrap: x and y must be 1D and of same size"
    assert isinstance(boot, int), "bootstrap: boot must be of type int"
    
    x_boot = np.ma.empty((xs, boot))
    for j in xrange(boot):
        ind = np.random.randint(xs,size=xs)
        x_boot[:,j] = x[ind]
            
    return x_boot

def gap2lai_bootstrap2(x, y, boot):
    """
    gap2lai: bootstrap two 1D arrays with same indices
    """
    xs, ys = x.size, y.size
    assert (xs==ys) & (x.ndim==1) & (y.ndim==1), "bootstrap: x and y must be 1D and of same size"
    assert isinstance(boot, int), "bootstrap: boot must be of type int"
    
    x_boot, y_boot = np.ma.empty(boot), np.ma.empty(boot)
    for j in xrange(boot):
        ind = np.random.randint(xs,size=xs)
        x_boot[j] = np.ma.mean(x[ind])
        y_boot[j] = np.ma.mean(y[ind])
    
    return x_boot, y_boot

if __name__ == '__main__':
     import doctest
     doctest.testmod()

     # from autostring import astr

     # # Create some data
     # tgap  = np.arange(0.,1.1,0.1)
     # lgap  = np.arange(0.,1.1,0.1)*0.1
     # G     = 0.5
     # alpha = 53.7
     
     # # compute each observation separately
     # Lt, Le, omega = gap2lai(tgap, lgap, G, alpha)
     # print(np.ma.round(Lt,2))
     # print(astr(Lt, 2, join=True))
     # # --   2.81 1.97 1.47 1.11 0.84 0.62 0.43 0.27 0.13 0.00
     # print(np.ma.round(Le,2))
     # print(astr(Le, 2, join=True))
     # # --   2.73 1.91 1.43 1.08 0.82 0.6 0.42 0.26 0.12 0.00
     # print(np.ma.round(omega,2))
     # print(astr(omega, 2, join=True))
     # # --   0.97 0.97 0.97 0.97 0.98 0.98 0.98 0.99 0.99 --
     
     # # compute mean and standard error of the mean
     # Lt, Le, omega = gap2lai(tgap, lgap, G, alpha, boot=10000)
     # print(np.ma.round(Lt,1))
     # print(astr(Lt, 1, join=True))
     # # 0.9 0.2
     # print(np.ma.round(Le,1))
     # print(astr(Le, 1, join=True))
     # # 0.8 0.2
     # print(np.ma.round(omega,1))
     # print(astr(omega, 1, join=True))
     # # 1.0 0.0
     # # if you don't have lgap, Le and Lt are equal
     # Lt, Le, omega = gap2lai(tgap, np.zeros_like(tgap), G, alpha)
     # print(np.ma.all(Lt==Le))        
     # # True
