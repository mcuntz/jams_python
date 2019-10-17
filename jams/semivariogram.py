#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def semivariogram(x, y, v, nL, di, td, stype='omnidirectional', negscat=0.,
                  model='exponential', graph=False, lunit='m', p0=(1.,1.,1.),
                  runtimediag=False):
    """
        Calculates single or multiple experimental semivariogram(s) for
        spatial distributed data. X and Y coordinates and corresponding
        values are required in separate numpy arrays. Different theoretical
        semivariograms can be fitted or the best fitting is chosen.

        Results can be plotted in various graphs and the fitted model parameters
        are given to the output.

        Small code parts are based on Alghalandis.com


        Definition
        ----------
        def semivariogram(x,y,v,nL,di,td,stype='omnidirectional',negscat=0.,
                          model='exponential',graph=False,lunit='m',
                          p0=(0.5,0.5,100),runtimediag=False):


        Input
        -----
        x         numpy array with longitude ('easting')
        y         numpy array with latitude ('northing')
        v         numpy array with values
        nL        int with number of lags to calculate
        di        list with angles to consider (in degree). di's > 180
                  or < -180 are not allowed.
        td        int with one sided angle tolerance, so td=30 means
                  an angular span of 60 (in degree). Choose td big
                  enough, so that no "empty" angular span with no
                  samples appear. Otherwise, an exeption is raised.
                  td > 180 is not allowed.


        Optional Input
        --------------
        stype      type of semivariogram

            'omnidirectional' semivariogram is calculated for the hole
                              circle from 0 to 180 and 0 to -180. So
                              every direction is considered. (default)
            'directional'     semivariogram is calculated only for the
                              given angles in di. Orientation is not
                              considered, so that e.g. 45 deg equals
                              -135 deg. Due to that, negative angles
                              in di make no sense (and are not allowed)
                              here.
            'directional+orientational'
                              semivariogram is calculated only for the
                              given angles in di. Orientation is
                              considered. Every angle between 0 to 180
                              and 0 to -180 are allowed.

        negscat   It is common in geostatistics to neglect large lags which
                  tend to scatter and contain outliers. By setting negscat to
                  a value between 0 and 1 (both excluded) you can cut off the
                  tail of the experimental semivariogram and the fitting will
                  be done only to the remaining part. If it's set to False, the
                  hole semivariogram is calculated. (default=False)

        model     type of theoretical semivariogram model to fit to the
                  experimental semivariogram.

            'exponential'     exponential semivariogram model (default)
            'spherical'       sperical semivariogram model
            'gaussian'        gaussian semivariogram model
            'power'           power semivariogram model
            'noidea'          fitts all three models and returns the one
                              with the smallest coefficient of variation
                              (That does not mean, the returning model
                              is the best model for your data, it only
                              has the smallest parameter variation, but
                              it may help you with the decision for the
                              right model)
            'nomodel'         model fitting is disabled

        graph     If True, plotting is enabled. If False, plotting is
                  disabled. (default=True)

        lunit     Unit of the longitude and latitude (default=meter,
                  only used for plot labeling)

        p0        Initial guess for parameter estimation of the
                  theoretical semivariogram.
                  1st entry: equals approximately the nugget
                  2nd entry: equals approximately nugget-sill
                  3rd entry: equals approximately the range
                  (default=(0.5,0.5,100))
                  Sometimes you have to play around with it if no
                  fitting result can be reached.

        runtimediag
                  If True, during run time of the function, diagnostics
                  are printed to the console. If False, runtime diagnostics
                  print is disabled. (default=True)


        Output
        ------
        if model is set to 'nomodel':
        h         list(s) of lags
        g         list(s) of variances
        c         list(s) of samples per lag

        otherwise:
        nugget    height of nugget(s) [(unit of v)**2], always >=0
        sill      height if sill(s) [(unit of v)**2], always >=0
        irange     distance of range [unit of x and y], always >=0
        vark      coefficient(s) of variation averaged over all three
                  model parameter (goodness-of-fit) [-]
        pcov      matrices of parameter variance and covariances
                  of the model parameters
        h         array(s) of lags
        g         array(s) of variances
        c         array(s) of samples per lag
        mod       model function, either the one you selected or the model with
                  the best fit if you selected model='noidea'
        popt      array of optimized parameters for the model

        graphs:
        Figure 1  shows a scatter plot of your original geodata
        Figure 2  shows the experimental and theoretical semivariogram
        Figure 3  shows the number of samples in each lag and angle
        Figure 4  visualize di and td, the angles and angle tolerances
                  you have chosen.


        References
        ----------
        Small code parts based on Alghalandis.com.


        Examples
        --------
        # provide you some sample data:
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)

        # easting
        >>> x = np.array([557509.27,557518.11,557526.95,557535.79,557544.63,
        ...               557553.47,557544.63,557535.79,557526.95,557518.11,
        ...               557526.95,557535.79,557544.63,557553.47,557562.31,
        ...               557571.15,557562.31,557553.47,557544.63,557535.79,
        ...               557544.63,557553.47,557562.31,557571.15,557579.99,
        ...               557597.66,557579.99,557562.31,557544.63,557526.7,
        ...               557509.27,557491.60,557473.92,557491.6,557509.27,
        ...               557526.95,557544.63,557562.31,557579.99,557597.66,
        ...               557615.34,557650.7,557686.05,557756.76,557827.47,
        ...               557898.18,557827.47,557756.76,557686.05,557615.34,
        ...               557650.70,557686.07,557721.41,557756.76,557721.41,
        ...               557686.05,557650.70,557686.05,557650.7,557579.99,
        ...               557615.34,557509.27,557579.99,557544.63,557509.27,
        ...               557473.92,557438.56,557403.21,557438.56,557473.92,
        ...               557509.27,557544.63,557579.99,557615.34,557615.34,
        ...               557579.99,557544.63,557509.27,557473.92,557438.56,
        ...               557403.21,557367.85,557332.5,557367.85,557403.21,
        ...               557438.56,557473.92,557332.50,557261.79,557191.08,
        ...               557261.79,557332.50,557403.21,557473.92,557544.63,
        ...               557615.34])

        # northing
        >>> y = np.array([4332422.55,4332413.71,4332404.87,4332396.03,
        ...               4332387.19,4332396.03,4332404.87,4332413.71,
        ...               4332422.55,4332431.39,4332440.23,4332431.39,
        ...               4332422.55,4332413.71,4332404.87,4332413.71,
        ...               4332422.55,4332431.39,4332440.23,4332449.07,
        ...               4332457.91,4332449.07,4332440.23,4332431.39,
        ...               4332422.55,4332404.87,4332387.19,4332369.52,
        ...               4332351.84,4332369.52,4332387.19,4332404.87,
        ...               4332422.55,4332440.23,4332457.91,4332475.58,
        ...               4332493.26,4332475.58,4332457.91,4332440.23,
        ...               4332422.55,4332316.48,4332210.42,4332281.13,
        ...               4332351.84,4332422.55,4332493.26,4332563.97,
        ...               4332634.68,4332563.97,4332528.62,4332493.25,
        ...               4332457.91,4332422.55,4332387.19,4332351.84,
        ...               4332387.19,4332422.55,4332457.91,4332599.33,
        ...               4332493.26,4332599.33,4332528.62,4332563.97,
        ...               4332528.62,4332493.26,4332457.91,4332422.55,
        ...               4332387.19,4332351.84,4332316.48,4332281.13,
        ...               4332316.48,4332351.84,4332281.13,4332245.77,
        ...               4332210.42,4332245.77,4332281.13,4332316.48,
        ...               4332351.84,4332387.19,4332422.55,4332457.91,
        ...               4332493.26,4332528.62,4332563.97,4332563.97,
        ...               4332493.26,4332422.55,4332351.84,4332281.13,
        ...               4332210.42,4332139.71,4332069.00,4332139.71])

        # value
        >>> v = np.array([9.94691161e-01,7.94158417e-02,0.00000000e+00,
        ...               1.75837990e+00,0.00000000e+00,0.00000000e+00,
        ...               0.00000000e+00,0.00000000e+00,1.02915310e+00,
        ...               2.69597379e+00,2.14552427e+00,2.18417112e+00,
        ...               0.00000000e+00,8.96101277e-01,1.14034753e+00,
        ...               3.46398689e-01,3.01418491e-01,0.00000000e+00,
        ...               1.17920343e+00,1.09682206e+00,4.79485665e-01,
        ...               0.00000000e+00,1.83183398e+00,0.00000000e+00,
        ...               0.00000000e+00,0.00000000e+00,9.86233407e-02,
        ...               7.68290376e-02,2.63911513e-01,0.00000000e+00,
        ...               2.10013460e+00,0.00000000e+00,2.47535521e+00,
        ...               1.47047869e+00,8.00371532e-01,2.39448347e+00,
        ...               0.00000000e+00,2.26426861e+00,0.00000000e+00,
        ...               0.00000000e+00,1.13769438e+00,1.01969271e+00,
        ...               2.26036007e+00,2.38991410e+00,1.82558084e-03,
        ...               0.00000000e+00,0.00000000e+00,2.52583544e+00,
        ...               6.35195403e-01,2.43778382e+00,0.00000000e+00,
        ...               2.47738704e+00,8.83280548e-01,2.42328547e+00,
        ...               0.00000000e+00,2.41534081e+00,2.45629467e+00,
        ...               0.00000000e+00,2.50770630e+00,1.30382267e+00,
        ...               2.06891940e+00,9.17384801e-02,0.00000000e+00,
        ...               1.10185544e-01,2.53460688e+00,2.15217780e+00,
        ...               1.16908154e+00,1.70072787e-01,1.60603658e-01,
        ...               2.15438377e+00,2.32464926e+00,3.26255002e-01,
        ...               0.00000000e+00,1.48404530e+00,2.10439439e+00,
        ...               0.00000000e+00,0.00000000e+00,0.00000000e+00,
        ...               2.34663663e-01,1.46993948e+00,2.67691613e+00,
        ...               2.13262460e-02,1.01551520e+00,1.10878523e+00,
        ...               1.80374874e+00,1.85571813e+00,2.93929948e+00,
        ...               4.43192829e-01,2.55962879e+00,0.00000000e+00,
        ...               1.46545683e+00,1.75659977e+00,0.00000000e+00,
        ...               2.37093751e+00,0.00000000e+00,0.00000000e+00])

        # omnidirectional semivariogram with exponential model and fifty lags
        >>> td = 180
        >>> di = [0]
        >>> nL = 50
        >>> h, g, c = semivariogram(x,y,v,nL,di,td,stype='omnidirectional',
        ...                         negscat=0., model='nomodel',graph=False,lunit='m',
        ...                         p0=(0.5,0.5,100),runtimediag=False)
        >>> from autostring import astr
        >>> print(astr(g[0][0:9], 3))
        ['0.550' '0.776' '0.619' '0.849' '0.991' '1.033' '1.067' '1.106' '1.079']
        >>> nugget,sill,irange,cork,h,g,c,semi_mod,semi_popt = semivariogram(
        ...      x,y,v,nL,di,td,stype='omnidirectional',
        ...      negscat=0.,model='exponential',graph=False,lunit='m',
        ...      p0=(0.5,0.5,1./100.), runtimediag=False)
        >>> print(astr(nugget, 1))
        ['0.5']
        >>> print(astr(sill, 1))
        ['1.1']
        >>> print(np.around(irange, -1))
        [80.]

        # directional semivariogram with spherical model and fifty lags
        >>> td = 45
        >>> di = [0,90]
        >>> nL = 50
        >>> nugget,sill,irange,cork,h,g,c,semi_mod,semi_popt = semivariogram(
        ...     x,y,v,nL,di,td,stype='directional',negscat=0.,
        ...     model='spherical',graph=False,lunit='m',p0=(0.5,0.5,100),
        ...     runtimediag=False)
        >>> print(astr(nugget, 1))
        ['0.5' '0.5']
        >>> print(astr(sill, 1))
        ['1.1' '1.1']

        # directional+orientational semivariogram with gaussian model and
        # fifty lags
        >>> td = 30
        >>> di = [0,45,90,135,180,-45,-90,-135]
        >>> nL = 50
        >>> nugget, sill, irange, cork, h, g, c, semi_mod, semi_popt = semivariogram(
        ...     x, y, v, nL, di, td, stype='directional+orientational',
        ...     negscat=0., model='spherical', graph=False, lunit='m',
        ...     p0=(0.5,0.5,100), runtimediag=False)
        >>> soll = [0.36, 0.42, 0.54, 0.50, 0.20, 0.57, 0.65, 0.45]
        >>> print(np.allclose(nugget, soll, atol=0.1))
        True
        >>> print(astr(sill, 1))
        ['1.0' '1.2' '1.1' '1.2' '1.0' '1.1' '1.1' '1.1']

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2011-2014 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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


        History
        -------
        Written,  AP, Feb 2011
        Modified, AP, Nov 2012 - Parameter bounds for fitting algr.
                               - Minor errors corrected
                  MC, Nov 2012 - default graph=False
                  AP, Dec 2012 - documentation change
                  MC, Feb 2013 - change from fmin_l_bgfs_n to fmin_tnc because the former printed
                                 "Positive dir derivative in projection"
                                 "Using the backtracking step"
                                 even with iprint=-1, disp=0.
                  MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """

    # check input data
    assert (np.shape(x) == np.shape(y)) & (np.shape(y) == np.shape(v)), 'SemivariogramError: x, y and v must have the same dimensions.'

    di = np.array(di)
    if (stype == 'directional') & np.any(di<0.):
        raise TypeError('SemivariogramError: if you choose stype=directional,'
                        ' no negative directions in di are allowed')

    assert not (np.any(di>180.) | np.any(di<(-180.))), 'SemivariogramError: elements of di > 180 or < -180 are not allowed.'

    assert td <= 180, 'SemivariogramError: td > 180 is not allowed.'

    if model != 'nomodel':
        try:
            # from scipy.optimize import fmin_l_bfgs_b
            from scipy.optimize import fmin_tnc
        except ImportError:
            model = 'nomodel'
            print('SemivariogramWarning: the package scipy.optimize.curve_fit can'
                  ' not be found. It requires at least scipy 0.8.0. No fitting of'
                  ' experimental semivariogram(s) is possible.')

    if runtimediag:
        import time
        start = time.time()

    #---------------------------------------
    #compute all distances and angles between vectors x and y
    r, t, xr, n = distang(x,y)

    t = np.where(t==-180, 180, t)

    if stype == 'omnidirectional':
        headtitle = 'Omnidirectional Semivariogram'
        di = np.array([0])
        td = 180
    elif stype == 'directional':
        headtitle = 'Directional Semivariograms'
    elif stype == 'directional+orientational':
        headtitle = 'Directional and Orientational Semivariograms'
    else:
        raise NameError('SemivariogramError: wrong semivariogram type! '
         'Choose omnidirectional, directional or directional+orientational')

    #---------------------------------------
    # compute semivariogram
    h = []
    g = []
    c = []
    for i in di:
        hi, gi, ci = semivario(r, t, xr, n, v, nL, i, td, stype=stype, negscat=negscat)
        h.append(hi)
        g.append(gi)
        c.append(ci)

    if runtimediag:
        stop1 = time.time()
        print('Experimental semivariogram took %0.3f seconds' %(stop1-start))

    if model == 'exponential':
        func  = [expvar]
        irange = [exprange]
        subtitles = ['Exponential Semivariogram Model']
    elif model == 'spherical':
        func  = [sphvar]
        irange = [sphrange]
        subtitles = ['Spherical Semivariogram Model']
    elif model == 'gaussian':
        func  = [gauvar]
        irange = [gaurange]
        subtitles = ['Gaussian Semivariogram Model' ]
    elif model == 'noidea':
        func  = [expvar, sphvar, gauvar]
        irange = [exprange,sphrange,gaurange]
        subtitles = ['Exponential Semivariogram Model',
                     'Spherical Semivariogram Model',
                     'Gaussian Semivariogram Model' ]
    elif model == 'nomodel':
        subtitles = ['No Semivariogram Model']
    else:
        raise NameError('SemivariogramError: wrong semivariogram model! '
         'Choose exponential, spherical or gaussian')

    #---------------------------------------
    # fit semivariogram model
    if model != 'nomodel':
        popt = np.zeros((len(func),di.size,3))
        cork = np.zeros((len(func),di.size))
        for j in range(len(func)):
            for i in range(di.size):
                try:
                    # popti, f, d = fmin_l_bfgs_b(absdif, p0,
                    #                             args=(h[i], g[i], func[j]),
                    #                             approx_grad=1,
                    #                             bounds=[(0.,None),(0.,None),
                    #                                     (0.,None)],
                    #                             iprint=-1, disp=0)
                    popti, f, d = fmin_tnc(absdif, p0,
                                                args=(h[i], g[i], func[j]),
                                                approx_grad=1,
                                                bounds=[(0.,None),(0.,None),
                                                        (0.,None)],
                                                disp=0)
                    popt[j,i,:] = popti
                    cork[j,i] = np.corrcoef(g[i], func[j](h[i],popti))[0,1]
                except AttributeError:
                    print(('SemivariogramWarning: %s can not be fitted to' +
                           ' %i deg angle') %(subtitles[j], di[i]))

        ind = np.argmax(np.mean(cork, axis=1))

        if runtimediag:
            stop2 = time.time()
            print('Theoretical semivariogram took %0.3f seconds' %(stop2-stop1))

        nugget  = popt[ind,:,0]
        sill    = popt[ind,:,0] + popt[ind,:,1]
        irange  = irange[ind](sill,popt[ind,:,:])
        cork    = cork[ind]

        if runtimediag:
            print('Model: %s'%subtitles[ind])
            print('Nugget(s): %s'%nugget)
            print('Sill(s): %s'%sill)
            print('Range(s): %s'%irange)
            print('Coefficient(s) of correlation: %s'%cork)

    #---------------------------------------
    # plot
    if graph:
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        mpl.rc('font', size=20)
        mpl.rc('lines', linewidth=2)
        mpl.rc('axes', linewidth=1.5)
        # mpl.rc('xtick.major', width=1.5)
        # mpl.rc('ytick.major', width=1.5)
        mpl.rcParams['lines.markersize']=6
        mpl.rcParams['lines.markeredgewidth']=1
        mpl.rcParams['grid.linewidth']=1.5
        # mpl.rcParams['legend.frameon']=False
        mpl.rcParams['legend.numpoints']=1
        mpl.rcParams['legend.handlelength']=1
        mpl.rcParams['mathtext.default']='regular'

        # scatterplot of input data
        # fig1 = plt.figure('scatter plot')
        fig1 = plt.figure(1)
        sub1 = fig1.add_subplot(111, aspect='equal')
        scat = sub1.scatter(x,y,c=v,s=40,cmap=plt.cm.jet)
        sub1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub1.xaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub1.yaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub1.grid('on')
        sub1.set_title('scatter plot')
        plt.xlabel('easting [%s]' %lunit)
        plt.ylabel('northing [%s]' %lunit)
        fig1.autofmt_xdate(rotation=45)
        # plt.tight_layout(pad=1, h_pad=0, w_pad=0)
        # cbar need to be below autofm_xdate !!!???
        cbar = fig1.colorbar(scat,orientation='vertical',pad=0.05,shrink = 0.7)
        cbar.set_label('value')

        # semivariogram plot with/without fitted model
        # fig2 = plt.figure('semivariogram')
        fig2 = plt.figure(2)
        sub2 = fig2.add_subplot(111)
        sub2.plot(0,0,'.',c='k',markersize=1)
        for i in range(di.size):
            sub2.plot(h[i], g[i], 'o--', c=cm.jet(256/di.size*i),
                      label='%i'%di[i])
            if model != 'nomodel':
                h_mod_step = (np.max(h[i])-np.min(h[i]))/100.
                h_mod = np.arange(0.,np.max(h[i])+h_mod_step,h_mod_step)
                sub2.plot(h_mod, func[ind](h_mod, popt[ind,i,:]), '-',
                          c=cm.jet(256/di.size*i))
                sub2.plot([irange[i],irange[i]],
                          [0.,func[ind](irange[i], popt[ind,i,:])], '-.',
                          c=cm.jet(256/di.size*i))
                sub2.plot(irange[i],func[ind](irange[i], popt[ind,i,:]), 'rs')
        plt.xlabel('h [%s]' %lunit)
        plt.ylabel('$\\gamma(h)$')
        if model != 'nomodel':
            sub2.set_title('%s \n %s' %(headtitle, subtitles[ind]))
        else:
            sub2.set_title('%s' %(headtitle))
        plt.axis('auto')
        plt.grid('on')
        plt.legend(loc='best')
        # plt.tight_layout(pad=1, h_pad=0, w_pad=0)

        # plot samples per lag
        # fig3 = plt.figure('samples per lag')
        fig3 = plt.figure(3)
        sub3 = fig3.add_subplot(111)
        for i in range(di.size):
            sub3.plot(h[i],c[i],'o--', c=cm.jet(256/di.size*i),
                      label='%i'%di[i])
        plt.xlabel('h [%s]' %lunit)
        plt.ylabel('Samples per lag')
        plt.title('%s' %(headtitle))
        plt.axis('auto')
        plt.grid('on')
        plt.legend(loc='best')
        # plt.tight_layout(pad=1, h_pad=0, w_pad=0)

        # plot of directions
        # fig4 = plt.figure('directions', figsize=(6,6))
        fig4 = plt.figure(4, figsize=(6,6))
        ax = fig4.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='black')
        if stype == 'directional':
            theta = np.deg2rad(np.append((di-td),(di-td)+180))
        else:
            theta = np.deg2rad(di-td)
        radii = np.ones_like(theta)
        radii[::2] += 0.1
        width = np.deg2rad(np.ones_like(theta)*td*2)
        xlab =[]
        for i in range(19):
            #xlab.append(r'$\sf{%i\degree}$' %(i*10))
            xlab.append('$%i\\degree$' %(i*10))
        for i in range(17):
            #xlab.append(r'$\sf{%i\degree}$' %(-170+i*10))
            xlab.append('$%i\\degree$' %(-170+i*10))
        bars = ax.bar(theta, radii, width=width, bottom=0.0)

        for r,bar in zip(radii, bars):
            bar.set_facecolor('#00ff04')
            bar.set_edgecolor('yellow')
            bar.set_alpha(0.5)
            ax.set_yticklabels(radii, alpha=0)
            ax.xaxis.set_major_locator(plt.LinearLocator(numticks=37,
                                                         presets=None))
            ax.set_xticklabels(xlab, fontsize=15)
            plt.grid(color='yellow')
        ax.set_title(headtitle, fontsize=15)
        plt.show()

    # output
    if model != 'nomodel':
        return nugget, sill, irange, cork, h, g, c, func[ind], popt[ind,:,:]
    else:
        return h, g, c

#---------------------------------------
# sub functions
#---------------------------------------
# function to compute all distances and angles between vectors x and y
def distang(x,y):
    n = x.size
    t = np.zeros(n*(n-1)//2)
    r = np.zeros(n*(n-1)//2)
    k = 0
    for o in range(n-1):
        for p in range(o+1,n):
            dx = x[p]-x[o]
            dy = y[p]-y[o]
            t[k] = np.arctan2(dy,dx)       # angle (theta)
            r[k] = np.sqrt(dx**2+dy**2)    # distance (ray)
            k += 1
    xr = max(r)
    return r, t, xr, n

#---------------------------------------
# function to compute semivariogram
def semivario(r, t, xr, n, z, nL, di, td, stype='omnidirectional', negscat=0.):
    L = xr/nL                           # lag size
    g = np.zeros(nL)
    c = np.zeros(nL)
    a = np.deg2rad(di)                     # angle = direction in radian
    ta = np.deg2rad(td)                    # tolerance of angle in radians
    b = np.deg2rad(0.01)                   # buffer for numerical problems
    for s in range(nL):
        k = 0
        q = 0
        g[s] = 0
        c[s] = 0
        for o in range(n-1):
            for p in range(o+1,n):
                if  (s*L < abs(r[k]) < (s+1)*L):
                    if stype == 'omnidirectional':
                        g[s] += (z[p]-z[o])**2
                        q += 1
                    if stype == 'directional':
                        if a<0:
                            m=np.deg2rad(di+180)
                        else:
                            m=np.deg2rad(di-180)

                        if a+ta+b > np.deg2rad(180):
                            if (((a-ta-b)<t[k]<a) or (a<t[k]<np.deg2rad(180)+b)
                                 or (-np.deg2rad(180)-b<t[k]<(-np.deg2rad(360)+
                                                              (a+ta+b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        elif a-ta-b < -np.deg2rad(180):
                            if (((a+ta+b)>t[k]>a) or (a>t[k]>-np.deg2rad(180)-b)
                                 or (np.deg2rad(180)+b>t[k]>(np.deg2rad(360)+
                                                             (a-ta-b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        else:
                            if (a-ta-b)<t[k]<(a+ta+b):
                                g[s] += (z[p]-z[o])**2
                                q += 1

                        if m+ta+b > np.deg2rad(180):
                            if (((m-ta-b)<t[k]<m) or (m<t[k]<np.deg2rad(180)+b)
                                 or (-np.deg2rad(180)-b<t[k]<(-np.deg2rad(360)+
                                                              (m+ta+b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        elif m-ta-b < -np.deg2rad(180):
                            if (((m+ta+b)>t[k]>m) or (m>t[k]>-np.deg2rad(180)-b)
                                 or (np.deg2rad(180)+b>t[k]>(np.deg2rad(360)+
                                                             (m-ta-b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        else:
                            if (m-ta-b)<t[k]<(m+ta+b):
                                g[s] += (z[p]-z[o])**2
                                q += 1

                    if stype == 'directional+orientational':
                        if a+ta+b > np.deg2rad(180):
                            if (((a-ta-b)<t[k]<a) or (a<t[k]<np.deg2rad(180)+b)
                                 or (-np.deg2rad(180)-b<t[k]<(-np.deg2rad(360)+
                                                              (a+ta+b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        elif a-ta-b < -np.deg2rad(180):
                            if (((a+ta+b)>t[k]>a) or (a>t[k]>-np.deg2rad(180)-b)
                                 or (np.deg2rad(180)+b>t[k]>(np.deg2rad(360)+
                                                             (a-ta-b)))):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                        else:
                            if (a-ta-b)<t[k]<(a+ta+b):
                                g[s] += (z[p]-z[o])**2
                                q += 1
                k += 1
        g[s] /= (q*2) if q>0 else np.NAN
        c[s] = q
    h = np.array(np.arange(nL))*L+L/2
    h = np.delete(h,np.where(np.isnan(g)))       # ranges
    c = np.delete(c,np.where(np.isnan(g)))       # number of samples per range
    g = np.delete(g,np.where(np.isnan(g)))       # variance per range

    if (negscat > 0.) and (negscat < 1.):
        h = h[0:int(np.size(h)*negscat)]
        c = c[0:int(np.size(c)*negscat)]
        g = g[0:int(np.size(g)*negscat)]

    return h,g,c

#---------------------------------------
# exponential semivariogram model
def expvar(h,p):
    # return p[0] + p[1]*(1.-np.exp(-abs(h)/p[2]))
    return p[0] + p[1]*(1.-np.exp(-abs(h)*p[2]))

# range of exponential model at a given sill
def exprange(sill,p):
    # return -np.log((-sill*0.95+p[:,0])/p[:,1]+1)*p[:,2]
    return -np.log((-sill*0.95+p[:,0])/p[:,1]+1)/p[:,2]

#---------------------------------------
# spherical semivariogram model
def sphvar(h,p):
    return np.where((0<=np.abs(h)) & (np.abs(h)<=p[2]) & (p[2]!=0.),
                   p[0]+p[1]*((3./2.)*(abs(h)/p[2])-(1./2.)*(abs(h)/p[2])**3),
                   p[0]+p[1])
    # out1 = p[0]+p[1]
    # out2 = out1
    # if p[2] != 0.: out2 = p[0]+p[1]*(1.5*(abs(h)/p[2]) - 0.5*(abs(h)/p[2])**3
    # return np.where((np.abs(h)>=0.) & (np.abs(h)<=p[2]), out1, out2)

# range of spherical model at a given sill
def sphrange(sill,p):
    return p[:,2]

#---------------------------------------
# gaussian semivariogram model
def gauvar(h,p):
    # return p[0] + p[1]*(1. - np.exp(-h**2/p[2]**2))
    return p[0] + p[1]*(1. - np.exp(-h**2 * p[2]**2))

# range of gaussian model at a given sill
def gaurange(sill,p):
    # return np.sqrt( -(p[:,2])**2 * np.log((-sill*0.95+p[:,0])/p[:,1]+1))
    return np.sqrt( -(1./p[:,2])**2 * np.log((-sill*0.95+p[:,0])/p[:,1]+1.))

#---------------------------------------
# objective function
def absdif(p,x,y,func):
        return np.sum(np.abs(y-func(x, p)))


# DOCTEST:
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # easting
    # x = np.array([557509.27,557518.11,557526.95,557535.79,557544.63,
    #               557553.47,557544.63,557535.79,557526.95,557518.11,
    #               557526.95,557535.79,557544.63,557553.47,557562.31,
    #               557571.15,557562.31,557553.47,557544.63,557535.79,
    #               557544.63,557553.47,557562.31,557571.15,557579.99,
    #               557597.66,557579.99,557562.31,557544.63,557526.7,
    #               557509.27,557491.60,557473.92,557491.6,557509.27,
    #               557526.95,557544.63,557562.31,557579.99,557597.66,
    #               557615.34,557650.7,557686.05,557756.76,557827.47,
    #               557898.18,557827.47,557756.76,557686.05,557615.34,
    #               557650.70,557686.07,557721.41,557756.76,557721.41,
    #               557686.05,557650.70,557686.05,557650.7,557579.99,
    #               557615.34,557509.27,557579.99,557544.63,557509.27,
    #               557473.92,557438.56,557403.21,557438.56,557473.92,
    #               557509.27,557544.63,557579.99,557615.34,557615.34,
    #               557579.99,557544.63,557509.27,557473.92,557438.56,
    #               557403.21,557367.85,557332.5,557367.85,557403.21,
    #               557438.56,557473.92,557332.50,557261.79,557191.08,
    #               557261.79,557332.50,557403.21,557473.92,557544.63,
    #               557615.34])

    # # northing
    # y = np.array([4332422.55,4332413.71,4332404.87,4332396.03,
    #               4332387.19,4332396.03,4332404.87,4332413.71,
    #               4332422.55,4332431.39,4332440.23,4332431.39,
    #               4332422.55,4332413.71,4332404.87,4332413.71,
    #               4332422.55,4332431.39,4332440.23,4332449.07,
    #               4332457.91,4332449.07,4332440.23,4332431.39,
    #               4332422.55,4332404.87,4332387.19,4332369.52,
    #               4332351.84,4332369.52,4332387.19,4332404.87,
    #               4332422.55,4332440.23,4332457.91,4332475.58,
    #               4332493.26,4332475.58,4332457.91,4332440.23,
    #               4332422.55,4332316.48,4332210.42,4332281.13,
    #               4332351.84,4332422.55,4332493.26,4332563.97,
    #               4332634.68,4332563.97,4332528.62,4332493.25,
    #               4332457.91,4332422.55,4332387.19,4332351.84,
    #               4332387.19,4332422.55,4332457.91,4332599.33,
    #               4332493.26,4332599.33,4332528.62,4332563.97,
    #               4332528.62,4332493.26,4332457.91,4332422.55,
    #               4332387.19,4332351.84,4332316.48,4332281.13,
    #               4332316.48,4332351.84,4332281.13,4332245.77,
    #               4332210.42,4332245.77,4332281.13,4332316.48,
    #               4332351.84,4332387.19,4332422.55,4332457.91,
    #               4332493.26,4332528.62,4332563.97,4332563.97,
    #               4332493.26,4332422.55,4332351.84,4332281.13,
    #               4332210.42,4332139.71,4332069.00,4332139.71])
    # # value
    # v = np.array([9.94691161e-01,7.94158417e-02,0.00000000e+00,
    #               1.75837990e+00,0.00000000e+00,0.00000000e+00,
    #               0.00000000e+00,0.00000000e+00,1.02915310e+00,
    #               2.69597379e+00,2.14552427e+00,2.18417112e+00,
    #               0.00000000e+00,8.96101277e-01,1.14034753e+00,
    #               3.46398689e-01,3.01418491e-01,0.00000000e+00,
    #               1.17920343e+00,1.09682206e+00,4.79485665e-01,
    #               0.00000000e+00,1.83183398e+00,0.00000000e+00,
    #               0.00000000e+00,0.00000000e+00,9.86233407e-02,
    #               7.68290376e-02,2.63911513e-01,0.00000000e+00,
    #               2.10013460e+00,0.00000000e+00,2.47535521e+00,
    #               1.47047869e+00,8.00371532e-01,2.39448347e+00,
    #               0.00000000e+00,2.26426861e+00,0.00000000e+00,
    #               0.00000000e+00,1.13769438e+00,1.01969271e+00,
    #               2.26036007e+00,2.38991410e+00,1.82558084e-03,
    #               0.00000000e+00,0.00000000e+00,2.52583544e+00,
    #               6.35195403e-01,2.43778382e+00,0.00000000e+00,
    #               2.47738704e+00,8.83280548e-01,2.42328547e+00,
    #               0.00000000e+00,2.41534081e+00,2.45629467e+00,
    #               0.00000000e+00,2.50770630e+00,1.30382267e+00,
    #               2.06891940e+00,9.17384801e-02,0.00000000e+00,
    #               1.10185544e-01,2.53460688e+00,2.15217780e+00,
    #               1.16908154e+00,1.70072787e-01,1.60603658e-01,
    #               2.15438377e+00,2.32464926e+00,3.26255002e-01,
    #               0.00000000e+00,1.48404530e+00,2.10439439e+00,
    #               0.00000000e+00,0.00000000e+00,0.00000000e+00,
    #               2.34663663e-01,1.46993948e+00,2.67691613e+00,
    #               2.13262460e-02,1.01551520e+00,1.10878523e+00,
    #               1.80374874e+00,1.85571813e+00,2.93929948e+00,
    #               4.43192829e-01,2.55962879e+00,0.00000000e+00,
    #               1.46545683e+00,1.75659977e+00,0.00000000e+00,
    #               2.37093751e+00,0.00000000e+00,0.00000000e+00])
    # # omnidirectional semivariogram with exponential model and fifty lags
    # td = 180
    # di = [0]
    # nL = 50
    # h, g, c = semivariogram(x, y, v, nL, di, td, stype='omnidirectional',
    #                         negscat=False, model='nomodel', graph=False, lunit='m',
    #                         p0=(0.5,0.5,100.), runtimediag=False)
    # print np.round(g,3)
    # # [[ 0.55   0.776  0.619  0.849  0.991  1.033  1.067  1.106  1.079  1.002
    # #    1.119  1.06   0.935  0.789  1.117  1.125  1.063  1.015  1.041  0.911
    # #    1.077  1.116  1.202  1.195  0.935  0.918  0.946  0.596  1.035  1.153
    # #    1.205  1.635  1.384  0.806  0.006  1.409  0.748  0.31   0.977  1.95
    # #    1.748  1.087]]
    # nugget, sill, irange, cork, h, g, c, semi_mod, semi_popt = semivariogram(
    #     x, y, v, nL, di, td, stype='omnidirectional',
    #     negscat=False, model='exponential', graph=False, lunit='m',
    #     p0=(0.5,0.5,1./100.), runtimediag=False)
    # print round(nugget, 2)
    # # 0.41
    # print round(sill, 2)
    # # 1.06
    # print round(irange, 0)
    # # 73.0

    # # directional semivariogram with spherical model and fifty lags
    # td = 45
    # di = [0, 90]
    # nL = 50
    # nugget, sill, irange, cork, h, g, c, semi_mod, semi_popt = semivariogram(
    #     x, y, v, nL, di, td, stype='directional',
    #     negscat=False, model='spherical', graph=False, lunit='m',
    #     p0=(0.5,0.5,100), runtimediag=False)
    # print np.round(nugget, 2)
    # # [ 0.49  0.5 ]
    # print np.round(sill, 2)
    # # [ 1.08  1.07]

    # # directional+orientational semivariogram with gaussian model and
    # # fifty lags
    # td = 30
    # di = [0, 45, 90, 135, 180, -45, -90, -135]
    # nL = 50
    # nugget, sill, irange, cork, h, g, c, semi_mod, semi_popt = semivariogram(
    #     x, y, v, nL, di, td, stype='directional+orientational',
    #     negscat=False, model='spherical', graph=True, lunit='m',
    #     p0=(0.5,0.5,100.), runtimediag=False)
    # print np.round(nugget, 2)
    # # [ 0.25  0.4   0.58  0.49  0.32  0.61  0.49  0.42]
    # print np.round(sill, 2)
    # # [ 0.98  1.12  1.15  1.18  1.01  1.11  0.99  1.09]

