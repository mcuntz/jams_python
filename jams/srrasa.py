#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def srrasa(xy, strata=5, n=3, plot=False):
    """
        Generates stratified random 2D points within a given rectangular area.


        Definition
        ----------
        def srrasa(xy, strata=5, n=3, plot=False):


        Input
        -----
        xy          list of floats (4), list with the x and y coordinates
                    enclosing the designated rectangle in the form [x1,x2,y1,y2]


        Optional Input
        --------------
        strata      int, number of strata per axis
        n           int, number of random points in each strata
        plot        bool, if True, stratas and points are plotted,
                    otherwise not


        Output
        ------
        rand_xy     ndarray (n,2), x and y coordinates of the stratified random
                    points in the given rectangular.


        Examples
        --------
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)
        >>> # gives within the rectangle of the given coordinates
        >>> # 16 (4**2) stratas with 3 random points in each one.
        >>> rand_xy = srrasa([652219.,652290.,5772970.,5773040.], strata=4, n=3, plot=False)
        >>> from autostring import astr
        >>> print(astr(rand_xy[0:4,0:2],6,pp=True))
        [['6.522264e+05' '5.772975e+06']
         ['6.522318e+05' '5.772973e+06']
         ['6.522190e+05' '5.772972e+06']
         ['6.522401e+05' '5.772979e+06']]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2013 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  AP, Nov 2012
        Modified, MC, Nov 2012 - default plot=False
                  AP, Dec 2012 - documentation change
                  MC, Feb 2013 - docstring
                  MC, Feb 2013 - ported to Python 3
    """

    # calculate strata steps
    sw = (xy[1]-xy[0])/strata
    sh = (xy[3]-xy[2])/strata
    xsteps = np.arange(xy[0],xy[1]+sw,sw)
    ysteps = np.arange(xy[2],xy[3]+sh,sh)

    # make output array
    rand_xy = np.empty((strata**2*n,2))

    # throw random points in each strata
    for j in range(strata):
        for i in range(strata):
            rand_xy[i*n+strata*n*j:(i+1)*n+strata*n*j,0] = (xsteps[i+1] - xsteps[i])*np.random.random(n) + xsteps[i]
            rand_xy[i*n+strata*n*j:(i+1)*n+strata*n*j,1] = (ysteps[j+1] - ysteps[j])*np.random.random(n) + ysteps[j]

    # plot stratas and random points within
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rc('font', size=20)
        mpl.rc('lines', linewidth=2)
        mpl.rc('axes', linewidth=1.5)
        mpl.rc('xtick.major', width=1.5)
        mpl.rc('ytick.major', width=1.5)
        mpl.rcParams['lines.markersize']=6

        fig = plt.figure('stratified random sampling')
        sub = fig.add_subplot(111, aspect='equal')
        sub.set_xlim(xy[0],xy[1])
        sub.set_ylim(xy[2],xy[3])
        for i in range(strata):
            sub.axhline(y=ysteps[i], color=(166/256., 206/256., 227/256.))
            sub.axvline(x=xsteps[i], color=(166/256., 206/256., 227/256.))
        sub.scatter(rand_xy[:,0],rand_xy[:,1],marker='+', s=60,
                    color=( 51/256., 160/256.,  44/256.))
        sub.set_xlabel('X')
        sub.set_ylabel('Y')
        sub.set_title('strata = %i, n = %i' %(strata,n))
        sub.xaxis.set_major_formatter(mpl.ticker.
                                      ScalarFormatter(useOffset=False))
        sub.yaxis.set_major_formatter(mpl.ticker.
                                      ScalarFormatter(useOffset=False))
        fig.autofmt_xdate(rotation=45)
        plt.tight_layout(pad=1, h_pad=0, w_pad=0)
        plt.show()

    return rand_xy



def srrasa_trans(xy,strata=5,n=3,num=3,rl=0.5,silent=True,plot=False):

    """
        Generates stratified random 2D transects within a given rectangular
        area.


        Definition
        ----------
        def srrasa(xy,strata=5,n=3,num=3,rl=0.5,silent=True,plot=False):


        Input
        -----
        xy          list of floats (4), list with the x and y coordinates
                    enclosing the designated rectangle in the form [x1,x2,y1,y2]


        Optional Input
        --------------
        strata      int, number of strata per axis
        n           int, number of random transects in each strata
        num         int, number of points in each transect
        rl          float [0. to 1.], relative length of transect with respect
                    to width of stratum
        silent      bool, if False, runtime diagnostics are printed to the
                    console, otherwise not
        plot        bool, if True, stratas and points are plotted,
                    otherwise not


        Output
        ------
        rand_xy     ndarray (n,2), x and y coordinates of the stratified random
                    transect points in the given rectangular.


        Examples
        --------
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)
        >>> # gives within the rectangle of the given coordinates
        >>> # 16 (4**2) stratas with 3 random transects in each one.
        >>> # Each transect is 0.5*width_of_strata long and contains 5 points logarithmical distributed.
        >>> rand_xy = srrasa_trans([652219.,652290.,5772970.,5773040.], strata=4,
        ...                        n=3, num=5, rl=0.5, silent=True, plot=False)
        >>> from autostring import astr
        >>> print(astr(rand_xy[0:4,0:2],6,pp=True))
        [['6.522264e+05' '5.772983e+06']
         ['6.522276e+05' '5.772983e+06']
         ['6.522292e+05' '5.772983e+06']
         ['6.522315e+05' '5.772983e+06']]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2012-2013 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  AP, Nov 2012
        Modified, AP, Dec 2012 - documentation change
                  MC, Feb 2013 - ported to Python 3
    """

    # calculate strata steps
    sw = (xy[1]-xy[0])/strata
    sh = (xy[3]-xy[2])/strata
    xsteps = np.arange(xy[0],xy[1]+sw,sw)
    ysteps = np.arange(xy[2],xy[3]+sh,sh)
    tl = sw*rl

    # make output array
    rand_xy = np.empty((strata**2*n*num,2))

    o = 0
    for j in range(strata):
        for i in range(strata):
            for k in range(n):

                goon = True
                it = 0
                while goon:
                    # random seed in strata
                    seedx=(xsteps[i+1]-xsteps[i])*np.random.random(1)+xsteps[i]
                    seedy=(ysteps[j+1]-ysteps[j])*np.random.random(1)+ysteps[j]

                    # make logarithmic transect
                    tx   =np.arange(1,num+1)
                    dis  =np.sort(tl-np.log(tx)/np.max(np.log(tx))*tl)
                    seedx=np.repeat(seedx,num)+dis
                    seedy=np.repeat(seedy,num)

                    # random angle in strata [deg]
                    angle = 360 * np.random.random(1)

                    # rotate transect to random angle
                    seedx_trans = (-(seedy-seedy[0])*np.sin(np.deg2rad(angle))+
                                    (seedx-seedx[0])*np.cos(np.deg2rad(angle))+
                                     seedx[0])
                    seedy_trans =  ((seedy-seedy[0])*np.cos(np.deg2rad(angle))+
                                    (seedx-seedx[0])*np.sin(np.deg2rad(angle))+
                                     seedy[0])

                    # test if transect is in strata
                    if (((seedx_trans>xsteps[i]).all()) &
                        ((seedx_trans<xsteps[i+1]).all()) &
                        ((seedy_trans>ysteps[j]).all()) &
                        ((seedy_trans<ysteps[j+1]).all())):
                        goon = False

                    if not silent:
                        print('strata= (', i, ',', j, ')', ' it= ', it)
                    it += 1

                rand_xy[o:o+num,0] = seedx_trans
                rand_xy[o:o+num,1] = seedy_trans
                o += num

    # plot stratas and random transect points within
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rc('font', size=20)
        mpl.rc('lines', linewidth=2)
        mpl.rc('axes', linewidth=1.5)
        mpl.rc('xtick.major', width=1.5)
        mpl.rc('ytick.major', width=1.5)
        mpl.rcParams['lines.markersize']=6

        fig = plt.figure('stratified random transect sampling')
        sub = fig.add_subplot(111, aspect='equal')
        sub.set_xlim(xy[0],xy[1])
        sub.set_ylim(xy[2],xy[3])
        for i in range(strata):
            sub.axhline(y=ysteps[i], color=(166/256., 206/256., 227/256.))
            sub.axvline(x=xsteps[i], color=(166/256., 206/256., 227/256.))
        sub.scatter(rand_xy[:,0],rand_xy[:,1],marker='+', s=60,
                    color=( 51/256., 160/256.,  44/256.))
        sub.set_xlabel('X')
        sub.set_ylabel('Y')
        sub.set_title('strata = %i, n = %i, num = %i' %(strata,n,num))
        sub.xaxis.set_major_formatter(mpl.ticker.
                                      ScalarFormatter(useOffset=False))
        sub.yaxis.set_major_formatter(mpl.ticker.
                                      ScalarFormatter(useOffset=False))
        fig.autofmt_xdate(rotation=45)
        plt.tight_layout(pad=1, h_pad=0, w_pad=0)
        plt.show()

    return rand_xy


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
