#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import pdist, squareform

def kriging(x, y, z, semi_mod, semi_popt, xnew=None, ynew=None, plot=False,
            masked=False, silent=True, eop=None, block=False):
    """
        Kriging a surface from a set of 2D points with a given semivariogram
        model and associated optimized parameters. Plot the surface and the
        corresponding kriging variance. The coordinates and values of the
        surface can be masked outside the convex hull of the given input points.
        Optional extraction of kriged values at distinct points within the
        surface is possible.


        Definition
        ----------
        def kriging(x,y,z,semi_mod,semi_popt,xnew=None,ynew=None,plot=False,
                    masked=False,silent=True,eop=None):


        Input
        -----
        x          array, x coordinates
        y          array, y coordinates
        z          array, values
        semi_mod   function, semivariogram model (e.g. output from the UFZ
                   semivariogram routine)
        semi_popt  array, parameters of the semivariogram model (e.g. output
                   from the UFZ semivariogram routine)
        xnew       array (n), x coordinates of the desired surface, they will be
                   used to generate a 2D mesh for the surface. If left None,
                   values will be kriged only for the points given in eop.
        ynew       array (m), y coordinates of the desired surface, they will be
                   used to generate a 2D mesh for the surface. If left None,
                   values will be kriged only for the points given in eop.
        eop        array (k,2), x and y coordinates of distinct points where
                   a kriged value is desired.


        Optional Input
        --------------
        plot       bool, plots will be generated if True, otherwise not.
        masked     bool, if True, the output arrays will be np.ma.masked_arrays
                   where coordinates and values outside of the convex hull of
                   the input data are masked. In the generated plots these
                   values will be hidden. If False, the output arrays will be
                   np.arrays and all values within the kriging rectangle are
                   visible in the plots.
        silent     bool, if True, no runtime diagnostics are printed to the
                   console.
        block      bool, if True, calculate block kriging
                   Note that kringing on a rectangular surface xnew,ynew with possible masking
                   and calculating the mean afterwards is almost always much faster,
                   except for very fine xnew,ynew grids.


        Output
        ------
        if eop is None:
        xnew     2D array (n,m), x coordinates of the surface grid
        ynew     2D array (n,m), y coordinates of the surface grid
        znew     2D array (n,m), values of the surface grid
        varnew   2D array (n,m), kriging variance of the surface grid

        if xnew is None and not block:
        eopz     array (k), kriged values at the desired distinct points of eop
        eopvar   array (k), kriging variance at the desired distinct points of
                 eop

        if block:
        bave     average over convex_hull of x,y or xnew,ynew.
        bvar     error of bave; not implemented, using dummy formular
       
        otherwise:
        xnew     2D array (n,m), x coordinates of the surface grid
        ynew     2D array (n,m), y coordinates of the surface grid
        znew     2D array (n,m), values of the surface grid
        varnew   2D array (n,m), kriging variance of the surface grid
        eopz     array (k), kriged values at the desired distinct points of eop
        eopvar   array (k), kriging variance at the desired distinct points of

        graphs:
        kriging surface, shows the kriged surface
        kriging variance, shows the kriging variance


        References
        ----------
        This routine is recoded and extended from a matlab script by Juliane
        Mai.


        Examples
        --------
        # provide you some sample data:
        >>> x = np.array([652225.,652175.,652205.,652235.,652265.,652165.,\
                          652195.,652225.,652255.,652285.,652175.,652205.,\
                          652235.,652265.,652175.,652205.,652235.,652265.,\
                          652195.,652225.,652255.,652285.,652235.,652265.,\
                          652225.,652255.,652285.,652195.,652200.,652200.,\
                          652240.,652230.,652260.,652260.,652265.])
        >>> y = np.array([5772960.,5772970.,5772970.,5772970.,5772970.,\
                          5772980.,5772980.,5772980.,5772980.,5772980.,\
                          5772990.,5772990.,5772990.,5772990.,5773000.,\
                          5773000.,5773000.,5773000.,5773010.,5773010.,\
                          5773010.,5773010.,5773020.,5773020.,5773030.,\
                          5773030.,5773030.,5772985.,5772990.,5772995.,\
                          5773015.,5773025.,5772985.,5772990.,5772995.])
        >>> z = np.array([2.16512767,4.97776467,4.2279204 ,0.        ,\
                          8.25658422,0.01238773,5.05858306,8.33503939,\
                          7.53470443,7.15304826,9.45150218,8.79359049,\
                          0.0536634 ,0.42101194,0.22721601,1.1458486 ,\
                          6.79183025,2.50622739,3.76725118,3.97934707,\
                          0.        ,0.24743279,1.4627512 ,0.38430722,\
                          5.30171261,0.        ,3.17667353,3.80908144,\
                          7.12445478,4.83891708,6.10898131,2.93801857,\
                          2.56170107,2.54503559,1.72767934])

        # make semivariogram
        >>> from semivariogram import semivariogram
        >>> nL = 40
        >>> di = [0]
        >>> td = 180
        >>> nugget,sill,range,vark,h,g,c,semi_mod,semi_popt=\
            semivariogram(x,y,z,nL,di,td,stype='omnidirectional',negscat=0.5,\
                          model='exponential',graph=False,lunit='m',\
                          p0=(0.,20.,1./8.),runtimediag=False)

        # x and y coordinates for the surface
        >>> xnew = np.arange(np.min(x),np.max(x),5.)
        >>> ynew = np.arange(np.min(y),np.max(y),5.)

        # krig the surface
        >>> xnew, ynew, znew, varnew = kriging(x,y,z,semi_mod,semi_popt,\
                                               xnew=xnew,ynew=ynew,silent=True,\
                                               plot=False,masked=False,eop=None)
        >>> print np.round(znew[0],3)
        [ 3.576  3.758  3.912  3.937  3.884  3.83   3.792  3.759  3.71   3.613
          3.407  2.981  2.165  2.366  2.458  2.797  3.304  3.817  4.298  4.717
          4.918  4.77   4.478  4.238]
        >>> print np.round(np.mean(znew),2)
        3.69

        # block krig the surface
        >>> bave, bvar = kriging(x, y, z, semi_mod, semi_popt,xnew=xnew, ynew=ynew,\
                                 silent=True,plot=False,masked=False,eop=None,block=True)
        >>> print np.round(bave,3)
        3.659

        # krig only at points of interest
        >>> poi = np.array([[652209.16,5772986.26],\
                            [652281.10,5773014.27],\
                            [652202.39,5772997.96],\
                            [652264.51,5772992.49],\
                            [652274.81,5772961.62],\
                            [652204.93,5772992.82],\
                            [652232.38,5773021.34],\
                            [652278.25,5773019.58],\
                            [652199.17,5773004.12],\
                            [652276.71,5773006.25]])
        >>> eopz, eopvar = kriging(x,y,z,semi_mod,semi_popt,xnew=None,\
                                   ynew=None,plot=False,masked=False,\
                                   silent=True,eop=poi)
        >>> print np.round(eopz,3)
        [ 6.409  1.677  3.168  1.262  4.636  6.534  2.244  2.255  2.996  2.111]

        # krig both, whole surface and on points of interest
        >>> xnew = np.arange(np.min(x),np.max(x),5.)
        >>> ynew = np.arange(np.min(y),np.max(y),5.)
        >>> xnew, ynew, znew, varnew, eopz, eopvar = kriging(x,y,z,semi_mod,\
                                                     semi_popt,xnew=xnew,\
                                                     ynew=ynew,plot=False,\
                                                     masked=False,silent=True,\
                                                     eop=poi)
        >>> print np.round(znew[0],3)
        [ 3.576  3.758  3.912  3.937  3.884  3.83   3.792  3.759  3.71   3.613
          3.407  2.981  2.165  2.366  2.458  2.797  3.304  3.817  4.298  4.717
          4.918  4.77   4.478  4.238]
        >>> print np.round(eopz,3)
        [ 6.409  1.677  3.168  1.262  4.636  6.534  2.244  2.255  2.996  2.111]


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public License as
        published by the Free Software Foundation, either version 3 of the License,
        or (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not,
        see <http://www.gnu.org/licenses/>.

        Copyright 2009-2012 Matthias Cuntz


        History
        -------
        Written, Arndt Piayda, Nov 2012
        Modified, Arndt Piayda, Dec 2012 - documentation change
        Modified, Matthias Cuntz & Juliane May, Feb 2013 - block
    """

    if not silent:
        import time
        print 'KRIG: prepare data...'
    # ironing :-)
    x, y, z   = x.flatten(), y.flatten(), z.flatten()
    semi_popt = semi_popt.flatten()
    if xnew is not None:
        xnew, ynew = xnew.flatten(), ynew.flatten()
    if eop is not None:
        eopx, eopy = eop[:,0].flatten(), eop[:,1].flatten()

    ###########################################################################
    # orignal x + y data
    # reshape and calculate lags and gammas
    if np.size(x)!=np.size(y):
            raise ValueError('kriging: x and y must have same dimensions')
    if np.size(x)!=np.size(z):
            raise ValueError('kriging: x and z must have same dimensions')

    xy  = np.vstack((x,y)).transpose()
    lag = squareform(pdist(xy, 'euclidean'))
    gamma = semi_mod(lag,semi_popt)

    # make A and append row and column of one's
    A = np.vstack((gamma,np.ones(np.shape(gamma)[1])))
    A = np.hstack((A,np.ones(np.shape(A)[0]).reshape(-1,1)))
    A[-1,-1] = 0.
    invA = np.linalg.inv(A)

    #######################################################################
    # calculate convex hull to hide outer areas
    if masked:
        from convex_hull import convex_hull
        from in_poly import in_poly
        if not silent:
            print 'KRIG: calculate hull...'
            start = time.time()
        hull_points = convex_hull(np.vstack((x,y)), graphic=False,
                                  smidgen=0.0075)
        if not silent:
            stop = time.time()
            print 'KRIG: calculating hull took %0.3f sec' %(stop-start)

    ###########################################################################
    # krig on grid
    if (xnew is not None) and (not block):
        if not silent:
            print 'KRIG: prepare mesh...'
        # make 2D mesh grid
        xnew, ynew = np.meshgrid(xnew, ynew)
        xnew_v = xnew.flatten()
        ynew_v = ynew.flatten()
        length = np.size(xnew_v)

        #######################################################################
        # calculate every znew of xnew and ynew
        if not silent:
            print 'KRIG: kriging...'
            start = time.time()

        znew = np.empty_like(xnew_v)
        varnew = np.empty_like(xnew_v)
        #lamnew = np.empty((np.size(xnew_v), np.size(x)))
        if masked:
            mask = np.empty_like(xnew_v, dtype=int)

        for grid in xrange(length):
            # make B
            b = np.sqrt((x-xnew_v[grid])**2 + (y-ynew_v[grid])**2)
            B = semi_mod(b,semi_popt)
            B = np.append(B, 1.)

            # calculate lambda
            lmd = np.dot(invA, B)

            # shorten it
            lmd = lmd[:-1]
            B   = B[:-1]

            znew[grid]   = np.dot(z,lmd)
            varnew[grid] = np.dot(lmd.transpose(), B)
            #lamnew[grid,:] = lmd

            ###################################################################
            # calculate convex hull to hide outer areas
            if masked:
                mask[grid] = in_poly([xnew_v[grid], ynew_v[grid]],
                                    hull_points[:,0],
                                    hull_points[:,1])

        znew   = znew.reshape(np.shape(xnew))
        varnew = varnew.reshape(np.shape(xnew))
        #lamnew = lamnew.reshape(np.shape(xnew))
        if masked:
            mask = mask.reshape(np.shape(xnew))
            mask = np.where(mask>0, 0, 1)

            xnew = np.ma.masked_array(xnew, mask)
            ynew = np.ma.masked_array(ynew, mask)
            znew = np.ma.masked_array(znew, mask)
            varnew = np.ma.masked_array(varnew, mask)

        if not silent:
            stop = time.time()
            print 'KRIG: kriging took %0.3f sec' %(stop-start)

    #######################################################################
    # krig on extraction points
    if eop is not None:
        length = np.size(eopx)

        #######################################################################
        # calculate every znew of xnew and ynew
        if not silent:
            print 'KRIG: kriging...'
            start = time.time()

        eopz = np.empty_like(eopx)
        eopvar = np.empty_like(eopx)
        #eoplam = np.empty((np.size(eopx), np.size(x)))
        if masked:
            mask = np.empty_like(eopx, dtype=int)

        for grid in xrange(length):
            # make B
            b = np.sqrt((x-eopx[grid])**2 + (y-eopy[grid])**2)
            B = semi_mod(b,semi_popt)
            B = np.append(B, 1.)

            # calculate lambda
            lmd = np.dot(invA, B)

            # shorten it
            lmd = lmd[:-1]
            B   = B[:-1]

            eopz[grid]   = np.dot(z,lmd)
            eopvar[grid] = np.dot(lmd.transpose(), B)
            #eoplam[grid,:] = lmd

            ###################################################################
            # calculate convex hull to hide outer areas
            if masked:
                mask[grid] = in_poly([eopx[grid], eopy[grid]],
                                    hull_points[:,0], hull_points[:,1])

        if masked:
            mask = np.where(mask>0, 0, 1)

            eopx = np.ma.masked_array(eopx, mask)
            eopy = np.ma.masked_array(eopy, mask)
            eopz = np.ma.masked_array(eopz, mask)
            eopvar = np.ma.masked_array(eopvar, mask)

        if not silent:
            stop = time.time()
            print 'KRIG: kriging took %0.3f sec' %(stop-start)


    ###########################################################################
    # block kriging
    if block:
        from scipy.spatial   import Delaunay    # for triangulation
        from scipy.integrate import dblquad     # area integral
        from convex_hull     import convex_hull # convex hull of data points
        from in_poly         import in_poly     # test if point is in polygon
        from area_poly       import area_poly   # the area of a polygon

        if not silent:
            print 'KRIG: block kriging...'
            start = time.time()

        # Functions for the three lines of a triangle
        # Use the global variable tria
        def a(xx): # line between x1 and x2
            if (tria[1,0]-tria[0,0]) == 0.:
                return b(xx)
            else:
                return (tria[1,1]-tria[0,1])/(tria[1,0]-tria[0,0]) * (xx-tria[0,0]) + tria[0,1]
        def b(xx): # line between x3 and x2
            return (tria[2,1]-tria[0,1])/(tria[2,0]-tria[0,0]) * (xx-tria[0,0]) + tria[0,1]
        def c(xx): # line between x2 and x3
            if (tria[2,0]-tria[1,0]) == 0.:
                return b(xx)
            else:
                return (tria[2,1]-tria[1,1])/(tria[2,0]-tria[1,0]) * (xx-tria[1,0]) + tria[1,1]
        # Case distinction if triangle points up or down
        def maxac(xx):
            return np.maximum(a(xx), c(xx))
        def minac(xx):
            return np.minimum(a(xx), c(xx))
        # Function to integrate
        def semi_modyx(yy,xx,shiftyx,obsyx,f,p):
            xx += shiftyx[1]
            yy += shiftyx[0]
            dis = np.sqrt((xx-obsyx[1])**2 + (yy-obsyx[0])**2)
            return f(dis,p) # semivariogram(distance, parameter)

        # x and y points
        if (xnew is not None) and (not masked):
            # Assume rectangle, construct 4 triangles
            xmin = np.amin(xnew)
            ymin = np.amin(ynew)
            xmax = np.amax(xnew)
            ymax = np.amax(ynew)
            xs   = 0.5*(xmax+xmin) # centre of gravity
            ys   = 0.5*(ymax+ymin)
            ntriangles = 4
            tri = np.empty((ntriangles,3,2), dtype=np.float)
            tri[0,0,:] = [xmin,ymin]
            tri[0,1,:] = [xmax,ymin]
            tri[0,2,:] = [xs,ys]
            tri[1,0,:] = [xmin,ymin]
            tri[1,1,:] = [xmin,ymax]
            tri[1,2,:] = [xs,ys]
            tri[2,0,:] = [xmin,ymax]
            tri[2,1,:] = [xmax,ymax]
            tri[2,2,:] = [xs,ys]
            tri[3,0,:] = [xmax,ymax]
            tri[3,1,:] = [xmax,ymin]
            tri[3,2,:] = [xs,ys]
            # Construct convex hull
            cxy = np.empty((ntriangles,2), dtype=np.float)
            cxy[0,:] = [xmin,ymin]
            cxy[1,:] = [xmax,ymin]
            cxy[2,:] = [xmax,ymax]
            cxy[3,:] = [xmin,ymax]
        else:
            # Get convex hull and vertices
            xy  = np.array(zip(x,y))
            d   = Delaunay(xy[:,:])
            cxy = convex_hull(xy.transpose())
            xs  = np.mean(cxy[:,0])
            ys  = np.mean(cxy[:,1])
            # # All triangles
            # tri = xy[d.vertices,:]
            # ntriangles = tri.shape[0]
            # Construct triangles from convex hull and centre of gravity
            ntriangles = d.convex_hull.shape[0]
            tri = np.empty((ntriangles,3,2), dtype=np.float)
            for i in xrange(ntriangles):
                tri[i,0,:] = xy[d.convex_hull[i,0],:]
                tri[i,1,:] = xy[d.convex_hull[i,1],:]
                tri[i,2,:] = [xs,ys]

        # Construct B-vector
        B     = np.empty(x.size+1, dtype=np.float)
        B[-1] = 1.
        for i in xrange(x.size):
            area = 0.
            tvol = 0.
            tvol_err = 0.
            # Calc mean semivariogramm over whole region
            for j in xrange(ntriangles):
                t  = tri[j,:,:2]
                ii = np.argsort(t[:,0])
                tria = t[ii,:]
                xs = np.mean(tria[:,0])
                ys = np.mean(tria[:,1])
                # Select only Delaunay triangles in polygon
                # this is always true now, because we construct the triangles
                # from the convex hull. If we were taking Delaunay output, then this
                # would allow concave polygons; but it needs much longer.
                if in_poly([xs,ys], cxy[:,0], cxy[:,1]) >= 0:
                    area += area_poly(tria[:,0],tria[:,1])
                    xmin  = np.amin(tria[:,0])
                    ymin  = np.amin(tria[:,1])
                    tria[:,0] -= xmin # shift for integral>0
                    tria[:,1] -= ymin
                    if (b(tria[1,0])>tria[1,1]): # triangle points down
                        vol, err = dblquad(semi_modyx, tria[0,0], tria[2,0],
                                           maxac, b,
                                           args=([ymin,xmin],[y[i],x[i]],semi_mod,semi_popt))
                    else:                        # triangle points up
                        vol, err = dblquad(semi_modyx, tria[0,0], tria[2,0],
                                           b, minac,
                                           args=([ymin,xmin],[y[i],x[i]],semi_mod,semi_popt))
                    tvol     += vol
                    tvol_err += err*err
                else:
                    print 'What the ...'
            B[i] = tvol / area

        # calculate lambda
        lmd = np.dot(invA, B)
        # shorten it
        lmd = lmd[:-1]
        B   = B[:-1]
        # averag
        baverage  = np.dot(z,lmd)
        # dummy equation for error - full equation not implemented yet.
        bvariance = np.dot(lmd.transpose(), B)

        if not silent:
            stop = time.time()
            print 'KRIG: block kriging took %0.3f sec' %(stop-start)

    ###########################################################################
    # plotting
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        if not silent:
            print 'KRIG: plotting...'

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

        # plotting contours of kriging
        # fig1 = plt.figure('kriging: surface', figsize=(15,10))
        fig1 = plt.figure(1, figsize=(15,10))
        sub1 = fig1.add_subplot(111, aspect='equal')#, aspect=1)
        if xnew is not None:
            lines    = sub1.contour(xnew,ynew,znew,10,linewidths=1.5,colors='k')
            fillings = sub1.contourf(xnew,ynew,znew,10,cmap=plt.cm.jet)
        if masked:
            hull = sub1.plot(np.hstack((hull_points[:,0],hull_points[0,0])),
                             np.hstack((hull_points[:,1],hull_points[0,1])),
                             color='k')
        if eop is not None:
            scat = sub1.scatter(eopx,eopy,marker='o',c='k',s=40)
        sub1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub1.xaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub1.yaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub1.grid('on')
        sub1.set_title('kriging')
        plt.xlabel('easting')
        plt.ylabel('northing')
        fig1.autofmt_xdate(rotation=45)

        # plt.tight_layout(pad=1, h_pad=0, w_pad=0)
        # cbar need to be below autofm_xdate !!!???
        if xnew is not None:
            cbar = fig1.colorbar(fillings, orientation='vertical', pad=0.05,
                                 shrink = 0.7)
            cbar.set_label('value')

        # plotting contours of variance
        # fig2 = plt.figure('kriging: variance', figsize=(15,10))
        fig2 = plt.figure(2, figsize=(15,10))
        sub2 = fig2.add_subplot(111, aspect='equal')#, aspect=1)
        if xnew is not None:
            lines    = sub2.contour(xnew,ynew,varnew,10,linewidths=1.5,
                                    colors='k')
            fillings = sub2.contourf(xnew,ynew,varnew,10,cmap=plt.cm.jet)
        if masked:
            hull = sub2.plot(np.hstack((hull_points[:,0],hull_points[0,0])),
                             np.hstack((hull_points[:,1],hull_points[0,1])),
                             color='k')
        if eop is not None:
            scat = sub2.scatter(eopx,eopy,marker='o',c='k',s=40)
        sub2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        sub2.xaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub2.yaxis.set_major_formatter(mpl.ticker.
                                       ScalarFormatter(useOffset=False))
        sub2.grid('on')
        sub2.set_title('variance')
        plt.xlabel('easting')
        plt.ylabel('northing')
        fig2.autofmt_xdate(rotation=45)

        # plt.tight_layout(pad=1, h_pad=0, w_pad=0)
        # cbar need to be below autofm_xdate !!!???
        if xnew is not None:
            cbar = fig2.colorbar(fillings, orientation='vertical', pad=0.05,
                                 shrink = 0.7)
            cbar.set_label('value')
        plt.show()

    if eop is None:
        if block:
            return baverage, bvariance
        else:
            return xnew, ynew, znew, varnew
    elif xnew is None:
        return eopz, eopvar
    else:
        return xnew, ynew, znew, varnew, eopz, eopvar

# DOCTEST:
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # x = np.array([652225.,652175.,652205.,652235.,652265.,652165.,\
    #               652195.,652225.,652255.,652285.,652175.,652205.,\
    #               652235.,652265.,652175.,652205.,652235.,652265.,\
    #               652195.,652225.,652255.,652285.,652235.,652265.,\
    #               652225.,652255.,652285.,652195.,652200.,652200.,\
    #               652240.,652230.,652260.,652260.,652265.])
    # y = np.array([5772960.,5772970.,5772970.,5772970.,5772970.,\
    #               5772980.,5772980.,5772980.,5772980.,5772980.,\
    #               5772990.,5772990.,5772990.,5772990.,5773000.,\
    #               5773000.,5773000.,5773000.,5773010.,5773010.,\
    #               5773010.,5773010.,5773020.,5773020.,5773030.,\
    #               5773030.,5773030.,5772985.,5772990.,5772995.,\
    #               5773015.,5773025.,5772985.,5772990.,5772995.])
    # z = np.array([2.16512767,4.97776467,4.2279204 ,0.        ,\
    #               8.25658422,0.01238773,5.05858306,8.33503939,\
    #               7.53470443,7.15304826,9.45150218,8.79359049,\
    #               0.0536634 ,0.42101194,0.22721601,1.1458486 ,\
    #               6.79183025,2.50622739,3.76725118,3.97934707,\
    #               0.        ,0.24743279,1.4627512 ,0.38430722,\
    #               5.30171261,0.        ,3.17667353,3.80908144,\
    #               7.12445478,4.83891708,6.10898131,2.93801857,\
    #               2.56170107,2.54503559,1.72767934])
    # # make semivariogram
    # from semivariogram import semivariogram
    # nL = 40
    # di = [0]
    # td = 180
    # nugget,sill,range,vark,h,g,c,semi_mod,semi_popt=\
    # semivariogram(x,y,z,nL,di,td,stype='omnidirectional',negscat=0.5,\
    #               model='exponential',graph=False,lunit='m',\
    #               p0=(0.,20.,1./8.),runtimediag=False)
    # # x and y coordinates for the surface
    # xnew = np.arange(np.amin(x),np.amax(x),5.)
    # ynew = np.arange(np.amin(y),np.amax(y),5.)
    # xnew, ynew, znew, varnew = kriging(x,y,z,semi_mod,semi_popt,\
    #                                    xnew=xnew,ynew=ynew,silent=True,\
    #                                    plot=False,masked=False,eop=None)
    # print np.round(znew[0],3)
    # # [ 3.576  3.758  3.912  3.937  3.884  3.83   3.792  3.759  3.71   3.613
    # #   3.407  2.981  2.165  2.366  2.458  2.797  3.304  3.817  4.298  4.717
    # #   4.918  4.77   4.478  4.238]
    # print np.round(np.mean(znew),2)
    # # 3.69

    # # block krig the surface
    # bave, bvar = kriging(x, y, z, semi_mod, semi_popt,
    #                      xnew=xnew, ynew=ynew,
    #                      silent=True,plot=False,masked=False,eop=None,block=True)
    # print np.round(bave,3)
    # # 3.659

    # # krig only at points of interest
    # poi = np.array([[652209.16,5772986.26],\
    #                 [652281.10,5773014.27],\
    #                 [652202.39,5772997.96],\
    #                 [652264.51,5772992.49],\
    #                 [652274.81,5772961.62],\
    #                 [652204.93,5772992.82],\
    #                 [652232.38,5773021.34],\
    #                 [652278.25,5773019.58],\
    #                 [652199.17,5773004.12],\
    #                 [652276.71,5773006.25]])
    # eopz, eopvar = kriging(x,y,z,semi_mod,semi_popt,xnew=None,\
    #                        ynew=None,plot=False,masked=False,\
    #                        silent=True,eop=poi)
    # print np.round(eopz,3)
    #     # [ 6.409  1.677  3.168  1.262  4.636  6.534  2.244  2.255  2.996  2.111]

    # # krig both, whole surface and on points of interest
    # xnew = np.arange(np.min(x),np.max(x),5.)
    # ynew = np.arange(np.min(y),np.max(y),5.)
    # xnew, ynew, znew, varnew, eopz, eopvar = kriging(x,y,z,semi_mod,\
    #                                          semi_popt,xnew=xnew,\
    #                                          ynew=ynew,plot=False,\
    #                                          masked=False,silent=True,\
    #                                          eop=poi)
    # print np.round(znew[0],3)
    #     # [ 3.576  3.758  3.912  3.937  3.884  3.83   3.792  3.759  3.71   3.613
    #     #   3.407  2.981  2.165  2.366  2.458  2.797  3.304  3.817  4.298  4.717
    #     #   4.918  4.77   4.478  4.238]
    # print np.round(eopz,3)
    #     # [ 6.409  1.677  3.168  1.262  4.636  6.534  2.244  2.255  2.996  2.111]
