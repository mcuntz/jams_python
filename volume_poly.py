#!/usr/bin/env python
import numpy as np
from scipy.spatial   import Delaunay    # for triangulation
from scipy.integrate import dblquad     # area integral
from convex_hull     import convex_hull # convex hull of data points
from in_poly         import in_poly     # test if point is in polygon
from area_poly       import area_poly   # the area of a polygon

def volume_poly(x, y, func, convexhull=False, area=False, **kwargs):
    """
        Volume of function above polygon. The polygon will be triangulated.
        Then the volume above each triangle is integrated, and alled summed up.
        If convexhull=True then the convex hull of the data points is taken
        as surface of the polygon.

        func must be func(y,x,*args) - note the y,x instead of x,y.
        This means that any unknown keyword in the call of volume_poly will be passed to func.


        Definition
        ----------
        def volume_poly(x, y, func, convexhull=False, area=False, **kwargs):


        Input
        -----
        x          1D array, x coordinates
        y          1D array, y coordinates
        func       function, func(y,x,*args)
                   any unknown keyword will be passed to func


        Optional Input
        --------------
        convexhull   bool, integrate convex hull of polygon instead of polygon itself
        area         bool, return also the area of the polygon (or the convex hull)


        Output
        ------
        volume, integration_error, polygon area if area


        References
        ----------
        Juliane Mai's notes.


        Examples
        --------
        >>> def f1(y,x):\
                return 1
        >>> def f2(y,x,r):\
                return np.sqrt(r*r-x*x-y*y)
        >>> r = 1.0
        >>> n = 1000
        >>> np.random.seed(1)
        >>> x = 2.*np.random.random(n) - 1.
        >>> y = 2.*np.random.random(n) - 1.
        >>> ii = np.where(x*x+y*y <= r*r)[0]
        >>> x  = x[ii]
        >>> y  = y[ii]
        >>> print np.round(np.pi * r*r,3)
        3.142
        >>> v = volume_poly(x, y, f1, convexhull=True)
        >>> print np.round(v[0],3)
        2.996
        >>> v = volume_poly(x, y, f1, convexhull=True, area=True)
        >>> print np.round(v[0],3), np.round(v[2],3)
        2.996 2.996

        >>> print np.round(0.5 * 4./3. * np.pi * r*r*r, 3)
        2.094
        >>> v = volume_poly(x, y, f2, r=r, convexhull=True)
        >>> print np.round(v[0], 3)
        2.072


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

        Copyright 2013 Matthias Cuntz, Juliane Mai


        History
        -------
        Written,  Matthias Cuntz & Juliane Mai, Feb 2013
    """

    # check
    if np.size(x) != np.size(y):
            raise ValueError('volume_poly: x and y must have same dimensions')

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
    def funcyx(yy,xx,shiftyx,*args):
        xx += shiftyx[1]
        yy += shiftyx[0]
        if args[0] == ():
            return func(yy,xx)
        else:
            kwargs = args[0]
            return func(yy,xx,**kwargs)

    # Get convex hull and vertices
    xy  = np.array(zip(x,y))
    d   = Delaunay(xy[:,:])
    cxy = convex_hull(xy.transpose())
    xs  = np.mean(cxy[:,0])
    ys  = np.mean(cxy[:,1])
    if convexhull:
        # Construct triangles from convex hull vertices and centre of gravity
        ntriangles = d.convex_hull.shape[0]
        tri = np.empty((ntriangles,3,2), dtype=np.float)
        for i in xrange(ntriangles):
            tri[i,0,:] = xy[d.convex_hull[i,0],:]
            tri[i,1,:] = xy[d.convex_hull[i,1],:]
            tri[i,2,:] = [xs,ys]
    else:
        # All triangles
        tri = xy[d.vertices,:]
        ntriangles = tri.shape[0]

    flaeche = 0.
    tvol = 0.
    tvol_err = 0.
    # Calc mean semivariogramm over whole region
    for j in xrange(ntriangles):
        t    = tri[j,:,:]
        ii   = np.argsort(t[:,0])
        tria = t[ii,:]
        xs   = np.mean(tria[:,0])
        ys   = np.mean(tria[:,1])
        # Select only Delaunay triangles that are inside the original polygon
        # i.e. exclude triangles in concave part of polygon
        # If convexhull=True then this is always true.
        if in_poly([xs,ys], cxy[:,0], cxy[:,1]) >= 0:
            flaeche   += area_poly(tria[:,0],tria[:,1])
            xmin       = np.amin(tria[:,0])
            ymin       = np.amin(tria[:,1])
            tria[:,0] -= xmin # shift for integral>0
            tria[:,1] -= ymin
            if (b(tria[1,0])>tria[1,1]): # triangle points down
                vol, err = dblquad(funcyx, tria[0,0], tria[2,0],
                                   maxac, b, args=([ymin,xmin],kwargs))
            else:                        # triangle points up
                vol, err = dblquad(funcyx, tria[0,0], tria[2,0],
                                   b, minac, args=([ymin,xmin],kwargs))
            tvol     += vol
            tvol_err += err*err
    tvol_err = np.sqrt(tvol_err)

    if area:
        return tvol, tvol_err, flaeche
    else:
        return tvol, tvol_err


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # def f1(y,x):
    #     return 1

    # def f2(y,x,r):
    #     return np.sqrt(r*r-x*x-y*y)

    # r = 1.0
    # n = 1000
    # np.random.seed(1)
    # x = 2.*np.random.random(n) - 1.
    # y = 2.*np.random.random(n) - 1.
    # ii = np.where(x*x+y*y <= r*r)[0]
    # x  = x[ii]
    # y  = y[ii]
    # print np.round(np.pi * r*r,3)

    # # v = volume_poly(x, y, f1)
    # # print v[0], v[1]
    # v = volume_poly(x, y, f1, convexhull=True)
    # print np.round(v[0],3)
    # v = volume_poly(x, y, f1, convexhull=True, area=True)
    # print np.round(v[0],3), np.round(v[2],3)

    # print np.round(0.5 * 4./3. * np.pi * r*r*r, 3)
    # v = volume_poly(x, y, f2, r=r, convexhull=True)
    # print np.round(v[0], 3)
