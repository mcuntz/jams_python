#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy.spatial    import Delaunay    # for triangulation
from scipy.integrate  import dblquad     # area integral
from jams.convex_hull import convex_hull # convex hull of data points
from jams.in_poly     import in_poly     # test if point is in polygon
from jams.area_poly   import area_poly   # the area of a polygon

def volume_poly(func, x=None, y=None, tri=None, convexhull=False, area=False, allvol=False, **kwargs):
    """
        Volume of function above polygon. The polygon will be triangulated.
        Then the volume above each triangle is integrated, and alled summed up.
        If convexhull=True then the convex hull of the data points is taken
        as surface of the polygon.

        func must be func(y,x,*args) - note the y,x instead of x,y.
        This means that any unknown keyword in the call of volume_poly will be passed to func.


        Definition
        ----------
        def volume_poly(func, x=None, y=None, tri=None, convexhull=False, area=False, **kwargs):


        Input
        -----
        func       function, func(y,x,*args)
                   any unknown keyword will be passed to func


        Optional Input
        --------------
        Either
          x          1D array, x coordinates
          y          1D array, y coordinates
        or
          tri        3D array (ntriangles,3,2), triangle coordinates (supercedes x, y)
                     (ntriangles, 3 corners, xy-coords for each corner)
        convexhull   bool, integrate convex hull of polygon instead of polygon itself (ignored if tri is given)
        area         bool, return also the area of the polygon (or the convex hull)
        allvol       bool, return also all volumes and areas of triangles


        Output
        ------
        volume, integration_error
        polygon area if area
        triangle volumes, areas of triangles if allvol


        References
        ----------
        Juliane Mai's notes.


        Examples
        --------
        >>> def f1(y,x):
        ...     return 1
        >>> def f2(y,x,r):
        ...     return np.sqrt(r*r-x*x-y*y)
        >>> r = 1.0
        >>> n = 1000
        >>> np.random.seed(1)
        >>> x = 2.*np.random.random(n) - 1.
        >>> y = 2.*np.random.random(n) - 1.
        >>> ii = np.where(x*x+y*y <= r*r)[0]
        >>> x  = x[ii]
        >>> y  = y[ii]
        >>> from autostring import astr
        >>> print(astr(np.pi * r*r,3))
        3.142
        >>> v = volume_poly(f1, x, y, convexhull=True)
        >>> print(astr(v[0],3))
        2.996
        >>> v = volume_poly(f1, x, y, convexhull=True, area=True)
        >>> print(astr(v[0],3))
        2.996
        >>> print(astr(v[2],3))
        2.996

        >>> print(astr(0.5 * 4./3. * np.pi * r*r*r, 3))
        2.094
        >>> v = volume_poly(f2, x, y, r=r, convexhull=True)
        >>> print(astr(v[0], 3))
        2.072

        # Triangles outside
        >>> xy  = np.array(list(zip(x,y)))
        >>> d   = Delaunay(xy[:,:])
        >>> cxy = convex_hull(xy.transpose())
        >>> xs  = np.mean(cxy[:,0])
        >>> ys  = np.mean(cxy[:,1])

        # Construct triangles from convex hull vertices and centre of gravity
        >>> ntriangles = d.convex_hull.shape[0]
        >>> tri = np.empty((ntriangles,3,2), dtype=np.float)
        >>> for i in range(ntriangles):
        ...     tri[i,0,:] = xy[d.convex_hull[i,0],:]
        ...     tri[i,1,:] = xy[d.convex_hull[i,1],:]
        ...     tri[i,2,:] = [xs,ys]
        >>> v = volume_poly(f1, tri=tri, convexhull=True)
        >>> print(astr(v[0],3))
        2.996
        >>> v = volume_poly(f1, tri=tri, convexhull=True, area=True)
        >>> print(astr(v[0],3))
        2.996
        >>> print(astr(v[2],3))
        2.996
        >>> v = volume_poly(f2, tri=tri, r=r, convexhull=True)
        >>> print(astr(v[0], 3))
        2.072

        # Return all volumes and areas
        >>> v, e, vv, aa = volume_poly(f2, tri=tri, r=r, convexhull=True, allvol=True)
        >>> print(astr(vv[0], 3))
        0.056


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2013 Matthias Cuntz, Juliane Mai - mc (at) macu (dot) de

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
        Written,  MC & JM, Feb 2013
        Modified, MC, Feb 2013 - tri
                  MC, Feb 2013 - ported to Python 3
    """

    # Functions for the three lines of a triangle
    # Use the global variable tria
    def a(xx): # line between x1 and x2
        if (tria[1,0]-tria[0,0]) == 0.:
            return b(xx)
        else:
            return (tria[1,1]-tria[0,1])/(tria[1,0]-tria[0,0]) * (xx-tria[0,0]) + tria[0,1]
    def b(xx): # line between x1 and x3
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

    if tri is not None:
        trigiven = True
        tri = np.array(tri) # assure
        ntriangles = tri.shape[0]
        if (tri.shape[1] != 3) | (tri.shape[2] != 2):
            raise ValueError('volume_poly: tri must be triangles with dimension (ntriangles, 3, 2)')
    else:
        trigiven = False
        # check
        if np.size(x) != np.size(y):
            raise ValueError('volume_poly: x and y must have same dimensions')
        # Get convex hull and vertices
        xy  = np.array(list(zip(x,y)))
        d   = Delaunay(xy[:,:])
        cxy = convex_hull(xy.transpose())
        xs  = np.mean(cxy[:,0])
        ys  = np.mean(cxy[:,1])
        if convexhull:
            # Construct triangles from convex hull vertices and centre of gravity
            ntriangles = d.convex_hull.shape[0]
            tri = np.empty((ntriangles,3,2), dtype=np.float)
            for i in range(ntriangles):
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
    volumes = np.zeros(ntriangles, dtype=np.float)
    areas   = np.zeros(ntriangles, dtype=np.float)
    # Calc mean semivariogramm over whole region
    thecount = 0
    for j in range(ntriangles):
        t    = tri[j,:,:]
        ii   = np.argsort(t[:,0])
        tria = t[ii,:]
        xs   = np.mean(tria[:,0])
        ys   = np.mean(tria[:,1])
        # Select only Delaunay triangles that are inside the original polygon
        # i.e. exclude triangles in concave part of polygon
        # If convexhull=True then this is always true.
        check_in_poly = False
        if not trigiven:
            check_in_poly = in_poly([xs,ys], cxy[:,0], cxy[:,1]) >= 0
        if trigiven | check_in_poly:
            thecount  += 1
            areas[j]   = area_poly(tria[:,0],tria[:,1])
            flaeche   += areas[j]
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
            volumes[j] = vol
            tvol      += vol
            tvol_err  += err*err
    tvol_err = np.sqrt(tvol_err)

    out  = [tvol]
    out += [tvol_err]
    if area: out += [flaeche]
    if allvol:
        out += [volumes[0:thecount]]
        out += [areas[0:thecount]]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

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
    # # 3.142
    # v = volume_poly(f1, x, y, convexhull=True)
    # print np.round(v[0],3)
    # # 2.996
    # v = volume_poly(f1, x, y, convexhull=True, area=True)
    # print np.round(v[0],3), np.round(v[2],3)
    # # 2.996 2.996

    # print np.round(0.5 * 4./3. * np.pi * r*r*r, 3)
    # # 2.094
    # v = volume_poly(f2, x, y, r=r, convexhull=True)
    # print np.round(v[0], 3)
    # # 2.072

    # # Triangles outside
    # xy  = np.array(zip(x,y))
    # d   = Delaunay(xy[:,:])
    # cxy = convex_hull(xy.transpose())
    # xs  = np.mean(cxy[:,0])
    # ys  = np.mean(cxy[:,1])
    # # Construct triangles from convex hull vertices and centre of gravity
    # ntriangles = d.convex_hull.shape[0]
    # tri = np.empty((ntriangles,3,2), dtype=np.float)
    # for i in range(ntriangles):
    #     tri[i,0,:] = xy[d.convex_hull[i,0],:]
    #     tri[i,1,:] = xy[d.convex_hull[i,1],:]
    #     tri[i,2,:] = [xs,ys]
    # v = volume_poly(f1, tri=tri, convexhull=True)
    # print np.round(v[0],3)
    # # 2.996
    # v = volume_poly(f1, tri=tri, convexhull=True, area=True)
    # print np.round(v[0],3), np.round(v[2],3)
    # # 2.996 2.996
    # v = volume_poly(f2, tri=tri, r=r, convexhull=True)
    # print np.round(v[0], 3)
    # # 2.072
    # v, e, vv, aa = volume_poly(f2, tri=tri, r=r, convexhull=True, allvol=True)
    # print np.round(vv[0], 3)
    # # 0.056

