#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def convex_hull(points, graphic=False, smidgen=0.0075):
    """
        Calculate subset of 2D points that make a convex hull around a set of
        2D points. Recursively eliminates points that lie inside two
        neighbouring points until only convex hull is remaining.


        Definition
        ----------
        def convex_hull(points, graphic=True, smidgen=0.0075)


        Input
        -----
        points       ndarray (2 x m), array of points for which to find hull


        Optional Input
        --------------
        graphic      bool, use pylab to show progress
        smidgen      float, offset for graphic number labels - useful
                     values depend on your data range


        Output
        ------
        hull_points  ndarray (2 x n), convex hull surrounding points


        References
        ----------
        This routine (with additional subroutines) was coded originally by
        Angus McMorland, 2007.
        It is a copy from the scipy cookbook page:
            http://www.scipy.org/Cookbook/Finding_Convex_Hull


        Examples
        --------
        # make some points
        >>> points = np.array([[2,3,2,4,5,5,7,5,5],[1,2,4,3,6,4,3,2,1]])
        >>> hull_xy = convex_hull(points, graphic=False, smidgen=0.075)
        >>> from autostring import astr
        >>> print(astr(hull_xy,pp=True))
        [['5' '1']
         ['7' '3']
         ['5' '6']
         ['2' '4']
         ['2' '1']]

        >>> hull_xy = convex_hull(points, graphic=False, smidgen=0.075)
        >>> print(astr(hull_xy,pp=True))
        [['5' '1']
         ['7' '3']
         ['5' '6']
         ['2' '4']
         ['2' '1']]

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2009-2013 Arndt Piayda

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
        Modified, MC, Feb 2013 - ported to Python 3
    """

    if graphic:
        import pylab as p
        p.clf()
        p.plot(points[0], points[1], 'ro')
    n_pts = points.shape[1]
    assert(n_pts > 5)
    centre = points.mean(1)
    if graphic: p.plot((centre[0],),(centre[1],),'bo')
    angles = np.apply_along_axis(_angle_to_point, 0, points, centre)
    pts_ord = points[:,angles.argsort()]
    if graphic:
        for i in range(n_pts):
            p.text(pts_ord[0,i] + smidgen, pts_ord[1,i] + smidgen, '%d' % i)
    pts = [x[0] for x in zip(pts_ord.transpose())]
    prev_pts = len(pts) + 1
    k = 0
    while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        if graphic: p.gca().patches = []
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if graphic:
                _draw_triangle(centre, pts[i], pts[(i + 1) % n_pts], graphic, facecolor='blue', alpha = 0.2)
                _draw_triangle(centre, pts[(i + 1) % n_pts], pts[(i + 2) % n_pts], graphic, facecolor='green', alpha = 0.2)
                _draw_triangle(centre, pts[i], pts[(i + 2) % n_pts], graphic, facecolor='red', alpha = 0.2)
            if Aij + Ajk < Aik:
                if graphic: p.plot((pts[i + 1][0],),(pts[i + 1][1],),'go')
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1
    if graphic: p.show()
    return np.asarray(pts)

def _angle_to_point(point, centre):
    '''calculate angle in 2-D between points and x axis'''
    delta = point - centre
    res = np.arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += np.pi
    return res

def _draw_triangle(p1, p2, p3, graphic, **kwargs):
    if graphic:
        import pylab as p
    tmp = np.vstack((p1,p2,p3))
    x,y = [x[0] for x in zip(tmp.transpose())]
    p.fill(x,y, **kwargs)

def area_of_triangle(p1, p2, p3):
    '''calculate area of any triangle given co-ordinates of the corners'''
    return np.linalg.norm(np.cross((p2 - p1), (p3 - p1)))/2.

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

