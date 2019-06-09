#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
"""
    Intersection of two curves from x,y coordinates.

    Inspired from Matlab code
        http://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections


    Example
    -------
        a, b = 1, 2
        phi  = np.linspace(3, 10, 100)
        x1   = a*phi - b*np.sin(phi)
        y1   = a - b*np.cos(phi)

        x2   = phi
        y2   = np.sin(phi)+2
        x, y = intersection(x1, y1, x2, y2)

        plt.plot(x1, y1, c="r")
        plt.plot(x2, y2, c="g")
        plt.plot(x, y, "*k")
        plt.show()


    Licence
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2017 Sukhbinder Singh

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

    Sukhbinder
    5 April 2017
"""

__all__ = ['intersection']

def _rect_inter_inner(x1, x2):
    n1 = x1.shape[0]-1
    n2 = x2.shape[0]-1
    X1 = np.c_[x1[:-1],x1[1:]]
    X2 = np.c_[x2[:-1],x2[1:]]
    
    S1 = np.tile(X1.min(axis=1),(n2,1)).T
    S2 = np.tile(X2.max(axis=1),(n1,1))
    S3 = np.tile(X1.max(axis=1),(n2,1)).T
    S4 = np.tile(X2.min(axis=1),(n1,1))

    return S1, S2, S3, S4


def _rectangle_intersection_(x1, y1, x2, y2):
    S1, S2, S3, S4 = _rect_inter_inner(x1,x2)
    S5, S6, S7, S8 = _rect_inter_inner(y1,y2)

    C1 = np.less_equal(S1,S2)
    C2 = np.greater_equal(S3,S4)
    C3 = np.less_equal(S5,S6)
    C4 = np.greater_equal(S7,S8)

    ii,jj = np.nonzero(C1 & C2 & C3 & C4)

    return ii, jj


def intersection(x1, y1, x2, y2):
    """
        Intersections of curves given by x,y coordinates.

        Computes the (x,y) locations where two curves intersect.
        The curves can be broken with NaNs or have vertical segments.


        Definition
        ----------
        def intersection(x1, y1, x2, y2):


        Input
        -----
        x1, y1     x,y coordinates of first curve
        x2, y2     x,y coordinates of second curve


        Output
        ------
        x,y = intersection(x1, y1, x2, y2)
        x,y coordinates of all intersection points


        Theory
        ------
        From the Matlab code:
        
        Given two line segments, L1 and L2,

          L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
          L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))

        we can write four equations with four unknowns and then solve them.  The
        four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
        L1 and L2, t1 is the distance from the starting point of L1 to the
        intersection relative to the length of L1 and t2 is the distance from the
        starting point of L2 to the intersection relative to the length of L2.

        So, the four equations are

           (x1(2) - x1(1))*t1 = x0 - x1(1)
           (x2(2) - x2(1))*t2 = x0 - x2(1)
           (y1(2) - y1(1))*t1 = y0 - y1(1)
           (y2(2) - y2(1))*t2 = y0 - y2(1)

        Rearranging and writing in matrix form,

         [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
               0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
          y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
               0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]

        Let's call that A*T = B.  We can solve for T with T = A\B.

        Once we have our solution we just have to look at t1 and t2 to determine
        whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
        line segments cross and we can include (x0,y0) in the output.

        In principle, we have to perform this computation on every pair of line
        segments in the input data.  This can be quite a large number of pairs so
        we will reduce it by doing a simple preliminary check to eliminate line
        segment pairs that could not possibly cross.  The check is to look at the
        smallest enclosing rectangles (with sides parallel to the axes) for each
        line segment pair and see if they overlap.  If they do then we have to
        compute t1 and t2 (via the A\B computation) to see if the line segments
        cross, but if they don't then the line segments cannot cross.  In a
        typical application, this technique will eliminate most of the potential
        line segment pairs.


        Examples
        --------
        >>> a, b = 1., 2.
        >>> phi  = np.linspace(3, 10, 100)
        >>> x1   = a*phi - b*np.sin(phi)
        >>> y1   = a - b*np.cos(phi)
        >>> x2   = phi
        >>> y2   = np.sin(phi) + 2.
        >>> x,y  = intersection(x1, y1, x2, y2)
        >>> print('{:.1f},{:.1f}'.format(x[0], y[0]))
        6.1,1.8
        >>> print('{:.1f},{:.1f}'.format(x[1], y[1]))
        8.4,2.9


        History
        -------
        Written,  Sukhbinder Singh, Jun 2017 - inspired by revision 1 of Matlab code of Douglas M. Schwarz
                  http://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
        Modified, MC, Nov 2018 - ported to JAMS
    """
    ii,jj = _rectangle_intersection_(x1,y1,x2,y2)
    n = len(ii)

    dxy1 = np.diff(np.c_[x1,y1],axis=0)
    dxy2 = np.diff(np.c_[x2,y2],axis=0)

    T  = np.zeros((4,n))
    AA = np.zeros((4,4,n))
    AA[0:2,2,:]  = -1.
    AA[2:4,3,:]  = -1.
    AA[0::2,0,:] = dxy1[ii,:].T
    AA[1::2,1,:] = dxy2[jj,:].T

    BB = np.zeros((4,n))
    BB[0,:] = -x1[ii].ravel()
    BB[1,:] = -x2[jj].ravel()
    BB[2,:] = -y1[ii].ravel()
    BB[3,:] = -y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i] = np.linalg.solve(AA[:,:,i], BB[:,i])
        except:
            T[:,i] = np.NaN

    in_range = (T[0,:]>=0.) & (T[1,:]>=0.) & (T[0,:]<=1.) & (T[1,:]<=1.)

    xy0 = T[2:,in_range]
    xy0 = xy0.T

    return xy0[:,0], xy0[:,1]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # a piece of a prolate cycloid, and am going to find
    # a, b = 1, 2
    # phi = np.linspace(3, 10, 100)
    # x1 = a*phi - b*np.sin(phi)
    # y1 = a - b*np.cos(phi)

    # x2=phi
    # y2=np.sin(phi)+2
    # x,y=intersection(x1,y1,x2,y2)
    # print(x,y)
    # import matplotlib.pyplot as plt
    # plt.plot(x1,y1,c='r')
    # plt.plot(x2,y2,c='g')
    # plt.plot(x,y,'*k')
    # plt.show()
