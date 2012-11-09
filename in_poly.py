#!/usr/bin/env python
import numpy as np

def in_poly(P, coord_x, coord_y):

    '''
        This routine is recoded from the UFZ Fortran library. 
        Copyright: Juliane Mai, 2012.
        The original version of the source code (pnpoly) was implemented by 
        W. Randolph Franklin. It was insufficiently assigning vertex/edge points.
        
        PURPOSE:
        Determines whether a 2D point falls in a polygon or is located outside or 
        lies on a vertex or an edge of the polygon. The polygon can be convex or
        not.
        
        DEFINITION:
        def in_poly(P, coord_x, coord_y):
        
        INPUT:
        P         : list, x and y coordinates of the point in question in the
                    form [x,y]
        coord_x   : np.array, x coordinates of the polygon
        coord_y   : np.array, y coordinates of the polygon
        
        RESTRICTIONS:
        The method is only applicable for 2D polygons and points.
        
        OUTPUT:
        erg       : integer, 1 = point inside polygon
                            -1 = point outside polygon
                             0 = point on vertex/edge
                             
        EXAMPLES:
        >>> coord_x = np.array([2.,2.,5.,7.,5.])
        >>> coord_y = np.array([1.,4.,6.,3.,1.])
        
        # point inside polygon
        >>> P = [4.,3.]
        >>> print in_poly(P, coord_x, coord_y)
        1
        
        # point outside polygon
        >>> P = [8.,6.]
        >>> print in_poly(P, coord_x, coord_y)
        -1
        
        # point on edge/vertex of polygon
        >>> P = [2.,2.]
        >>> print in_poly(P, coord_x, coord_y)
        0
        
        LICENSE:
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
    
        HISTORY:
        Written, Arndt Piayda, Nov 2012

    '''
    ########################################################################### 
    # ironing :-)
    coord_x, coord_y = coord_x.flatten(), coord_y.flatten()
    
    # test input sizes
    if np.size(coord_x)!=np.size(coord_y):
        raise ValueError('in_poly: coord_x and coord_y must have same size')
    
    # make working arrays
    N = np.size(coord_x)
    X = np.empty_like(coord_x)
    Y = np.empty_like(coord_y)

    # result is outside as long as no other test works
    erg = -1

    for I in xrange(N):
        # Edge test
        X[I] = coord_x[I] - P[0]
        Y[I] = coord_y[I] - P[1]
        
        if (X[I] == 0.) and (Y[I] == 0.):
            erg = 0

    for I in xrange(N):
        #vertical Vertex test
        J = 0 + (I+1)%N 
        if (coord_x[I] == coord_x[J]) and (coord_x[I] == P[0]):
            ly = (P[1]-coord_y[J]) / (coord_y[I]-coord_y[J])
            if (ly >= 0.) and (ly <= 1.):
                erg=0
        
        # horizontal Vertex test
        if (coord_y[I] == coord_y[J]) and (coord_y[I] == P[1]):
            lx = (P[0]-coord_x[J]) / (coord_x[I]-coord_x[J])
            if (lx >= 0.) and (lx <= 1.):
                erg=0
         
        # Inside test
        MX = X[I] >= 0. 
        NX = X[J] >= 0. 
        MY = Y[I] >= 0. 
        NY = Y[J] >= 0. 

        test1 = not((MY or NY) and (MX or NX)) or (MX and NX)
        test2 = not(MY and NY and (MX or NX) and not(MX and NX))

        if (not test1):
            if test2:
                if ((Y[I]*X[J] - X[I]*Y[J]) / (X[J] - X[I]) < 0.):
                    continue
                else:
                    if ((Y[I]*X[J] - X[I]*Y[J]) / (X[J]-X[I]) > 0.):
                        erg=-erg
                        continue
                    else:
                        erg=0
            else:
                erg=-erg
    return erg

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    