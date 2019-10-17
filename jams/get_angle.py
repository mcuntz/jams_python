#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def get_angle(xy1, xy2):
    """
        Returns the angle in radiant from each point in xy1 to each point in 
        xy2. Angles range from 0 (north) turning clockwise to 6.28... rad
        (360 degree).


        Definition
        ----------
        get_angle(xy1, xy2):


        Input
        -----
        xy1       2D np.array (n,2),
                  x and y coordinates of the points from where angles are desired 
        xy2       2D np.array (m,2), x and y coordinates of the points to where
                  angles are desired 


        Output
        ------
        angle     2D np.array (n,m), angles from each point in xy1 to all points
                  in xy2 in radiant.


        Restrictions
        ------------
        No support for 1D arrays.


        Example
        --------
        >>> ang = np.rad2deg(get_angle(np.array([[0.,0.]]), np.array([[0.,1.]])))
        >>> print(np.ma.round(ang).astype(int))
        [[0]]
        
        >>> ang = np.rad2deg(get_angle(np.array([[0.,0.]]), np.array([[1.,0.]])))
        >>> print(np.ma.round(ang).astype(int))
        [[90]]
        
        >>> ang = np.rad2deg(get_angle(np.array([[0.,0.]]), np.array([[-1.,0.]])))
        >>> print(np.ma.round(ang).astype(int))
        [[270]]
        
        >>> ang = np.rad2deg(get_angle(np.array([[0.,0.],[1.,1.]]), np.array([[1.,1.],[1.,0.]])))
        >>> print(np.ma.round(ang).astype(int))
        [[45 90]
         [-- 180]]
         

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014 Arndt Piayda

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
        Written,  AP, Jun 2014
    """
    
    oldsettings = np.geterr()
    np.seterr(divide='ignore', invalid='ignore')
    
    rad180 = np.deg2rad(180.)
    rad360 = np.deg2rad(360.)
    angle = np.ma.empty((xy1.shape[0],xy2.shape[0]))
    for i, item in enumerate(xy1):
        adjacent_leg = xy2[:,1]-item[1]
        opposite_leg = xy2[:,0]-item[0]
        hypothenuse  = np.ma.sqrt((xy2[:,0]-item[0])**2 + (xy2[:,1]-item[1])**2)
        rad          = np.ma.arcsin(opposite_leg/hypothenuse)
        
        angle[i,:] = np.ma.where(
                                 # 1. quadrant
                                 (adjacent_leg >= 0.) & (opposite_leg >= 0.), rad,
                                 # 2. quadrant
                                 np.where((adjacent_leg <  0.) & (opposite_leg >= 0.), rad180-rad,
                                 # 3. quadrant
                                 np.where((adjacent_leg <  0.) & (opposite_leg <  0.), np.abs(rad)+rad180,
                                 # 4. quadrant
                                 rad360+rad)))
    
    np.seterr(**oldsettings)
    
    return angle

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
