#!/usr/bin/env python
from __future__ import print_function
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
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Arndt Piayda


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
