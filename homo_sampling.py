#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import cdist

def homo_sampling(p1, p2, n, func=None, plot=False, maxit=1000):
    """       
        Generates homogeneously, randomly distributed sample points within a
        given rectangular area spanned by p1 (lower left) and p2 (upper right)
        corner. The area can be masked by a user given function to exclude
        certain areas of the rectangle from sampling.
        
        
        Definition
        ----------
        homo_sampling(p1, p2, n, func=None, plot=False, maxit=1000):
        
        
        Input
        ----- 
        p1          np.array((1,2)), coordinates of lower left corner (x,y) 
        p2          np.array((1,2)), coordinates of upper right corner (x,y)
        n           int, number of points to be generated
        
        
        Optional Input
        --------------
        func        callable function of the form func(xy) defined by user with
                    xy beeing a np.array((m,2)), returning a boolean
                    np.array((m)) with True when points in xy are NOT in the
                    desired area and False when points are in the desired area
        plot        bool, if True results are plotted (default=False)
        maxit       int, maximum number of iterations for maximizing the minimum
                    distance between points (default=1000)
                    
        
        Output
        ------
        xy_out      np.array((n,2)), x and y coordinates of points
        
        
        Restrictions
        ------------
        Does not deal with masked arrays.
        TODO:
            - a priory estimation of masked area fraction to increase n
              internally to reduce necessary iterations
            - include running window of max(minimum distance between points)
              to stop iteration when no further increase occurs within
              window-width=maxit

                             
        Examples
        --------
        >>> # Define the corner points of a rectangle
        >>> p1 = np.array([644242.398, 5766315.198])
        >>> p2 = np.array([644555.611, 5766578.860])
         
        >>> # Generate 10 random points within the rectangle
        >>> xy = homo_sampling(p1, p2, 10)
         
        >>> # Increase maximum number of iterations to distribute points
        >>> # even more homogeneous
        >>> xy = homo_sampling(p1, p2, 10, maxit=10000)
         
        >>> # Define a function that excludes the area outside a given radius r
        >>> # around a central point cp
        >>> def func(xy):
        ...     cp = np.array([[644394.412, 5766429.022]])
        ...     r = 100.
        ...     from scipy.spatial.distance import cdist
        ...     return cdist(cp, xy, metric='euclidean')[0] > r
         
        >>> # Generate points only within the radius around cp        
        >>> xy = homo_sampling(p1, p2, 10, func=func)
         
        >>> # Since the result are random numbers and can not be restricted to
        >>> # a certain seed, no result comparison within doctest can be performed.
        
        
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
    
        Copyright 2014 Arndt Piayda
    
    
        History
        -------
        Written,  AP, Jul 2014
        
    """
    # check
    if (len(p1)!=2) or (len(p2)!=2):
        raise ValueError('homo_sampling: p1 and p2 must be of length 2')
    elif np.any(p2-p1<0) or np.any(p2==p1):
        raise ValueError('homo_sampling: you mixed up the corners')
    elif func:
        if not callable(func):
            raise ValueError('homo_sampling: func must be a callable function')
        if np.ndim(func(np.array([[1.,1.]])))!=1:
            raise ValueError('homo_sampling: func must return 1D np.array')
        if func(np.array([[1.,1.]])).dtype!=np.bool:
            raise ValueError('homo_sampling: func must return boolean np.array')
    
    # generate 1st random point set
    xy  = np.empty((n,2))
    rem = np.ones((n), dtype=np.bool)
    while np.any(rem):
        xy[rem,0] = (p2[0] - p1[0])*np.random.random(np.sum(rem)) + p1[0]
        xy[rem,1] = (p2[1] - p1[1])*np.random.random(np.sum(rem)) + p1[1]
        
        rem = func(xy) if func else False 
        
    # maximize the distance between points
    conv =[]
    tempxy = []
    
    for i in xrange(maxit):
        d   = cdist(xy, xy, metric='euclidean')
        d   = np.ma.array(d, mask=d==0.)
        rep = np.argmin(np.ma.amin(d, axis=0))
        
        conv   += [np.ma.amin(d)]
        tempxy += [np.copy(xy)]
        
        rem = np.ones((1), dtype=np.bool)
        while np.any(rem):
            xnew = (p2[0] - p1[0])*np.random.random(1) + p1[0]
            ynew = (p2[1] - p1[1])*np.random.random(1) + p1[1]
    
            rem = func(np.array([[xnew[0],ynew[0]]])) if func else False
            
        xy[rep,:] = [xnew, ynew]
    
    xy_out = tempxy[np.argmax(conv)]
    
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        fig = plt.figure('homo_sampling: minimum distance=%f'%(np.amax(conv)))
        sub1 = fig.add_subplot(122, aspect='equal')
        sub1.set_xlim(p1[0],p2[0])
        sub1.set_ylim(p1[1],p2[1])
        sub1.scatter(xy_out[:,0],xy_out[:,1])
        sub1.set_xlabel('X')
        sub1.set_ylabel('Y')
        
        sub2 = fig.add_subplot(121)
        sub2.plot(np.arange(maxit),conv)
        sub2.plot(np.argmax(conv),np.amax(conv),'o')
        sub2.set_xlabel('iterations')
        sub2.set_ylabel('minimum distance between points')
        plt.show()
    
        print np.argmax(conv)
        print np.amax(conv)
    
    return xy_out

if __name__ == '__main__':
    import doctest
    doctest.testmod()
