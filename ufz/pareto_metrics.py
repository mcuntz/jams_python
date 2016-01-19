#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['sn']

def sn(front, reference_front):
    """
        The run-time search quality of a multi-objective optimization algorithm can be measured by the number of solutions
        that form an approxiamte Pareto front at searching time t that are nondominated by the Pareto front Z* in objective space.
        This number of points is called success number SN(t) of the algorithm. SN(t) is basically counting the number of points
        of Z* detected by the algorithm after t function/model evaluations.

        For large problems the reference Pareto front Z* might be hard to determine and is hence often unknown.
        In this case SN can not be determined.

        Comparing two algorithms a larger value of SN(t) indicates the algorithm's greater ability to find solutions along
        the Pareto front Z*.
        

        Definition
        ----------
        def sn(front, reference_front):


        Input
        -----
        front			    members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
	    reference_front		all members of the best Pareto front (often unknown)
                            array with n2 rows (number of points) and m cols (number of objectives)


        Optional Input
        --------------
        None


        Output
        ------
        Scalar integer giving the number of points of the current front which are nondominated by the
        reference front. Scalar is less equal n2.


        Literature
        ----------
        F Zheng, AC Zecchin, HR Maier, and AR Simpson
        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water
        distribution system design problems
        [Draft]


        Restrictions
        ------------
        Objectives are assumed to be minimized.


        Examples
        --------
        >>> from autostring import astr
        >>> reference_front = np.array([[3,13],[3,11],[3,9],[4,8],[5,7],[6,6],[7,5],[8,4],[10,3],[12,2],[14,1],[16,1],[18,1],[20,1],[22,1]])
        >>> front = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> success_number = sn(front, reference_front)
        >>> print(astr(success_number))
        2
        

        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2016 Juliane Mai


        History
        -------
        Written,  JM, Jan 2016 
    """

    # checks if elements of current front are members of reference front
    return np.sum(np.array([ np.where(np.all(reference_front==ff,axis=1))[0].shape[0] for ff in list(front) ]))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # area = ellipse_area(1, 1)
    # print(astr(area, 3))

    # area = ellipse_area(2)
    # print(astr(area, 3))





def NDcheck(nobj, new, old_set):
    # checks if new_set is dominating some points of old_set
    # returns set which only consists of non-dominated points

    dominanceFlag = 0
    
    if len(old_set) == 0:
        new_set = np.array(new)
        return new_set

    if (len(np.shape(old_set)) == 1):
        old_set = np.expand_dims(old_set, axis=0)

    nsets = np.shape(old_set)[0]

    ii = -1 
    while (ii < nsets-1):
        ii += 1
        numeql = 0
        numimp = 0
        numdeg = 0

        for jj in range(nobj):
            if (new[jj] == old_set[ii,jj]):
                numeql += 1
            if (new[jj] <  old_set[ii,jj]):
                numimp += 1
            if (new[jj] >  old_set[ii,jj]):
                numdeg += 1
                
        if (numimp == 0) and (numdeg > 0):
            # new is dominated --> return old set
            dominanceFlag = -1
            return old_set
        else:
            if (numeql == nobj):
                # Objective functions are the same for new and archived solution ii
                # Replace solution ii in old_set with new
                old_set[ii] = new
                dominanceFlag = 0
                return old_set
            else:
                if (numimp > 0) and (numdeg == 0):
                    # new dominates ii-th solution in the old_set
                    if (nsets > 1):
                        # Remove solution ii from PF_set
                        old_set = np.delete(old_set,ii,axis=0)
                        nsets  -= 1
                    else:
                        # Remove last point in old_set
                        # This means all old points are dominated by new candidate
                        old_set = np.array([])
                        nsets -= 1
                    ii -= 1
                    dominanceFlag = 1

    if (nsets == 0):
        new_set = np.array(new)
    else:
        new_set = np.vstack((old_set,np.array(new)))

    return new_set
