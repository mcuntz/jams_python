#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['sn','cz','hi']

def is_dominated(nobj, point, front):
    # checks whether a new canidate "point" is dominated by the "front" or not
    # does not manipulate "front"
    #
    # returns: 1 --> point is not dominated and hence should added to the "front"
    #          0 --> point is already part of the "front"
    #         -1 --> point is dominated
    
    nsets = np.shape(front)[0]

    dominanceFlag = 0
    ii = -1
    # print(' ')
    # print('point: ',point)
    while (ii < nsets-1):
        ii += 1
        numeql = 0
        numimp = 0
        numdeg = 0

        for jj in range(nobj):
            if (point[jj] == front[ii,jj]):
                numeql += 1
            if (point[jj] <  front[ii,jj]):
                numimp += 1
            if (point[jj] >  front[ii,jj]):
                numdeg += 1
        # print('numeql = ',numeql,'   numimp = ',numimp, '   numdeg = ',numdeg)
                
        if (numimp == 0) and (numdeg > 0):
            # "point" is dominated
            dominanceFlag = -1
            # print('return 1: ',dominanceFlag)
            return dominanceFlag
        else:
            if (numeql == nobj):
                # Objective functions are the same for "point" and archived solution ii
                dominanceFlag = 0
                # print('return 2: ',dominanceFlag)
                return dominanceFlag
            else:
                if (numimp > 0) and (numdeg == 0):
                    # "point" dominates ii-th solution in the front
                    dominanceFlag = 1
                    # print('return 3: ',dominanceFlag)
                    return dominanceFlag
    # print('return 4: ',dominanceFlag)
    return dominanceFlag

def point_to_front(nobj, point, front):
    # It checks if "point" is dominating some points of "front".
    # If "point" dominates others, it adds "point" to "front" and deletes dominated ones. At the end the updated "front"
    # is returned. It hence is doing the same as the routine "is_dominated", but returns updated front instead of only
    # information about the domination status.
    #
    # (originally routine NDcheck in Matlab version by Masoud)
    #
    # returns: set which only consists of non-dominated points

    dominanceFlag = 0
    
    if len(front) == 0:
        new_set = np.array(point)
        return new_set

    if (len(np.shape(front)) == 1):
        front = np.expand_dims(front, axis=0)

    nsets = np.shape(front)[0]

    ii = -1 
    while (ii < nsets-1):
        ii += 1
        numeql = 0
        numimp = 0
        numdeg = 0

        for jj in range(nobj):
            if (point[jj] == front[ii,jj]):
                numeql += 1
            if (point[jj] <  front[ii,jj]):
                numimp += 1
            if (point[jj] >  front[ii,jj]):
                numdeg += 1
                
        if (numimp == 0) and (numdeg > 0):
            # "point" is dominated --> return "front" unchanged
            dominanceFlag = -1
            return front
        else:
            if (numeql == nobj):
                # Objective functions are the same for "point" and archived solution ii
                # Replace solution ii in "front" with "point"
                front[ii] = point
                dominanceFlag = 0
                return front
            else:
                if (numimp > 0) and (numdeg == 0):
                    # "point" dominates ii-th solution in the "front"
                    if (nsets > 1):
                        # Remove solution ii from "front"
                        front = np.delete(front,ii,axis=0)
                        nsets  -= 1
                    else:
                        # Remove last point in "front"
                        # This means all points on "front" are dominated by new candidate "point"
                        front = np.array([])
                        nsets -= 1
                    ii -= 1
                    dominanceFlag = 1

    if (dominanceFlag == 0):
        # "point" is a new edge of the "front"
        # hence "point" will be added and extent of "front" increases
        np.vstack([front,np.array(point)])

    if (nsets == 0):
        new_set = np.array(point)
    else:
        new_set = np.vstack((front,np.array(point)))

    return new_set

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

def cz(front, all_fronts, all_cz=False):
    """
        To charcterize and compare the performance of different Multi-objective algorithms (MOAs), the contribution of
        each algorithm to the aggregated best front from all MOEAs can be used.
        
        CZ = Contribution to aggregated best front Z
        

        Definition
        ----------
        def cz(front, all_fronts, all_cz=False):


        Input
        -----
        front			    members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
	    all_fronts   		fronts of all fronts under consideration
                            these will be aggregated and the amount of points from "front" in the aggreagtion counted
                            list of k arrays with n2(k) rows (number of points) and m cols (number of objectives)


        Optional Input
        --------------
        all_cz              if True the CZ of every member of all_fronts is returned as a list of integers
                            "front" is then ignored


        Output
        ------
        Scalar integer giving the number of points of the current front which are nondominated by any other point of
        all given fronts "all_fronts". Scalar is less equal sum( n2(k) ).
        If "all_cz" True, then output is list of integers.


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
        >>> front_1 = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> front_2 = np.array([[5,10],[5,9],[6,8],[7,7],[8,6],[9,5],[10,4],[12,3],[14,3],[15,2],[16,1],[18,1]])
        >>> front_3 = np.array([[7,10],[7,8],[8,7],[10,6],[12,5],[13,4],[16,4],[18,3],[20,3]])
        >>> all_fronts = [ front_1, front_2, front_3 ]
        >>> contrib_z = cz(front_1, [ front_1, front_2, front_3 ])
        >>> print(astr(contrib_z))
        5
        >>> contrib_z = cz(front_2, all_fronts)
        >>> print(astr(contrib_z))
        3
        >>> contrib_z = cz(front_3, all_fronts)
        >>> print(astr(contrib_z))
        0
        >>> contrib_z = cz(front_3, all_fronts, all_cz=True)
        >>> print(astr(contrib_z))
        ['5' '3' '0']
        

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
        Written,  JM, Feb 2016 
    """
    # aggregate the fronts
    nfronts = len(all_fronts)
    nobj    = all_fronts[0].shape[1]
    aggretated_front = []
    for ifronts in range(nfronts):
        nsets = all_fronts[ifronts].shape[0]
        for isets in range(nsets):
            aggretated_front = point_to_front(nobj, all_fronts[ifronts][isets], aggretated_front) 
        
    # counts elements of current front which are members of aggregated front
    if not(all_cz):
        return np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(front) ]))
    else:
        return [ np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(kk) ])) for kk in all_fronts ]


def hi(front, reference_point, nsamples=None, reference_front=None):
    """
        (Convergence metric)
        
        The hypervolume indicator (HI) is a measure for soultion distribution and uniformity in objective space (Zitzler and Thiele, 1999).
        It measures the hypervolume of the multi-dimensional region enclosed by a front with respect to a reference point.
        It represents the overall searching performance of solution quality, solution diversity and the uniformity of the solutions on the
        front (Hadka and Reed, 2012).

        By default, the routine returns an estimation of the hypervolume dominated by the approximated Pareto front
        set "front" and bounded by the reference point "reference_point". This means the volume covered by the "front" relative to the
        overall feasible hypervolume.

                ^
                |
             10 |---------------------------------------* (reference point)
                |                 |                     |
                |                 |       A ~ 20        |
                |                 *                     |              ==> HI = 20 / 50 = 0.4
                |                  *                    |
                |            (front) *                  |
                |                        *--------------|
                |                                       |
                -------------------------------------------------->
                                                  5

        Optionally, HI can be returned as the ratio of the hypervolume relative to that of the known best Pareto front Z* when
        argument "reference front" is given. Hence, HI lies within [0,1] which a large value representing a hypervolume closer to Z*.

                ^
                |
             10 |---------------------------------------* (reference point)
                |       |         |                     |
                |       | B ~ 30  |      A ~ 20         |
                |       *         *                     |              ==> HI = A / B = 20 / 30 = 0.667
                |        *         *   (front)          |
                |         *          *                  |
                |   (best)  *            *--------------|
                |               *-----------------------|
                -------------------------------------------------->
                                                  5        

        Definition
        ----------
        def hi(front, reference_point, nsamples=None, reference_front=None):


        Input
        -----
        front			    members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
	    reference_point		worst possible point in objective space


        Optional Input
        --------------
        nsamples            number of random points used to approximate the hypervolume
                            (default: 1000)
        reference_front     best known pareto front
                            if given the HI is returned as the ratio between HI of tested front and HI of the reference front
                            value is hence bounded between 0 and 1 where large values indicate a front closer to the reference front

                            
        Output
        ------
        if( reference_front == None):
            Scalar floating number of the ratio between hypervolume covered by the "front" and
            hypervolume bounded by "reference point"
        else:
            Scalar floating number of the ratio between hypervolume covered by the "front" and
            hypervolume covered by "reference front"

            
        Literature
        ----------
        F Zheng, AC Zecchin, HR Maier, and AR Simpson
        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water distribution system design problems
        [Draft]

        Hadka, D., and Reed, P. (2012).
        "Diagnostic Assessment of Search Controls and Failure Modes in Many-Objective Evolutionary Optimization."
        Evolutionary Computation, 20(3), 423-452.

        Zitzler, E., and Thiele, L. (1999).
        "Multiobjective evolutionary algorithms: a comparative case study and the strength Pareto approach."
        Evolutionary Computation, IEEE Transactions on, 3(4), 257-271.

        
        Restrictions
        ------------
        Objectives are assumed to be minimized.


        Examples
        --------
        >>> from autostring import astr
        >>> front           = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> reference_point = np.array([25,15])
        >>> reference_front = np.array([[3,13],[3,11],[3,9],[4,8],[5,7],[6,6],[7,5],[8,4],[10,3],[12,2],[14,1],[16,1],[18,1],[20,1],[22,1]])
        >>> hypervolume_indicator = hi(front, reference_point)
        >>> print(astr(hypervolume_indicator,prec=4))
        0.6060
        >>> hypervolume_indicator = hi(front, reference_point, reference_front=reference_front)
        >>> print(astr(hypervolume_indicator,prec=4))
        0.8511
        

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
        Written,  JM, Feb 2016 
    """

    import random

    # check for optionals
    if (nsamples == None):
        nsamp = 1000
    else:
        nsamp = nsamples
    

    # hypervolume of "front"
    nsets_front = front.shape[0]
    nobj        = front.shape[1]
    random.seed(1000)

    # add edge points to "front"
    best_all_directions = np.min(front,axis=0)
    for iobj in range(nobj):
        edge = best_all_directions.copy()
        edge[iobj] = reference_point[iobj]
        # print(edge)
        # Or is it ??
        # edge = reference_point
        # edge[iobj] = best_all_directions[iobj]
        front = np.vstack([front,np.array(edge)])

    # hypervolume of "front"
    mc_sample = [ [ random.random()*reference_point[ii] for ii in range(nobj) ] for jj in range(nsamp) ]
    dom_flags = [ is_dominated(nobj, sample, front) for sample in mc_sample ]
    n_below_front = np.where(np.array(dom_flags)== 1)[0].shape[0]
    n_above_front = np.where(np.array(dom_flags)==-1)[0].shape[0]  # this will be the hypervolume
    n_on_front    = np.where(np.array(dom_flags)== 0)[0].shape[0]   # should never happen, otherwise adding of edges is wrong or it is an integer problem

    if (reference_front != None):
        # hypervolume of "reference_front"
        dom_flags = [ is_dominated(nobj, sample, reference_front) for sample in mc_sample ]
        n_below_reffront = np.where(np.array(dom_flags)== 1)[0].shape[0]
        n_above_reffront = np.where(np.array(dom_flags)==-1)[0].shape[0]  # this will be the hypervolume
        n_on_reffront    = np.where(np.array(dom_flags)== 0)[0].shape[0]   # should never happen, otherwise adding of edges is wrong or it is an integer problem
    
    if (reference_front == None):
        hi = (1.0 * n_above_front) / (1.0 * nsamp)
    else:
        hi = (1.0 * n_above_front) / (1.0 * n_above_reffront)
        
    return hi


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # front_1 = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
    # front_2 = np.array([[5,10],[5,9],[6,8],[7,7],[8,6],[9,5],[10,4],[12,3],[14,3],[15,2],[16,1],[18,1]])
    # front_3 = np.array([[7,10],[7,8],[8,7],[10,6],[12,5],[13,4],[16,4],[18,3],[20,3]])
    # all_fronts = [ front_1, front_2, front_3 ]

    # # aggregate the fronts
    # nfronts = len(all_fronts)
    # nobj    = all_fronts[0].shape[1]
    # aggretated_front = []
    # for ifronts in range(nfronts):
    #     nsets = all_fronts[ifronts].shape[0]
    #     for isets in range(nsets):
    #         aggretated_front = point_to_front(nobj, all_fronts[ifronts][isets], aggretated_front)

    # print(aggretated_front)
        
    # # counts elements of current front which are members of aggregated front
    # print( np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(front_1) ])) )
    # print( np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(front_2) ])) )
    # print( np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(front_3) ])) )
    # print( [ np.sum(np.array([ np.where(np.all(aggretated_front==ff,axis=1))[0].shape[0] for ff in list(kk) ])) for kk in all_fronts ] )


    


