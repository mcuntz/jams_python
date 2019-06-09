#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['sn','cz','hi','ef','aed']

# Overview of the metrics (Table 1 in Zheng F et al., [Draft])
#
#   -------------------------------------------------------------------------------------------------------------------
#   Metric       Symbol      Metric range              Run-time search behavior   Comments
#   Category                 Value interpretation      characteristic
#   -------------------------------------------------------------------------------------------------------------------
#   Solution     SN          [0,Infty]                 MOA's ability to find      Wang et al. (2015), but 
#   quality                  larger value is better    solutions that are on      current study is the first to 
#   metric                                             Pareto front               consider its run-time statistics
#   -------------------------------------------------------------------------------------------------------------------
#   Solution     CZ          [0,Infty]                 Relative performance of    Developed by Zheng et al. [Draft]
#   quality                  larger value is better    MOAs in providing non-     
#   metric                                             dominated solutions
#   -------------------------------------------------------------------------------------------------------------------
#   Spacing      EF          [0,1]                     Front extend within the    Developed by Zheng et al. [Draft]
#   metric                   larger value is better    objective space     
#                                                     
#   -------------------------------------------------------------------------------------------------------------------
#   Convergence  HI          [0,1]                     Solution quality,          Hadka and Reed (2012), but
#   metric                   larger value is better    diversity and uniformity   current study is the first to
#                                                                                 consider its run-time variation
#   -------------------------------------------------------------------------------------------------------------------
#   Convergence  AED         [0,1]                     Convergence in objective   Kollat and Reed (2006)
#   metric                   smaller value is better   space   
#                                                      
#   -------------------------------------------------------------------------------------------------------------------
#   Convergence  Var         [0,1]                     Convergence in decision    Introduced in Zharie (2002), but
#   metric                                             space                      current study is the first to
#                                                                                 consider it in multi-objective space
#   -------------------------------------------------------------------------------------------------------------------
#
#
# Literature
# ----------
#        Hadka, D., & Reed, P. M. (2013).
#        Borg: An Auto-Adaptive Many-Objective Evolutionary Computing Framework.
#        Evolutionary Computation, 21(2), 231-259.
#
#        Kollat, J. B., & Reed, P. M. (2006).
#        Comparing state-of-the-art evolutionary multi-objective algorithms for long-term groundwater monitoring design.
#        Advances in Water Resources, 29(6), 792-807. doi:10.1016/j.advwatres.2005.07.010
#
#        Wang, Q., Guidolin, M., Savic, D., & Kapelan, Z. (2015).
#        Two-Objective Design of Benchmark Problems of a Water Distribution System via MOEAs: Towards the Best-Known Approximation of the True Pareto Front.
#        Journal of Water Resources Planning and Management, 141(3), 04014060. doi:10.1061/(ASCE)WR.1943-5452.0000460
#
#        Zaharie, D. (2002).
#        Critical values for the control parameters of differential evolution algorithms (pp. 62-67).
#        Presented at the 8th International Mendel Conference on Soft Computation.
#
#        F Zheng, AC Zecchin, HR Maier, and AR Simpson
#        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water
#        distribution system design problems
#        [Draft]

def is_dominated(nobj, point, front):
    # checks whether a new canidate 'point' is dominated by the 'front' or not
    # does not manipulate 'front'
    #
    # returns: 1 --> point is dominating other points on front and hence should be added to the 'front' and others removed
    #          0 --> point is neither dominating others nor dominated by others, ie. it is either part of the 'front' or
    #                is a new point filling gaps or extends front and should hence be added
    #         -1 --> point is dominated
    
    nsets = np.shape(front)[0]

    dominanceFlag = 0
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
            # 'point' is dominated
            dominanceFlag = -1
            return dominanceFlag
        else:
            if (numeql == nobj):
                # Objective functions are the same for 'point' and archived solution ii
                dominanceFlag = 0
                return dominanceFlag
            else:
                if (numimp > 0) and (numdeg == 0):
                    # 'point' dominates ii-th solution in the front
                    dominanceFlag = 1
                    return dominanceFlag
    
    return dominanceFlag

def point_to_front(nobj, point, front):
    # It checks if 'point' is dominating some points of 'front'.
    # If 'point' dominates others, it adds 'point' to 'front' and deletes dominated ones. At the end the updated 'front'
    # is returned. It hence is doing the same as the routine 'is_dominated', but returns updated front instead of only
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
            # 'point' is dominated --> return 'front' unchanged
            dominanceFlag = -1
            return front
        else:
            if (numeql == nobj):
                # Objective functions are the same for 'point' and archived solution ii
                # Replace solution ii in 'front' with 'point'
                front[ii] = point
                dominanceFlag = 0
                return front
            else:
                if (numimp > 0) and (numdeg == 0):
                    # 'point' dominates ii-th solution in the 'front'
                    if (nsets > 1):
                        # Remove solution ii from 'front'
                        front = np.delete(front,ii,axis=0)
                        nsets  -= 1
                    else:
                        # Remove last point in 'front'
                        # This means all points on 'front' are dominated by new candidate 'point'
                        front = np.array([])
                        nsets -= 1
                    ii -= 1
                    dominanceFlag = 1

    if (dominanceFlag == 0):
        # 'point' is a new point of the 'front' but does not dominate any current member
        # hence 'point' will be added
        np.vstack([front,np.array(point)])

    if (nsets == 0):
        new_set = np.array(point)
    else:
        new_set = np.vstack((front,np.array(point)))

    return new_set

def sn(front, reference_front):
    """
        (Solution quality metric)
        
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
        front               members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
        reference_front     all members of the best Pareto front (often unknown)
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
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2016 Juliane Mai

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
        Written,  JM, Jan 2016 
    """

    # checks if elements of current front are members of reference front
    return np.sum(np.array([ np.where(np.all(reference_front==ff,axis=1))[0].shape[0] for ff in list(front) ]))

def cz(front, all_fronts, all_cz=False):
    """
        (Solution quality metric)
        
        To charcterize and compare the performance of different Multi-objective algorithms (MOAs), the contribution of
        each algorithm to the aggregated best front from all MOEAs can be used.
        
        CZ = Contribution to aggregated best front Z
        

        Definition
        ----------
        def cz(front, all_fronts, all_cz=False):


        Input
        -----
        front               members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
        all_fronts          fronts of all fronts under consideration
                            these will be aggregated and the amount of points from 'front' in the aggreagtion counted
                            list of k arrays with n2(k) rows (number of points) and m cols (number of objectives)


        Optional Input
        --------------
        all_cz              if True the CZ of every member of all_fronts is returned as a list of integers
                            'front' is then ignored


        Output
        ------
        Scalar integer giving the number of points of the current front which are nondominated by any other point of
        all given fronts 'all_fronts'. Scalar is less equal sum( n2(k) ).
        If 'all_cz' True, then output is list of integers.


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
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2016 Juliane Mai

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


def hi(front, reference_point, nsamples=None, reference_front=None, hi_range=False):
    """
        (Convergence metric)
        
        The hypervolume indicator (HI) is a measure for soultion distribution and uniformity in objective space (Zitzler and Thiele, 1999).
        It measures the hypervolume of the multi-dimensional region enclosed by a front with respect to a reference point.
        It represents the overall searching performance of solution quality, solution diversity and the uniformity of the solutions on the
        front (Hadka and Reed, 2012).

        By default, the routine returns an estimation of the hypervolume dominated by the approximated Pareto front
        set 'front' and bounded by the reference point 'reference_point'. This means the volume covered by the 'front' relative to the
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
        argument 'reference front' is given. Hence, HI lies within [0,1] which a large value representing a hypervolume closer to Z*.

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
        def hi(front, reference_point, nsamples=None, reference_front=None, hi_range=False):


        Input
        -----
        front               members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
        reference_point     worst possible point in objective space


        Optional Input
        --------------
        nsamples            number of random points used to approximate the hypervolume
                            (default: 1000)
        reference_front     best known pareto front
                            if given the HI is returned as the ratio between HI of tested front and HI of the reference front
                            value is hence bounded between 0 and 1 where large values indicate a front closer to the reference front
        hi_range            if True an lower and upper value for HI will be returned
                            this range is due to pointwise given 'front', i.e. the coarser the 'front' the more randomly sampled points
                            have a dominance flag of zero and the more uncertain is the HI

                            
        Output
        ------
        if( reference_front == None):
            Scalar floating number of the ratio between hypervolume covered by the 'front' and
            hypervolume bounded by 'reference point'
        else:
            Scalar floating number of the ratio between hypervolume covered by the 'front' and
            hypervolume covered by 'reference front'

            
        Literature
        ----------
        F Zheng, AC Zecchin, HR Maier, and AR Simpson
        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water distribution system design problems
        [Draft]

        Hadka, D., and Reed, P. (2012).
        'Diagnostic Assessment of Search Controls and Failure Modes in Many-Objective Evolutionary Optimization.'
        Evolutionary Computation, 20(3), 423-452.

        Zitzler, E., and Thiele, L. (1999).
        'Multiobjective evolutionary algorithms: a comparative case study and the strength Pareto approach.'
        Evolutionary Computation, IEEE Transactions on, 3(4), 257-271.

        
        Restrictions
        ------------
        Objectives are assumed to be minimized.


        Examples
        --------
        >>> from autostring import astr
        >>> front           = np.array([[10,12],[12.5,11],[15,10],[17.5,9],[20,8]])
        >>> reference_point = np.array([25,15])
        >>> reference_front = np.array([[8,12],[11,10],[14,8],[17,6],[20,4]])
        >>> hypervolume_indicator = hi(front, reference_point, nsamples=10000)
        >>> # theoretical value: 85/375 = 0.226667
        >>> print(astr(hypervolume_indicator,prec=4))
        0.2321
        >>> hypervolume_indicator = hi(front, reference_point, nsamples=10000, hi_range=True)
        >>> # theoretical value: 85/375 = 0.226667
        >>> print(astr(hypervolume_indicator,prec=4))
        ['0.2195' '0.2448']
        >>> hypervolume_indicator = hi(front, reference_point, reference_front=reference_front, nsamples=10000)
        >>> # theoretical value: 85/139 = 0.6115
        >>> print(astr(hypervolume_indicator,prec=4))
        0.6258
        >>> hypervolume_indicator = hi(front, reference_point, reference_front=reference_front, nsamples=10000, hi_range=True)
        >>> # theoretical value: 85/139 = 0.6115
        >>> print(astr(hypervolume_indicator,prec=4))
        ['0.5440' '0.7234']
        >>> front           = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> reference_point = np.array([25,15])
        >>> reference_front = np.array([[3,13],[3,11],[3,9],[4,8],[5,7],[6,6],[7,5],[8,4],[10,3],[12,2],[14,1],[16,1],[18,1],[20,1],[22,1]])
        >>> hypervolume_indicator = hi(front, reference_point)
        >>> print(astr(hypervolume_indicator,prec=4))
        0.6180
        >>> hypervolume_indicator = hi(front, reference_point, reference_front=reference_front)
        >>> print(astr(hypervolume_indicator,prec=4))
        0.8477
        

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2016 Juliane Mai

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
        Written,  JM, Feb 2016 
    """

    import random
    from autostring import astr

    # check for optionals
    if (nsamples is None):
        nsamp = 1000
    else:
        nsamp = nsamples
    

    # initialization
    nsets_front = front.shape[0]
    nobj        = front.shape[1]
    random.seed(1000)

    # add edge points to 'front'
    best_all_directions = np.min(front,axis=0)
    for iobj in range(nobj):
        edge = best_all_directions.copy()
        edge[iobj] = reference_point[iobj]
        # print(edge)
        # Or is it ??
        # edge = reference_point
        # edge[iobj] = best_all_directions[iobj]
        front = np.vstack([front,np.array(edge)])

    # add edge points to 'reference_front'
    if ( not(reference_front is None) ):
        best_all_directions = np.min(reference_front,axis=0)
        for iobj in range(nobj):
            edge = best_all_directions.copy()
            edge[iobj] = reference_point[iobj]
            # print(edge)
            # Or is it ??
            # edge = reference_point
            # edge[iobj] = best_all_directions[iobj]
            reference_front = np.vstack([reference_front,np.array(edge)])

    # draw random points in feasible domain bounded by 0 and 'reference_point'
    mc_sample = [ [ random.random()*reference_point[ii] for ii in range(nobj) ] for jj in range(nsamp) ]

    # hypervolume of 'front'
    dom_flags = [ is_dominated(nobj, sample, front) for sample in mc_sample ]
    n_below_front = np.where(np.array(dom_flags)== 1)[0].shape[0]
    n_above_front = np.where(np.array(dom_flags)==-1)[0].shape[0]   # this will be the hypervolume
    n_on_front    = np.where(np.array(dom_flags)== 0)[0].shape[0]   # should never happen, otherwise adding of edges is wrong or it is an integer problem

    # front_file = open('hypervolume_front.out','w')
    # for isamp in range(nsamp):
    #     front_file.write(astr(dom_flags[isamp]) + '  ' + astr(mc_sample[isamp][0],prec=4) + '  ' + astr(mc_sample[isamp][1],prec=4) + '  ' +  '\n') 
    # front_file.close()

    # hypervolume of 'reference_front'
    if ( not(reference_front is None) ):
        dom_flags = [ is_dominated(nobj, sample, reference_front) for sample in mc_sample ]
        n_below_reffront = np.where(np.array(dom_flags)== 1)[0].shape[0]
        n_above_reffront = np.where(np.array(dom_flags)==-1)[0].shape[0]   # this will be the hypervolume
        n_on_reffront    = np.where(np.array(dom_flags)== 0)[0].shape[0]   # should never happen, otherwise adding of edges is wrong or it is an integer problem

        # reffront_file = open('hypervolume_reffront.out','w')
        # for isamp in range(nsamp):
        #     reffront_file.write(astr(dom_flags[isamp]) + '  ' + astr(mc_sample[isamp][0],prec=4) + '  ' + astr(mc_sample[isamp][1],prec=4) + '  ' + '\n') 
        # reffront_file.close()
    
    if (reference_front is None):
        if hi_range:
            hi = np.array([(1.0*n_above_front + 0.0*n_on_front) / (1.0*nsamp),
                           (1.0*n_above_front + 1.0*n_on_front) / (1.0*nsamp) ])
        else:
            hi = (1.0*n_above_front + 0.5*n_on_front) / (1.0*nsamp)
    else:
        if hi_range:
            hi = np.array([(1.0*n_above_front + 0.0*n_on_front) / (1.0*n_above_reffront + 1.0*n_on_reffront),
                           (1.0*n_above_front + 1.0*n_on_front) / (1.0*n_above_reffront + 0.0*n_on_reffront) ])
        else:
            hi = (1.0*n_above_front + 0.5*n_on_front) / (1.0*n_above_reffront + 0.5*n_on_reffront)
        
    return hi


def ef(front, reference_front):
    """
        (Spacing metric)
        
        This metric measures the approximate extend of the front (EF) relative to the known best Pareto front Z* given by

             EF = max_(i,j = 1,...,N) ( dist(F_i, F_j) ) / max_(F,G in Z*) ( dist(F, G) )

        where the nominator indicates the maximum distance between two points on an algorithm's front and the denominator
        is the maximum distance between two points on the known Pareto front Z*. dist(X,Y) is the Euclidean distance between
        two vectors X and Y. To calculate EF, each dimension of the solution vectors in Z* is normalized to [0,1] initially,
        followed by the normalization on each dimension of the 'front' Z using the data ranges from Z*.

        A larger value of EF suggests that the approximate front is closer to the Pareto front in terms of solution extend.
        It should be noted that the proposed EF only measures the coverage of the two-objective Pareto approximate fronts, and
        that the overall quality of the fronts has to be assessed using other metrics.
   

        Definition
        ----------
        def ef(front, reference_front)


        Input
        -----
        front               members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
        reference_front     all members of the best Pareto front (often unknown)
                            array with n2 rows (number of points) and m cols (number of objectives)


        Optional Input
        --------------
        None

                            
        Output
        ------
        Scalar floating number determining the extend of the 'front' relative to the extend of the 'reference_front'.

            
        Literature
        ----------
        F Zheng, AC Zecchin, HR Maier, and AR Simpson
        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water distribution system design problems
        [Draft]

        
        Restrictions
        ------------
        Objectives are assumed to be minimized.


        Examples
        --------
        >>> from autostring import astr
        >>> reference_front = np.array([[3,13],[3,11],[3,9],[4,8],[5,7],[6,6],[7,5],[8,4],[10,3],[12,2],[14,1],[16,1],[18,1],[20,1],[22,1]])
        >>> front_1 = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> front_2 = np.array([[5,10],[5,9],[6,8],[7,7],[8,6],[9,5],[10,4],[12,3],[14,3],[15,2],[16,1],[18,1]])
        >>> front_3 = np.array([[7,10],[7,8],[8,7],[10,6],[12,5],[13,4],[16,4],[18,3],[20,3]])
        >>> extend_front = ef(front_1, reference_front)
        >>> print(astr(extend_front,prec=4))
        0.6933
        >>> extend_front = ef(front_2, reference_front)
        >>> print(astr(extend_front,prec=4))
        0.7179
        >>> extend_front = ef(front_3, reference_front)
        >>> print(astr(extend_front,prec=4))
        0.6358
        

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2016 Juliane Mai

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
        Written,  JM, Feb 2016 
    """
    from scipy.spatial.distance import pdist

    nobj = reference_front.shape[1]
    if ( nobj != 2 ):
        raise ValueError('pareto_metrics: This metric is only applicable for 2-dimensional fronts!')

    # assure that input is floating and not integer
    reference_front = reference_front.astype(float)
    front           = front.astype(float)

    # column-wise range of reference-front used to scale 'front' and 'reference_front'
    col_min = np.min(reference_front,axis=0)
    col_max = np.max(reference_front,axis=0)
    
    reference_front     = (reference_front - col_min ) / (col_max - col_min)
    front               = (          front - col_min ) / (col_max - col_min)
    extend_front = np.max(pdist(front)) / np.max(pdist(reference_front))
    return extend_front


def aed(front, reference_front):
    """
        (Convergence metric)
        
        The average Euclidean distance (AED) between an approximate front and the known best Pareto front Z* is, besides HI,
        typically used to represent an MOA's convergence status in the objective space (Kollat and Reed, 2006).
        The average Euclidean distance represents the average of all minimal Euclidean distances between a point of the 
        'front' to all points of the 'reference front'.

            AED = 1/N * sum_{i=1}^N min( dist(F_i, F) ),  where dist(F_i,F) is Euclidean distance of point F_i of 'front'
                                                                            to the points of 'reference front' F

        To allow for a comparison of AED values of different 'fronts', a normalization of the AED is applied. Therefore, all
        dimensions of the reference front Z* are normalized using the minimal and maxiaml values of each dimension. All
        dimensions are hence normaized to [0,1]. The same normalization is subsequently applied to the dimensions 'front' using the
        minimal and maximal values of the 'reference front'. The AED is thus within the range [0,1] where a lower value indicates a
        better front in general, as it possesses an overall shorter distance to the Pareto front in objective space.

        It should be noted that this metric does not take the extend of the 'front' into account as it measures the mean of the
        Euclidean distance between each solution F_i of the 'front' and its corresponding nearest member in the refence front Z*.
   

        Definition
        ----------
        def aed(front, reference_front)


        Input
        -----
        front               members of current front
                            array with n1 rows (number of points) and m cols (number of objectives)
        reference_front     all members of the best Pareto front (often unknown)
                            array with n2 rows (number of points) and m cols (number of objectives)


        Optional Input
        --------------
        None

                            
        Output
        ------
        Scalar floating number determining the average euclidean distance of the 'front' to the 'reference_front'

            
        Literature
        ----------
        F Zheng, AC Zecchin, HR Maier, and AR Simpson
        Comparison of the searching behavior of NSGA-II, SAMODE and Borg MOEAs applied to water distribution system design problems
        [Draft]

        
        Restrictions
        ------------
        Objectives are assumed to be minimized.


        Examples
        --------
        >>> from autostring import astr
        >>> reference_front = np.array([[3,13],[3,11],[3,9],[4,8],[5,7],[6,6],[7,5],[8,4],[10,3],[12,2],[14,1],[16,1],[18,1],[20,1],[22,1]])
        >>> front_1 = np.array([[6,11],[6,8],[7,6],[8,5],[10,3],[12,2],[15,2],[18,2]])
        >>> front_2 = np.array([[5,10],[5,9],[6,8],[7,7],[8,6],[9,5],[10,4],[12,3],[14,3],[15,2],[16,1],[18,1]])
        >>> front_3 = np.array([[7,10],[7,8],[8,7],[10,6],[12,5],[13,4],[16,4],[18,3],[20,3]])
        >>> avergae_euclid_dist = aed(front_1, reference_front)
        >>> print(astr(avergae_euclid_dist,prec=4))
        0.0071
        >>> avergae_euclid_dist = aed(front_2, reference_front)
        >>> print(astr(avergae_euclid_dist,prec=4))
        0.0090
        >>> avergae_euclid_dist = aed(front_3, reference_front)
        >>> print(astr(avergae_euclid_dist,prec=4))
        0.0341
        

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2016 Juliane Mai

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
        Written,  JM, Feb 2016 
    """

    # assure that input is floating and not integer
    reference_front = reference_front.astype(float)
    front           = front.astype(float)

    # column-wise range of reference-front used to scale 'front' and 'reference_front'
    col_min = np.min(reference_front,axis=0)
    col_max = np.max(reference_front,axis=0)
    
    reference_front     = (reference_front - col_min ) / (col_max - col_min)
    front               = (          front - col_min ) / (col_max - col_min)
    average_euclid_dist = np.mean([ np.min(np.sum((reference_front-ifront)**2,axis=1)) for ifront in front ])
    return average_euclid_dist



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)



    


