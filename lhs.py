#!/usr/bin/env python
import numpy as np
import scipy.stats as stats

def lhs(dist, param, nsample):
    """
        Latin Hypercube Sampling of any distribution without correlations after Stein (1987).

        Definition
        ----------
        def lhs(dist, param, nsample):


        Input
        -----
        dist       random number generator (list) from scipy.stats such as stats.norm, stats.beta, etc.
        param      tuple of parameters as required for dist
        nsample    number of samples per parameter


        Output
        ------
        Latin hypercube Sample array of [size(nsample),nsample]


        Restrictions
        ------------
        No correlations between parameters possible.


        References
        ----------
        Stein, M. 1987. Large Sample Properties of Simulations Using Latin Hypercube Sampling. 
            Technometrics 29:143-151


        Examples
        --------
        >>> import numpy as np
        >>> import scipy.stats as stats
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)
        >>> dist = [stats.norm, stats.uniform]
        >>> pars = [(50,2),(1,5)]
        >>> c    = lhs(dist, pars, 20)
        >>> print c
        [[ 52.8216393   51.95643181  46.71040364  50.58536004  52.13686967
           48.67993776  50.29845536  50.86481869  48.0903917   49.10324748
           49.60112225  49.92107575  50.05125419  51.31086834  49.2367479
           47.78022149  53.78732955  48.50866727  45.08745643  51.41199449]
         [  4.95018614   2.49206539   2.07835604   4.67308065   3.96909729
            5.72365167   1.27126105   5.2597637    5.0424576    1.96953563
            4.02458671   3.35527691   4.48947238   2.88329132   3.67296928
            1.57887891   1.17162523   3.20865642   2.50457207   5.93753608]]


        History
        -------
        Written, MC, May 2012 - combination of Matlab routines of Budiman (2003)
                                and Python routines of Flavio Codeco Coelho (2008)
    """
    #
    # Check input
    if not isinstance(dist,(list,tuple)):
        dist  = [dist]
        param = [param]
    else:
        assert len(dist) == len(param)
    ndist = len(dist)
    
    # LHS
    ran    = np.random.uniform(0.,1.,(ndist,nsample))
    lhsout = np.empty((ndist,nsample))
    for j,d in enumerate(dist):
        if not isinstance(d, (stats.rv_discrete,stats.rv_continuous)):
            raise TypeError('dist is not a scipy.stats distribution object.')
        # force type to float for sage compatibility
        pars = tuple([np.float(k) for k in param[j]])
        idx = np.array(np.random.permutation(nsample), dtype=np.float)
        p   = (idx+ran[j,:])/np.float(nsample) # probability of cdf
        lhsout[j,:] = d(*pars).ppf(p)          # inverse of cdf

    return np.squeeze(lhsout)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # import matplotlib.pyplot as plt
    # dist = [stats.norm, stats.uniform]
    # pars = [(50,2),(1,5)]
    # c    = lhs(dist, pars, 20000)
    # plt.figure()
    # plt.hist(c[0,:])
    # plt.figure()
    # plt.hist(c[1,:])
    
    # dist = [stats.uniform, stats.uniform]
    # pars = [(50,2),(1,5)]
    # c    = lhs(dist, pars, 20000)
    # plt.figure()
    # plt.plot(c[0,:],c[1,:],'ko',markersize=1.0)
    # plt.show()
