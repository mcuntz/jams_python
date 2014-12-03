#!/usr/bin/env python
from __future__ import print_function
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
        >>> dist = [stats.norm, stats.uniform] # for uniform (min, max-min)
        >>> pars = [(50,2),(1,5)]
        >>> c    = lhs(dist, pars, 20)
        >>> from autostring import astr
        >>> print(astr(c[0:2,0:4],3,pp=True))
        [['52.822' '51.956' '46.710' '50.585']
         [' 4.950' ' 2.492' ' 2.078' ' 4.673']]


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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012 - combination of Matlab routines of Budiman (2003)
                                 and Python routines of Flavio Codeco Coelho (2008)
        Modified, MC, Feb 2013 - ported to Python 3
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
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

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

