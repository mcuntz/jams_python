#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
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

        >>> np.random.seed(1)
        >>> dist = [stats.norm]
        >>> pars = [(50,2)]
        >>> c    = lhs(dist, pars, 20)
        >>> print(c.shape)
        (1, 20)
        >>> print(astr(c[0,0:4],3,pp=True))
        ['51.171' '48.562' '51.683' '50.585']

        >>> np.random.seed(1)
        >>> dist = stats.norm
        >>> pars = (50,2)
        >>> c    = lhs(dist, pars, 20)
        >>> print(c.shape)
        (20,)
        >>> print(astr(c[0:4],3,pp=True))
        ['51.171' '48.562' '51.683' '50.585']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2015 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012 - combination of Matlab routines of Budiman (2003)
                                 and Python routines of Flavio Codeco Coelho (2008)
        Modified, MC, Feb 2013 - ported to Python 3
                      Nov 2016 - preserve shape <- nodim
    """
    #
    # Check input
    if not isinstance(dist,(list,tuple)):
        nodim = True
        dist  = [dist]
        param = [param]
    else:
        nodim = False
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

    if nodim:
        return lhsout[0,:]
    else:
        return lhsout


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

