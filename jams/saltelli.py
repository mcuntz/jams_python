#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np

def saltelli(params, nbase, lhs=False, nskip=1):
    """
        Samples parameters in given ranges so they can be used
        for the calculation of Sobol indices as described by
        Saltelli (2002) and Saltelli et al. (2007).


        Definition
        ----------
        def saltelli(params, nbase, lhs=False, nskip=1):


        Input
        -----
        params    (nparams,2) array of parameter ranges
        nbase     base sample size


        Optional Input
        --------------
        lhs       if True: use latin hypercube sampling in uniform distributions
                  if False: use sobol sequence
                            This uses the package sobol of Corrado Chisari
        nskip     number of sobol sequences to skip at the beginning


        Output
        ------
        array with (nparams,nbase*(nparam+2)) with sample values in the right order.


        Restrictions
        ------------
        1. Only uniform latin hypercube sampling.
        2. No parameter dependencies possible.


        References
        ----------
        Saltelli, A. (2002). Making best use of model evaluations to compute sensitivity indices.
            Computer Physics Communications, 145(2), 280-297.
        Saltelli, A. et al. (2008). Global sensitivity analysis. The primer.
            John Wiley & Sons Inc., NJ, USA, ISBN 978-0-470-05997-5 (pp. 1-292)


        Examples
        --------
        >>> import numpy as np
        >>> params = np.array([[  1.00000000e+00,   1.00000000e+06],
        ...                    [  1.00000000e+00,   2.00000000e+01]])
        >>> nbase = 10
        >>> # seed for reproducible results in doctest
        >>> np.random.seed(1)
        >>> out = saltelli(params, nbase)
        >>> from autostring import astr
        >>> print(astr(out[0:2,0:4],3,pp=True))
        [['     1.000' '500000.500' '750000.250' '250000.750']
         ['     1.000' '    10.500' '     5.750' '    15.250']]

        >>> out = saltelli(params, nbase, nskip=2)
        >>> print(astr(out[0:2,0:4],3,pp=True))
        [['500000.500' '750000.250' '250000.750' '375000.625']
         ['    10.500' '     5.750' '    15.250' '     8.125']]

        >>> out = saltelli(params, nbase, nskip=2, lhs=True)
        >>> print(astr(out[0:2,0:4],3,pp=True))
        [['341702.859' '872032.577' '500011.937' '930233.327']
         ['     1.796' '     6.102' '    16.588' '     4.568']]


        License
        -------
        This file is part of the JAMS Python package.

        It is NOT released under the GNU Lesser General Public License, yet.

        If you use this routine, please contact Matthias Cuntz.

        Copyright 2012-2014 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    #
    # Check input
    assert np.size(params[0,:]) == 2, 'parameter ranges must be given in form (nparams,2).'
    nparams = np.size(params[:,0])
    #
    nso = nbase*(nparams+2)
    # A=0:nbase, B=nbase:nbase*2, C=nbase*2:nso
    zoff = params[:,0]
    zmul = params[:,1] - params[:,0]
    # Two random samples A and B
    pA   = np.empty((nparams,nbase))
    pB   = np.empty((nparams,nbase))
    if lhs:
        import scipy.stats as stats
        from jams.lhs import lhs
        dist = [stats.uniform     for i in range(nparams)]
        pars = [(zoff[i],zmul[i]) for i in range(nparams)]
        dist = dist + dist # 2*nparams
        pars = pars + pars
        lat  = lhs(dist, pars, nbase)
        for i in range(nparams):
            pA[i,:] = lat[i,:]
            pB[i,:] = lat[i+nparams,:]
    else:
        import sobol
        sob = sobol.i4_sobol_generate(2*nparams,nbase,nskip)
        for i in range(nparams):
            pA[i,:] = zoff[i] + zmul[i]*sob[i,:]
            pB[i,:] = zoff[i] + zmul[i]*sob[i+nparams,:]
    # The C sample is nparams the B sammple
    pC = np.array([pB for i in range(nparams)])
    # where on each repeat one column is replaced by the column of A
    for i in range(nparams):
        pC[i,i,:] = pA[i,:]

    # Reshape so that one can do runs over 2nd dim and then use jams.sobol_index
    pout = np.empty((nparams,nso))
    for i in range(nparams):
        pout[i,:] = np.concatenate((pA[i,:],pB[i,:],np.ravel(pC[:,i,:])))

    return pout


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    # params = np.array([[  1.00000000e+00,   1.00000000e+06],
    #                    [  1.00000000e+00,   2.00000000e+01]])
    # nbase = 10
    # out = saltelli(params, nbase)
    # print out
    # out = saltelli(params, nbase, nskip=2)
    # print out
    # out = saltelli(params, nbase, nskip=2,lhs=True)
    # print out

