#!/usr/bin/env python
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
                            This uses the library sobol of Corrado Chisari
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
        >>> params = np.array([[  1.00000000e+00,   1.00000000e+06], \
                               [  1.00000000e+00,   2.00000000e+01]])
        >>> nbase = 10
        >>> out = saltelli(params, nbase)
        >>> print out
        [[  1.00000000e+00   5.00000500e+05   7.50000250e+05   2.50000750e+05
            3.75000625e+05   8.75000125e+05   6.25000375e+05   1.25000875e+05
            1.87500812e+05   6.87500312e+05   1.00000000e+00   5.00000500e+05
            7.50000250e+05   2.50000750e+05   6.25000375e+05   1.25000875e+05
            3.75000625e+05   8.75000125e+05   3.12500688e+05   8.12500188e+05
            1.00000000e+00   5.00000500e+05   7.50000250e+05   2.50000750e+05
            3.75000625e+05   8.75000125e+05   6.25000375e+05   1.25000875e+05
            1.87500812e+05   6.87500312e+05   1.00000000e+00   5.00000500e+05
            7.50000250e+05   2.50000750e+05   6.25000375e+05   1.25000875e+05
            3.75000625e+05   8.75000125e+05   3.12500688e+05   8.12500188e+05]
         [  1.00000000e+00   1.05000000e+01   5.75000000e+00   1.52500000e+01
            8.12500000e+00   1.76250000e+01   3.37500000e+00   1.28750000e+01
            6.93750000e+00   1.64375000e+01   1.00000000e+00   1.05000000e+01
            5.75000000e+00   1.52500000e+01   3.37500000e+00   1.28750000e+01
            8.12500000e+00   1.76250000e+01   1.40625000e+01   4.56250000e+00
            1.00000000e+00   1.05000000e+01   5.75000000e+00   1.52500000e+01
            3.37500000e+00   1.28750000e+01   8.12500000e+00   1.76250000e+01
            1.40625000e+01   4.56250000e+00   1.00000000e+00   1.05000000e+01
            5.75000000e+00   1.52500000e+01   8.12500000e+00   1.76250000e+01
            3.37500000e+00   1.28750000e+01   6.93750000e+00   1.64375000e+01]]

        >>> out = saltelli(params, nbase, nskip=2)
        >>> print out
        [[  5.00000500e+05   7.50000250e+05   2.50000750e+05   3.75000625e+05
            8.75000125e+05   6.25000375e+05   1.25000875e+05   1.87500812e+05
            6.87500312e+05   9.37500062e+05   5.00000500e+05   7.50000250e+05
            2.50000750e+05   6.25000375e+05   1.25000875e+05   3.75000625e+05
            8.75000125e+05   3.12500688e+05   8.12500188e+05   5.62500438e+05
            5.00000500e+05   7.50000250e+05   2.50000750e+05   3.75000625e+05
            8.75000125e+05   6.25000375e+05   1.25000875e+05   1.87500812e+05
            6.87500312e+05   9.37500062e+05   5.00000500e+05   7.50000250e+05
            2.50000750e+05   6.25000375e+05   1.25000875e+05   3.75000625e+05
            8.75000125e+05   3.12500688e+05   8.12500188e+05   5.62500438e+05]
         [  1.05000000e+01   5.75000000e+00   1.52500000e+01   8.12500000e+00
            1.76250000e+01   3.37500000e+00   1.28750000e+01   6.93750000e+00
            1.64375000e+01   2.18750000e+00   1.05000000e+01   5.75000000e+00
            1.52500000e+01   3.37500000e+00   1.28750000e+01   8.12500000e+00
            1.76250000e+01   1.40625000e+01   4.56250000e+00   1.88125000e+01
            1.05000000e+01   5.75000000e+00   1.52500000e+01   3.37500000e+00
            1.28750000e+01   8.12500000e+00   1.76250000e+01   1.40625000e+01
            4.56250000e+00   1.88125000e+01   1.05000000e+01   5.75000000e+00
            1.52500000e+01   8.12500000e+00   1.76250000e+01   3.37500000e+00
            1.28750000e+01   6.93750000e+00   1.64375000e+01   2.18750000e+00]]

        >>> out = saltelli(params, nbase, nskip=2, lhs=True)


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012 Matthias Cuntz


        History
        -------
        Written, MC, May 2012
    """
    #
    # Check input
    if np.size(params[0,:]) != 2:
        raise ValueError('parameter ranges must be given in form (nparams,2).')
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
        from ufz import lhs as latin
        dist = [stats.uniform     for i in xrange(nparams)]
        pars = [(zoff[i],zmul[i]) for i in xrange(nparams)]
        dist = dist + dist # 2*nparams
        pars = pars + pars
        lat  = latin(dist, pars, nbase)
        for i in xrange(nparams):
            pA[i,:] = lat[i,:]
            pB[i,:] = lat[i+nparams,:]
    else:
        import sobol
        sob = sobol.i4_sobol_generate(2*nparams,nbase,nskip)
        for i in xrange(nparams):
            pA[i,:] = zoff[i] + zmul[i]*sob[i,:]
            pB[i,:] = zoff[i] + zmul[i]*sob[i+nparams,:]
    # The C sample is nparams the B sammple
    pC = np.array([pB for i in xrange(nparams)])
    # where on each repeat one column is replaced by the column of A
    for i in xrange(nparams):
        pC[i,i,:] = pA[i,:]

    # Reshape so that one can do runs over 2nd dim and then use ufz.sobol_index
    pout = np.empty((nparams,nso))
    for i in xrange(nparams):
        pout[i,:] = np.concatenate((pA[i,:],pB[i,:],np.ravel(pC[:,i,:])))

    return pout


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # params = np.array([[  1.00000000e+00,   1.00000000e+06],
    #                    [  1.00000000e+00,   2.00000000e+01]])
    # nbase = 10
    # out = saltelli(params, nbase)
    # print out
    # out = saltelli(params, nbase, nskip=2)
    # print out
    # out = saltelli(params, nbase, nskip=2,lhs=True)
    # print out
