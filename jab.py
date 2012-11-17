#!/usr/bin/env python
import numpy as np

def jab(arr, ind, nind=None, weight=False, nsteps=1):
    """
        Jackknife-after-Bootstrap error, an approximation of the Jackknife error.


        Definition
        ----------
        def jab(arr, ind, nind=None, weight=False, nsteps=1):


        Input
        -----
        arr       1D/2D-array of bootstrap standard error outputs.
                  1st dim is number of bootstraps. Possible 2nd dim is number of bootstrap outputs.
        
        ind       2D-array of bootstrap indices.
                  1st dim is number of bootstraps. 2nd dim is number of samples per bootstrap.

        Optional Input
        --------------
        nind      # of data points, i.e. maximal jackknife index.
                  Default: ind.shape[1], which assumes that each bootstraps took N samples from N data points.
        weight    if True: do weighted estimate after Wang et al. (1997)
                  Default: False
        nsteps    Give nsteps estimates of JAB error, i.e. after each nboot/nsteps bootstraps.
                  Default: 1

        
        Output
        ------
        Array of Jacknife-after-boostrap errors.
        This should converge towards the full Jackknif error after a large number of bootstraps.


        Restrictions
        ------------
        None


        References
        ----------
        Efron B. (1992)
            Jackknife-after-bootstrap standard errors and influence functions
            J. Royal. Statist. Soc. B 54, 83-127
        Wang et al. (1997)
            Weighted jackknife-after-bootstrap: A heuristic approach
            In Andradottir et al. (Eds.) Proceedings of the 29th conference on Winter simulation,
            7-10 Dec 1997, Atlanta, Georgia, IEEE Computer Society Washington, DC, USA.,
            pp. 240-245, doi:10.1145/268437.268486


        Example
        --------
        >>> import numpy as np
        >>> # 5 Bootstrap with 2 standard error outputs, e.g. for 2 parameters
        >>> Bstat = np.array([[ 3.625,  5.375],\
                              [ 3.625,  4.125],\
                              [ 4.625,  3.75 ],\
                              [ 4.   ,  4.5  ],\
                              [ 4.25 ,  3.625]])
        >>> # Indices of the bootstraps
        >>> indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],\
                                [1, 2, 3, 0, 5, 5, 4, 1],\
                                [3, 0, 3, 5, 4, 4, 1, 7],\
                                [0, 3, 5, 2, 2, 1, 4, 4],\
                                [4, 4, 5, 5, 3, 3, 2, 1]])

        >>> # Normal JAB
        >>> print jab(Bstat, indices)
        [ 0.41639433  0.683304  ]

        >>> # Weighted JAB
        >>> print jab(Bstat, indices, weight=True)
        [ 0.17766158  0.29154304]

        >>> # Normal JAB but bootstrap was done by taking 6 samples out of 10 data points
        >>> print jab(Bstat, indices, nind=10)
        [ 0.52922539  0.87936986]

        >>> # Weighted JAB with 6 out of 10 bootstrap
        >>> print jab(Bstat[:,0], indices, nind=10)
        0.529225385911

        >>> # Normal JAB, 2 steps
        >>> print jab(Bstat, indices, nsteps=5)
        [[ 0.          0.          0.57282196  0.44779553  0.41639433]
         [ 0.          0.77951196  0.67549895  0.58545341  0.683304  ]]


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

        Copyright 2012 Maren Goehler, Matthias Cuntz


        History
        -------
        Written,  MG, Aug 2012
        Modified, MC, Nov 2012 - rewrite
    """

    # Check input
    if len(ind.shape) != 2:
        raise ValueError('Index array must have 2 dimensions.')
    onlyone = False
    if len(arr.shape) != 2:
        arr     = arr[:,np.newaxis]
        onlyone = True
    if ind.shape[0] != arr.shape[0]:
        raise ValueError('First dimension of boostrap standard error array and index array have to match.')
    nboot = ind.shape[0]
    if nind != None:
        if nind < ind.shape[1]:
            print "Warning: nind < ind.shape[1]"
    else:
        nind = ind.shape[1]
    nout = arr.shape[1]

    # Produce mask with True where data point was used in bootstrap
    mask = np.zeros((nboot,nind), dtype=np.bool)
    for i in xrange(nboot):
        mask[i,ind[i,:]] = True

    seb_i      = np.ma.empty((nind,nout,nsteps))
    seb_i.mask = np.zeros(seb_i.shape, np.bool)  # make mask array
    b_i        = np.empty((nind,nsteps))         # number of indices
    step       = nboot//nsteps
    for k in xrange(nsteps):
        iarr    = arr[0:(k+1)*step,:]
        imask   = mask[0:(k+1)*step,:]
        for i in xrange(nind):
            # Search all rows where i-th column is False, i.e. were not used in bootstrap
            isnoti = np.where(imask[:,i]==False)[0]
            # std and number of rows for possible weighting
            b_i[i,k] = isnoti.size
            if isnoti.size == 0:
                seb_i.mask[i,:,k] = True
            elif isnoti.size == 1:
                seb_i[i,:,k] = 0.
            else:
                seb_i[i,:,k] = np.ma.std(iarr[isnoti,:],axis=0)

    # proposed empirical weight of Wang et al. (1997)
    if weight:
        w_i    = b_i.astype(np.float) / (np.float(np.sum(b_i,axis=0))/np.float(nind) + np.float(nind))
        seb_i *= w_i[:,np.newaxis,:]
    # jackknife-after-bootstrap error
    se_jab = np.sqrt(nind-1) * np.ma.std(seb_i,axis=0)

    if onlyone:
        if nsteps==1:
            return se_jab[0,0]
        else:
            return se_jab[0,:]
    else:
        if nsteps==1:
            return se_jab[:,0]
        else:
            return se_jab

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # # 5 Bootstrap with 2 standard error outputs, e.g. for 2 parameters
    # Bstat = np.array([[ 3.625,  5.375],\
    #                   [ 3.625,  4.125],\
    #                   [ 4.625,  3.75 ],\
    #                   [ 4.   ,  4.5  ],\
    #                   [ 4.25 ,  3.625]])
    # # Indices of the bootstraps
    # indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],\
    #                     [1, 2, 3, 0, 5, 5, 4, 1],\
    #                     [3, 0, 3, 5, 4, 4, 1, 7],\
    #                     [0, 3, 5, 2, 2, 1, 4, 4],\
    #                     [4, 4, 5, 5, 3, 3, 2, 1]])

    # # Normal JAB
    # print jab(Bstat, indices)
    # #[ 0.41639433  0.683304  ]

    # # Weighted JAB
    # print jab(Bstat, indices, weight=True)
    # #[ 0.17766158  0.29154304]

    # # Normal JAB but bootstrap was done by taking 6 samples out of 10 data points
    # print jab(Bstat, indices, nind=10)
    # #[ 0.52922539  0.87936986]

    # # Weighted JAB with 6 out of 10 bootstrap
    # print jab(Bstat[:,0], indices, nind=10)
    # #0.529225385911

    # # Normal JAB, 2 steps
    # print jab(Bstat, indices, nsteps=5)
    # #[[ 0.          0.          0.57282196  0.44779553  0.41639433]
    # # [ 0.          0.77951196  0.67549895  0.58545341  0.683304  ]]
