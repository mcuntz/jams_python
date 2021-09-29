#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['jab']


def jab(arr, ind=None, nind=None, mask=None, weight=False, nsteps=1):
    """
    Jackknife-after-Bootstrap error, an approximation of the Jackknife error.

    If you one gives the indices of the bootstrap, a mask will be produced that
    tells if this sample was used in the current bootstrap:
        mask = np.zeros((nboot,nind), dtype=bool)
        for i in range(nboot):
            mask[i,ind[i,:]] = True
    One can also give this mask directly, reducing the memory amount
    by a factor of 64.


    Definition
    ----------
    def jab(arr, ind=None, nind=None, mask=None, weight=False, nsteps=1):


    Input
    -----
    arr       1D/2D-array of bootstrap standard error outputs.
              1st dim is number of bootstraps. Possible 2nd dim is number
              of bootstrap outputs.


    Optional Input
    --------------
    ind       2D-array of bootstrap indices. mask has priority over ind.
              1st dim is number of bootstraps. 2nd dim is number of samples
              per bootstrap.
    nind      # of data points, i.e. maximal jackknife index.
              Default: ind.shape[1], which assumes that each bootstraps took
              N samples from N data points.
              Ignored if mask is given.
    mask      2D-array of bollean if index used in this bootstrap.
              mask has priority over ind.
              1st dim is number of bootstraps. 2nd dim is number of samples
              per bootstrap.
    weight    if True: do weighted estimate after Wang et al. (1997)
              Default: False
    nsteps    Give nsteps estimates of JAB error, i.e. after each
              nboot/nsteps bootstraps.
              Default: 1


    Output
    ------
    Array of Jacknife-after-boostrap errors.
    This should converge towards the full Jackknif error after a large number
    of bootstraps.


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
        In Andradottir et al. (Eds.) Proceedings of the 29th conference on
        Winter simulation, 7-10 Dec 1997, Atlanta, Georgia, IEEE Computer
        Society Washington, DC, USA., pp. 240-245, doi:10.1145/268437.268486


    Example
    --------
    # 5 Bootstrap with 2 standard error outputs, e.g. for 2 parameters
    >>> import numpy as np
    >>> Bstat = np.array([[ 3.625,  5.375],
    ...                   [ 3.625,  4.125],
    ...                   [ 4.625,  3.75 ],
    ...                   [ 4.   ,  4.5  ],
    ...                   [ 4.25 ,  3.625]])

    # Indices of the bootstraps
    >>> indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],
    ...                     [1, 2, 3, 0, 5, 5, 4, 1],
    ...                     [3, 0, 3, 5, 4, 4, 1, 7],
    ...                     [0, 3, 5, 2, 2, 1, 4, 4],
    ...                     [4, 4, 5, 5, 3, 3, 2, 1]])

    # Normal JAB
    >>> from autostring import astr
    >>> print(astr(jab(Bstat, indices),3,pp=True))
    ['0.131' '0.393']

    # Normal JAB giving weights
    >>> mask = np.zeros((Bstat.shape[0],indices.shape[1]), dtype=bool)
    >>> for i in range(Bstat.shape[0]):
    ...     mask[i,indices[i,:]] = True
    >>> print(astr(jab(Bstat, mask=mask),3,pp=True))
    ['0.131' '0.393']

    # Weighted JAB
    >>> print(astr(jab(Bstat, indices, weight=True),3,pp=True))
    ['0.056' '0.168']

    # Normal JAB but bootstrap was done by taking 6 samples out of
    # 10 data points
    >>> print(astr(jab(Bstat, indices, nind=10),3,pp=True))
    ['0.146' '0.378']

    # Weighted JAB with 6 out of 10 bootstrap
    >>> print(astr(jab(Bstat[:,0], indices, nind=10),3,pp=True))
    0.146

    # Normal JAB, 2 steps
    >>> print(astr(jab(Bstat, indices, nsteps=5),3,pp=True))
    [['--   ' '0.000' '0.661' '0.312' '0.131']
     ['--   ' '0.000' '0.579' '0.288' '0.393']]


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python
    library, Department of Computational Hydrosystems, Helmholtz Centre for
    Environmental Research - UFZ, Leipzig, Germany.

    Copyright (c) 2012-2021 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Maren Goehler,  Aug 2012
    Modified, Matthias Cuntz, Nov 2012 - complete rewrite
              Matthias Cuntz, Dec 2012 - mask
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Apr 2014 - assert
              Matthias Cuntz, Sep 2021 - code refactoring
    """
    # Shapes: make nboot bootstraps with nind samples per bootstrap
    # for nout parameters
    # arr[nboot, nout]
    # ind[nboot, nind]
    # mask[nboot, nind]

    # Check input
    onlyone = False
    if len(arr.shape) != 2:
        arr     = arr[:, np.newaxis]
        onlyone = True
    if mask is not None:
        assert len(mask.shape) == 2, 'Mask array must have 2 dimensions.'
        assert mask.shape[0] == arr.shape[0], 'First dimension of boostrap standard error array and mask array have to match.'
        nboot = mask.shape[0]
        nind  = mask.shape[1]
    else:
        assert len(ind.shape) == 2, 'Index array must have 2 dimensions.'
        assert ind.shape[0] == arr.shape[0], 'First dimension of boostrap standard error array and index array have to match.'
        nboot = ind.shape[0]
        if nind is not None:
            if nind < ind.shape[1]:
                print("Warning: nind < ind.shape[1]")
        else:
            nind = ind.shape[1]
    nout = arr.shape[1]

    if mask is None:
        # Produce mask with True where data point was used in bootstrap
        mask = np.zeros((nboot, nind), dtype=bool)
        for i in range(nboot):
            mask[i, ind[i, :]] = True

    seb_i      = np.ma.empty((nind, nout, nsteps))     # std jab_i
    seb_i.mask = np.zeros(seb_i.shape, dtype=bool)  # make mask an array
    b_i        = np.empty((nind, nsteps))              # number of indices
    step       = nboot // nsteps
    isnoti     = [np.empty((0), dtype=int) for i in range(nind)]
    for k in range(nsteps):
        # Searching the mask only for the current step and appending
        # the indices seems to be marginally faster than always searching
        # the mask from 0 to the current step
        #   imask   = mask[0:(k+1)*step,:]
        # and then later
        #   isnoti = np.where(imask[:,i]==False)[0]
        imask = mask[k*step:(k+1)*step, :]  # current step
        iarr  = arr[0:(k+1)*step, :]        # 1st to current step
        for i in range(nind):
            # Search all rows (in current step) where i-th column is False,
            # i.e. were not used in bootstrap
            iisnoti = k*step + np.where(imask[:, i] == False)[0]
            # Concat with all steps before
            isnoti[i] = np.concatenate((isnoti[i], iisnoti))
            # std and number of rows for possible weighting
            isize     = isnoti[i].size
            b_i[i, k] = isize
            if isize == 0:
                seb_i.mask[i, :, k] = True
            elif isize == 1:
                seb_i.mask[i, :, k] = True
            else:  # Wang et al. (1997), Eq. 3
                seb_i[i, :, k] = np.ma.std(iarr[isnoti[i], :], axis=0, ddof=0)

    # proposed empirical weight of Wang et al. (1997)
    if weight:  # Wang et al. (1997), Eq. 9
        w_i    = b_i.astype(float) / (
            float(np.sum(b_i, axis=0))/float(nind) + float(nind))
        seb_i *= w_i[:, np.newaxis, :]
    # jackknife-after-bootstrap error, Wang et al. (1997), Eq. 4
    se_jab = np.sqrt(float(nind-1)) * np.ma.std(seb_i, axis=0, ddof=0)

    if not np.any(se_jab.mask):
        se_jab = se_jab.data
    if onlyone:
        if nsteps == 1:
            return se_jab[0, 0]
        else:
            return se_jab[0, :]
    else:
        if nsteps == 1:
            return se_jab[:, 0]
        else:
            return se_jab


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # 5 Bootstrap with 2 standard error outputs, e.g. for 2 parameters
    # Bstat = np.array([[ 3.625,  5.375],
    #                   [ 3.625,  4.125],
    #                   [ 4.625,  3.75 ],
    #                   [ 4.   ,  4.5  ],
    #                   [ 4.25 ,  3.625]])
    # # Indices of the bootstraps
    # indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],
    #                     [1, 2, 3, 0, 5, 5, 4, 1],
    #                     [3, 0, 3, 5, 4, 4, 1, 7],
    #                     [0, 3, 5, 2, 2, 1, 4, 4],
    #                     [4, 4, 5, 5, 3, 3, 2, 1]])

    # # Normal JAB
    # print jab(Bstat, indices)
    # #[ 0.13132114  0.39335039]

    # # Normal JAB giving weights
    # mask = np.zeros((Bstat.shape[0],indices.shape[1]), dtype=bool)
    # for i in range(Bstat.shape[0]):
    #     mask[i,indices[i,:]] = True
    # print jab(Bstat, mask=mask)
    # #[ 0.13132114  0.39335039]

    # # Weighted JAB
    # print jab(Bstat, indices, weight=True)
    # #[ 0.05603035  0.1678295 ]

    # # Normal JAB but bootstrap was done by taking 6 samples out of 10 data points
    # print jab(Bstat, indices, nind=10)
    # #[ 0.14620639  0.37764585]

    # # Weighted JAB with 6 out of 10 bootstrap
    # print jab(Bstat[:,0], indices, nind=10)
    # #0.146206392124

    # # Normal JAB, 5 steps
    # print jab(Bstat, indices, nsteps=5)
    # #[[-- 0.0 0.661437827766 0.311804782231 0.131321139382]
    # # [-- 0.0 0.578758099295 0.287799087529 0.393350393185]]
