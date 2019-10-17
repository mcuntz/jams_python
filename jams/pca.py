#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['pca']

def pca(mat, corr=False, ndim=None, rvar=None):
    """
        Principal component analysis (PCA) upon the first dimension of an 2D-array.

        If you have an Nd-array or an array with masked or NaN values,
        you may want to use Andrew Dawson's eofs package:
            https://github.com/ajdawson/eofs


        Definition
        ----------
        def pca(mat, corr=False, ndim=None, rvar=None):


        Input
        -----
        mat       2D-array, where the first dimension is assumed to represent time.


        Optional Input
        --------------
        corr      True: PCA on correlation matrix
                  False: PCA on covariance matrix
        ndim      # of principal components to retain (default: None, i.e. all)
        rvar      Percent (0-1) of cumulative variance to retain (default: None, i.e. > 1, i.e. all)
                  This uses numpy.searchsorted, which searches for lesser than (<) rvar
                  and not lesser equal (<=) rvar.
                  So rvar=1 retains all but the last principal component, and
                  rvar<np.amin(eigenvalues) retains no principal component at all.


        Output
        ------
        principal_components, eigenvalues, eigenvectors

        If the input is nxk array and there is no dimensionality reduction (ndim or rvar), then
        principal_components is nxk, eigenvalues k, eigenvectors kxk arrays.

        If the input is nxk array and dimensionality reduction selects m principal components, then
        principal_components is nxm, eigenvalues k, eigenvectors kxm arrays.


        Examples
        --------
        >>> import numpy as np
        >>> A = np.array([ [2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9],
        ...                [2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1] ]).T

        >>> comps, evals, evecs = pca(A)
        >>> Anew = np.dot(comps, evecs.T) + np.mean(A, axis=0)
        >>> print(np.allclose(A, Anew))
        True

        >>> comps, evals, evecs = pca(A, corr=True)
        >>> Anew = np.dot(comps, evecs.T)*np.std(A, axis=0, ddof=1) + np.mean(A, axis=0)
        >>> print(np.allclose(A, Anew))
        True

        >>> comps, evals, evecs = pca(A, corr=True, ndim=1)
        >>> Anew = np.dot(comps, evecs.T)*np.std(A, axis=0, ddof=1) + np.mean(A, axis=0)
        >>> print(np.allclose(A, Anew))
        False
        >>> print(np.allclose(A, Anew, rtol=0.2))
        True

        >>> comps, evals, evecs = pca(A, corr=True, rvar=0.97)
        >>> Anew = np.dot(comps, evecs.T)*np.std(A, axis=0, ddof=1) + np.mean(A, axis=0)
        >>> print(np.allclose(A, Anew, rtol=0.2))
        True


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Nov 2014
    """
    from scipy import linalg
    imat = mat.copy()
    #
    # Check input
    nd = imat.ndim  # save for later
    assert nd > 1, 'array must be at 2D-ND.'
    n, k = imat.shape
    #
    # Covariance or corelation matrix
    n1 = 1./(float(n)-1.)
    if corr:
        imat /= np.std(imat, axis=0, ddof=1)
    imat -= imat.mean(axis=0)
    S  = n1 * np.dot(imat.T, imat)
    #
    # Eigenvectors and eigenvalues with linalg.eigh
    # rather than linlag.eig since S is symmetric.
    # The performance gain is substantial.
    evals, evecs = linalg.eigh(S)
    #
    # Sort eigenvalues in decreasing order
    # Eigenvalues from linalg.eigh are order increasingly
    # while linalg.eig eigenvalues do not have to be sorted.
    # Use argsort to be on the save side.
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    # sort eigenvectors accordingly
    evecs = evecs[:,idx]
    #
    # Select the first m eigenvectors either by number
    # or by percent of explained variance
    if ndim is not None:
        evecs = evecs[:, :ndim]
    if rvar is not None:
        cumevals  = np.cumsum(evals)
        cumevals /= cumevals[-1]
        cut   = np.searchsorted(cumevals, rvar)
        evecs = evecs[:, :cut]
    #
    # Calculate principal components and return
    # (reduced) components, eigenvalues, (reduced) eigenvectors
    return np.dot(evecs.T, imat.T).T, evals, evecs


def check_pca(data, comps, eigenvectors, corr=False, rtol=1e-05, atol=1e-08):
    '''
        Data recovered from PCA against original data.


        Definition
        ----------
        def check_pca(data, comps, eigenvectors, corr=False, rtol=1e-05, atol=1e-08):


        Input
        -----
        def check_pca(data, comps, eigenvectors, corr=False, rtol=1e-05, atol=1e-08):
        data           original 2D data (nxk)
        comps          2D (reduced) principal components (nxm)
        eigenvectors   (Reduced) Eigenvectors (kxm)


        Optional Input
        --------------
        corr      True: PCA on correlation matrix
                  False: PCA on covariance matrix
        rtol      numpy.allclose rtol keyword: relative tolerance.
        atol      numpy.allclose atol keyword: absolute tolerance.


        Output
        ------
        True or False: numpy.allclose(data, data_recovered, rtol=rtol, atol=atol)


        Examples
        --------
        >>> import numpy as np
        >>> A = np.array([ [2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9],
        ...                [2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1] ]).T

        >>> comps, evals, evecs = pca(A)
        >>> check_pca(A, comps, evecs)

        >>> comps, evals, evecs = pca(A, corr=True)
        >>> check_pca(A, comps, evecs, corr=True)

        >>> comps, evals, evecs = pca(A, corr=True, ndim=1)
        >>> check_pca(A, comps, evecs, corr=True, rtol=0.2)

        >>> comps, evals, evecs = pca(A, corr=True, rvar=0.97)
        >>> check_pca(A, comps, evecs, corr=True, rtol=0.2)


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Nov 2014
    '''
    # reconstruction
    if corr:
        data_check = np.dot(comps, eigenvectors.T)*np.std(data, axis=0, ddof=1) + np.mean(data, axis=0)
    else:
        data_check = np.dot(comps, eigenvectors.T) + np.mean(data, axis=0)
    assert np.allclose(data, data_check, rtol=rtol, atol=atol), 'data != data_check'


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from scipy import linalg
    # from scipy.misc import lena
    # import matplotlib.pyplot as plt

    # # http://glowingpython.blogspot.sg/2011/07/principal-component-analysis-with-numpy.html
    # A = np.array([ [2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9],
    #                [2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1] ]).T

    # fig = plt.figure()

    # comps, evals, evecs = pca(A)
    # check_pca(A, comps, evecs)

    # sub = plt.subplot(221)
    # # every eigenvector describes the direction of a principal component
    # m = np.mean(A, axis=0)
    # sub.plot([0, -evecs[0,0]*1]+m[0], [0, -evecs[0,1]*1]+m[1], 'k--')
    # sub.plot([0,  evecs[1,0]*1]+m[0], [0,  evecs[1,1]*1]+m[1], 'k--')
    # sub.plot(A[:,0], A[:,1], 'bo') # the data
    # sub.axis('equal')
    # # transformed data
    # sub = plt.subplot(222)
    # sub.plot(comps[:,0], comps[:,1], 'ro')
    # sub.axis('equal')

    # comps, evals, evecs = pca(A, corr=True)
    # check_pca(A, comps, evecs, corr=True)

    # sub = plt.subplot(223)
    # # every eigenvector describes the direction of a principal component
    # m = np.mean(A, axis=0)
    # sub.plot([0, -evecs[0,0]*1]+m[0], [0, -evecs[0,1]*1]+m[1], 'k--')
    # sub.plot([0,  evecs[1,0]*1]+m[0], [0,  evecs[1,1]*1]+m[1], 'k--')
    # sub.plot(A[:,0], A[:,1], 'bo') # the data
    # sub.axis('equal')
    # # transformed data
    # sub = plt.subplot(224)
    # sub.plot(comps[:,0], comps[:,1], 'ro')
    # sub.axis('equal')

    
    # # http://glowingpython.blogspot.it/2011/07/pca-and-image-compression-with-numpy.html
    # fig = plt.figure()
    # A = lena().astype(np.float)

    # full_pc = np.size(A, axis=1) # numbers of all the principal components
    # i=1
    # dist = []
    # nn = 50
    # for numpc in range(1,full_pc+nn,nn): # 1, 51, 101 ... full_pc
    #     comps, evals, evecs = pca(A, ndim=numpc)
    #     Ar = np.dot(comps, evecs.T) + np.mean(A, axis=0) # image reconstruction
    #     # Difference in Frobenius norm
    #     dist.append(linalg.norm(A-Ar,'fro'))
    #     # Reconstructed fig with less PCs
    #     ax = plt.subplot(((full_pc+nn)//nn//3)+1, 3, i, frame_on=False)
    #     ax.xaxis.set_major_locator(plt.NullLocator()) # remove ticks
    #     ax.yaxis.set_major_locator(plt.NullLocator())
    #     i += 1
    #     ax.imshow(Ar)
    #     ax.set_title('PCs # '+str(numpc)+'/'+str(full_pc))
    #     plt.gray()

    # check_pca(A, comps, evecs)

    # plt.figure()
    # perc = np.cumsum(evals)/np.sum(evals)
    # dist = np.array(dist)/np.amax(dist)
    # plt.plot(range(perc.size), perc, 'b')
    # plt.plot(range(1,full_pc+50,50), dist, 'r')
    # plt.axis([1,full_pc,0,1.1])

    # plt.figure()
    # plt.imshow(Ar)
    # plt.title('last numpc')
    # plt.gray()

    # plt.figure()
    # plt.imshow(A)
    # plt.title('numpc FULL')
    # plt.gray()

    
    # # other example
    # # http://sebastianraschka.com/Articles/2014_pca_step_by_step.html
    
    # plt.show()
