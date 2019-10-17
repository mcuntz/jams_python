#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def pi(s=None, m=None, norm=None, b=False, evalues=False, ematrix=False):
    """
        Calculates the parameter importance index PI or alternatively the B index.
        The indices can be normalised and the normalisation can be output as well,
        e.g. to calculate PI-tilde aferwards.


        Definition
        ----------
        def pi(s=s, m=m, norm=None, b=False, evalues=False, ematrix=False):


        Optional Input
        --------------
        Either
          s      sensitivity matrix[nparams,n]
        or
          m      transpose(s).s matrix[nparams,nparams]
        norm     normalised indexes; also output normalisation
                 'sum':   index /= sum(index)
                 'ev':    index /= sum(eigenvalues)
                 'evsum': both index /= sum(eigenvalues) then index /= sum(index)
        b        True:  output B index
                 False: output PI index
        evalues  True:  output eigenvalues
        ematrix  True:  output eigenmatrix


        Output
        ------
        Indices and normalisation, optionally eigenvalues and eigenmatrix


        Restrictions
        ------------
        Be aware that s is (Nparams x N) matrix and Python is column-major so that
        m is s.transpose(s) (not transpose(s).s).


        References
        ----------
        Goehler M, J Mai, and M Cuntz (2013)
            Use of eigendecomposition in a parameter sensitivity analysis of the Community Land Model,
            J Geophys Res 188, 904-921, doi:10.1002/jgrg.20072


        Examples
        --------
        >>> import numpy as np
        >>> s = np.array([[ 0.000630519,  0.000699618,  0.001640544,  0.000429059,  0.000460696],
        ...           [ 0.000181144,  0.000200996,  0.000477661,  0.000109410,  0.000117478],
        ...           [-0.000048173, -0.000053453, -0.000126809, -0.000031822, -0.000034169],
        ...           [ 0.000032940,  0.000036550, -0.000170568,  0.000021916,  0.000023532],
        ...           [ 0.000000000,  0.000101740,  0.000000000,  0.000000000,  0.000068367],
        ...           [-0.000020657, -0.000136825, -0.000054175, -0.000014864, -0.000092501],
        ...           [ 0.000025603,  0.000722743,  0.000067146,  0.000018423,  0.000486356],
        ...           [-0.000006285, -0.000191060, -0.000016484, -0.000004523, -0.000128557],
        ...           [-0.000099740, -0.000110671, -0.000003137, -0.000081818, -0.000087852],
        ...           [ 0.000000015,  0.000000017, -0.000000186,  0.000000006,  0.000000006]])
        >>> pi1 = pi(s)
        >>> from autostring import astr
        >>> print(astr(pi1[0:8],3,pp=True))
        ['4.407e-06' '1.270e-06' '3.390e-07' '3.135e-07' '2.170e-07' '3.793e-07' '1.650e-06' '4.342e-07']

        >>> pi2 = pi(s=s)
        >>> print(astr(pi2[0:8],3,pp=True))
        ['4.407e-06' '1.270e-06' '3.390e-07' '3.135e-07' '2.170e-07' '3.793e-07' '1.650e-06' '4.342e-07']

        >>> pi3 = pi(m=np.dot(s,np.transpose(s)))
        >>> print(astr(pi3[0:8],3,pp=True))
        ['4.407e-06' '1.270e-06' '3.390e-07' '3.135e-07' '2.170e-07' '3.793e-07' '1.650e-06' '4.342e-07']

        >>> pi4, ii4 = pi(s, norm='ev')
        >>> print(astr(pi4[0:8],3,pp=True))
        ['0.838' '0.242' '0.064' '0.060' '0.041' '0.072' '0.314' '0.083']
        >>> print(astr(ii4,3,pp=True))
        5.258e-06

        >>> pi5, ii5 = pi(s, norm='sum')
        >>> print(astr(pi5[0:8],3,pp=True))
        ['0.471' '0.136' '0.036' '0.033' '0.023' '0.041' '0.176' '0.046']
        >>> print(astr(ii5,3,pp=True))
        9.364e-06

        >>> pi6, ii6 = pi(s, norm='evsum')
        >>> print(astr(pi6[0:8],3,pp=True))
        ['0.471' '0.136' '0.036' '0.033' '0.023' '0.041' '0.176' '0.046']
        >>> print(astr(ii6,3,pp=True))
        9.364e-06

        >>> b1 = pi(s, b=True)
        >>> print(astr(b1[0:8],3,pp=True))
        ['3.975e-06' '3.271e-07' '2.344e-08' '3.255e-08' '1.503e-08' '3.086e-08' '7.644e-07' '5.336e-08']

        >>> b2, ii2 = pi(s, b=True, norm='ev')
        >>> print(astr(b2[0:8],3,pp=True))
        ['0.756' '0.062' '0.004' '0.006' '0.003' '0.006' '0.145' '0.010']

        >>> b3, ii3 = pi(s, b=True, norm='sum')
        >>> print(astr(b3[0:8],3,pp=True))
        ['0.756' '0.062' '0.004' '0.006' '0.003' '0.006' '0.145' '0.010']

        >>> b4, ii4 = pi(s, b=True, norm='evsum')
        >>> print(astr(b4[0:8],3,pp=True))
        ['0.756' '0.062' '0.004' '0.006' '0.003' '0.006' '0.145' '0.010']


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    # Check input
    assert not ((s is None) and (m is None)), 'No input matrix given.'
    assert not ((s is not None) and (m is not None)), 'Only s or m matrix possible.'
    if (norm is not None):
        if norm.lower() not in ['sum','ev','evsum']:
            raise ValueError("norm must be None, 'sum', 'ev', or 'evsum'.")
    #
    if (s is not None):
      # check if s is masked array
      if type(s) == type(np.ones(1)):
         m  = np.dot(s,np.transpose(s))
      else:
         m  = np.ma.dot(s,np.transpose(s))
    mm  = m.shape
    assert mm[0] == mm[1], 'm matrix not NxN.'
    inn = mm[0]
    #
    ind = np.empty(inn)
    if b:
        for i in range(inn):
            ind[i] = m[i,i]
    else:
        import scipy.linalg as la
        ev, em = la.eigh(m)

        for i in range(inn):
            ind[i] = np.sum(np.abs(ev[:]) * np.abs(em[i,:]))

    if norm is not None:
        if norm.lower() == 'sum':
            inorm = np.sum(ind)
            ind *= 1./inorm
        elif norm.lower() == 'ev':
            if b:
                inorm = np.trace(m)
            else:
                inorm = np.sum(np.abs(ev[:]))
            ind *= 1./inorm
        elif norm.lower() == 'evsum':
            if b:
                inorm = np.trace(m)
            else:
                inorm = np.sum(np.abs(ev[:]))
            ind *= 1./inorm
            inorm2 = np.sum(ind)
            ind *= 1./inorm2
            inorm *= inorm2
        out = [ind, inorm]
        if evalues:
            out += [ev]
        if ematrix:
            out += [em]
    else:
        if evalues | ematrix:
            out = [ind]
            if evalues:
                out += [ev]
            if ematrix:
                out += [em]
        else:
            out = ind
    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    # s = np.array([[ 0.000630519,  0.000699618,  0.001640544,  0.000429059,  0.000460696],
    #               [ 0.000181144,  0.000200996,  0.000477661,  0.000109410,  0.000117478],
    #               [-0.000048173, -0.000053453, -0.000126809, -0.000031822, -0.000034169],
    #               [ 0.000032940,  0.000036550, -0.000170568,  0.000021916,  0.000023532],
    #               [ 0.000000000,  0.000101740,  0.000000000,  0.000000000,  0.000068367],
    #               [-0.000020657, -0.000136825, -0.000054175, -0.000014864, -0.000092501],
    #               [ 0.000025603,  0.000722743,  0.000067146,  0.000018423,  0.000486356],
    #               [-0.000006285, -0.000191060, -0.000016484, -0.000004523, -0.000128557],
    #               [-0.000099740, -0.000110671, -0.000003137, -0.000081818, -0.000087852],
    #               [ 0.000000015,  0.000000017, -0.000000186,  0.000000006,  0.000000006]])
    # pi1 = pi(s)
    # print pi1
    # pi2 = pi(s=s)
    # print pi2
    # pi3 = pi(m=np.dot(s,np.transpose(s)))
    # print pi3
    # pi4, ii4 = pi(s, norm='ev')
    # print pi4, ii4
    # pi5, ii5 = pi(s, norm='sum')
    # print pi5, ii5
    # pi6, ii6 = pi(s, norm='evsum')
    # print pi6, ii6
    # b1 = pi(s, b=True)
    # print b1
    # b2, ii2 = pi(s, b=True, norm='ev')
    # print b2
    # b3, ii3 = pi(s, b=True, norm='sum')
    # print b3
    # b4, ii4 = pi(s, b=True, norm='evsum')
    # print b4

