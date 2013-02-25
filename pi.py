#!/usr/bin/env python
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
          s       sensitivity matrix[nparams,n]
        or
          m       transpose(s).s matrix[nparams,nparams]
        norm     normalise incides; also output normalisation
                 'sum':   index /= sum(index)
                 'ev':    index /= sum(eigenvectors)
                 'evsum': both index /= sum(eigenvectors) then index /= sum(index)
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
        Goehler et al. (2013). Use of eigendecomposition in sensitivity analsyis of a land surface model,
            J Geophys Res.


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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
    """
    # Check input
    if (s==None) & (m==None):
        raise ValueError('No input matrix given.')
    if (s!=None) & (m!=None):
        raise ValueError('Only s or m matrix possible.')
    if (norm!=None):
        if norm.lower() not in ['sum','ev','evsum']:
            raise ValueError("norm must be None, 'sum', 'ev', or 'evsum'.")
    #
    if (s!=None):
      # check if s is masked array
      if type(s) == type(np.ones(1)):
         m  = np.dot(s,np.transpose(s))
      else:
         m  = np.ma.dot(s,np.transpose(s))
    mm  = m.shape
    if (mm[0] != mm[1]):
        raise ValueError('m matrix not NxN.')
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

    if norm != None:
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
    doctest.testmod()
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

