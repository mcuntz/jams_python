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
        Eigenvalues are trasnformed to positive and the correpondings eigenvectors are mirrored.


        References
        ----------
        Goehler et al. (2013). Use of eigendecomposition in sensitivity analsyis of a land surface model,
            J Geophys Res.


        Examples
        --------
        >>> import numpy as np
        >>> s = np.array([[ 0.000630519,  0.000699618,  0.001640544,  0.000429059,  0.000460696], \
                      [ 0.000181144,  0.000200996,  0.000477661,  0.000109410,  0.000117478], \
                      [-0.000048173, -0.000053453, -0.000126809, -0.000031822, -0.000034169], \
                      [ 0.000032940,  0.000036550, -0.000170568,  0.000021916,  0.000023532], \
                      [ 0.000000000,  0.000101740,  0.000000000,  0.000000000,  0.000068367], \
                      [-0.000020657, -0.000136825, -0.000054175, -0.000014864, -0.000092501], \
                      [ 0.000025603,  0.000722743,  0.000067146,  0.000018423,  0.000486356], \
                      [-0.000006285, -0.000191060, -0.000016484, -0.000004523, -0.000128557], \
                      [-0.000099740, -0.000110671, -0.000003137, -0.000081818, -0.000087852], \
                      [ 0.000000015,  0.000000017, -0.000000186,  0.000000006,  0.000000006]])
        >>> pi1 = pi(s)
        >>> print pi1
        [  4.40722045e-06   1.27021008e-06   3.39027139e-07   3.13508021e-07
           2.17042558e-07   3.79277689e-07   1.65009950e-06   4.34168308e-07
           3.52585373e-07   3.75045831e-10]

        >>> pi2 = pi(s=s)
        >>> print pi2
        [  4.40722045e-06   1.27021008e-06   3.39027139e-07   3.13508021e-07
           2.17042558e-07   3.79277689e-07   1.65009950e-06   4.34168308e-07
           3.52585373e-07   3.75045831e-10]

        >>> pi3 = pi(m=np.dot(s,np.transpose(s)))
        >>> print pi3
        [  4.40722045e-06   1.27021008e-06   3.39027139e-07   3.13508021e-07
           2.17042558e-07   3.79277689e-07   1.65009950e-06   4.34168308e-07
           3.52585373e-07   3.75045831e-10]

        >>> pi4, ii4 = pi(s, norm='ev')
        >>> print pi4, ii4
        [  8.38171515e-01   2.41570378e-01   6.44766681e-02   5.96234056e-02
           4.12774653e-02   7.21315756e-02   3.13818293e-01   8.25707523e-02
           6.70551926e-02   7.13267547e-05] 5.25813675905e-06

        >>> pi5, ii5 = pi(s, norm='sum')
        >>> print pi5, ii5
        [  4.70680171e-01   1.35655274e-01   3.62072543e-02   3.34818760e-02
           2.31796048e-02   4.05059129e-02   1.76226518e-01   4.63680943e-02
           3.76552400e-02   4.00539610e-05] 9.36351416965e-06

        >>> pi6, ii6 = pi(s, norm='evsum')
        >>> print pi6, ii6
        [  4.70680171e-01   1.35655274e-01   3.62072543e-02   3.34818760e-02
           2.31796048e-02   4.05059129e-02   1.76226518e-01   4.63680943e-02
           3.76552400e-02   4.00539610e-05] 9.36351416965e-06

        >>> b1 = pi(s, b=True)
        >>> print b1
        [  3.97473660e-06   3.27144200e-07   2.34385439e-08   3.25484548e-08
           1.50250743e-08   3.08600964e-08   7.64403109e-07   5.33625069e-08
           3.66181376e-08   3.51820000e-14]

        >>> b2, ii2 = pi(s, b=True, norm='ev')
        >>> print b2
        [  7.55921115e-01   6.22167538e-02   4.45757593e-03   6.19011188e-03
           2.85749021e-03   5.86901745e-03   1.45375281e-01   1.01485582e-02
           6.96409000e-03   6.69096328e-09]

        >>> b3, ii3 = pi(s, b=True, norm='sum')
        >>> print b3
        [  7.55921115e-01   6.22167538e-02   4.45757593e-03   6.19011188e-03
           2.85749021e-03   5.86901745e-03   1.45375281e-01   1.01485582e-02
           6.96409000e-03   6.69096328e-09]

        >>> b4, ii4 = pi(s, b=True, norm='evsum')
        >>> print b4
        [  7.55921115e-01   6.22167538e-02   4.45757593e-03   6.19011188e-03
           2.85749021e-03   5.86901745e-03   1.45375281e-01   1.01485582e-02
           6.96409000e-03   6.69096328e-09]


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
        for i in xrange(inn):
            ind[i] = m[i,i]
    else:
        import scipy.linalg as la
        ev, em = la.eigh(m)

        if np.any(ev < 0.):
            ii = np.squeeze(np.where(ev < 0.))
            ev[ii]   *= -1.
            em[:,ii] *= -1.
            ii = np.argsort(ev)
            ev = ev[ii]
            em = em[:,ii]
        
        for i in xrange(inn):
            ind[i] = np.sum(ev[:] * np.abs(em[i,:]))

    if norm != None:
        if norm.lower() == 'sum':
            inorm = np.sum(ind)
            ind *= 1./inorm
        elif norm.lower() == 'ev':
            inorm = np.trace(m)
            ind *= 1./inorm
        elif norm.lower() == 'evsum':
            inorm = np.trace(m)
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
    # s = np.array([[ 0.000630519,  0.000699618,  0.001640544,  0.000429059,  0.000460696], \
    #               [ 0.000181144,  0.000200996,  0.000477661,  0.000109410,  0.000117478], \
    #               [-0.000048173, -0.000053453, -0.000126809, -0.000031822, -0.000034169], \
    #               [ 0.000032940,  0.000036550, -0.000170568,  0.000021916,  0.000023532], \
    #               [ 0.000000000,  0.000101740,  0.000000000,  0.000000000,  0.000068367], \
    #               [-0.000020657, -0.000136825, -0.000054175, -0.000014864, -0.000092501], \
    #               [ 0.000025603,  0.000722743,  0.000067146,  0.000018423,  0.000486356], \
    #               [-0.000006285, -0.000191060, -0.000016484, -0.000004523, -0.000128557], \
    #               [-0.000099740, -0.000110671, -0.000003137, -0.000081818, -0.000087852], \
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
