#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None, si=True, sti=True, mean=False, wmean=False):
    """
        Calculates the first-order Si and total STi variance-based sensitivity indices
        after Saltelli et al. (2008) with improvements of Saltelli (2002).

        Definition
        ----------
        def sobol_index(s=None, ns=None, ya=None, yb=None, yc=None, si=True, sti=True,
                        mean=False, wmean=False):


        Optional Input
        --------------
        Either
          s        ns*(k+2) model outputs
          ns       base sample number
        or
          ya       ns model outputs
          yb       second ns model outputs
          yc       (ns,k) model outputs
        with k: # or parameters. See Saltelli et al. (2008, p. 164ff) for definitions.
        
        si         if True (default): output Si
        sti        if True (default): output STi
        mean       if True: output mean Si and/or STi (if 2D/3D input)
        wmean      if True: output variance weighted mean Si and/or STi (if 2D/3D input)

        If 2D input (3D for yc), then Si/STi will be calculated for each element in first dimension.
        Assumes that first dimension is time, i.e. calculates indices per time step.
        mean and wmean give then the mean and variance weighted mean of the time steps, resp.


        Output
        ------
        First-order sensitivity indices Si and Total sensitivity indices STi.


        Restrictions
        ------------
        The right order and size is assumed if s and ns are given.
        If three positional arguments are give, i.e. would be s, ns and ya, then it is assumed
        that the user gave actually ya, yb and yc.


        References
        ----------
        Saltelli, A. (2002). Making best use of model evaluations to compute sensitivity indices.
            Computer Physics Communications, 145(2), 280-297.
        Saltelli, A. et al. (2008). Global sensitivity analysis. The primer.
            John Wiley & Sons Inc., NJ, USA, ISBN 978-0-470-05997-5 (pp. 1-292)


        Examples
        --------
        >>> import numpy as np
        >>> iya = np.array([-257.0921, -152.8434,  -68.7691, -129.8188,  -89.9667])
        >>> iyb = np.array([-174.3796, -147.3429, -119.2169, -102.7477, -160.2185])
        >>> iyc = np.array([[-101.2691, -139.5885, -127.6342, -106.9594, -115.8792],
        ...                [-130.0901, -154.8311,  -24.3989, -133.4675, -169.0250],
        ...                [-180.9130, -140.7924, -127.3886,  -93.2011, -119.7757],
        ...                [-171.2342, -146.7013, -118.9833, -102.5973, -163.5437],
        ...                [-255.2163, -179.8984, -131.5042, -103.7205, -167.4869],
        ...                [-258.7999, -125.7803, -116.5700, -129.3884, -157.6105],
        ...                [-169.9223, -163.9877, -120.3642,  -89.8831, -160.3463],
        ...                [-173.9842, -144.7646, -119.3019,  -95.3019, -182.1460],
        ...                [-193.2946, -155.0912, -141.4387, -107.3154, -171.4209],
        ...                [-174.3416, -147.3517, -119.2144, -102.7560, -160.2369]])

        # Give 3 positional arguments
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc)
        >>> from autostring import astr
        >>> print(astr(isi1,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isti1,3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

        # Give 3 optional arguments
        >>> isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc)
        >>> print(astr(isi2,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

        # Give 2 positional arguments
        >>> s = np.concatenate((iya, iyb, np.ravel(iyc)))
        >>> ns = iya.size
        >>> isi3, isti3 = sobol_index(s, ns)
        >>> print(astr(isi3,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

        # Give 2 optional arguments
        >>> isi4, isti4 = sobol_index(s=s, ns=ns)
        >>> print(astr(isi4,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

        # 2 optional arguments and no STi output
        >>> isi5, = sobol_index(s=s, ns=ns, si=True, sti=False)
        >>> print(astr(isi5,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

        # 2 optional arguments and no Si output
        >>> isti5, = sobol_index(s=s, ns=ns, si=False, sti=True)
        >>> print(astr(isti5,3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

        # Indices per time step
        # 3 times the same
        >>> iya = np.asarray([iya, iya, iya])
        >>> iyb = np.asarray([iyb, iyb, iyb])
        >>> iyc = np.asarray([iyc, iyc, iyc])

        # standard
        >>> isi1, isti1 = sobol_index(iya,iyb,iyc)
        >>> print(astr(isi1[0,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isi1[2,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isti1[0,:],3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
        >>> print(astr(isti1[2,:],3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

        # Incl. standard mean of all time steps
        >>> isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc, mean=True)
        >>> print(astr(isi2[0,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isi2[2,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(msi2,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(msti2,3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

        # Incl. weighted mean of all time steps
        >>> s  = np.asarray([s,s,s])
        >>> ns = iya.shape[1]
        >>> isi3, isti3, wsi2, wsti2 = sobol_index(s, ns, wmean=True)
        >>> print(astr(isi3[0,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isi3[1,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(wsi2,3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(wsti2,3,pp=True))
        [' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
        >>> isi4, isti4 = sobol_index(s=s, ns=ns)
        >>> print(astr(isi4[0,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
        >>> print(astr(isi4[2,:],3,pp=True))
        ['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']



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
                  MC, Aug 2013 - Si per time step
    """
    # Check input
    if (si==False) & (sti==False):
        raise ValueError('No output chosen: si=False and sti=False.')
    # s, ns or ya, yb, yc
    isone = False
    if ((s != None) & (ns != None) & (ya == None)):
        if s.ndim == 1:
            isone = True
            s = s[np.newaxis,:]
        nsa   = ns
        ntime = s.shape[0]
        ss    = s.shape[1]
        nn    = (ss//nsa)-2
        iyA   = s[:,0:nsa]
        iyB   = s[:,nsa:2*nsa]
        iyC   = np.reshape(s[:,2*nsa:],(ntime,nn,nsa))
    else:
        if not (((s != None) & (ns != None) & (ya != None)) | ((ya != None) & (yb !=None) & (yc !=None))):
            raise ValueError('Either s and ns has to be given or ya, yb and yc.')
        if ((s != None) & (ns != None) & (ya != None)):
            if s.ndim == 1:
                isone = True
                s  = s[np.newaxis,:]
                ns = ns[np.newaxis,:]
                ya = ya[np.newaxis,:,:]
            iyA     = s
            iyB     = ns
            iyC     = ya
        else:
            if ya.ndim == 1:
                isone = True
                ya = ya[np.newaxis,:]
                yb = yb[np.newaxis,:]
                yc = yc[np.newaxis,:,:]
            iyA     = ya
            iyB     = yb
            iyC     = yc
        ntime = iyA.shape[0]
        nsa   = iyA.shape[1]
        if not ((nsa == iyB.shape[1]) & (nsa == iyC.shape[2])):
            raise ValueError('ya and yb must have same size as yc[1].')
        nn = iyC.shape[1]
    fnsa  = 1. / np.float(nsa)
    fnsa1 = 1. / np.float(nsa-1)
    f02    = fnsa*np.sum(iyA*iyB, axis=1)
    iyAiyA = np.mean(iyA**2, axis=1)
    varA   = iyAiyA - f02
    isi  = np.empty((ntime,nn))
    isti = np.empty((ntime,nn))
    for i in range(nn):
        iyCi = iyC[:,i,:]
        iyAiyC = fnsa1*np.sum(iyA*iyCi, axis=1)
        isi[:,i] = (iyAiyC - f02) / varA
        iyBiyC  = fnsa1*np.sum(iyB*iyCi, axis=1)
        isti[:,i] = 1. - (iyBiyC - f02) / varA

    # simple mean
    msi  = np.mean(isi,  axis=0)
    msti = np.mean(isti, axis=0)
    # weighted mean
    varA = varA[:,np.newaxis]
    denom = 1./np.sum(varA, axis=0)
    wsi   = np.sum(isi*varA,  axis=0) * denom
    wsti  = np.sum(isti*varA, axis=0) * denom

    if isone:
        isi  = isi[0,:]
        isti = isti[0,:]
        
    out = []
    if si:  out = out + [isi]
    if sti: out = out + [isti]
    if not isone:
        if mean:
            if si:  out = out + [msi]
            if sti: out = out + [msti]
        if wmean:
            if si:  out = out + [wsi]
            if sti: out = out + [wsti]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # import numpy as np
    # print('1D')
    # iya = np.array([-257.0921, -152.8434,  -68.7691, -129.8188,  -89.9667])
    # iyb = np.array([-174.3796, -147.3429, -119.2169, -102.7477, -160.2185])
    # iyc = np.array([[-101.2691, -139.5885, -127.6342, -106.9594, -115.8792],
    #                [-130.0901, -154.8311,  -24.3989, -133.4675, -169.0250],
    #                [-180.9130, -140.7924, -127.3886,  -93.2011, -119.7757],
    #                [-171.2342, -146.7013, -118.9833, -102.5973, -163.5437],
    #                [-255.2163, -179.8984, -131.5042, -103.7205, -167.4869],
    #                [-258.7999, -125.7803, -116.5700, -129.3884, -157.6105],
    #                [-169.9223, -163.9877, -120.3642,  -89.8831, -160.3463],
    #                [-173.9842, -144.7646, -119.3019,  -95.3019, -182.1460],
    #                [-193.2946, -155.0912, -141.4387, -107.3154, -171.4209],
    #                [-174.3416, -147.3517, -119.2144, -102.7560, -160.2369]])
    # # Give 3 positional arguments
    # isi1, isti1 = sobol_index(iya,iyb,iyc)
    # from autostring import astr
    # print(astr(isi1,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isti1,3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

    # # Give 3 optional arguments
    # isi2, isti2 = sobol_index(ya=iya,yb=iyb,yc=iyc)
    # print(astr(isi2,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

    # # Give 2 positional arguments
    # s = np.concatenate((iya, iyb, np.ravel(iyc)))
    # ns = iya.size
    # isi3, isti3 = sobol_index(s, ns)
    # print(astr(isi3,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

    # # Give 2 optional arguments
    # isi4, isti4 = sobol_index(s=s, ns=ns)
    # print(astr(isi4,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

    # # 2 optional arguments and no STi output
    # isi5, = sobol_index(s=s, ns=ns, si=True, sti=False)
    # print(astr(isi5,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']

    # # 2 optional arguments and no Si output
    # isti5, = sobol_index(s=s, ns=ns, si=False, sti=True)
    # print(astr(isti5,3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']

    # print('2D')
    # # 3 time steps the same
    # iya = np.asarray([iya, iya, iya])
    # iyb = np.asarray([iyb, iyb, iyb])
    # iyc = np.asarray([iyc, iyc, iyc])
    # isi1, isti1 = sobol_index(iya,iyb,iyc)
    # print(astr(isi1[0,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isi1[2,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isti1[0,:],3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
    # print(astr(isti1[2,:],3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
    # isi2, isti2, msi2, msti2 = sobol_index(ya=iya,yb=iyb,yc=iyc, mean=True)
    # print(astr(isi2[0,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isi2[2,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(msi2,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(msti2,3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
    # s  = np.asarray([s,s,s])
    # ns = iya.shape[1]
    # isi3, isti3, wsi2, wsti2 = sobol_index(s, ns, wmean=True)
    # print(astr(isi3[0,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isi3[1,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(wsi2,3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(wsti2,3,pp=True))
    # #[' 0.972' ' 0.482' '-0.074' '-0.560' '-2.280' '-1.642' '-0.613' '-0.755' '-1.311' '-0.572']
    # isi4, isti4 = sobol_index(s=s, ns=ns)
    # print(astr(isi4[0,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
    # print(astr(isi4[2,:],3,pp=True))
    # #['-0.172' ' 0.685' ' 1.344' ' 1.581' ' 3.794' ' 3.325' ' 1.617' ' 1.672' ' 2.356' ' 1.631']
