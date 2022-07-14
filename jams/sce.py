#!/usr/bin/env python
"""
Shuffled-Complex-Evolution (SCE) algorithm for function minimization

References
----------
Duan, Sorooshian and Gupta (1992)
  Effective and efficient global optimization for conceptual rainfall-runoff
  models, Water Resour Res 28, 1015-1031, https://doi.org/10.1029/91WR02985

This code is based on a Fortran program of Qingyun Duan (2004), ported to
Python by Stijn Van Hoey (2011). It was taken up, debugged, enhanced and is
maintained by Matthias Cuntz while at Department of Computational Hydrosystems
(CHS), Helmholtz Centre for Environmental Research - UFZ, Leipzig, Germany, and
continued while at Institut National de Recherche pour l'Agriculture,
l'Alimentation et l'Environnement (INRAE), Nancy, France.

Copyright (c) 2004-2020 Qingyun Duan, Stijn Van Hoey, Matthias Cuntz
- mc (at) macu (dot) de
Released under the MIT License.

* Written in Fortran by Q Duan, Sep 2004
* Ported to Python by Stijn Van Hoey, 2011
  https://github.com/stijnvanhoey/Optimization_SCE
* Synchronised with enhanced Fortran version of CHS, Oct 2013, Matthias Cuntz
* Added functionality to call external executable, Nov 2016, Matthias Cuntz
* Treat NaN and Inf in function output, Nov 2016, Matthias Cuntz
* Possibility to exclude (mask) parameters from optimisation,
  Nov 2016, Matthias Cuntz
* Added restart possibility, Nov 2016, Matthias Cuntz
* Return also function value of best parameter set if maxit==True,
  Nov 2016, Matthias Cuntz
* Removed functionality to call external executable, Dec 2017, Matthias Cuntz
* Print out number of function evaluations with printit=1,
  Mar 2018, Matthias Cuntz
* Mask parameters with degenerated ranges, e.g. upper<lower bound,
  Mar 2018, Matthias Cuntz
* Use only masked parameters in calculation of geometric range,
  Mar 2018, Matthias Cuntz
* Removed bug that calculated the size of the complexes using all parameters
  and not only the masked parameters, Mar 2018, Matthias Cuntz
* Fixed bug where masked parameters were always out of bounds,
  Mar 2018, Matthias Cuntz
* Allow scalar bounds, which will be taken for all parameters,
  Mar 2018, Matthias Cuntz
* Define number of parameters by inquiring x0 in case of restart,
  May 2018, Matthias Cuntz
* Removed multiplication with hundred in criter_change,
  regarded a bug compared to Fortran code, May 2018, Matthias Cuntz
* Removed exec command to make restart work with Python 3,
  May 2020, Matthias Cuntz
* Code refactoring, Sep 2021, Matthias Cuntz
* Added keywords args and kwargs to pass to function,
  Apr 2022, Matthias Cuntz
* Different sampling of input parameters, Jul 2022, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

.. autosummary::
    sce
"""
from __future__ import division, absolute_import, print_function
from distutils.util import strtobool
from numpy import savez_compressed
import numpy as np

# ToDo: write tmp/population files (of Fortran code)

__all__ = ['sce']


def _SampleInputMatrix(nrows, bl, bu, rnd, sampling='half-open'):
    """
    Create input parameter matrix (nrows,npars) for
    nrows simulations and npars parameters with bounds bl and bu


    Input
    -----
    nrows
        # of simulations
    bl
        (npars) lower bounds of parameters
    bu
        (npars) upper bounds of parameters
    rnd
        RandomState instance


    Optional Input
    --------------
    sampling
        string or (npars) list of strings of options for sampling
        random numbers. Options can be
        'half-open': same as 'right-half-open' (default)
        'left-half-open': sample random floats in half-open interval (bl, bu]
        'right-half-open': sample random floats in half-open interval [bl, bu)
        'open': sample random floats in open interval (bl, bu)
        'log': sample half-open interval [log(bl), log(bu)),
               which samples better intervals spanning orders or magnitude
               such as bl=10^-9 and bu=10^-4


    Output
    ------
    parameter matrix (nrows, npars) with new parameter samples


    References
    ----------
    Duan, Q., S. Sorooshian, and V. Gupta,
      Effective and efficient global optimization for conceptual
      rainfall-runoff models,
      Water Resour. Res., 28, 1015-1031, 1992.


    History
    -------
    Written, Q. Duan, Sep 2004
    Modified,
        S. Van Hoey 2011 - ported to Python
        Matthias Cuntz, Oct 2013
            - adapted to JAMS package and sync with JAMS Fortran version
        Matthias Cuntz, Mar 2018
            - removed npars from call, get it from lower boundary
        Matthias Cuntz, May 2020 - underscore before function name
        Matthias Cuntz, Jul 2022 - pass RandomState
        Matthias Cuntz, Jul 2022 - sampling

    """
    npars = len(bl)
    if isinstance(sampling, str):
        isampling = [sampling] * npars
    else:
        isampling = sampling
    assert len(isampling) == npars, (
        f'sampling must be string or list of strings'
        f' with {npars} entries')
    x = np.zeros((nrows, npars))
    bound = bu - bl
    for i in range(nrows):
        irnd = rnd.random_sample(npars)
        for j in range(npars):
            opt = isampling[j].lower()
            if opt == 'half-open':
                x[i, j] = bl[j] + irnd[j] * bound[j]
            elif opt == 'left-half-open':
                irnd[j] = 1. - irnd[j]
                x[i, j] = bl[j] + irnd[j] * bound[j]
            elif opt == 'right-half-open':
                x[i, j] = bl[j] + irnd[j] * bound[j]
            elif opt == 'open':
                iirnd = irnd[j]
                while not (iirnd > 0.):
                    iirnd = rnd.random_sample(1)
                x[i, j] = bl[j] + iirnd * bound[j]
            elif opt == 'log':
                # x must be > 0. for ln(x)
                xshift = 0.
                xmul = 1.
                if (bl[j] * bu[j]) < 0.:
                    # bl < 0, bu > 0 -> shift to > 0
                    xshift = np.maximum(np.abs(bl[j]), np.abs(bu[j]))
                elif (bl[j] < 0.) or (bu[j] < 0.):
                    # (bl =< 0 and bu =< 0) -> shift to > 0
                    xmul = -1.
                if (bl[j] * bu[j]) == 0.:
                    if bl[j] == 0.:
                        # (bl == 0 and bu > 0) -> shift to [bu, 2*bu)
                        xshift = bu[j]
                    if bu[j] == 0.:
                        # (bl < 0 and bu == 0) -> shift to [2*bl, bl) < 0.
                        # but xmul==-1
                        xshift = bl[j]  # < 0
                        assert xmul == -1., (
                            'wrong logic in sampling option log')
                lnbl = np.log((bl[j] + xshift) * xmul)
                lnbu = np.log((bu[j] + xshift) * xmul)
                x[i, j] = (np.exp(lnbl +
                                  irnd[j] * (lnbu - lnbl)) -
                           xshift) * xmul
            else:
                raise ValueError(f'unknown sampling option {isampling[j]}.\n'
                                 f'Known samplings are: half-open,'
                                 f' left-half-open, right-half-open, open,'
                                 f' log')

    return x


def _cce(func, s, sf, bl, bu, mask, sampling, icall, maxn, alpha, beta,
         maxit, printit, rnd, args=(), kwargs={}):
    """
    Generate a new point in a simplex


    Input
    -----
    func
        function to minimize
    s                 2D-array
        the sorted simplex in order of increasing function values
    sf                1D-array
        function values in increasing order
    bl                (npars)
        lower bounds of parameters
    bu                (npars)
        upper bounds of parameters
    mask              (npars)
        mask to include (1) or exclude (0) parameter from optimisation
    sampling
        string or (npars) list of strings of options for sampling
        random numbers. Options can be
        'half-open': same as 'right-half-open' (default)
        'left-half-open': sample random floats in half-open interval (bl, bu]
        'right-half-open': sample random floats in half-open interval [bl, bu)
        'open': sample random floats in open interval (bl, bu)
        'log': sample half-open interval [log(bl), log(bu)),
               which samples better intervals spanning orders or magnitude
               such as bl=10^-9 and bu=10^-4
    icall
        counter of function calls
    maxn
        maximum number of function evaluations allowed during optimization
    alpha
        parameter for reflection  of points in complex
    beta
        parameter for contraction of points in complex
    maxit
        if True: maximise instead of minimize func
    printit
        if 1: print each function evaluation
    rnd
        RandomState instance


    Optional Input
    --------------
    args : tuple, optional
        Extra arguments passed to the function *func*
    kwargs : dict, optional
        Extra keyword arguments passed to the function *func*


    Output
    ------
    new parameter set, function value of new set, number of function calls


    References
    ----------
    Duan, Q., S. Sorooshian, and V. Gupta,
      Effective and efficient global optimization for conceptual
      rainfall-runoff models, Water Resour. Res., 28, 1015-1031, 1992.


    History
    -------
    Written,  Q. Duan, Sep 2004
    Modified,
        S. Van Hoey 2011
            - ported to Python
        Matthias Cuntz, Oct 2013
            - adapted to JAMS package and sync with JAMS Fortran version
        Matthias Cuntz, Mar 2018
            - changed indentation of random point last resort
            - fixed masked values that were always out of bounds
        Matthias Cuntz, May 2020 - underscore before function name
        Matthias Cuntz, Apr 2022 - args and kwargs passed to function
        Matthias Cuntz, Jul 2022 - pass RandomState
        Matthias Cuntz, Jul 2022 - sampling

    """

    """
    List of local variables
      sb(.) = the best point of the simplex
      sw(.) = the worst point of the simplex
      w2(.) = the second worst point of the simplex
      fw = function value of the worst point
      ce(.) = the centroid of the simplex excluding wo
      snew(.) = new point generated from the simplex
      iviol = flag indicating if constraints are violated
            = 1 , yes
            = 0 , no
    """
    # Assign the best and worst points:
    sb = s[0, :]
    # fb = sf[0]
    sw = s[-1, :]
    fw = sf[-1]

    # Compute the centroid of the simplex excluding the worst point:
    ce = np.mean(s[:-1, :], axis=0)

    # Attempt a reflection point
    snew = ce + alpha * (ce - sw)
    # sb should have initial params at mask==False
    snew = np.where(mask, snew, sb)

    # Check if is outside the bounds:
    ibound = 0
    # s1 = snew-bl
    # idx = (s1<0).nonzero()
    # if idx[0].size != 0: ibound = 1
    if np.ma.any(np.ma.array(snew - bl, mask=~mask) < 0.):
        ibound = 1

    # s1 = bu-snew
    # idx = (s1<0).nonzero()
    # if idx[0].size != 0: ibound = 2
    if np.ma.any(np.ma.array(bu - snew, mask=~mask) < 0.):
        ibound = 2

    if ibound >= 1:
        snew = _SampleInputMatrix(1, bl, bu, rnd, sampling=sampling)[0]
        snew = np.where(mask, snew, sb)

    fuc = func(snew, *args, **kwargs)
    fnew = -fuc if maxit else fuc
    icall += 1
    if printit == 1:
        print('  i, f, X: ', icall, fnew, snew)

    # Reflection failed; now attempt a contraction point:
    if fnew > fw:
        snew = sw + beta * (ce - sw)
        snew = np.where(mask, snew, sb)
        fuc = func(snew, *args, **kwargs)
        fnew = -fuc if maxit else fuc
        icall += 1
        if printit == 1:
            print('  i, f, X: ', icall, fnew, snew)

    # Both reflection and contraction have failed, attempt a random point;
    if fnew > fw:
        snew = _SampleInputMatrix(1, bl, bu, rnd, sampling=sampling)[0]
        snew = np.where(mask, snew, sb)
        fuc = func(snew, *args, **kwargs)
        fnew = -fuc if maxit else fuc
        icall += 1
        if printit == 1:
            print('  i, f, X: ', icall, fnew, snew)

    # end of _cce
    return snew, fnew, icall


def sce(func, x0, bl, bu,
        mask=None,
        sampling='half-open',
        maxn=1000, kstop=10, pcento=0.0001,
        ngs=2, npg=None, nps=None, nspl=None, mings=None,
        peps=0.001, seed=None, iniflg=True,
        alpha=0.8, beta=0.45, maxit=False, printit=2,
        outf=False, outhist=False, outcall=False,
        restart=False, restartfile1='sce.restart.npz',
        restartfile2='sce.restart.txt',
        args=(), kwargs={}):
    """
    Shuffled-Complex-Evolution algorithm for function minimization

    Parameters
    ----------
    func : callable
        Python function to minimize, callable as `func(x)` with `x` the
        function parameters.
    x0 : array_like
        Parameter values used in initial complex and with `mask==0`.
    bl : array_like
        Lower bounds of parameters.
    bu : array_like
        Upper bounds of parameters.
    mask : array_like, optional
        Include (1,True) or exclude (0,False) parameters in minimization
        (default: include all parameters).

        'nopt = sum(mask)`
    sampling : string or array_like of strings, optional
        Options for sampling random numbers. Options can be:
        'half-open': same as 'right-half-open' (default)

        'left-half-open': sample random floats in half-open interval (bl, bu]

        'right-half-open': sample random floats in half-open interval [bl, bu)

        'open': sample random floats in open interval (bl, bu)

        'log': sample half-open interval [log(bl), log(bu)), which samples
        better intervals spanning orders or magnitude such as
        bl=10^-9 and bu=10^-4
    maxn : int, optional
        Maximum number of function evaluations allowed during minimization
        (default: 1000).
    kstop : int, optional
        Maximum number of evolution loops before convergency (default: 10).
    pcento : float, optional
        Percentage change allowed in kstop loops before convergency
        (default: 0.0001).
    peps : float, optional
        Value of normalised geometric range needed for convergence
        (default: 0.001).
    ngs : int, optional
        Number of complexes (default: 2).
    npg : int, optional
        Number of points in each complex (default: `2*nopt+1`).
    nps : int, optional
        Number of points in each sub-complex (default: `nopt+1`).
    mings : int, optional
        Minimum number of complexes required if the number of complexes is
        allowed to reduce as the optimization proceeds (default: `ngs`).
    nspl : int, optional
        Number of evolution steps allowed for each complex before complex
        shuffling (default: `2*nopt+1`).
    seed : int, optional
        Random number generator seed (default: None).
    iniflg : bool, optional
        If True: include initial parameter in initial population
        (default: True).
    alpha : float, optional
        Parameter for reflection of points in complex (default: 0.8).
    beta : float, optional
        Parameter for contraction of points in complex (default: 0.45).
    maxit : bool, optional
        If True: maximize instead of minimize func (default: False).
    printit : int, optional
        Controlling print-out (default: 2)

        0: print information on the best point of the population

        1: print information on each function evaluation

        2: no printing.
    outf : bool, optional
        If True: return best function value (default: False).
    outhist : bool, optional
        If True: return parameters and function values of each evolution loop
        (default: False).
    outcall : bool, optional
        If True: return number of function evaluations (default: False).
    restart : bool, optional
        If True, continue from saved state in `restartfile1` and `restartfile2`
        (default: False).
    restartfile1 : str, optional
        Filename for saving state of SCE, array variables
        (default: `sce.restart.npz`)

        Note: state will be always written, except if `restartfile1=None`.
    restartfile2 : int, optional
        Filename for saving state of SCE, non-array variables
        (default: `sce.restart.txt`)
    args : tuple, optional
        Extra arguments passed to the function *func*
    kwargs : dict, optional
        Extra keyword arguments passed to the function *func*

    Returns
    -------
    list
       list[0] is array with best parameter set

       list[1:] depends on `outf`, `outhist` and `outcall`.

       If `outf=True`, `outhist=True` and `outcall=True`

        list[1] : function value at best parameter set

        list[2] : vector of best parameters sets per evolution, function values
                  of the best parameter vector

        list[3] : number of function evaluations

    References
    ----------
    Duan, Sorooshian and Gupta (1992)
        Effective and efficient global optimization for conceptual
        rainfall-runoff models, Water Resources Research 28, 1015-1031,
        https://doi.org/10.1029/91WR02985
    Duan, Sorooshian, and Gupta (1994)
        Optimal use of the SCE-UA global optimization method for calibrating
        watershed models, Journal of Hydrology 58, 265–284,
        https://doi.org/10.1016/0022-1694(94)90057-4
    Behrangi, Khakbaz, Vrugt, Duan, and Sorooshian (2008)
        Comment on “Dynamically dimensioned search algorithm for
        computationally efficient watershed model calibration”
        by Bryan A. Tolson and Christine A. Shoemaker, Water Resources
        Research 44, W12603, http://doi.org/10.1029/2007WR006429

    Example
    -------
    >>> import numpy as np
    >>> from scipy.optimize import rosen

    >>> bl = np.array([-5.,-2.])
    >>> bu = np.array([5.,8.])
    >>> x0 = np.array([-2.,7.])
    >>> bestx, bestf, icall = sce(
    ...     rosen, x0, bl, bu, seed=1, maxn=1000, outf=True, outcall=True,
    ...     printit=2, restartfile1=None)
    >>> print(icall)
    298
    >>> print('{:.3f}'.format(bestf))
    0.001

    >>> nopt = 10
    >>> bl = np.ones((10))*(-5.)
    >>> bu = np.ones((10))*(5.)
    >>> x0 = np.ones((10))*(0.5)
    >>> bestx, bestf, icall = sce(
    ...     rosen, x0, bl, bu, maxn=30000, kstop=10, pcento=0.0001, seed=12358,
    ...     ngs=5, npg=5*nopt+1, nps=nopt+1, nspl=5*nopt+1, mings=2,
    ...     iniflg=True, printit=2, alpha=0.8, beta=0.45,
    ...     outf=True, outcall=True, restartfile1=None)
    >>> print(icall)
    30228
    >>> print('{:.3g}'.format(bestf))
    3.38e-12


    History
    -------
    Written,  Qingyun Duan,   Sep 2004
    Modified, Stijn Van Hoey,     2011 - ported to Python
              Matthias Cuntz, Oct 2013
                  - adapted to JAMS package and sync with JAMS Fortran version
              Matthias Cuntz, Nov 2016 - call external programs
              Matthias Cuntz, Nov 2016 - NaN and Inf
              Matthias Cuntz, Nov 2016 - mask
              Matthias Cuntz, Nov 2016 - restart - only Python 2
              Matthias Cuntz, Nov 2016 - restartfile1=None
              Matthias Cuntz, Nov 2016 - return -bestf if maxit
              Matthias Cuntz, Dec 2017 - rm call of external programs
              Matthias Cuntz, Mar 2018
                  - print out also number of function evaluations with
                    printit==1
                  - mask parameters with degenerated ranges, e.g. upper<lower
                  - use only masked in normalized geometric range of parameters
                  - bug: calc of size of complexes used all parameters not only
                    masked
                  - bug: masked points were always out of bounds
                  - allow scalar bounds, which will be taken for all parameters
              Matthias Cuntz, May 2018
                  - define nn=len(x0) in case of restart
                  - removed *100 from criter_change, which was a bug compared
                    to Fortran code
              Matthias Cuntz, May 2020
                  - removed exec commands for read/write of restart files
                  - use numpy.savez_compressed to be independent of JAMS
              Matthias Cuntz, Apr 2022 - args and kwargs passed to function
              Matthias Cuntz, Jul 2022 - pass RandomState between routines
              Matthias Cuntz, Jul 2022 - sampling

    """

    """
       List of local variables
       x(.,.)    = coordinates of points in the population
       xf(.)     = function values of x(.,.)
       xx(.)     = coordinates of a single point in x
       cx(.,.)   = coordinates of points in a complex
       cf(.)     = function values of cx(.,.)
       s(.,.)    = coordinates of points in the current simplex
       sf(.)     = function values of s(.,.)
       bestx(.)  = best point at current shuffling loop
       bestf     = function value of bestx(.)
       worstx(.) = worst point at current shuffling loop
       worstf    = function value of worstx(.)
       xnstd(.)  = standard deviation of parameters in the population
       gnrng     = normalized geometric mean of parameter ranges
       lcs(.)    = indices locating position of s(.,.) in x(.,.)
       bound(.)  = bound on ith variable being optimized
       ngs1      = number of complexes in current population
       ngs2      = number of complexes in last population
       iseed1    = current random seed
       criter(.) = vector containing the best criterion values of the last
                   10 shuffling loops
    """
    if not restart:
        # Initialize SCE parameters:
        nn = len(x0)
        if mask is None:
            mask = np.ones(nn, dtype=bool)
        nopt = np.sum(mask)
        if npg is None:
            npg   = 2 * nopt + 1
        if nps is None:
            nps   = nopt + 1
        if nspl is None:
            nspl  = 2 * nopt + 1
        if mings is None:
            mings = ngs
        npt = npg * ngs

        # assure numpy array with size(x0)=nn
        bl = np.array(bl)
        bu = np.array(bu)
        # same bounds for all parameters
        if len(bl) == 1:
            bl = np.full(nn, bl[0])
        if len(bu) == 1:
            bu = np.full(nn, bu[0])

        bound = bu - bl

        # Seed random number generator
        rnd = np.random.RandomState(seed=seed)

        # degenerated bounds
        if np.any(bound <= 0.):
            mask[np.where(bound <= 0.)] = False

        large = 0.5 * np.finfo(float).max

        # Create an initial population to fill array x(npt,nn):
        x = _SampleInputMatrix(npt, bl, bu, rnd, sampling=sampling)
        for i in range(npt):
            x[i, :] = np.where(mask, x[i, :], x0)
        if iniflg == 1:
            x[0, :] = x0

        icall = 0
        xf = np.zeros(npt)
        for i in range(npt):
            fuc = func(x[i, :], *args, **kwargs)
            xf[i] = -fuc if maxit else fuc
            icall += 1
            if printit == 1:
                print('  i, f, X: ', icall, xf[i], x[i, :])

        # remember largest for treating of NaNs
        large = xf[np.isfinite(xf)].max()
        large = 1.1 * large if large > 0. else 0.9 * large

        # Sort the population in order of increasing function values;
        xf = np.where(np.isfinite(xf), xf, large)
        idx = np.argsort(xf)
        xf  = np.sort(xf)
        x   = x[idx, :]

        # Record the best and worst points;
        bestx  = x[0, :]
        bestf  = xf[0]
        worstx = x[-1, :]
        worstf = xf[-1]

        allbestf = bestf
        allbestx = bestx

        # # Compute the standard deviation for each parameter
        # xnstd = np.std(x, axis=0)

        # Computes the normalized geometric range of the parameters
        # gnrng = np.exp(
        #     np.mean(np.log((np.max(x,axis=0)-np.min(x,axis=0))/bound)))
        rrange = (np.ma.array(x.max(axis=0) - x.min(axis=0), mask=~mask) /
                  np.ma.array(bound, mask=~mask))
        gnrng  = np.ma.exp(np.ma.mean(np.ma.log(rrange)))

        if printit < 2:
            print('')
            print('Initial complex: Best f: {:f}, worst f {:f}'.format(
                bestf, worstf))
            print('  best X:')
            print(bestx)
            print('')

        # Check for convergence
        if icall >= maxn:
            if printit < 2:
                print('Optimisation terminated because trial number {:d}'
                      ' reached maximum number of trials {:d} at the'
                      ' initial loop.'.format(maxn, icall))

        if gnrng < peps:
            if printit < 2:
                print('The population has converged to a small parameter'
                      ' space {:f} (<{:f}).'.format(gnrng, peps))

        # Begin evolution loops:
        nloop  = 0
        criter = []
        criter_change = 1e+5

        # save restart
        if restartfile1 is not None:
            rs1, rs2, rs3, rs4, rs5 = rnd.get_state()
            # Only arrays with savez_compressed - restartfile1
            savez_compressed(restartfile1,
                             bl=bl, bu=bu, bound=bound, mask=mask,
                             criter=criter,
                             x=x, xf=xf,
                             bestx=bestx, worstx=worstx, allbestf=allbestf,
                             allbestx=allbestx,
                             rs2=rs2)
            p = open(restartfile2, 'w')
            iout  = [nopt, npg, nps, nspl, mings, npt, nloop, icall, rs3, rs4]
            fform = '{:d},' * (len(iout) - 1) + '{:d}'
            print(fform.format(*iout), file=p)
            fout = [gnrng, criter_change, bestf, worstf, rs5]
            fform = '{:.15g},' * (len(fout) - 1) + '{:.15g}'
            print(fform.format(*fout), file=p)
            print(maxit, file=p)
            print(rs1, file=p)
            p.close()

    else:  # if not restart
        nn = len(x0)

        # load restart
        p1 = open(restartfile1, 'rb')
        pp = np.load(p1)
        bl       = pp['bl']
        bu       = pp['bu']
        bound    = pp['bound']
        mask     = pp['mask']
        criter   = pp['criter']
        x        = pp['x']
        xf       = pp['xf']
        bestx    = pp['bestx']
        worstx   = pp['worstx']
        allbestf = pp['allbestf']
        allbestx = pp['allbestx']
        rs2      = pp['rs2']
        p1.close()
        p2 = open(restartfile2, 'r')
        nopt, npg, nps, nspl, mings, npt, nloop, icall, rs3, rs4 = [
            int(inn) for inn in p2.readline().rstrip().split(',') ]
        gnrng, criter_change, bestf, worstf, rs5 = [
            float(inn) for inn in p2.readline().rstrip().split(',') ]
        maxit = bool(strtobool(p2.readline().rstrip()))
        rs1 = p2.readline().rstrip()
        p2.close()
        rnd = np.random.RandomState(seed=seed)
        rnd.set_state((rs1, rs2, rs3, rs4, rs5))

    # Outer Loop
    while icall < maxn and gnrng > peps and criter_change > pcento:
        nloop += 1

        # Loop on complexes (sub-populations);
        for igs in range(ngs):
            # Partition the population into complexes (sub-populations);
            cx = np.zeros((npg, nn))
            cf = np.zeros((npg))

            k1 = np.array(range(npg))
            k2 = k1 * ngs + igs
            cx[k1, :] = x[k2, :]
            cf[k1]    = xf[k2]

            # Evolve sub-population igs for nspl steps:
            for loop in range(nspl):

                # Select simplex by sampling the complex according to a linear
                # probability distribution
                lcs    = np.zeros(nps, dtype=int)
                lcs[0] = 1
                for k3 in range(1, nps):
                    for i in range(1000):
                        # lpos = 1 + int(np.floor(
                        #    npg + 0.5 -
                        # np.sqrt((npg+0.5)**2 - npg*(npg+1)*rnd.random())))
                        lpos = int(np.floor(npg + 0.5 -
                                            np.sqrt((npg + 0.5)**2 -
                                                    npg * (npg + 1) *
                                                    rnd.rand())))
                        # check of element al eens gekozen
                        # idx=find(lcs(1:k3-1)==lpos)
                        idx = (lcs[0:k3] == lpos).nonzero()
                        if idx[0].size == 0:
                            break
                    lcs[k3] = lpos
                lcs.sort()

                # Construct the simplex:
                s  = np.zeros((nps, nn))
                s  = cx[lcs, :]
                sf = cf[lcs]

                # remember largest for treating of NaNs
                large = cf[np.isfinite(cf)].max()
                large = 1.1 * large if large > 0. else 0.9 * large

                snew, fnew, icall = _cce(func, s, sf, bl, bu, mask,
                                         sampling, icall, maxn, alpha,
                                         beta, maxit, printit, rnd,
                                         args=args, kwargs=kwargs)
                # Replace the worst point in Simplex with the new point:
                s[-1, :] = snew
                sf[-1]   = fnew

                # Replace the simplex into the complex;
                cx[lcs, :] = s
                cf[lcs]    = sf

                # Sort the complex;
                cf  = np.where(np.isfinite(cf), cf, large)
                idx = np.argsort(cf)
                cf  = np.sort(cf)
                cx  = cx[idx, :]
            # End of Inner Loop for Competitive Evolution of Simplexes:
            # for loop in range(nspl):
            # i.e. end of evolve sub-population igs for nspl steps

            # Replace the complex back into the population;
            x[k2, :] = cx[k1, :]
            xf[k2]   = cf[k1]
        # End of Loop on Complex Evolution: for igs in range(ngs):

        # Shuffled the complexes;
        idx = np.argsort(xf)
        xf  = np.sort(xf)
        x   = x[idx, :]

        # Record the best and worst points;
        bestx  = x[0, :]
        bestf  = xf[0]
        worstx = x[-1, :]
        worstf = xf[-1]

        # appenden en op einde reshapen!!
        allbestx = np.append(allbestx, bestx, axis=0)
        allbestf = np.append(allbestf, bestf)

        # # Compute the standard deviation for each parameter
        # xnstd = np.std(x, axis=0)

        # Computes the normalized geometric range of the parameters
        # gnrng=np.exp(np.mean(np.log((np.max(x,axis=0)-np.min(x,axis=0))/bound)))
        rrange = (np.ma.array(x.max(axis=0) - x.min(axis=0), mask=~mask) /
                  np.ma.array(bound, mask=~mask))
        gnrng  = np.ma.exp(np.ma.mean(np.ma.log(rrange)))

        if printit < 2:
            print('')
            print('Evolution loop {0:d}, trials {1:d}. Best f: {2:f},'
                  ' worst f {3:f}'.format(nloop, icall, bestf, worstf))
            print('  best X:')
            print(bestx)
            print('')

        # Check for convergency;
        if icall >= maxn:
            if printit < 2:
                print('Optimisation terminated because trial number {:d} '
                      'reached maximum number of trials'
                      ' {:d}.'.format(maxn, icall))

        if gnrng < peps:
            if printit < 2:
                print('The population has converged to a small parameter '
                      ' space {:f} (<{:f}).'.format(gnrng, peps))

        criter = np.append(criter, bestf)

        if nloop >= kstop:  # nodig zodat minimum zoveel doorlopen worden
            criter_change = np.abs(criter[nloop - 1] - criter[nloop - kstop])
            criter_change = criter_change / np.maximum(
                np.mean(np.abs(criter[nloop - kstop:nloop])), 1e-15)
            if criter_change < pcento:
                if printit < 2:
                    print('The best point has improved by less then {:f} in'
                          ' the last {:d} loops.'.format(pcento, kstop))

        # save restart
        if restartfile1 is not None:
            rs1, rs2, rs3, rs4, rs5 = rnd.get_state()
            savez_compressed(restartfile1,
                             bl=bl, bu=bu, bound=bound, mask=mask,
                             criter=criter,
                             x=x, xf=xf,
                             bestx=bestx, worstx=worstx, allbestf=allbestf,
                             allbestx=allbestx, rs2=rs2)
            p = open(restartfile2, 'w')
            iout  = [nopt, npg, nps, nspl, mings, npt, nloop, icall, rs3, rs4]
            fform = '{:d},' * (len(iout) - 1) + '{:d}'
            print(fform.format(*iout), file=p)
            fout = [gnrng, criter_change, bestf, worstf, rs5]
            fform = '{:.15g},' * (len(fout) - 1) + '{:.15g}'
            print(fform.format(*fout), file=p)
            print(maxit, file=p)
            print(rs1, file=p)
            p.close()
    # End of the Outer Loop: while icall<maxn and gnrng>peps
    # and criter_change>pcento

    if printit < 2:
        print('Search stopped at trial number {0:d} with normalized'
              ' geometric range {1:f}. '.format(icall, gnrng))
        print('The best point has improved by {:f} in the last'
              ' {:d} loops.'.format(criter_change, kstop))

    # reshape allbestx
    allbestx = allbestx.reshape(allbestx.size // nn, nn)

    # end of subroutine sce
    if maxit:
        bestf    *= -1.
        allbestf *= -1.
    out = [bestx]
    if outf:
        out += [bestf]
    if outhist:
        out += [allbestx, allbestf]
    if outcall:
        out += [icall]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # maxn = 10000
    # from pyjams.functions import ackley, griewank, goldstein_price
    # from pyjams.functions import rastrigin, rosenbrock, six_hump_camelback
    # seed = 12345

    # """
    # This is the Ackley Function
    # Global Optimum (n>=2): 0.0 at origin
    # """
    # npara = 10
    # bl = -10 * np.ones(npara)
    # bu = 10 * np.ones(npara)
    # x0 = np.random.rand(npara) * 10.
    # bestx, bestf = sce(ackley, x0, bl, bu, maxn=maxn, outf=True,
    #                    restartfile1=None, seed=seed)
    # print('Ackley ', bestx, bestf)

    # """
    # This is the Griewank Function (2-D or 10-D)
    # Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
    #    Global Optimum: 0, at origin
    # """
    # npara = 10
    # bl = -600 * np.ones(npara)
    # bu = 600 * np.ones(npara)
    # x0 = np.random.rand(npara) * 600.
    # bestx, bestf = sce(griewank, x0, bl, bu, maxn=maxn, outf=True,
    #                    restartfile1=None, seed=seed)
    # print('Griewank ', bestx, bestf)

    # """
    # This is the Goldstein-Price Function
    # Bound X1=[-2,2], X2=[-2,2]
    # Global Optimum: 3.0,(0.0,-1.0)
    # """
    # npara = 2
    # bl = -2 * np.ones(npara)
    # bu = 2 * np.ones(npara)
    # x0 = np.random.rand(npara) * 2.
    # bestx, bestf = sce(goldstein_price, x0, bl, bu, maxn=maxn,
    #                    outf=True, restartfile1=None, seed=seed)
    # print('Goldstein ', bestx, bestf)

    # """
    # This is the Rastrigin Function
    # Bound: X1=[-1,1], X2=[-1,1]
    # Global Optimum: -2, (0,0)
    # """
    # npara = 2
    # bl = -1 * np.ones(npara)
    # bu = 1 * np.ones(npara)
    # x0 = np.random.rand(npara) * 1.
    # bestx, bestf = sce(rastrigin, x0, bl, bu, maxn=maxn,
    #                    outf=True, restartfile1=None, seed=seed)
    # print('Rastrigin ', bestx, bestf)

    # """
    # This is the Rosenbrock Function
    # Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
    #        bl=[-5 -5]; bu=[5 5]; x0=[1 1];
    # """
    # npara = 2
    # bl = -2 * np.ones(npara)
    # bu = 5 * np.ones(npara)
    # x0 = np.random.rand(npara) * 2.
    # bestx, bestf = sce(rosenbrock, x0, bl, bu, maxn=maxn,
    #                    outf=True, restartfile1=None, seed=seed)
    # print('Rosenbrock ', bestx, bestf)

    # """
    # This is the Six-hump Camelback Function.
    # Bound: X1=[-5,5], X2=[-5,5]
    # True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    # """
    # npara = 2
    # bl = -5 * np.ones(npara)
    # bu = 5 * np.ones(npara)
    # x0 = np.random.rand(npara) * 5.
    # bestx, bestf = sce(six_hump_camelback, x0, bl, bu, maxn=maxn,
    #                    outf=True, restartfile1=None, seed=seed)
    # print('Six_hump_camelback ', bestx, bestf)

    # # Restart
    # maxn = 500
    # npara = 2
    # lb = -2 * np.ones(npara)
    # ub = 5 * np.ones(npara)
    # x0 = np.zeros(npara)
    # bestx, bestf = sce(rosenbrock, x0, lb, ub, maxn=maxn,
    #                    outf=True, restartfile1=None,
    #                    seed=seed, restart=False)
    # print('Rosenbrock Reference - ', bestx, bestf)
    # bestx, bestf = sce(rosenbrock, x0, lb, ub, maxn=maxn // 2, outf=True,
    #                    seed=seed, restart=False)
    # print('Rosenbrock Restart 1 - ', bestx, bestf)
    # bestx, bestf = sce(rosenbrock, x0, lb, ub, maxn=maxn, outf=True,
    #                    seed=seed, restart=True)
    # print('Rosenbrock Restart 2 - ', bestx, bestf)
    # # Clean restart
    # import os
    # os.remove('sce.restart.npz')
    # os.remove('sce.restart.txt')

    # # different sampling
    # npara = 2
    # bl = -2 * np.ones(npara)
    # bu = 5 * np.ones(npara)
    # x0 = np.random.rand(npara) * 2.
    # # sampling is str
    # for sampling in ['half-open', 'left-half-open',
    #                  'right-half-open', 'open', 'log']:
    #     bestx, bestf = sce(rosenbrock, x0, bl, bu, sampling=sampling,
    #                        maxn=maxn, outf=True, restartfile1=None, seed=seed)
    #     print('Rosenbrock0 ', sampling, bestx, bestf)
    # # sampling is list
    # sampling = ['half-open', 'log']
    # bestx, bestf = sce(rosenbrock, x0, bl, bu, sampling=sampling,
    #                    maxn=maxn, outf=True, restartfile1=None, seed=seed)
    # print('Rosenbrock1 ', sampling, bestx, bestf)
    # # bounds include 0
    # bl = [0., -2.]
    # bu = [5., 8.]
    # x0 = np.random.rand(npara) * 2.
    # for sampling in ['half-open', 'left-half-open',
    #                  'right-half-open', 'open', 'log']:
    #     bestx, bestf = sce(rosenbrock, x0, bl, bu, sampling=sampling,
    #                        maxn=maxn, outf=True, restartfile1=None, seed=seed)
    #     print('Rosenbrock2 ', sampling, bestx, bestf)
    # # bounds all > 0
    # bl = [0.1, 0.1]
    # bu = [5., 8.]
    # x0 = np.random.rand(npara) * 2.
    # for sampling in ['half-open', 'left-half-open',
    #                  'right-half-open', 'Open', 'log']:
    #     bestx, bestf = sce(rosenbrock, x0, bl, bu, sampling=sampling,
    #                        maxn=maxn, outf=True, restartfile1=None, seed=seed)
    #     print('Rosenbrock3 ', sampling, bestx, bestf)
