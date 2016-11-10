#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import subprocess
from jams.const import huge
# ToDo: restart
# ToDo: tmp/population files

def SampleInputMatrix(nrows, npars, bl, bu, distname='randomUniform'):
    '''
        Create input parameter matrix (nrows,npars) for
        nrows simulations and npars parameters with bounds bl and bu


        Definition
        ----------
        def SampleInputMatrix(nrows, npars, bl, bu, distname='randomUniform'):


        Input
        -----
        nrows             # of simulations
        npars             # of parameters
        bl                (npars) lower bounds of parameters
        bu                (npars) upper bounds of parameters


        Optional Input
        --------------
        distname          initial sampling ditribution (not implemented yet, takes uniform distribution)


        Output
        ------
        parameter matrix (nrows, npars) with new parameter samples


        References
        ----------
        Duan, Q., S. Sorooshian, and V. Gupta,
            Effective and efficient global optimization for conceptual rainfall-runoff models,
            Water Resour. Res., 28, 1015-1031, 1992.


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2013 Matthias Cuntz


        History
        -------
        Written,  Q. Duan, Sep 2004
        Modified, S. Van Hoey 2011 - ported to Python
                  MC, Oct 2013     - adapted to JAMS package and sync with JAMS Fortran version
    '''
    x = np.zeros((nrows,npars))
    bound = bu-bl
    for i in range(nrows):
        # x[i,:]= bl + DistSelector([0.0,1.0,npars],distname='randomUniform')*bound # only used in full Vhoeys-framework
        x[i,:]= bl + np.random.rand(1,npars)*bound
    return x


def cce(functn, s, sf, bl, bu, mask, icall, maxn, alpha, beta, maxit, printit,
        parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug):
    '''
        Generate a new point in a simplex


        Definition
        ----------
        def cce(functn, s, sf, bl, bu, mask, icall, maxn, alpha, beta, maxit, printit):


        Input
        -----
        functn            function to minimise (python function or string for external executable)
        s                 2D-array, the sorted simplex in order of increasing function values
        sf                1D-array, function values in increasing order
        bl                (npars) lower bounds of parameters
        bu                (npars) upper bounds of parameters
        mask              (npars) mask to include (1) or exclude (0) parameter from optimisation
        icall             counter of function calls
        maxn              maximum number of function evaluations allowed during optimization
        alpha             parameter for reflection  of points in complex
        beta              parameter for contraction of points in complex
        maxit             if True: maximise instead of minimise functn
        printit           if ==1: print each function evaluation
        parameterfile     Parameter file for executable; must be given if functn is name of executable
        parameterwriter   Python function for writing parameter file if functn is name of executable
        objectivefile     File with objective value from executable; must be given if functn is name of executable
        objectivereader   Python function for reading objective value if functn is name of executable
        shell             If True, the specified command will be executed through the shell.
        debug             If True, model output is displayed for executable.


        Optional Input
        --------------
        None


        Output
        ------
        new parameter set, function value of new set, number of function calls


        References
        ----------
        Duan, Q., S. Sorooshian, and V. Gupta,
            Effective and efficient global optimization for conceptual rainfall-runoff models,
            Water Resour. Res., 28, 1015-1031, 1992.


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2013 Matthias Cuntz


        History
        -------
        Written,  Q. Duan, Sep 2004
        Modified, S. Van Hoey 2011 - ported to Python
                  MC, Oct 2013 - adapted to JAMS package and sync with JAMS Fortran version
    '''

    '''
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
    '''
    nps, nopt = s.shape
    n = nps
    m = nopt

    # Assign the best and worst points:
    sb = s[0,:]
    fb = sf[0]
    sw = s[-1,:]
    fw = sf[-1]

    # Compute the centroid of the simplex excluding the worst point:
    ce = np.mean(s[:-1,:],axis=0)

    # Attempt a reflection point
    snew = ce + alpha*(ce-sw)
    snew = np.where(mask, snew, sb) # sb should have initial params at mask==False

    # Check if is outside the bounds:
    ibound = 0
    s1 = snew-bl
    idx = (s1<0).nonzero()
    if idx[0].size != 0:
        ibound = 1

    s1 = bu-snew
    idx = (s1<0).nonzero()
    if idx[0].size != 0: ibound = 2

    if ibound >= 1: 
        snew = SampleInputMatrix(1,nopt,bl,bu,distname='randomUniform')[0]  #checken!!
        snew = np.where(mask, snew, sb)

    fuc = call_function(functn, snew, bl, bu, mask,
                        parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
    fnew = -fuc if maxit else fuc
    if printit==1: print('  f, X: ', fnew, snew)
    icall += 1

    # Reflection failed; now attempt a contraction point:
    if fnew > fw:
        snew = sw + beta*(ce-sw)
        snew = np.where(mask, snew, sb)
        fuc = call_function(functn, snew, bl, bu, mask,
                            parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
        fnew = -fuc if maxit else fuc
        if printit==1: print('  f, X: ', fnew, snew)
        icall += 1

    # Both reflection and contraction have failed, attempt a random point;
        if fnew > fw:
            snew = SampleInputMatrix(1,nopt,bl,bu,distname='randomUniform')[0]  #checken!!
            snew = np.where(mask, snew, sb)
            fuc = call_function(functn, snew, bl, bu, mask,
                                parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
            fnew = -fuc if maxit else fuc
            if printit==1: print('  f, X: ', fnew, snew)
            icall += 1

    # END OF CCE
    return snew, fnew, icall


def call_function(functn, params, bl, bu, mask,
                  parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug):
    '''
        Call python function or external executable


        Definition
        ----------
        def call_function(functn, params, bl, bu, mask,
                          parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug):


        Input
        -----
        functn            function to minimise (python function or string for external executable)
        params            (npars) parameter set
        bl                (npars) lower bounds of parameters
        bu                (npars) upper bounds of parameters
        mask              (npars) mask to include (1) or exclude (0) parameter from optimisation
        parameterfile     Parameter file for executable; must be given if functn is name of executable
        parameterwriter   Python function for writing parameter file if functn is name of executable
        objectivefile     File with objective value from executable; must be given if functn is name of executable
        objectivereader   Python function for reading objective value if functn is name of executable
        shell             If True, the specified command will be executed through the shell.
        debug             If True, model output is displayed for executable.

        Output
        ------
        Function value


        History
        -------
        Written,  MC, Nov 2016
    '''
    if isinstance(functn, (str,list)):
        parameterwriter(parameterfile, params, bl, bu, mask)
        if debug:
            err = subprocess.call(functn, shell=shell)
        else:
            err = subprocess.check_output(functn, shell=shell)
        obj = objectivereader(objectivefile)
        return obj
    else:
        return functn(params)
    

def sce(functn, x0, bl, bu,
        mask=None,
        maxn=1000, kstop=10, pcento=0.0001,
        ngs=2, npg=None, nps=None, nspl=None, mings=None,
        peps=0.001, seed=0, iniflg=True,
        alpha=0.8, beta=0.45, maxit=False, printit=2,
        outf=False, outhist=False, outcall=False,
        parameterfile=None, parameterwriter=None,
        objectivefile=None, objectivereader=None,
        shell=False, debug=False):
    '''
        Shuffled-Complex-Evolution algorithm for function minimalisation


        Definition
        ----------
        def sce(functn, x0, bl, bu,
                mask=None,
                maxn=1000, kstop=10, pcento=0.0001,
                ngs=2, npg=None, nps=None, nspl=None, mings=None,
                peps=0.001, seed=0, iniflg=True,
                alpha=0.8, beta=0.45, maxit=False, printit=0,
                outf=False, outhist=False, outcall=False,
                parameterfile=None, parameterwriter=None,
                objectivefile=None, objectivereader=None,
                shell=False, debug=False):


        Input
        -----
        functn            function to minimise (python function or string for external executable)
        x0                1D-array of size nopt, initial parameter set
        bl                (npars) lower bounds of parameters
        bu                (npars) upper bounds of parameters


        Optional Input
        --------------
        mask              (npars) include (1,True) or exclude (0,False) parameters in optimisation
        maxn              maximum number of function evaluations allowed during optimization (default: 1000)
        kstop             maximum number of evolution loops before convergency (default: 10)
        pcento            the percentage change allowed in kstop loops before convergency (default: 0.0001)
        peps              value of normalised geometric range needed for convergence (default: 0.001)
        ngs               number of complexes (default: 2)
        npg               number of points in each complex (default: 2*nopt+1)
        nps               number of points in each sub-complex (default: nopt+1)
        mings             minimum number of complexes required if the number of complexes is allowed to reduce as the
                          optimization proceeds (default: ngs)
        nspl              number of evolution steps allowed for each complex before complex shuffling (default: 2*nopt+1)
        seed              if >0, the random number seed (default: 0)
        iniflg            if True: include initial parameter in initial population (default: True)
        alpha             parameter for reflection  of points in complex (default: 0.8)
        beta              parameter for contraction of points in complex (default: 0.45)
        maxit             if True: maximise instead of minimise functn (default: false)
        printit           controlling print-out (default: 2)
                          0: print information on the best point of the population
                          1: print information on each function evaluation
                          2: no printing
        outf              if True: return best function value (default: False)
        outhist           if True: return parameters and function values of each evolution loop (default: False)
        outcall           if True: return number of function evaluations (default: False)
        parameterfile     Parameter file for executable; must be given if functn is name of executable
        parameterwriter   Python function for writing parameter file if functn is name of executable
        objectivefile     File with objective value from executable; must be given if functn is name of executable
        objectivereader   Python function for reading objective value if functn is name of executable
        shell             If True, the specified executable will be executed through the shell (default: False).
        debug             If True, model output is displayed for executable (default: False).


        Output
        ------
        best parameter set
        depending on outf, outhist and outcall also
            function value at best parameter set
            vector of best parameters sets per evolution, function values of the best parameter vector
            number of function evaluations
        If functn is name of executable, parameterwriter will be called at the end,
        i.e. parameter file contains best parameter set.


        References
        ----------
        Duan, Q., S. Sorooshian, and V. Gupta,
            Effective and efficient global optimization for conceptual rainfall-runoff models,
            Water Resour. Res., 28, 1015-1031, 1992.


        Example
        -------
        >>> import numpy as np
        >>> from functions import ackley, rosenbrock
        >>> from autostring import astr

        >>> bl = np.array([-5.,-2.])
        >>> bu = np.array([5.,8.])
        >>> x0 = np.array([-2.,7.])
        >>> bestx, bestf, icall = sce(rosenbrock, x0, bl, bu, seed=1, maxn=1000, outf=True, outcall=True, printit=2)
        >>> print(astr(icall))
        298
        >>> print(astr(bestf,3))
        0.001

        >>> nopt = 30
        >>> bl = np.ones((30))*(-10.)
        >>> bu = np.ones((30))*(10.)
        >>> x0 = np.ones((30))*(0.5)
        >>> bestx, bestf, icall = sce(ackley, x0, bl, bu,
        ...                           maxn=30000, kstop=10, pcento=0.0001, seed=10987,
        ...                           ngs=2, npg=2*nopt+1, nps=nopt+1, nspl=2*nopt+1, mings=2,
        ...                           iniflg=True, printit=2, alpha=0.8, beta=0.45,
        ...                           outf=True, outcall=True)
        >>> print(astr(icall))
        4608
        >>> print(astr(bestf,3))
        0.014


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2013 Matthias Cuntz


        History
        -------
        Written,  Q. Duan, Sep 2004
        Modified, S. Van Hoey 2011 - ported to Python
                  MC, Oct 2013 - adapted to JAMS package and sync with JAMS Fortran version
                  MC, Nov 2016 - call external programs
                  MC, Nov 2016 - NaN and Inf
                  MC, Nov 2016 - mask
    '''

    '''
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
    '''
    # Initialize SCE parameters:
    nopt = len(x0)
    if npg   is None: npg   = 2*nopt+1
    if nps   is None: nps   = nopt+1
    if nspl  is None: nspl  = 2*nopt+1
    if mings is None: mings = ngs
    npt = npg*ngs

    # assure numpy array
    bl = np.array(bl)
    bu = np.array(bu)

    bound = bu-bl

    if seed>0: np.random.seed(seed)

    if mask is None: mask = np.ones(nopt, dtype=np.bool)

    large = -0.5*huge if maxit else 0.5*huge

    # Check parameterfile etc. if functn is name of executable
    if isinstance(functn, (str,list)):
        if parameterfile is None:
            raise IOError('parameterfile must be given if functn is name of executable.')
        if parameterwriter is None:
            raise IOError('parameterwrite must be given if functn is name of executable.')
        if objectivefile is None:
            raise IOError('objectivefile must be given if functn is name of executable.')
        if objectivereader is None:
            raise IOError('objectivereader must be given if functn is name of executable.')

    # Create an initial population to fill array x(npt,nopt):
    # ToDo: include mask
    x = SampleInputMatrix(npt, nopt, bl, bu, distname='randomUniform')
    for i in range(npt):
        x[i,:] = np.where(mask, x[i,:], x0)
    if iniflg==1: x[0,:] = x0

    nloop = 0
    icall = 0
    xf = np.zeros(npt)
    for i in range (npt):
        fuc = call_function(functn, x[i,:], bl, bu, mask,
                            parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
        xf[i] = -fuc if maxit else fuc
        if printit==1: print('  f, X: ', xf[i], x[i,:])
        icall += 1
    f0 = xf[0]

    # remember largest for treating of NaNs
    if maxit:
        large = xf[np.isfinite(xf)].min()
        large = 0.9*large if large>0. else 1.1*large
    else:
        large = xf[np.isfinite(xf)].max()
        large = 1.1*large if large>0. else 0.9*large

    # Sort the population in order of increasing function values;
    xf = np.where(np.isfinite(xf), xf, large)
    idx = np.argsort(xf)
    xf  = np.sort(xf)
    x   = x[idx,:]

    # Record the best and worst points;
    bestx  = x[0,:]
    bestf  = xf[0]
    worstx = x[-1,:]
    worstf = xf[-1]

    allbestf = bestf
    allbestx = bestx

    # Compute the standard deviation for each parameter
    xnstd = np.std(x,axis=0)

    # Computes the normalized geometric range of the parameters
    gnrng = np.exp(np.mean(np.log((np.max(x,axis=0)-np.min(x,axis=0))/bound)))

    if printit<2:
        print('The Initial Loop: 0. Best f: {:f}, worst f {:f}'.format(bestf, worstf))
        print('  best X: ', bestx)
        print('')

    # Check for convergence
    if icall >= maxn:
        if printit<2:
            print('Optimisation terminated because trial number {:d} '
                  'reached maximum number of trials {:d} at the initial loop.'.format(maxn,icall))

    if gnrng < peps:
        if printit<2:
            print('The population has converged to a small parameter space {:f} (<{:f}).'.format(gnrng, peps))

    # Begin evolution loops:
    nloop  = 0
    criter = []
    criter_change = 1e+5

    while icall<maxn and gnrng>peps and criter_change>pcento:
        nloop += 1

        # Loop on complexes (sub-populations);
        for igs in range(ngs):
            # Partition the population into complexes (sub-populations);
            cx = np.zeros((npg,nopt))
            cf = np.zeros((npg))

            k1 = np.array(range(npg))
            k2 = k1*ngs+igs
            cx[k1,:] = x[k2,:]
            cf[k1]   = xf[k2]

            # Evolve sub-population igs for nspl steps:
            for loop in range(nspl):

                # Select simplex by sampling the complex according to a linear
                # probability distribution
                lcs    = np.array([0]*nps)
                lcs[0] = 1
                for k3 in range(1,nps):
                    for i in range(1000):
                        # lpos = 1 + int(np.floor(npg+0.5-np.sqrt((npg+0.5)**2 - npg*(npg+1)*random.random())))
                        lpos = int(np.floor(npg+0.5-np.sqrt((npg+0.5)**2 - npg*(npg+1)*np.random.rand())))
                        # idx=find(lcs(1:k3-1)==lpos)
                        idx = (lcs[0:k3]==lpos).nonzero()  # check of element al eens gekozen
                        if idx[0].size == 0:
                            break
                    lcs[k3] = lpos
                lcs.sort()

                # Construct the simplex:
                s  = np.zeros((nps,nopt))
                s  = cx[lcs,:]
                sf = cf[lcs]

                # remember largest for treating of NaNs
                if maxit:
                    large = cf[np.isfinite(cf)].min()
                    large = 0.9*large if large>0. else 1.1*large
                else:
                    large = cf[np.isfinite(cf)].max()
                    large = 1.1*large if large>0. else 0.9*large

                snew, fnew, icall = cce(functn, s, sf, bl, bu, mask, icall, maxn, alpha, beta, maxit, printit,
                                        parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
                # Replace the worst point in Simplex with the new point:
                s[-1,:] = snew
                sf[-1]  = fnew

                # Replace the simplex into the complex;
                cx[lcs,:] = s
                cf[lcs]   = sf

                # Sort the complex;
                cf  = np.where(np.isfinite(cf), cf, large)
                idx = np.argsort(cf)
                cf  = np.sort(cf)
                cx  = cx[idx,:]
            # End of Inner Loop for Competitive Evolution of Simplexes: for loop in range(nspl):
            # i.e. end of evolve sub-population igs for nspl steps

            # Replace the complex back into the population;
            x[k2,:] = cx[k1,:]
            xf[k2]  = cf[k1]
        # End of Loop on Complex Evolution: for igs in range(ngs):

        # Shuffled the complexes;
        idx = np.argsort(xf)
        xf  = np.sort(xf)
        x   = x[idx,:]

        PX = x
        PF = xf

        # Record the best and worst points;
        bestx  = x[0,:]
        bestf  = xf[0]
        worstx = x[-1,:]
        worstf = xf[-1]

        allbestx = np.append(allbestx,bestx, axis=0) #appenden en op einde reshapen!!
        allbestf = np.append(allbestf,bestf)

        # Compute the standard deviation for each parameter
        xnstd=np.std(x,axis=0)

        # Computes the normalized geometric range of the parameters
        gnrng=np.exp(np.mean(np.log((np.max(x,axis=0)-np.min(x,axis=0))/bound)))

        if printit<2:
            print('Evolution loop {:d}, trials {:d}. Best f: {:f}, worst f {:f}'.format(nloop, icall, bestf, worstf))
            print('  best X: ', bestx)
            print('')

        # Check for convergency;
        if icall >= maxn:
            if printit<2:
                print('Optimisation terminated because trial number {:d} '
                      'reached maximum number of trials {:d}.'.format(maxn,icall))

        if gnrng < peps:
            if printit<2:
                print('The population has converged to a small parameter space {:f} (<{:f}).'.format(gnrng, peps))

        criter = np.append(criter,bestf)

        if nloop >= kstop: # nodig zodat minimum zoveel doorlopen worden
            criter_change = np.abs(criter[nloop-1]-criter[nloop-kstop])*100
            criter_change = criter_change/np.mean(np.abs(criter[nloop-kstop:nloop]))
            if criter_change < pcento:
                if printit<2:
                    print('The best point has improved by less then {:f} in the last {:d} loops.'.format(pcento, kstop))
    # End of the Outer Loop: while icall<maxn and gnrng>peps and criter_change>pcento

    if printit<2:
        print('Search stopped at trial number {:d} with normalized geometric range {:f}. '.format(icall, gnrng))
        print('The best point has improved by {:f} in the last {:d} loops.'.format(criter_change, kstop))

    # reshape allbestx
    allbestx = allbestx.reshape(allbestx.size/nopt,nopt)

    # end of subroutine sce
    out = [bestx]
    if outf:    out += [bestf]
    if outhist: out += [allbestx, allbestf]
    if outcall: out += [icall]

    # write parameter file with best parameters
    if isinstance(functn, (str,list)):
        parameterwriter(parameterfile, bestx, bl, bu, mask)
        

    return out


def griewank(x):
    '''
    This is the Griewank Function (2-D or 10-D)
    Bound: X(i)=[-600,600], for i=1,2,...,10
    Global Optimum: 0, at origin
    '''
    nopt = np.size(x)
    #if (nopt == 2) | (nopt == 10):
    xx = x
    if nopt==2:
        d = 200.0
    else:
        d = 4000.0

    u1 = 0.0
    u2 = 1.0
    for j in range(nopt):
        u1 = u1 + xx[j]**2/d
        u2 = u2 * np.cos(xx[j]/np.sqrt(float(j+1)))

    f = u1 - u2 + 1
    return f


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # #-------------------------------------------------------------------------------
    # # Name:        Test for Shuffled Complex Evolution Algorithm implementation
    # # Purpose:      Call one of the example functions to run and find the optimum
    # #               Visualize the 2D (or 3D) objective function +trace of the allbestx
    # #
    # # Author:      VHOEYS
    # #
    # # Created:     11/10/2011
    # # Copyright:   (c) VHOEYS 2011
    # # Licence:     <your licence>
    # #-------------------------------------------------------------------------------

    # import os
    # import sys
    # import time
    # from functions import goldstein_price, griewank, rastrigin, rosenbrock, six_hump_camelback

    # # PARAMETERS TO TUNE THE ALGORITHM
    # # Definition:import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # #  seed = the random seed number (for repetetive testing purpose;pos integers)
    # #  iniflg = flag for initial parameter array (=1, included it in initial
    # #           population; otherwise, not included)
    # #  ngs = number of complexes (sub-populations)
    # #  peps = value of NORMALIZED GEOMETRIC RANGE needed for convergence
    # #  maxn = maximum number of function evaluations allowed during optimization
    # #  kstop = maximum number of evolution loops before convergency
    # #  pcento = the percentage change allowed in kstop loops before convergency
    # maxn=10000

    # # PARAMETERS FOR OPTIMIZATION PROBLEM
    # # Definition:
    # #  x0 = the initial parameter array at the start; np.array
    # #     = the optimized parameter array at the end;
    # #  bl = the lower bound of the parameters; np.array
    # #  bu = the upper bound of the parameters; np.array

    # foo=input('Please enter an Example number (1-5) for example: ')
    # start = time.clock()

    # if foo==1:
    #     '''1
    #      This is the Goldstein-Price Function
    #      Bound X1=[-2,2], X2=[-2,2]; Global Optimum: 3.0,(0.0,-1.0)
    #     '''
    #     bl=np.array([-2,-2])
    #     bu=np.array([2,2])
    #     x0=np.array([2,2])
    #     bestx, allbestx, allbestf = sce(goldstein_price, x0, bl, bu, maxn=maxn, outhist=True)

    # elif foo==2:
    #     '''2
    #       This is the Rosenbrock Function
    #       Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
    #         bl=[-5 -5]; bu=[5 5]; x0=[1 1];
    #     '''
    #     bl=np.array([-5.,-2.])
    #     bu=np.array([5.,8.])
    #     x0=np.array([-2.,7.])
    #     bestx, allbestx, allbestf = sce(rosenbrock, x0, bl, bu, maxn=maxn, outhist=True)

    # elif foo==3:
    #     '''3
    #     %  This is the Six-hump Camelback Function.
    #     %  Bound: X1=[-5,5], X2=[-5,5]
    #     %  True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    #     '''

    #     bl=np.array([-2.,-2.])
    #     bu=np.array([2.,2.])
    #     x0=np.array([-1.,1.])
    #     bestx, allbestx, allbestf = sce(six_hump_camelback, x0, bl, bu, maxn=maxn, outhist=True)

    # elif foo==4:
    #     '''4
    #     %  This is the Rastrigin Function
    #     %  Bound: X1=[-1,1], X2=[-1,1]
    #     %  Global Optimum: -2, (0,0)
    #     '''
    #     bl=np.array([-5.,-5.])
    #     bu=np.array([5.,5.])
    #     x0=np.array([-1.,1.])
    #     bestx, allbestx, allbestf = sce(rastrigin, x0, bl, bu, maxn=maxn, outhist=True, ngs=3, kstop=30, pcento=0.001)

    # elif foo==5:
    #     '''5
    #       This is the Griewank Function (2-D or 10-D)
    #       Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
    #       Global Optimum: 0, at origin
    #     '''
    #     bl=-600*np.ones(2)
    #     bu=600*np.ones(2)
    #     x0=np.ones(2)
    #     bestx, allbestx, allbestf = sce(griewank, x0, bl, bu, maxn=maxn, outhist=True, ngs=3, kstop=30, pcento=0.001)

    # else:
    #     raise ValueError('Example number > 5')

    # elapsed = (time.clock() - start)
    # print('The calculation of the SCE algorithm took %f seconds' %elapsed)

    # '''
    # plot the trace of the parametersvalue
    # '''
    # import matplotlib.pyplot as plt

    # fig=plt.figure()

    # ax1 = plt.subplot(121)
    # ax1.plot(allbestx)
    # plt.title('Trace of the different parameters')
    # plt.xlabel('Evolution Loop')
    # plt.ylabel('Parvalue')

    # '''
    # Plot the parmaeterspace in 2D with trace of the allbestx
    # '''
    # #   make these smaller to increase the resolution
    # dx, dy = 0.05, 0.05
    # if foo == 5: dx, dy = 5., 5.

    # x = np.arange(bl[0], bu[0], dx)
    # y = np.arange(bl[1], bu[1], dy)
    # X, Y = np.meshgrid(x, y)

    # if foo==1:
    #     parspace = goldstein_price([X,Y])
    # elif foo==2:
    #     parspace = rosenbrock([X,Y])
    # elif foo==3:
    #     parspace = six_hump_camelback([X,Y])
    # elif foo==4:
    #     parspace = rastrigin([X,Y])
    # elif foo==5:
    #     parspace = np.empty((x.size,y.size))
    #     for i in range(x.size):
    #         for j in range(y.size):
    #             parspace[i,j] = griewank([x[i],y[j]])

    # ax2 = plt.subplot(122)
    # ax2.pcolor(X, Y, parspace)
    # ##plt.colorbar()
    # ax2.plot(allbestx[:,0],allbestx[:,1],'*')
    # plt.title('Trace of the allbestx parameter combinations')
    # plt.xlabel('PAR 1')
    # plt.ylabel('PAR 2')

    # # # Plot the parmaeterspace in 3D - commented out
    # # from mpl_toolkits.mplot3d import Axes3D
    # # from matplotlib import cm
    # # from matplotlib.ticker import LinearLocator, FormatStrFormatter
    # # fig = plt.figure()
    # # ax = fig.gca(projection='3d')
    # # surf = ax.plot_surface(X, Y, parspace, rstride=8, cstride=8, cmap=cm.jet,linewidth=0, antialiased=False)
    # # cset = ax.contourf(X, Y, parspace, zdir='z', offset=-100)
    # # fig.colorbar(surf, shrink=0.5, aspect=5)

    # plt.show()



    # from functions import goldstein_price, griewank, rastrigin, rosenbrock, six_hump_camelback
    # '''
    #     This is the Griewank Function (2-D or 10-D)
    #     Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
    #        Global Optimum: 0, at origin
    # '''
    # maxn  = 30000 # maximum number of function evaluations allowed during optimization (default: 1000)
    # kstop = 20    # maximum number of evolution loops before convergency (default: 10)
    # bl = -600*np.ones(10)
    # bu = 600*np.ones(10)
    # x0 = np.random.rand(10)*600.
    # bestx, bestf, numfunccalls = sce(griewank, x0, bl, bu, maxn=maxn, outf=True, outcall=True, kstop=kstop)
    # print(bestx, bestf, numfunccalls)
