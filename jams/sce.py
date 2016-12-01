#!/usr/bin/env python
from __future__ import print_function
import subprocess
from distutils.util import strtobool
import numpy as np
from jams.const import huge
from jams import savez_compressed
# ToDo: write tmp/population files (of Fortran)
# ToDo: write out also in logfile if not None (use jams.tee as in joptimise)

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
        peps=0.001, seed=None, iniflg=True,
        alpha=0.8, beta=0.45, maxit=False, printit=2,
        outf=False, outhist=False, outcall=False,
        restart=False, restartfile1='sce.restart.npz', restartfile2='sce.restart.txt',
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
                peps=0.001, seed=None, iniflg=True,
                alpha=0.8, beta=0.45, maxit=False, printit=0,
                outf=False, outhist=False, outcall=False,
                restart=False, restartfile1='sce.restart.npz', restartfile2='sce.restart.txt',
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
        seed              Random number generator seed (default: None)
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
        restart           if True, continue from saved state in restartfile1/2
        restartfile1/2    File names for saving state of SCE (default: sce.restart.npz and sce.restart.txt)
                          State will be always written, except if restartfile1=None.
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
        >>> bestx, bestf, icall = sce(rosenbrock, x0, bl, bu, seed=1, maxn=1000, outf=True, outcall=True, printit=2, restartfile1=None)
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
        ...                           outf=True, outcall=True, restartfile1=None)
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
                  MC, Nov 2016 - restart - only Python 2
                  MC, Nov 2016 - restartfile1=None
                  MC, Nov 2016 - return -bestf if maxit
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
    if restartfile1 is not None:
        # Only arrays with savez_compressed - restartfile1
        restartarray  = ['bl', 'bu', 'bound', 'mask',
                         'criter',
                         'x', 'xf',
                         'bestx', 'worstx', 'allbestf', 'allbestx',
                         'rs2']
        # Save scalars in simple text file - restartfile2
        restartint    = ['nopt', 'npg', 'nps', 'nspl', 'mings', 'npt',
                         'nloop', 'icall', 'rs3', 'rs4']
        restartfloat  = ['gnrng', 'criter_change', 'bestf', 'worstf', 'rs5']
        restartbool   = ['maxit']
        restartstring = ['rs1']
        saveargarray = '"'+restartfile1+'"'
        for j in restartarray: saveargarray = saveargarray + ', '+j+'='+j
        saveargint    = ','.join(restartint)
        saveargfloat  = ','.join(restartfloat)
        saveargbool   = ','.join(restartbool)
        saveargstring = ','.join(restartstring)

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

    if not restart:
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

        # Seed random number generator
        np.random.seed(seed=seed)

        if mask is None: mask = np.ones(nopt, dtype=np.bool)

        large = 0.5*huge

        # Create an initial population to fill array x(npt,nopt):
        x = SampleInputMatrix(npt, nopt, bl, bu, distname='randomUniform')
        for i in range(npt):
            x[i,:] = np.where(mask, x[i,:], x0)
        if iniflg==1: x[0,:] = x0

        icall = 0
        xf = np.zeros(npt)
        for i in range (npt):
            fuc = call_function(functn, x[i,:], bl, bu, mask,
                                parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
            xf[i] = -fuc if maxit else fuc
            if printit==1: print('  f, X: ', xf[i], x[i,:])
            icall += 1

        # remember largest for treating of NaNs
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

        # save restart
        if restartfile1 is not None:
            rs1, rs2, rs3, rs4, rs5 = np.random.get_state()
            exec("savez_compressed("+saveargarray+")")
            p = open(restartfile2, 'w')
            exec("print("+saveargint+", file=p)")
            exec("print("+saveargfloat+", file=p)")
            exec("print("+saveargbool+", file=p)")
            exec("print("+saveargstring+", file=p)")
            p.close()

    else: # if no restart

        # load restart
        p1 = open(restartfile1, 'rb')
        pp = np.load(p1)
        for i in pp.files: exec(i+" = pp['"+i+"']")
        p1.close()
        p2 = open(restartfile2, 'r')
        for i, pp in enumerate(p2.readline().rstrip().split()): exec(restartint[i]+" = int(pp)")
        for i, pp in enumerate(p2.readline().rstrip().split()): exec(restartfloat[i]+" = float(pp)")
        for i, pp in enumerate(p2.readline().rstrip().split()): exec(restartbool[i]+" = bool(strtobool(pp))")
        for i, pp in enumerate(p2.readline().rstrip().split()): exec(restartstring[i]+" = pp")
        p2.close()
        np.random.set_state((rs1, rs2, rs3, rs4, rs5))

    # Outer Loop
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
                lcs    = np.zeros(nps, dtype=np.int)
                lcs[0] = 1
                for k3 in range(1,nps):
                    for i in range(1000):
                        # lpos = 1 + int(np.floor(npg+0.5-np.sqrt((npg+0.5)**2 - npg*(npg+1)*np.random.random())))
                        lpos = int(np.floor(npg+0.5-np.sqrt((npg+0.5)**2 - npg*(npg+1)*np.random.rand())))
                        # idx=find(lcs(1:k3-1)==lpos)
                        idx = (lcs[0:k3]==lpos).nonzero()  # check of element al eens gekozen
                        if idx[0].size == 0: break
                    lcs[k3] = lpos
                lcs.sort()

                # Construct the simplex:
                s  = np.zeros((nps,nopt))
                s  = cx[lcs,:]
                sf = cf[lcs]

                # remember largest for treating of NaNs
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
            print('Evolution loop {0:d}, trials {1:d}. Best f: {2:f}, worst f {3:f}'.format(nloop, icall, bestf, worstf))
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

        # save restart
        if restartfile1 is not None:
            rs1, rs2, rs3, rs4, rs5 = np.random.get_state()
            exec("savez_compressed("+saveargarray+")")
            p = open(restartfile2, 'w')
            exec("print("+saveargint+", file=p)")
            exec("print("+saveargfloat+", file=p)")
            exec("print("+saveargbool+", file=p)")
            exec("print("+saveargstring+", file=p)")
            p.close()
    # End of the Outer Loop: while icall<maxn and gnrng>peps and criter_change>pcento

    if printit<2:
        print('Search stopped at trial number {0:d} with normalized geometric range {1:f}. '.format(icall, gnrng))
        print('The best point has improved by {:f} in the last {:d} loops.'.format(criter_change, kstop))

    # reshape allbestx
    allbestx = allbestx.reshape(allbestx.size/nopt,nopt)

    # end of subroutine sce
    if maxit:
        bestf    *= -1.
        allbestf *= -1.
    out = [bestx]
    if outf:    out += [bestf]
    if outhist: out += [allbestx, allbestf]
    if outcall: out += [icall]

    # write parameter file with best parameters
    if isinstance(functn, (str,list)):
        parameterwriter(parameterfile, bestx, bl, bu, mask)

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    # maxn = 10000
    # from jams.functions import ackley, griewank, goldstein_price, rastrigin, rosenbrock, six_hump_camelback
    # '''
    # This is the Ackley Function
    # Global Optimum (n>=2): 0.0 at origin
    # '''
    # npara = 10
    # bl = -10*np.ones(npara)
    # bu = 10*np.ones(npara)
    # x0 = np.random.rand(npara)*10.
    # bestx, bestf = sce(ackley, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Ackley ', bestx, bestf)
    # '''
    #     This is the Griewank Function (2-D or 10-D)
    #     Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
    #        Global Optimum: 0, at origin
    # '''
    # npara = 10
    # bl = -600*np.ones(npara)
    # bu = 600*np.ones(npara)
    # x0 = np.random.rand(npara)*600.
    # bestx, bestf = sce(griewank, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Griewank ', bestx, bestf)
    # '''
    # This is the Goldstein-Price Function
    # Bound X1=[-2,2], X2=[-2,2]
    # Global Optimum: 3.0,(0.0,-1.0)
    # '''
    # npara = 2
    # bl = -2*np.ones(npara)
    # bu = 2*np.ones(npara)
    # x0 = np.random.rand(npara)*2.
    # bestx, bestf = sce(goldstein_price, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Goldstein ', bestx, bestf)
    # '''
    # This is the Rastrigin Function
    # Bound: X1=[-1,1], X2=[-1,1]
    # Global Optimum: -2, (0,0)
    # '''
    # npara = 2
    # bl = -1*np.ones(npara)
    # bu = 1*np.ones(npara)
    # x0 = np.random.rand(npara)*1.
    # bestx, bestf = sce(rastrigin, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Rastrigin ', bestx, bestf)
    # '''
    # This is the Rosenbrock Function
    # Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
    #        bl=[-5 -5]; bu=[5 5]; x0=[1 1];
    # '''
    # npara = 2
    # bl = -2*np.ones(npara)
    # bu = 5*np.ones(npara)
    # x0 = np.random.rand(npara)*2.
    # bestx, bestf = sce(rosenbrock, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Rosenbrock ', bestx, bestf)
    # '''
    # This is the Six-hump Camelback Function.
    # Bound: X1=[-5,5], X2=[-5,5]
    # True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    # '''
    # npara = 2
    # bl = -5*np.ones(npara)
    # bu = 5*np.ones(npara)
    # x0 = np.random.rand(npara)*5.
    # bestx, bestf = sce(six_hump_camelback, x0, bl, bu, maxn=maxn, outf=True, restartfile1=None)
    # print('Six_hump_camelback ', bestx, bestf)
