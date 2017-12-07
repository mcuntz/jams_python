#!/usr/bin/env python
from __future__ import print_function
from functools import partial
import subprocess
import numpy as np
from jams import morris_sampling, elementary_effects

def _ext_obj_wrapper(func, lb, ub, mask,
                     parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug, rank,
                     params):
    '''
        Wrapper function for external program to be optimised
        to be used with partial:
            obj = partial(_ext_obj_wrapper, func, lb, ub, mask,
                          parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug, rank)
        This allows then calling obj with only the argument params:
            fx = obj(params)


        Definition
        ----------
        def _ext_obj_wrapper(func, lb, ub, mask,
                             parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug, rank,
                             params):


        Input
        -----
        func              function to minimise (python function or string for external executable)
        lb                (npars) lower bounds of parameters
        ub                (npars) upper bounds of parameters
        mask              (npars) mask to include (1) or exclude (0) parameter from optimisation
        parameterfile     Parameter file for executable; must be given if func is name of executable
        parameterwriter   Python function for writing parameter file if func is name of executable
        objectivefile     File with objective value from executable; must be given if func is name of executable
        objectivereader   Python function for reading objective value if func is name of executable
        shell             If True, the specified command will be executed through the shell.
        debug             If True, model output is displayed for executable.
        rank              Process rank will be given as last argument to external program.
        params            (npars) parameter set


        Output
        ------
        Function value


        History
        -------
        Written,  MC, Nov 2016
    '''
    if isinstance(func, (str,list)):
        parameterwriter(parameterfile, params, lb, ub, mask)
        if rank is None:
            func1 = func
        else:
            if isinstance(func, str):
                func1 = [func, str(rank)]
            else:
                func1 = func+[str(rank)]
        if debug:
            err = subprocess.call(func1, shell=shell)
        else:
            err = subprocess.check_output(func1, shell=shell)
        obj = objectivereader(objectivefile)
        return obj
    else:
        return func(params)


# Function wrappers for objective and constraints
# used with functools.partial
def _obj_wrapper(func, arg, kwarg, x):
    '''
        Wrapper function for function to be optimised
        to be used with partial:
            obj = partial(_obj_wrapper, func, arg, kwarg)
        This allows then calling obj with only the argument x:
            fx = obj(x)
        which translates to:
            fx = func(x, *arg, **kwarg)
    '''
    return func(x, *arg, **kwarg)


# Screening
def screening(func, x0, lb, ub, mask=None,
              arg=(), kwarg={},
              nt=-1, nsteps=6, ntotal=-1,
              processes=1, pool=None,
              verbose=0,
              parameterfile=None, parameterwriter=None,
              objectivefile=None, objectivereader=None,
              shell=False, debug=False):
    """
        Morris Screening


        Definition
        ----------
        def screening(func, x0, lb, ub, mask=None,
                      arg=(), kwarg={},
                      nt=-1, nsteps=6, ntotal=-1,
                      processes=1,
                      verbose=0,
                      parameterfile=None, parameterwriter=None,
                      objectivefile=None, objectivereader=None,
                      shell=False, debug=False):


        Input
        -----
        func        python function or string for external executable
                    The objective function to screen.
        x0          1D-array
                    Will be taken at dimensions with mask==False.
        lb          1D-array
                    Lower bounds of parameters.
        ub          1D-array
                    Upper bounds of parameters.


        Optional Input
        --------------
        mask        1D-array (Default: include all parameters)
                    Include (1,True) or exclude (0,False) parameters in screening.
        arg         tuple (Default: empty tuple)
                    Additional arguments passed to objective and constraint functions (only for Python functions).
        kwarg       dict (Default: empty dict)
                    Additional keyword arguments passed to objective and constraint functions (only for Python functions).
        nt          int (Default: len(x0))
                    The number of trajectories to return.
        ntotal      int (Default: max(nt*10,nt**2))
                    Total number of trajectories for selecting optimal nt trajectories.
        nsteps      int (Default: 6)
                    Number of steps between lower and upper bounds.
        processes   int (Default: 1)
                    The number of processes to use to evaluate objective function and constraints.
        pool        schwimmbad pool object (Default: None)
                    Generic map function used from schwimmbad, which provides, serial, multiprocessor,
                    and MPI mapping functions. The pool is chosen with
                    schwimmbad.choose_pool(mpi=True/False, processes=nummultiprocessors).
                    The pool will be chosen automatically if pool is None.
                    MPI pools can only be opened and closed once. So if you want to use screening several
                    times in one program, then you have to choose the pool and later close the pool in the
                    wrapping progam and pass it to screening.                    
        parameterfile     string (Default: None)
                          Parameter file for executable; must be given if func is name of executable
        parameterwriter   function (Default: None)
                          Python function for writing parameter file if func is name of executable
        objectivefile     string (Default: None)
                          File with objective value from executable; must be given if func is name of executable
        objectivereader   function (Default: None)
                          Python function for reading objective value if func is name of executable
        shell             boolean (Default: False)
                          If True, the specified executable will be executed through the shell.
        debug             boolean (Default: False)
                          If True, model output is displayed for executable.


        Output
        ------
        if nt==1:
            1D-array - (nparameter,) with elementary effects of each parameter

        if nt>1:
            2D-array - (nparameter,3) with per parameter
                       1. mean of absolute elementary effects over all nt trajectories (mu*)
                       2. mean of elementary effects over all nt trajectories (mu)
                       3. standard deviation of elementary effects over all nt trajectories (sigma)


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

        Copyright 2017-2018 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2017
        Modified, MC, Jan 2018
    """
    # Get MPI communicator
    try:
        from mpi4py import MPI
        comm  = MPI.COMM_WORLD
        csize = comm.Get_size()
        crank = comm.Get_rank()
    except ImportError:
        comm  = None
        csize = 1
        crank = 0
    if csize > 1:
        passrank = crank
    else:
        passrank = None

    # Checks
    # Function
    assert hasattr(func, '__call__') or isinstance(func, (str,list)), 'Invalid function handle or external call.'
    # Bounds
    assert len(x0)==len(lb), 'Initial values and lower bounds must have the same lengths.'
    assert len(lb)==len(ub), 'Lower- and upper-bounds must have the same lengths.'
    x0 = np.array(x0)
    lb = np.array(lb)
    ub = np.array(ub)
    # Mask
    if mask is None:
        assert np.all(ub >= lb), 'All upper-bounds must be greater or equal than lower-bounds.'
    else:
        assert len(mask)==len(ub), 'Mask and bounds must have the same lengths.'
        if not np.all(mask):
            assert len(mask)==len(x0), 'Mask and x0 must have the same lengths.'
        assert np.all(ub[mask] >= lb[mask]), 'All unmasked upper-bounds must be greater or equal than lower-bounds.'
    # Parameterfile etc. keywords if func is name of executable
    if isinstance(func, (str,list)):
        if parameterfile is None:
            raise IOError('parameterfile must be given if func is name of executable.')
        else:
            if csize > 1:
                parameterfile1 = parameterfile + '.' + str(crank)
            else:
                parameterfile1 = parameterfile
        if parameterwriter is None:
            raise IOError('parameterwrite must be given if func is name of executable.')
        if objectivefile is None:
            raise IOError('objectivefile must be given if func is name of executable.')
        else:
            if csize > 1:
                objectivefile1 = objectivefile + '.' + str(crank)
            else:
                objectivefile1 = objectivefile
        if objectivereader is None:
            raise IOError('objectivereader must be given if func is name of executable.')

    # Set defaults
    npara = len(lb)
    if mask is None:
        imask  = np.ones(npara, dtype=np.bool)
    else:
        imask  = mask
    nmask = np.sum(imask)
    if nt <= 0:
        nt = npara
    if ntotal <= 0:
        nt = max(nt**2, 10*nt)

    # Partialise objective function
    if isinstance(func, (str,list)):
        obj = partial(_ext_obj_wrapper, func, lb, ub, imask,
                      parameterfile1, parameterwriter, objectivefile1, objectivereader, shell, debug, passrank)
    else:
        obj = partial(_obj_wrapper, func, arg, kwarg)

    # Sample trajectories
    if (crank==0) and (verbose > 0):
        print('Sample trajectories')
    tmatrix, tvec = morris_sampling(nmask, lb[imask], ub[imask], N=ntotal, p=nsteps, r=nt)

    # Set input vector to trajectories and masked elements = x0
    x = np.tile(x0,tvec.size).reshape(tvec.size,npara) # default to x0
    x[:,imask] = tmatrix                               # replaced unmasked with trajectorie values

    # Choose the right mapping function: single, multiprocessor or mpi
    if pool is None:
        import schwimmbad
        ipool = schwimmbad.choose_pool(mpi=False if csize==1 else True, processes=processes)
    else:
        ipool = pool

    # Calculate all model runs
    if (crank==0) and (verbose > 0):
        print('Calculate objective functions')
    fx = np.array(ipool.map(obj, x))
    if pool is None:
        ipool.close()

    # Calc elementary effects
    if (crank==0) and (verbose > 0):
        print('Calculate elementary effects')
    sa, res = elementary_effects(nmask, tmatrix, tvec, fx, p=nsteps)

    # Output with zero for all masked parameters
    if nt == 1:
        out = np.zeros(npara)
        out[imask] = sa[:,0]
    else:
        out = np.zeros((npara,3))
        out[imask,:] = res

    return out

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import schwimmbad
    # import jams

    # import time as ptime
    # t1 = ptime.time()
    
    # try:
    #     from mpi4py import MPI
    #     comm  = MPI.COMM_WORLD
    #     csize = comm.Get_size()
    #     crank = comm.Get_rank()
    # except ImportError:
    #     comm  = None
    #     csize = 1
    #     crank = 0


    # seed = 1025

    # if seed is not None:
    #     np.random.seed(seed=seed)

    # if False:
    #     # work with global parameters
    #     def morris(X):
    #         return jams.functions.morris(X, beta0, beta1, beta2, beta3, beta4)
    #     func = morris
    #     npars = 20
    #     beta0              = 0.
    #     beta1              = np.random.standard_normal(npars)
    #     beta1[:10]         = 20.
    #     beta2              = np.random.standard_normal((npars,npars))
    #     beta2[:6,:6]       = -15.
    #     beta3              = np.zeros((npars,npars,npars))
    #     beta3[:5,:5,:5]    = -10.
    #     beta4              = np.zeros((npars,npars,npars,npars))
    #     beta4[:4,:4,:4,:4] = 5.
    #     x0    = np.ones(npars)
    #     lb    = np.zeros(npars)
    #     ub    = np.ones(npars)
    #     nt      = 20
    #     ntotal  = 10*nt
    #     nsteps  = 6
    #     verbose = 1

    #     out = screening(func, x0, lb, ub, mask=None, nt=nt, nsteps=nsteps, ntotal=ntotal, processes=4, verbose=1)

    #     t2    = ptime.time()
        
    #     if crank == 0:
    #         strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    #         print('Time elapsed: ', strin)
    #         print('mu*, mu, std: ', out)
    #         out1 = out / out.max(axis=0)
    #         mustar = out1[:,0]
    #         mustar = np.sort(mustar)
    #         import matplotlib.pyplot as plt
    #         plt.plot(mustar, 'ro')
    #         plt.show()

    # if False:
    #     # pass parameters
    #     func = jams.functions.morris
    #     npars = 20
    #     x0    = np.ones(npars)
    #     lb    = np.zeros(npars)
    #     ub    = np.ones(npars)
    #     beta0              = 0.
    #     beta1              = np.random.standard_normal(npars)
    #     beta1[:10]         = 20.
    #     beta2              = np.random.standard_normal((npars,npars))
    #     beta2[:6,:6]       = -15.
    #     beta3              = np.zeros((npars,npars,npars))
    #     beta3[:5,:5,:5]    = -10.
    #     beta4              = np.zeros((npars,npars,npars,npars))
    #     beta4[:4,:4,:4,:4] = 5.
    #     args = [beta0, beta1, beta2, beta3, beta4] # Morris
    #     nt = 19
    #     ntotal  = 10*nt
    #     nsteps = 6
    #     verbose = 1

    #     out = screening(func, x0, lb, ub, arg=args, mask=None, nt=nt, nsteps=nsteps, ntotal=ntotal, processes=4, verbose=1)
    #     t2    = ptime.time()
        
    #     if crank == 0:
    #         strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    #         print('Time elapsed: ', strin)
    #         print('mu*, mu, std: ', out)
    #         out1 = out / out.max(axis=0)
    #         mustar = out1[:,0]
    #         mustar = np.sort(mustar)
    #         import matplotlib.pyplot as plt
    #         plt.plot(mustar, 'ro')
    #         plt.show()

    # if True:
    #     # pass parameters and pool
    #     processes = 4
    #     pool = schwimmbad.choose_pool(mpi=False if csize==1 else True, processes=processes)
    #     func = jams.functions.morris
    #     npars = 20
    #     x0    = np.ones(npars)
    #     lb    = np.zeros(npars)
    #     ub    = np.ones(npars)
    #     beta0              = 0.
    #     beta1              = np.random.standard_normal(npars)
    #     beta1[:10]         = 20.
    #     beta2              = np.random.standard_normal((npars,npars))
    #     beta2[:6,:6]       = -15.
    #     beta3              = np.zeros((npars,npars,npars))
    #     beta3[:5,:5,:5]    = -10.
    #     beta4              = np.zeros((npars,npars,npars,npars))
    #     beta4[:4,:4,:4,:4] = 5.
    #     args = [beta0, beta1, beta2, beta3, beta4] # Morris
    #     nt = 19
    #     ntotal  = 10*nt
    #     nsteps = 6
    #     verbose = 1

    #     out = screening(func, x0, lb, ub, arg=args, mask=None, nt=nt, nsteps=nsteps, ntotal=ntotal,
    #                     processes=processes, pool=pool, verbose=1)
    #     t2    = ptime.time()
        
    #     if crank == 0:
    #         strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    #         print('Time elapsed: ', strin)
    #         print('mu*, mu, std: ', out)
    #         out1 = out / out.max(axis=0)
    #         mustar = out1[:,0]
    #         mustar = np.sort(mustar)
    #         import matplotlib.pyplot as plt
    #         plt.plot(mustar, 'ro')
    #         plt.show()
    #     pool.close()

