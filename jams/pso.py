#!/usr/bin/env python
from __future__ import print_function
from functools import partial
import subprocess
from distutils.util import strtobool
import numpy as np
from jams.const import huge
from jams import savez_compressed
from mpi4py import MPI
# ToDo:
#   Handling constraints
#   write tmp/population files (as in SCE of Fortran)
#   write out also in logfile if not None (use jams.tee as in joptimise)
#   crossover with quadratic function (QIPSO) ?

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
        parameterfile     Parameter file for executable; must be given if functn is name of executable
        parameterwriter   Python function for writing parameter file if functn is name of executable
        objectivefile     File with objective value from executable; must be given if functn is name of executable
        objectivereader   Python function for reading objective value if functn is name of executable
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


def _is_feasible_wrapper(func, x):
    '''
        Checking that all constraints are met.

        func is one of the 'partialised' constraint functions
            _cons_none_wrapper, _cons_ieqcons_wrapper, and _cons_f_ieqcons_wrapper.

        This function is partialised:
            is_feasible = partial(_is_feasible_wrapper, cons)
        which is unnecessary for now. It could still simply be:
            is_feasible = cons
    '''
    return np.all(func(x)>=0.)


def _cons_none_wrapper(x):
    '''
        Dummy constraints function if no constraints.
        Returns 0.
    '''
    return np.array([0])


def _cons_ieqcons_wrapper(ieqcons, arg, kwarg, x):
    '''
        Wrapper function for constraints to be used with partial in the case
        that ieqcons is given, i.e. a list of constraint functions, returning
        >=0. if constraints are met.
    '''
    return np.array([y(x, *arg, **kwarg) for y in ieqcons])


def _cons_f_ieqcons_wrapper(f_ieqcons, arg, kwarg, x):
    '''
        Wrapper function for constraints to be used with partial in the case
        that f_ieqcons is given, i.e. a single function returning an 1-D array
        values >=0. where constraints are met.
    '''
    return np.array(f_ieqcons(x, *arg, **kwarg))


def diversity(x, lb, ub, mask):
    '''
        Diversity in swarm positions.
        From [Pant et al., 2009] http://www.iaeng.org/IJCS/issues_v36/issue_2/IJCS_36_2_02.pdf
        changed to RMSE.

        RMSE is generally about 1/2 or 1/3 of normalized geometric range (SCE).
    '''
    # all dimensions [0,1]
    x01 = np.where(mask, (x-lb)/(ub-lb), x)

    S = float(x01.shape[0])
    D = float(x01.shape[1])

    # average particle position
    pmean = np.sum(x01, axis=0)/S

    # # Eq. 5 of Pant et al. (2009)
    # div = np.sum(np.sqrt(np.sum((x01-pmean)**2,axis=1)),axis=0) / S

    # average RMSE
    div = np.sum(np.sqrt(np.sum((x01-pmean)**2,axis=1)/D),axis=0) / S

    # # Normalized geometric range of SCE
    # imax = x01.max(axis=0)
    # imin = x01.min(axis=0)
    # imaxmin = np.ma.array(imax-imin, mask=(imax==imin))
    # div = np.ma.exp(np.ma.mean(np.ma.log(imaxmin)))

    return div


def get_neighbor_indeces(n, S, topology, kl=1):
    '''
        Get the indices of the neighbors for the current particle given the topology.

        Input
        -----
        n           integer
                    Particle index.
        S           integer
                    Swarm size.
        topology    string
                    Neighborhood topologies. These are rather social than geographical topologies.
                    All neighborhoods comprise the current particle as well.
                    [Kennedy & Mendes, 2002] http://dx.doi.org/10.1109/CEC.2002.1004493
                    'gbest'    Neighborhood is entire swarm.
                    'lbest'    Partciles aranged in a ring, in which each particle communicates with
                               kl particles on each side, i.e. particle i has the neighborhood
                               i-kl, i-kl+1, ..., i, i+1, ..., i+kl-1, i+kl
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80
                    'ring'     'lbest' with kl=1
                    'neumann'  Neighborhood of a point including all points at a Hamming distance of 1.
                               Particles are arranges in a lattice, where each particle interacts with
                               its immediate 4 neighbors to the N, S, E, and W.
                               [Kennedy and Mendes, 2006] http://dx.doi.org/10.1109/TSMCC.2006.875410
                               The von Neumann neighborhood is configured into r rows and c columns,
                               where r is the highest integer less than or equal to sqrt(n) that evenly
                               divides n and c = n / r
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80


        Optional Input
        --------------
        kl          integer
                    Neighborhood distance in topology 'lbest'.
                    (Default: 1 = ring)
    '''
    if topology.lower() == 'ring': kl=1

    if topology.lower() == 'gbest':
        ii = list(np.arange(S))
    elif (topology.lower() == 'lbest') or (topology.lower() == 'ring'):
        ii = [n]                  # self
        for kk in range(1,kl+1):
            ii.append((n-kk) % S) # left
            ii.append((n+kk) % S) # right
    elif topology.lower() == 'neumann':
        rows = int(np.floor(np.sqrt(S)))
        while (S % rows) != 0: rows -= 1
        cols = S // rows
        left  = lambda x, c: (x - 1) % c + x//c * c
        right = lambda x, c: (x + 1) % c + x//c * c
        above = (n - rows) % S
        below = (n + rows) % S
        ii = [left(above, cols),  above, right(above, cols),
              left(n, cols),      n,     right(n, cols),
              left(below, cols),  below, right(below, cols)]

    return ii


def get_best_neighbor(p, fp, topology, kl=1):
    '''
        Get the best neighbor for a given topology

        Input
        -----
        p           ND-array
                    The best known position of each particle.
        fp          1D-array
                    The objective values at each position in p.
        topology    string
                    Neighborhood topologies. These are rather social than geographical topologies.
                    All neighborhoods comprise the current particle as well.
                    [Kennedy & Mendes, 2002] http://dx.doi.org/10.1109/CEC.2002.1004493
                    'gbest'    Neighborhood is entire swarm.
                    'lbest'    Partciles aranged in a ring, in which each particle communicates with
                               kl particles on each side, i.e. particle i has the neighborhood
                               i-kl, i-kl+1, ..., i, i+1, ..., i+kl-1, i+kl
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80
                    'ring'     'lbest' with kl=1
                    'neumann'  Neighborhood of a point including all points at a Hamming distance of 1.
                               Particles are arranges in a lattice, where each particle interacts with
                               its immediate 4 neighbors to the N, S, E, and W.
                               [Kennedy and Mendes, 2006] http://dx.doi.org/10.1109/TSMCC.2006.875410
                               The von Neumann neighborhood is configured into r rows and c columns,
                               where r is the highest integer less than or equal to sqrt(n) that evenly
                               divides n and c = n / r
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80


        Optional Input
        --------------
        kl          integer
                    Neighborhood distance in topology 'lbest'.
                    (Default: 1 = ring)
    '''
    if topology.lower() == 'ring': kl=1
    S = p.shape[0]
    D = p.shape[1]

    g  = np.ones((S,D))*np.inf
    fg = np.ones(S)*np.inf
    if topology.lower() == 'gbest':
        i_min = np.argmin(fp)
        ig  = p[i_min,:] # overall best
        ifg = fp[i_min]
        for ss in range(S):
            g[ss,:] = ig
            fg[ss]  = ifg
    elif (topology.lower() == 'lbest') or (topology.lower() == 'ring') or (topology.lower() == 'neumann'):
        g  = np.ones((S,D))*np.inf
        fg = np.ones(S)*np.inf
        for ss in range(S):
            ii = get_neighbor_indeces(ss, S, topology, kl=kl)
            pp  = p[ii,:]
            fpp = fp[ii]
            i_min = np.argmin(fpp)
            g[ss,:] = pp[i_min,:].copy()
            fg[ss]  = fpp[i_min]

    return [g, fg]


# Particle Swarm Optimisation
def pso(func, x0, lb, ub,
        mask=None,
        ieqcons=None, f_ieqcons=None,
        arg=(), kwarg={},
        swarmsize=None, inertia=None, phip=None, phig=None, maxn=250,
        minstep=1e-8, minobj=1e-8, maxit=False,
        init='lhs', strategy='fips', topology='gbest', kl=1,
        includex0=False, seed=None,
        processes=1,
        verbose=0, pout=False, cout=False,
        restart=False, restartfile1='pso.restart.npz', restartfile2='pso.restart.txt',
        parameterfile=None, parameterwriter=None,
        objectivefile=None, objectivereader=None,
        shell=False, debug=False):
    """
        Particle Swarm Optimization (PSO)


        Definition
        ----------
        def pso(func, x0, lb, ub,
                mask=None,
                ieqcons=None, f_ieqcons=None,
                arg=(), kwarg={},
                swarmsize=None, inertia=None, phip=None, phig=None, maxn=250,
                minstep=1e-8, minobj=1e-8, maxit=False,
                init='lhs', strategy='fips', topology='gbest', kl=1,
                includex0=False, seed=None,
                processes=1,
                verbose=0, pout=False, cout=False,
                restart=False, restartfile1='pso.restart.npz', restartfile2='pso.restart.txt',
                parameterfile=None, parameterwriter=None,
                objectivefile=None, objectivereader=None,
                shell=False, debug=False):


        Input
        -----
        func        python function or string for external executable
                    The function to be minimized.
        x0          1D-array
                    Will be taken at dimensions with mask==False.
        lb          1D-array
                    The lower bounds of the design variable(s).
        ub          1D-array
                    The upper bounds of the design variable(s).


        Optional Input
        --------------
        mask        1D-array
                    include (1,True) or exclude (0,False) parameters in optimisation.
                    (Default: include all dimensions)
        ieqcons     list
                    A list of functions of length n such that ieqcons[j](x,*arg) >= 0.0 in
                    a successfully optimized problem.
                    (Default: None)
        f_ieqcons   function
                    Returns a 1-D array in which each element must be greater or equal
                    to 0.0 in a successfully optimized problem. If f_ieqcons is specified,
                    ieqcons is ignored.
                    (Default: None)
        arg         tuple
                    Additional arguments passed to objective and constraint functions.
                    (Default: empty tuple)
        kwarg       dict
                    Additional keyword arguments passed to objective and constraint functions.
                    (Default: empty dict)
        swarmsize   int
                    The number of particles in the swarm.
                    (Default: max(min(3*len(lb),40),10))
        inertia     scalar
                    Particle velocity scaling factor.
                    Default depends on algorithm (strategy).
        phip        scalar
                    Scaling factor to search away from the particle's best known position.
                    Default depends on algorithm (strategy).
        phig        scalar
                    Scaling factor to search away from the swarm's (neighbor's) best known position.
                    Default depends on algorithm (strategy).
        maxn        int
                    The maximum number of iterations for the swarm to search.
                    (Default: 250)
        minstep     scalar
                    The minimum stepsize of swarm's best position before the search terminates.
                    (Default: 1e-8)
        minobj      scalar
                    Objective function defining convergence. Attention at maxit=True.
                    (Default: 1e-8)
        maxit       boolean
                    Minimise or maximise func.
                    (Default: False)
                    False: minimise objective function down to minobj
                    True:  maximise objective function, i.e. minimise -objective function down to minobj
        init        string
                    How to sample the initial swarm positions and velocities.
                    (Default: 'lhs')
                    'random': random sampling from uniform distribution
                    'lhs':    latin hypercube sampling from uniform distributions
                    'sobol':  quasirandom Sobol sequence (only up to 40 dimensions)
        strategy    string
                    PSO variants.
                    (Default: 'fips')
                    'original':   Textbook particle swarm algorithm with inertia weight
                                  x = current position
                                  p = particles best position
                                  g = neighborhood best position
                                  rg, rp = np.random.uniform(size=(S,D))
                                  inertia=0.5, phip=2., phig=2.
                                  v = inertia*v + rp*phip*(p-x) + rg*phig*(g-x)
                                  x = x + v
                    'inertia':    Same as 'original' but with inertia weight decreasing from 0.9 to 0.4
                                  over time, i.e. over iterations (it)
                                  imax = 0.9
                                  imin = 0.4
                                  inertia = imax - float(it)/float(maxn-1) * (imax-imin)
                    'canonical':  Clerk & Kennedy (2000) with fixed constriction factor
                                  From PaGMO (esa.github.io/pagmo):
                                  Clerc's analysis of the iterative system led him to propose a strategy for the
                                  placement of "constriction coefficients" on the terms of the formulas; these
                                  coefficients controlled the convergence of the particle and allowed an elegant and
                                  well-explained method for preventing explosion, ensuring convergence, and
                                  eliminating the arbitrary Vmax parameter. The analysis also takes the guesswork
                                  out of setting the values of phi_1 and phi_2.
                                  "This is the canonical particle swarm algorithm of today."
                                  [Poli et al., 2007] http://dx.doi.org/10.1007/s11721-007-0002-0
                                  [Clerc & Kennedy, 2002] http://dx.doi.org/10.1109/4235.985692
                                  inertia=0.7289, phip=2.05, phig=2.05
                                  v = inertia * (v + phip*rp*(p - x) + phig*rg*(g - x))
                                  x = x + v
                    'fips':       Fully Informed Particle Swarm
                                  From PaGMO (esa.github.io/pagmo):
                                  Whereas in the traditional algorithm each particle is affected by its own
                                  previous performance and the single best success found in its neighborhood, in
                                  Mendes' fully informed particle swarm (FIPS), the particle is affected by all its
                                  neighbors, sometimes with no influence from its own previous success.
                                  "With good parameters, FIPS appears to find better solutions in fewer iterations
                                  than the canonical algorithm, but it is much more dependent on the population
                                  topology."
                                  [Poli et al., 2007] http://dx.doi.org/10.1007/s11721-007-0002-0
                                  [Mendes et al., 2004] http://dx.doi.org/10.1109/TEVC.2004.826074
                                  ri = np.random.uniform(size=nneighbor)
                                  inertia = 0.7289
                                  acc_coeff = phip+phig = 4.1
                                  v = inertia * (v + np.sum(ri*acc_coeff/nneighbor*(p[:]-x)))
                    'nips':       Neighborhood Informed Particle Swarm
                                  'fips' but particles are not informed by all other particles
                                  but only the particles in its neighborhood given by topology.
                                  ri = np.random.uniform(size=nneighbor)
                                  inertia = 0.7289, acc_coeff = phip+phig = 4.1
                                  v = inertia * (v + np.sum(ri*acc_coeff/nneighbor*(p[neighbor[:]]-x)))
        topology    string
                    Neighborhood topologies. These are rather social than geographical topologies.
                    All neighborhoods comprise the current particle as well.
                    [Kennedy & Mendes, 2002] http://dx.doi.org/10.1109/CEC.2002.1004493
                    (Default: 'gbest')
                    'gbest'    Neighborhood is entire swarm.
                    'lbest'    Partciles aranged in a ring, in which each particle communicates with
                               kl particles on each side, i.e. particle i has the neighborhood
                               i-kl, i-kl+1, ..., i, i+1, ..., i+kl-1, i+kl
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80
                    'ring'     'lbest' with kl=1
                    'neumann'  Neighborhood of a point including all points at a Hamming distance of 1.
                               Particles are arranges in a lattice, where each particle interacts with
                               its immediate 4 neighbors to the N, S, E, and W.
                               [Kennedy and Mendes, 2006] http://dx.doi.org/10.1109/TSMCC.2006.875410
                               The von Neumann neighborhood is configured into r rows and c columns,
                               where r is the highest integer less than or equal to sqrt(n) that evenly
                               divides n and c = n / r
                               [Mohais et al., 2005] http://dx.doi.org/10.1007/11589990_80
        kl          integer
                    Neighborhood distance in topology 'lbest'.
                    (Default: 1 = ring)
        verbose     integer
                    Controlling amount of print-out.
                    (Default: 0)
                    0: No print-out.
                    1: Printing convergence criteria, etc.
                    2: Printing after each step.
        includex0   boolean
                    True: include x0 in initial swarm
                    (Default: False)
        seed        int or array_like
                    Seed for numpy's random number generator.
                    (Default: None)
        processes   int
                    The number of processes to use to evaluate objective function and constraints.
                    (Default: 1)
        pout        boolean
                    True: include best per-particle positions and their objective values in output.
                    (Default: False)
        cout        boolean
                    True: include number of function calls in output.
                    (Default: False)
        restart           boolean
                          if True, continue from saved state in restartfile1/2.
                          (Default: False)
        restartfile1/2    string
                          File names for saving current state of PSO.
                          (Default: pso.restart.npz and pso.restart.txt)
                          State will be always written, except if restartfile1=None.
        parameterfile     string
                          Parameter file for executable; must be given if func is name of executable
                          (Default: None)
        parameterwriter   function
                          Python function for writing parameter file if func is name of executable
                          (Default: None)
        parameterfile     string
                          Parameter file for executable; must be given if func is name of executable
                          (Default: None)
        parameterwriter   function
                          Python function for writing parameter file if func is name of executable
                          (Default: None)
        objectivefile     string
                          File with objective value from executable; must be given if func is name of executable
                          (Default: None)
        objectivereader   function
                          Python function for reading objective value if func is name of executable
                          (Default: None)
        shell             boolean
                          If True, the specified executable will be executed through the shell.
                          (Default: False)
        debug             boolean
                          If True, model output is displayed for executable.
                          (Default: False)


        Output
        ------
        g           1D-array
                    The swarm's best known position.
        fg          scalar
                    The objective value at g.
        p           ND-array
                    The best known position of each particle.
        fp          1D-array
                    The objective values at each position in p.


        License
        -------
        The original code of Abraham Lee (pyswarm) was published under the BSD license.

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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  Abraham Lee, 2013 - https://github.com/tisimst/pyswarm
        Modified, MC, Nov 2016 - adapted to JAMS package
                  MC, Nov 2016 - Changed defaults from swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxn=100
                               - Include vmax for original: clip(v, vmin, vmax)
                               - Sobol sequences and latin hypercube sampling for initial swarm positions
                               - Different PSO algorithms: original, decreasing inertia weights, constricted
                                 fully informed
                               - Stop if function below minobj
                               - debug -> verbose, maxit, cout
                               - neighborhoods
                               - external function - mask, x0, parameterfile, parameterwriter,
                                                     objectivefile, objectivereader, shell, debug
                  MC, Dec 2016 - includex0, restart, mpi
    """
    # Get MPI communicator
    comm = MPI.COMM_WORLD
    csize = comm.Get_size()
    crank = comm.Get_rank()
    if csize > 1:
        passrank = crank
    else:
        passrank = None

    dodiv = False

    # Different variabels types (array, float, int, ...) for restart
    if restartfile1 is not None:
        # Only arrays with savez_compressed - restartfile1
        restartarray  = ['lb', 'ub', 'mask1', 'mask2', 'x02',
                         'v', 'x', 'fx', 'fs', 'p', 'fp', 'gp', 'fgp',
                         'rs2']
        if dodiv: restartarray.extend(['gx'])
        # Save scalars in simple text file - restartfile2
        restartint    = ['D', 'S', 'it', 'iS', 'crank',
                         'rs3', 'rs4']
        restartfloat  = ['rs5']
        restartbool   = ['maxit']
        restartstring = ['rs1']
        if csize > 1:
            restartfile1 = restartfile1[0:restartfile1.rfind(".")] + '.' + str(crank) + restartfile1[restartfile1.rfind("."):]
            restartfile2 = restartfile2[0:restartfile2.rfind(".")] + '.' + str(crank) + restartfile2[restartfile2.rfind("."):]
        saveargarray = '"'+restartfile1+'"'
        for j in restartarray: saveargarray = saveargarray + ', '+j+'='+j
        saveargint    = ','.join(restartint)
        saveargfloat  = ','.join(restartfloat)
        saveargbool   = ','.join(restartbool)
        saveargstring = ','.join(restartstring)

    # Checks
    # Function
    assert hasattr(func, '__call__') or isinstance(func, (str,list)), 'Invalid function handle or external call.'
    # Bounds
    assert len(lb)==len(ub), 'Lower- and upper-bounds must have the same lengths.'
    lb = np.array(lb)
    ub = np.array(ub)
    assert np.all(ub >= lb), 'All upper-bounds must be greater than lower-bounds.'
    # Mask
    if mask is not None:
        assert len(mask)==len(ub), 'Mask and bounds must have the same lengths.'
        if not np.all(mask):
            assert len(mask)==len(x0), 'Mask and x0 must have the same lengths.'
    # Initialisation keyword
    inits = ['random', 'lhs', 'sobol']
    assert init.lower() in inits, 'Initialisation {:} not in {:}'.format(init, inits)
    if init.lower() == 'sobol':
        assert len(lb) <= 40, "Sobol' sequences only work up to 40 dimensions."
    # Strategy keyword
    ptypes = ['original', 'inertia', 'canonical', 'fips', 'nips']
    assert strategy.lower() in ptypes, 'PSO implementation {:} not in {:}'.format(strategy, ptypes)
    # Topology keyword
    ttypes = ['gbest', 'lbest', 'ring', 'neumann']
    assert topology.lower() in ttypes, 'Topology {:} not in {:}'.format(topology, ttypes)
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

    # Set defaults per strategy
    if strategy.lower() == 'original':    # Kennedy & Eberhart, 2001
        if inertia is None: inertia=0.5
        if phip  is None: phip=2.
        if phig  is None: phig=2.
    elif strategy.lower() == 'inertia':   # Shi & Eberhart (1998)
        imax = 0.9
        imin = 0.4
        if phip is None: phip=2.
        if phig is None: phig=2.
    elif strategy.lower() == 'canonical': # Clerc & Kennedy (2000)
        if inertia is None: inertia=0.7289
        if phip  is None: phip=2.05
        if phig  is None: phig=2.05
    elif strategy.lower() == 'fips':      # Mendes & Kennedy (2004)
        if inertia is None: inertia=0.7289
        if phip  is None: phip=2.05
        if phig  is None: phig=2.05
    elif strategy.lower() == 'nips':
        if inertia is None: inertia=0.7289
        if phip  is None: phip=2.05
        if phig  is None: phig=2.05

    if not restart:
        # Problem sizes
        D = len(lb) # dimension of each particle
        if swarmsize is None:
            S = max(min(3*D,40),10)
            S = max(min(S//csize*csize, (40//csize+1)*csize), (10//csize+1)*csize)
        else:
            S = swarmsize

        # Local swarmsize
        if S % csize != 0:
            raise ValueError("Swarmsize "+str(S)+" must be multiple of number of processes "+str(csize)+".")
        iS = S//csize # local swarmsize

        # Initialise 1D mask
        if mask is not None:
            mask1 = mask
        else:
            mask1 = np.ones(D, dtype=np.bool)

        # Partialise objective function
        if isinstance(func, (str,list)):
            obj = partial(_ext_obj_wrapper, func, lb, ub, mask1,
                          parameterfile1, parameterwriter, objectivefile1, objectivereader, shell, debug, passrank)
        else:
            obj = partial(_obj_wrapper, func, arg, kwarg)

        # Check for constraint function(s) and partialise them
        if f_ieqcons is None:
            if ieqcons is None:
                if (verbose>=1) and (crank == 0):
                    print('No constraints given.')
                cons = _cons_none_wrapper
            else:
                if (verbose>=1) and (crank == 0):
                    print('Converting ieqcons to a single constraint function.')
                cons = partial(_cons_ieqcons_wrapper, ieqcons, arg, kwarg)
        else:
            if (verbose>=1) and (crank == 0):
                print('Single constraint function given in f_ieqcons.')
            cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, arg, kwarg)
        is_feasible = partial(_is_feasible_wrapper, cons)

        # Initialize the multiprocessing module if necessary
        if processes > 1:
            import multiprocessing
            mp_pool = multiprocessing.Pool(processes)

        # 2D mask
        mask2 = np.tile(mask1,iS).reshape((iS,D))
        x02   = np.tile(x0,iS).reshape((iS,D))

        # Deal with NaN and Inf
        large = huge / S

        # Seed random number generator
        if crank == 0:
            np.random.seed(seed=seed)

        # Initialize the particle swarm
        # current particle positions
        # current particle velocities
        if init.lower() == 'random':
            # Random numbers only on rank 0 for reproducible results
            if crank == 0:
                rand = np.random.uniform(size=(2*S,D))
            else:
                rand = np.empty((2*S,D), dtype=np.float)
            # Scatter has different ordering than needed for reproducible results
            # Do it manually
            comm.Bcast(rand, root=0)
            x = rand[crank*iS:crank*iS+iS,:]
            v = rand[S+crank*iS:S+crank*iS+iS,:]
        elif init.lower() == 'sobol':
            import sobol
            nskip = D*S + crank*D*iS
            x = sobol.i4_sobol_generate(D,iS,nskip).transpose()
            nskip = 2*D*S + crank*D*iS
            v = sobol.i4_sobol_generate(D,iS,nskip).transpose()
        elif init.lower() == 'lhs':
            x = np.empty((iS,D), dtype=np.float64)
            v = np.empty((iS,D), dtype=np.float64)
            if crank == 0:
                import scipy.stats as stats
                from jams.lhs import lhs
                dist = [stats.uniform for i in range(D)]
                pars = [(0,1) for i in range(D)]
                gx   = lhs(dist, pars, S).transpose()
                gv   = lhs(dist, pars, S).transpose()
            else:
                gx = np.empty((S,D), dtype=np.float64)
                gv = np.empty((S,D), dtype=np.float64)
            comm.Scatter([gx, MPI.DOUBLE], [x, MPI.DOUBLE])
            comm.Scatter([gv, MPI.DOUBLE], [v, MPI.DOUBLE])
        fx  = np.ones(iS)*large                            # local current particles function values
        if dodiv: gx  = np.empty((iS,D), dtype=np.float64) # global current individual particle positions
        fs  = np.zeros(iS, dtype=bool)                     # current combined feasibility for each local particle
        p   = np.ones((iS,D), dtype=np.float64)*large      # local particles individual best positions
        fp  = np.ones(iS, dtype=np.float64)*large          # local particles individual best function values
        gp  = np.ones((S,D), dtype=np.float64)*large       # global particles individual best positions
        fgp = np.ones(S, dtype=np.float64)*large           # global particles individual best function values

        # Maximum velocity
        vmax = np.abs(ub - lb)
        vmin = -vmax

        # Initialize particle positions and velocities
        v = vmin + v*(vmax-vmin)
        x = lb + x*(ub - lb)
        if (crank == 0) and includex0: x[-1,:] = x0
        x = np.where(mask2, x, x02)

        # Calculate first objective and constraints for each particle
        if processes > 1:
            fs = np.array(mp_pool.map(is_feasible, x))
            ii = np.where(fs)[0]
            if ii.size > 0:
                fx[ii] = np.array(mp_pool.map(obj, x[ii,:]))
        else:
            for i in range(iS):
                fs[i] = is_feasible(x[i,:])
                if fs[i]:
                    print(i, x[i:])
                    fx[i] = obj(x[i,:])
        # maximise
        if maxit: fx *= -1.

        # NaN/Inf - ToDo: Check
        large = max(fp.max(), fx[np.isfinite(fx)].max())
        large = 1.1*large if large>0. else 0.9*large
        fx = np.where(np.isfinite(fx), fx, large)

        # print(1, x[np.where(fs)[0],:])
        # Store particles best positions (if constraints are satisfied)
        i_update = (fx < fp) & fs
        if np.any(i_update):
            p[i_update,:] = x[i_update,:].copy()
            fp[i_update]  = fx[i_update]

        # gather local best particles into global best particles
        comm.Allgather([p, MPI.DOUBLE], [gp, MPI.DOUBLE])
        comm.Allgather([fp, MPI.DOUBLE], [fgp, MPI.DOUBLE])

        if dodiv:
            comm.Allgather([x, MPI.DOUBLE], [gx, MPI.DOUBLE])
            gmask = np.tile(mask1,S).reshape((S,D))
            divers = diversity(gx, lb, ub, gmask)
            if crank == 0: print(1, divers, fgp.min())

        # Iterate until termination criterion met
        it = 1

        # save restart
        if restartfile1 is not None:
            if crank == 0:
                rs1, rs2, rs3, rs4, rs5 = np.random.get_state()
            else:
                rs1, rs2, rs3, rs4, rs5 = 'MT19937', np.array(624, dtype=np.uint), 0, 0, 0.
            exec("savez_compressed("+saveargarray+")")
            p2 = open(restartfile2, 'w')
            exec("print("+saveargint+", file=p2)")
            exec("print("+saveargfloat+", file=p2)")
            exec("print("+saveargbool+", file=p2)")
            exec("print("+saveargstring+", file=p2)")
            p2.close()

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
        if crank == 0:
            np.random.set_state((rs1, rs2, rs3, rs4, rs5))

        # Partialise objective function
        if isinstance(func, (str,list)):
            obj = partial(_ext_obj_wrapper, func, lb, ub, mask1,
                          parameterfile1, parameterwriter, objectivefile1, objectivereader, shell, debug, passrank)
        else:
            obj = partial(_obj_wrapper, func, arg, kwarg)

        # Check for constraint function(s) and partialise them
        if f_ieqcons is None:
            if ieqcons is None:
                if (verbose>=1) and (crank == 0):
                    print('No constraints given.')
                cons = _cons_none_wrapper
            else:
                if (verbose>=1) and (crank == 0):
                    print('Converting ieqcons to a single constraint function.')
                cons = partial(_cons_ieqcons_wrapper, ieqcons, arg, kwarg)
        else:
            if (verbose>=1) and (crank == 0):
                print('Single constraint function given in f_ieqcons.')
            cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, arg, kwarg)
        is_feasible = partial(_is_feasible_wrapper, cons)

        if processes > 1:
            import multiprocessing
            mp_pool = multiprocessing.Pool(processes)

    # Iterate swarm
    while it < maxn:
        # Stop if minimum found
        if fgp.min() < minobj:
            if (verbose>=1) and (crank == 0): print('minobj found.')
            break

        # Update neighbors best positions
        g  = np.empty((iS,D), dtype=np.float64) # local best neighbors
        fg = np.empty(iS, dtype=np.float64)
        if crank == 0:
            gg, fgg = get_best_neighbor(gp, fgp, topology, kl=kl)  # global best neighbors
        else:
            gg  = np.empty((S,D), dtype=np.float64)
            fgg = np.empty(S, dtype=np.float64)
        # Scatter global best neighbors into local best neighbors
        comm.Scatter([gg, MPI.DOUBLE], [g, MPI.DOUBLE])
        comm.Scatter([fgg, MPI.DOUBLE], [fg, MPI.DOUBLE])

        # Update the particles velocities
        if crank == 0:
            rand = np.random.uniform(size=(2*S+S*S,D))
        else:
            rand = np.empty((2*S+S*S,D), dtype=np.float)
        # Scatter has different ordering than needed for reproducible results
        # Do it manually
        comm.Bcast(rand, root=0)

        rp = rand[crank*iS:crank*iS+iS,:]
        rg = rand[S+crank*iS:S+crank*iS+iS,:]
        if strategy.lower() == 'original':    # Kennedy & Eberhart, 2001
            v = inertia*v + phip*rp*(p - x) + phig*rg*(g - x)
            v = np.clip(v, vmin, vmax)
            v = np.clip(v, vmin, vmax)
        elif strategy.lower() == 'inertia':   # Shi & Eberhart (1998)
            inertia = imax - float(it)/float(maxn-1) * (imax-imin)
            v = inertia*v + phip*rp*(p - x) + phig*rg*(g - x)
            v = np.clip(v, vmin, vmax)
        elif strategy.lower() == 'canonical': # Clerc & Kennedy (2000)
            v = inertia * (v + phip*rp*(p - x) + phig*rg*(g - x))
        elif strategy.lower() == 'fips':      # Mendes & Kennedy (2004)
            acc_coeff = (phip + phig) / float(S)
            for i in range(iS):
                ri = rand[2*S+crank*iS*S+i*S:2*S+crank*iS*S+(i+1)*S,:]
                v[i,:] = inertia * (v[i,:] + np.sum(ri[:,:]*acc_coeff*(gp[:,:]-x[i,:]), axis=0))
        elif strategy.lower() == 'nips':      # Mendes & Kennedy (2004)
            acc_coeff = (phip + phig) / float(S)
            for i in range(iS):
                ri = rand[2*S+crank*iS*S+i*S:2*S+crank*iS*S+(i+1)*S,:]
                ii = get_neighbor_indeces(crank*iS+i, S, topology, kl=kl)
                v[i,:] = inertia * (v[i,:] + np.sum(ri[ii,:]*acc_coeff*(gp[ii,:]-x[i,:]), axis=0))

        # Update the particles positions
        x = x + v
        x = np.where(mask2, x, x02)

        # Limit to bounds
        x = np.clip(x, lb, ub)

        # Update objectives and constraints
        if processes > 1:
            fs = np.array(mp_pool.map(is_feasible, x))
            ii = np.where(fs)[0]
            if ii.size > 0:
                fx[ii] = np.array(mp_pool.map(obj, x[ii,:]))
                if maxit: fx[ii] *= -1.
        else:
            for i in range(iS):
                fs[i] = is_feasible(x[i,:])
                if fs[i]:
                    print(i, x[i:])
                    fx[i] = obj(x[i,:])
                    if maxit: fx[i] *= -1.

        # NaN/Inf
        large = max(fp.max(), fx[np.isfinite(fx)].max())
        large = 1.1*large if large>0. else 0.9*large
        fx = np.where(np.isfinite(fx), fx, large)

        # print(2, x[np.where(fs)[0],:])
        # Store particles best positions (if constraints are satisfied)
        i_update = (fx < fp) & fs
        if np.any(i_update):
            p[i_update,:] = x[i_update,:].copy()
            fp[i_update]  = fx[i_update]

        # gather local best particles into global best particles
        comm.Allgather([p, MPI.DOUBLE], [gp, MPI.DOUBLE] )
        comm.Allgather([fp, MPI.DOUBLE], [fgp, MPI.DOUBLE] )

        if dodiv:
            comm.Allgather([x, MPI.DOUBLE], [gx, MPI.DOUBLE])
            gmask = np.tile(mask1,S).reshape((S,D))
            divers = diversity(gx, lb, ub, gmask)
            if crank == 0: print(2, divers, fgp.min())

        it += 1

        # save restart
        if restartfile1 is not None:
            if crank == 0:
                rs1, rs2, rs3, rs4, rs5 = np.random.get_state()
            else:
                rs1, rs2, rs3, rs4, rs5 = 'MT19937', np.array(624, dtype=np.uint), 0, 0, 0.
            exec("savez_compressed("+saveargarray+")")
            p2 = open(restartfile2, 'w')
            exec("print("+saveargint+", file=p2)")
            exec("print("+saveargfloat+", file=p2)")
            exec("print("+saveargbool+", file=p2)")
            exec("print("+saveargstring+", file=p2)")
            p2.close()
    # end of swarm iteration: while it < maxn:

    if (it == maxn) and (verbose>=1) and (crank == 0): print('Maximum iterations reached --> {:}.'.format(maxn))

    # global best
    i_min = np.argmin(fgp)
    bestx = gp[i_min,:].copy()
    bestf = fgp[i_min]

    if maxit:
        bestf *= -1.
        fgp   *= -1.

    if not any(map(is_feasible, gp)) and (crank == 0): print("PSO could not find any feasible point in the search space.")

    # write parameter file with best parameters
    if crank == 0:
        if isinstance(func, (str,list)):
            parameterwriter(parameterfile, bestx, lb, ub, mask1)

    out = [bestx, bestf]
    if pout:
        out += [p, fp]
    if cout:
        out += [it*S]

    return out


if __name__ == '__main__':
    # import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    comm = MPI.COMM_WORLD
    csize = comm.Get_size()
    crank = comm.Get_rank()

    algo = 'fips'
    init = 'lhs'
    swarmsize = 40
    maxn = 250
    topology = 'neumann'

    if swarmsize % csize != 0:
        raise ValueError("Swarmsize "+str(swarmsize)+" must be multiple of number of processes "+str(csize)+".")

    from jams.functions import ackley, griewank, goldstein_price, rastrigin, rosenbrock, six_hump_camelback
    '''
    This is the Ackley Function
    Global Optimum (n>=2): 0.0 at origin
    '''
    npara = 10
    lb = -10*np.ones(npara)
    ub = 10*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(ackley, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Ackley 0 at origin ', bestx, bestf)
    '''
        This is the Griewank Function (2-D or 10-D)
        Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
           Global Optimum: 0, at origin
    '''
    npara = 10
    lb = -600*np.ones(npara)
    ub = 600*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(griewank, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Griewank 0 at origin ', bestx, bestf)
    '''
    This is the Goldstein-Price Function
    Bound X1=[-2,2], X2=[-2,2]
    Global Optimum: 3.0,(0.0,-1.0)
    '''
    npara = 2
    lb = -2*np.ones(npara)
    ub = 2*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(goldstein_price, x0, lb, ub, processes=4, minobj=3.+1e-8,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Goldstein 3 at (0,-1) ', bestx, bestf)
    '''
    This is the Rastrigin Function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, (0,0)
    '''
    npara = 2
    lb = -1*np.ones(npara)
    ub = 1*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(rastrigin, x0, lb, ub, processes=4, minobj=-2.+1e-8,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Rastrigin -2 at origin ', bestx, bestf)
    '''
    This is the Rosenbrock Function
    Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
           lb=[-5 -5]; ub=[5 5]; x0=[1 1];
    '''
    npara = 2
    lb = -2*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Rosenbrock 0 at (1,1) ', bestx, bestf)
    '''
    This is the Six-hump Camelback Function.
    Bound: X1=[-5,5], X2=[-5,5]
    True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    npara = 2
    lb = -5*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(six_hump_camelback, x0, lb, ub, processes=4, minobj=-1.031628453489877+1e-8,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None)
    if crank == 0: print('Six_hump_camelback -1.03... +-(-0.08983,0.7126) ', bestx, bestf)

    # Restart
    algo = 'fips'
    init = 'lhs'
    swarmsize = 40
    maxn = 250
    topology = 'neumann'
    if swarmsize % csize != 0:
        raise ValueError("Swarmsize "+str(swarmsize)+" must be multiple of number of processes "+str(csize)+".")
    from jams.functions import rosenbrock
    npara = 2
    lb = -2*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    seed = 1234
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None,
                       seed=seed, restart=False)
    if crank == 0: print('Rosenbrock Reference - ', bestx, bestf)
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn//2,
                       seed=seed, restart=False)
    if crank == 0: print('Rosenbrock Restart 1 - ', bestx, bestf)
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4,
                       init=init, strategy=algo, topology=topology, verbose=0,
                       swarmsize=swarmsize, maxn=maxn,
                       seed=seed, restart=True)
    if crank == 0: print('Rosenbrock Restart 2 - ', bestx, bestf)

    # Constraints
    algo = 'fips'
    init = 'lhs'
    swarmsize = 40
    maxn = 250
    topology = 'neumann'
    if swarmsize % csize != 0:
        raise ValueError("Swarmsize "+str(swarmsize)+" must be multiple of number of processes "+str(csize)+".")
    from jams.functions import rosenbrock
    npara = 2
    lb = -2*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    seed = 1234
    def fcon(x): # feasible = np.all(fcon(x)>=0.)
        return -np.sign(x)
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4, f_ieqcons=fcon,
                       init=init, strategy=algo, topology=topology, verbose=2,
                       swarmsize=swarmsize, maxn=maxn, restartfile1=None,
                       seed=seed, restart=False)
    if crank == 0: print('Rosenbrock constraint - ', bestx, bestf)
