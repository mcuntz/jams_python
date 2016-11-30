#!/usr/bin/env python
from __future__ import print_function
from functools import partial
import subprocess
import numpy as np
from jams.const import huge
# ToDo:
#   include x0
#   seed
#   restart
#   crossover with quadratic function (QIPSO)
#   iPython parallel
#   MPI

def _ext_obj_wrapper(func, lb, ub, mask,
                     parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug,
                     params):
    '''
        Wrapper function for external program to be optimised
        to be used with partial:
            obj = partial(_ext_obj_wrapper, func, lb, ub, mask,
                          parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
        This allows then calling obj with only the argument params:
            fx = obj(params)


        Definition
        ----------
        def _ext_obj_wrapper(func, lb, ub, mask,
                             parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug,
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
        if debug:
            err = subprocess.call(func, shell=shell)
        else:
            err = subprocess.check_output(func, shell=shell)
        obj = objectivereader(objectivefile)
        return obj
    else:
        return func(params)


# Function wrappers for objective and constraints
# used with functools.partial
def _obj_wrapper(func, args, kwargs, x):
    '''
        Wrapper function for function to be optimised
        to be used with partial:
            obj = partial(_obj_wrapper, func, args, kwargs)
        This allows then calling obj with only the argument x:
            fx = obj(x)
        which translates to:
            fx = func(x, *args, **kwargs)
    '''
    return func(x, *args, **kwargs)


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
    return np.all(func(x)>=0)


def _cons_none_wrapper(x):
    '''
        Dummy constraints function if no constraints.
        Returns 0.
    '''
    return np.array([0])


def _cons_ieqcons_wrapper(ieqcons, args, kwargs, x):
    '''
        Wrapper function for constraints to be used with partial in the case
        that ieqcons is given, i.e. a list of constraint functions, returning
        >=0. if constraints are met.
    '''
    return np.array([y(x, *args, **kwargs) for y in ieqcons])


def _cons_f_ieqcons_wrapper(f_ieqcons, args, kwargs, x):
    '''
        Wrapper function for constraints to be used with partial in the case
        that f_ieqcons is given, i.e. a single function returning an 1-D array
        values >=0. where constraints are met.
    '''
    return np.array(f_ieqcons(x, *args, **kwargs))


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

    if topology.lower() == 'gbest':
        i_min = np.argmin(fp)
        g  = p[i_min,:].copy()         # overall best
        fg = fp[i_min]
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


# Particle Swarm Optimisation
def pso(func, x0, lb, ub,
        mask=None,
        ieqcons=[], f_ieqcons=None,
        args=(), kwargs={},
        swarmsize=None, inertia=None, phip=None, phig=None, maxn=250,
        minstep=1e-8, minobj=1e-8, maxit=False,
        init='lhs', strategy='fips', topology='gbest', kl=1,
        processes=1,
        verbose=0, pout=False, cout=False,
        parameterfile=None, parameterwriter=None,
        objectivefile=None, objectivereader=None,
        shell=False, debug=False):
    """
        Particle Swarm Optimization (PSO)


        Definition
        ----------
        def pso(func, x0, lb, ub,
                mask=None,
                ieqcons=[], f_ieqcons=None,
                args=(), kwargs={},
                swarmsize=None, inertia=None, phip=None, phig=None, maxn=250,
                minstep=1e-8, minobj=1e-8, maxit=False,
                init='lhs', strategy='fips', topology='gbest', kl=1,
                processes=1,
                verbose=0, pout=False, cout=False,
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
                    A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in
                    a successfully optimized problem.
                    (Default: None)
        f_ieqcons   function
                    Returns a 1-D array in which each element must be greater or equal
                    to 0.0 in a successfully optimized problem. If f_ieqcons is specified,
                    ieqcons is ignored.
                    (Default: None)
        args        tuple
                    Additional arguments passed to objective and constraint functions.
                    (Default: empty tuple)
        kwargs      dict
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
        processes   int
                    The number of processes to use to evaluate objective function and constraints.
                    (Default: 1)
        pout        boolean
                    True: include best per-particle positions and their objective values in output.
                    (Default: False)
        cout        boolean
                    True: include number of function calls in output.
                    (Default: False)
        parameterfile     Parameter file for executable; must be given if func is name of executable
        parameterwriter   Python function for writing parameter file if func is name of executable
        objectivefile     File with objective value from executable; must be given if func is name of executable
        objectivereader   Python function for reading objective value if func is name of executable
        shell             If True, the specified executable will be executed through the shell (default: False).
        debug             If True, model output is displayed for executable (default: False).


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
    """
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
        if parameterwriter is None:
            raise IOError('parameterwrite must be given if func is name of executable.')
        if objectivefile is None:
            raise IOError('objectivefile must be given if func is name of executable.')
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

    # Problem sizes
    D = len(lb) # dimension of each particle
    if swarmsize is None:
        S = max(min(3*D,40),10)
    else:
        S = swarmsize

    # Initialise 1D mask
    if mask is not None:
        mask1 = mask
    else:
        mask1 = np.ones(D, dtype=np.bool)
    # 2D mask
    mask2 = np.tile(mask1,S).reshape((S,D))
    x02   = np.tile(x0,S).reshape((S,D))

    # Partialise objective function
    if isinstance(func, (str,list)):
        obj = partial(_ext_obj_wrapper, func, lb, ub, mask1,
                      parameterfile, parameterwriter, objectivefile, objectivereader, shell, debug)
    else:
        obj = partial(_obj_wrapper, func, args, kwargs)

    # Check for constraint function(s) and partialise them
    if f_ieqcons is None:
        if ieqcons is None:
            if verbose>=1:
                print('No constraints given.')
            cons = _cons_none_wrapper
        else:
            if verbose>=1:
                print('Converting ieqcons to a single constraint function.')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        if verbose>=1:
            print('Single constraint function given in f_ieqcons.')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the multiprocessing module if necessary
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)

    # Deal with NaN and Inf
    large = 0.5*huge

    # Initialize the particle swarm
    # current particle positions
    # current particle velocities
    if init.lower() == 'random':
        x  = np.random.uniform(size=(S,D))
        v  = np.random.uniform(size=(S,D))
    elif init.lower() == 'sobol':
        import sobol
        nskip = D*S
        x = sobol.i4_sobol_generate(D,S,nskip).transpose()
        nskip = (D*S)**2
        v = sobol.i4_sobol_generate(D,S,nskip).transpose()
    elif init.lower() == 'lhs':
        import scipy.stats as stats
        from jams.lhs import lhs
        dist = [stats.uniform for i in range(D)]
        pars = [(0,1) for i in range(D)]
        x    = lhs(dist, pars, S).transpose()
        v    = lhs(dist, pars, S).transpose()
    fx = np.ones(S)*large              # current particle function values
    fs = np.zeros(S, dtype=bool)       # current combined feasibility for each particle
    p  = np.ones((S,D))*large          # particles individual best positions
    fp = np.ones(S)*large              # particles individual best function values

    # Maximum velocity
    vmax = np.abs(ub - lb)
    vmin = -vmax

    # Initialize particle positions and velocities
    v = vmin + v*(vmax-vmin)
    x = lb + x*(ub - lb)
    x = np.where(mask2, x, x02)

    # Calculate first objective and constraints for each particle
    if processes > 1:
        fx = np.array(mp_pool.map(obj, x))
        fs = np.array(mp_pool.map(is_feasible, x))
    else:
        for i in range(S):
            fx[i] = obj(x[i,:])
            fs[i] = is_feasible(x[i,:])
    # maximise
    if maxit: fx *= -1.

    # NaN/Inf
    large = max(fp.max(), fx[np.isfinite(fx)].max())
    large = 1.1*large if large>0. else 0.9*large
    fx = np.where(np.isfinite(fx), fx, large)

    # Store particles best positions (if constraints are satisfied)
    i_update = (fx < fp) & fs
    if np.any(i_update):
        p[i_update,:] = x[i_update,:].copy()
        fp[i_update]  = fx[i_update]
    
    # Iterate until termination criterion met
    it = 1
    while (it < maxn):
        # Stop if minimum found
        if fp.min() < minobj:
            if verbose>=1: print('minobj found.')
            break

        # Update neighbors best positions
        g, fg = get_best_neighbor(p, fp, topology, kl=kl)

        # Update the particles velocities
        rp = np.random.uniform(size=(S,D))
        rg = np.random.uniform(size=(S,D))
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
            for i in range(S):
                ri = np.random.uniform(size=(S,D))
                v[i,:] = inertia * (v[i,:] + np.sum(ri[:,:]*acc_coeff*(p[:,:]-x[i,:]), axis=0))
        elif strategy.lower() == 'nips':      # Mendes & Kennedy (2004)
            acc_coeff = (phip + phig) / float(S)
            for i in range(S):
                ri = np.random.uniform(size=(S,D))
                ii = get_neighbor_indeces(i, S, topology, kl=kl)
                v[i,:] = inertia * (v[i,:] + np.sum(ri[ii,:]*acc_coeff*(p[ii,:]-x[i,:]), axis=0))

        # Update the particles positions
        x = x + v
        x = np.where(mask2, x, x02)

        # Limit to bounds
        x = np.clip(x, lb, ub)

        # Update objectives and constraints
        if processes > 1:
            fx = np.array(mp_pool.map(obj, x))
            fs = np.array(mp_pool.map(is_feasible, x))
        else:
            for i in range(S):
                fx[i] = obj(x[i,:])
                fs[i] = is_feasible(x[i,:])
        if maxit: fx *= -1.

        # NaN/Inf
        large = max(fp.max(), fx[np.isfinite(fx)].max())
        large = 1.1*large if large>0. else 0.9*large
        fx = np.where(np.isfinite(fx), fx, large)

        # Store particles best positions (if constraints are satisfied)
        i_update = (fx < fp) & fs
        if np.any(i_update):
            p[i_update,:] = x[i_update,:].copy()
            fp[i_update]  = fx[i_update]

        it += 1

    if (it == maxn) and (verbose>=1): print('Maximum iterations reached --> {:}.'.format(maxn))

    # global best
    i_min = np.argmin(fp)
    bestx = p[i_min,:].copy()
    bestf = fp[i_min]

    if maxit:
        bestf *= -1.
        fp    *= -1.

    if not any(map(is_feasible, p)): print("PSO could not find any feasible point in the search space.")

    # write parameter file with best parameters
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

    algo = 'nips'
    init = 'lhs'
    swarmsize = None
    maxn = 250
    topology = 'gbest'

    from jams.functions import ackley, griewank, goldstein_price, rastrigin, rosenbrock, six_hump_camelback
    '''
    This is the Ackley Function
    Global Optimum (n>=2): 0.0 at origin
    '''
    npara = 10
    lb = -10*np.ones(npara)
    ub = 10*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(ackley, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Ackley ', bestx, bestf)
    '''
        This is the Griewank Function (2-D or 10-D)
        Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
           Global Optimum: 0, at origin
    '''
    npara = 10
    lb = -600*np.ones(npara)
    ub = 600*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(griewank, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Griewank ', bestx, bestf)
    '''
    This is the Goldstein-Price Function
    Bound X1=[-2,2], X2=[-2,2]
    Global Optimum: 3.0,(0.0,-1.0)
    '''
    npara = 2
    lb = -2*np.ones(npara)
    ub = 2*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(goldstein_price, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Goldstein ', bestx, bestf)
    '''
    This is the Rastrigin Function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, (0,0)
    '''
    npara = 2
    lb = -1*np.ones(npara)
    ub = 1*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(rastrigin, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Rastrigin ', bestx, bestf)
    '''
    This is the Rosenbrock Function
    Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
           lb=[-5 -5]; ub=[5 5]; x0=[1 1];
    '''
    npara = 2
    lb = -2*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(rosenbrock, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Rosenbrock ', bestx, bestf)
    '''
    This is the Six-hump Camelback Function.
    Bound: X1=[-5,5], X2=[-5,5]
    True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    npara = 2
    lb = -5*np.ones(npara)
    ub = 5*np.ones(npara)
    x0 = np.zeros(npara)
    bestx, bestf = pso(six_hump_camelback, x0, lb, ub, processes=4, init=init, strategy=algo, topology=topology, verbose=0, swarmsize=swarmsize, maxn=maxn)
    print('Six_hump_camelback ', bestx, bestf)
