#!/usr/bin/env python
from __future__ import print_function
from functools import partial
import numpy as np
# ToDo:
#   local best / topologies
#   crossover with quadratic function
#   external function
#   iPython parallel
#   MPI

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


# Particle Swarm optimisation
def pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={},
        swarmsize=40, omega=None, phip=None, phig=None, maxiter=250,
        minstep=1e-8, minobj=1e-8,
        init='random', psotype='fips',
        verbose=0, processes=1, particle_output=False):
    """
        Particle Swarm Optimization (PSO)


        Definition
        ----------
        def pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={},
                swarmsize=40, omega=None, phip=None, phig=None, maxiter=250,
                minstep=1e-8, minobj=1e-8,
                init='random', psotype='fips',
                verbose=0, processes=1, particle_output=False):


        Input
        -----
        func        function
                    The function to be minimized
        lb          1D-array
                    The lower bounds of the design variable(s)
        ub          1D-array
                    The upper bounds of the design variable(s)


        Optional Input
        --------------
        ieqcons     list
                    A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in
                    a successfully optimized problem (Default: [])
        f_ieqcons   function
                    Returns a 1-D array in which each element must be greater or equal
                    to 0.0 in a successfully optimized problem. If f_ieqcons is specified,
                    ieqcons is ignored (Default: None)
        args        tuple
                    Additional arguments passed to objective and constraint functions
                    (Default: empty tuple)
        kwargs      dict
                    Additional keyword arguments passed to objective and constraint
                    functions (Default: empty dict)
        swarmsize   int
                    The number of particles in the swarm (Default: 40)
        omega       scalar
                    Particle velocity scaling factor.
                    Default depends on algorithm (psotype).
        phip        scalar
                    Scaling factor to search away from the particle's best known position.
                    Default depends on algorithm (psotype).
        phig        scalar
                    Scaling factor to search away from the swarm's best known position.
                    Default depends on algorithm (psotype).
        maxiter     int
                    The maximum number of iterations for the swarm to search (Default: 250)
        minstep     scalar
                    The minimum stepsize of swarm's best position before the search
                    terminates (Default: 1e-8)
        minobj     scalar
                    Objective function defining convergence (Default: 1e-8)
        init        string
                    How to sample the initial swarm positions and velocities (Default: 'lhs')
                    'random': random sampling from uniform distribution
                    'lhs':    latin hypercube sampling from uniform distributions
                    'sobol':  quasirandom Sobol sequence (only up to 40 dimensions)
        psotype     string
                    PSO algorithm (Default: 'fips')
                    'original':   Textbook particle swarm algorithm with inertia weight
                                  x = current position
                                  p = particles best position
                                  g = neighborhood best position
                                  rg, rp = np.random.uniform(size=(S,D))
                                  omega=0.5, phip=2., phig=2.
                                  v = omega*v + rp*phip*(p-x) + rg*phig*(g-x)
                                  x = x + v
                    'inertia':    Same as 'original' but with inertia weight decreasing from 0.9 to 0.4
                                  over time, i.e. over iterations (it)
                                  omax = 0.9
                                  omin = 0.4
                                  omega = omax - float(it)/float(maxiter-1) * (omax-omin)
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
                                  omega=0.7289, phip=2.05, phig=2.05
                                  v = omega * (v + phip*rp*(p - x) + phig*rg*(g - x))
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
                                  omega = 0.7289
                                  acc_coeff = phip+phig = 4.1
                                  v = omega * (v + np.sum(ri*acc_coeff/nneighbor*(p[neighbor[:]]-x)))
        verbose     integer
                    Controlling amount of print-out (default: 0)
                    0: No print-out
                    1: Printing convergence criteria, etc.
                    2: Printing after each step.
        processes   int
                    The number of processes to use to evaluate objective function and
                    constraints (default: 1)
        particle_output   boolean
                    Whether to include the best per-particle position and the objective
                    values at those.

        Output
        ------
        g           1D-array
                    The swarm's best known position (optimal design)
        f           scalar
                    The objective value at ``g``
        p           ND-array
                    The best known position per particle
        pf          1D-array
                    The objective values at each position in p


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
                  MC, Nov 2016 - swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxiter=100
                                 -> swarmsize=40, omega=0.5, phip=2., phig=2., maxiter=250
                               - include vmax: v.clip(vmin,vmax)
                               - Sobol sequences and latin hypercube sampling for initial swarm positions
                               - Different PSO algorithms
                               - minfunc -> minobj
                               - debug -> verbose
    """
    # Check arguments
    assert len(lb)==len(ub), 'Lower- and upper-bounds must be the same length'
    assert hasattr(func, '__call__'), 'Invalid function handle'
    lb = np.array(lb) # assert ND-array
    ub = np.array(ub)
    assert np.all(ub > lb), 'All upper-bound values must be greater than lower-bound values'
    inits = ['random', 'lhs', 'sobol']
    assert init.lower() in inits, 'Initialisation not in {:}'.format(inits)
    if init.lower() == 'sobol':
        assert len(lb) <= 40, "Sobol' sequences only work up to 40 dimensions."
    ptypes = ['original', 'inertia', 'canonical', 'fips']
    assert psotype.lower() in ptypes, 'PSO type not in {:}'.format(ptypes)

    # Defaults
    if psotype.lower() == 'original':   # Kennedy & Eberhart, 2001
        if omega is None: omega=0.5
        if phip  is None: phip=2.
        if phig  is None: phig=2.
    elif psotype.lower() == 'inertia':  # Shi & Eberhart (1998)
        omax = 0.9
        omin = 0.4
        if phip is None: phip=2.
        if phig is None: phig=2.
    elif psotype.lower() == 'canonical': # Clerc & Kennedy (2000)
        if omega is None: omega=0.7289
        if phip  is None: phip=2.05
        if phig  is None: phig=2.05
    elif psotype.lower() == 'fips': # Mendes & Kennedy (2004)
        if omega is None: omega=0.7289
        if phip  is None: phip=2.05
        if phig  is None: phig=2.05

    # Maximum velocity
    vmax = np.abs(ub - lb)
    vmin = -vmax

    # Initialize objective function
    obj = partial(_obj_wrapper, func, args, kwargs)

    # Check for constraint function(s)
    if f_ieqcons is None:
        if not len(ieqcons):
            if verbose>=1:
                print('No constraints given.')
            cons = _cons_none_wrapper
        else:
            if verbose>=1:
                print('Converting ieqcons to a single constraint function')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        if verbose>=1:
            print('Single constraint function given in f_ieqcons')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the multiprocessing module if necessary
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)

    # Initialize the particle swarm
    S  = swarmsize
    D  = len(lb)                       # dimension of each particle
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
    fx = np.ones(S)*np.inf             # current particle function values
    fs = np.zeros(S, dtype=bool)       # current combined feasibility for each particle
    p  = np.ones((S,D))*np.inf         # particle's individual best positions
    fp = np.ones(S)*np.inf             # particle's individual best function values
    g  = np.ones(D)*np.inf             # swarm's best position
    fg = np.inf                        # swarm's best function value

    # Initialize the particles positions and velocities
    x = lb + x*(ub - lb)
    v = vmin + v*(vmax-vmin)

    # Calculate objective and constraints (may be dummy) for each particle
    if processes > 1:
        fx = np.array(mp_pool.map(obj, x))
        fs = np.array(mp_pool.map(is_feasible, x))
    else:
        for i in range(S):
            fx[i] = obj(x[i,:])
            fs[i] = is_feasible(x[i,:])

    # Store particle's best position (if constraints are satisfied)
    i_update = (fx < fp) & fs
    if np.any(i_update):
        p[i_update,:] = x[i_update,:].copy()
        fp[i_update]  = fx[i_update]

    # Update swarm's best position
    i_min = np.argmin(fp)
    if fp[i_min] < fg:
        g  = p[i_min,:].copy()
        fg = fp[i_min]
    else:
        # At the start, there may not be any feasible starting point, so just
        # give it a temporary "best" point since it's likely to change
        g  = x[0, :].copy()
        fg = fp[0]

    # Iterate until termination criterion met
    it = 1
    while (it <= maxiter):
        # Update the particles velocities
        rp = np.random.uniform(size=(S,D))
        rg = np.random.uniform(size=(S,D))
        if psotype.lower() == 'original':    # Kennedy & Eberhart, 2001
            v = omega*v + phip*rp*(p - x) + phig*rg*(g - x)
            v.clip(vmin, vmax)
        elif psotype.lower() == 'inertia':   # Shi & Eberhart (1998)
            omega = omax - float(it)/float(maxiter-1) * (omax-omin)
            v = omega*v + phip*rp*(p - x) + phig*rg*(g - x)
            v.clip(vmin, vmax)
        elif psotype.lower() == 'canonical': # Clerc & Kennedy (2000)
            v = omega * (v + phip*rp*(p - x) + phig*rg*(g - x))
        elif psotype.lower() == 'fips':      # Mendes & Kennedy (2004)
            acc_coeff = (phip + phig) / float(S)
            for i in range(S):
                ri = np.random.uniform(size=(S,D))
                v[i,:] = omega * (v[i,:] + np.sum(ri[:,:]*acc_coeff*(p[:,:]-x[i,:]), axis=0))

        # Update the particles positions
        x = x + v

        # Limit to bounds
        x.clip(lb, ub)

        # Update objectives and constraints
        if processes > 1:
            fx = np.array(mp_pool.map(obj, x))
            fs = np.array(mp_pool.map(is_feasible, x))
        else:
            for i in range(S):
                fx[i] = obj(x[i,:])
                fs[i] = is_feasible(x[i,:])

        # Store particle's best position (if constraints are satisfied)
        i_update = (fx < fp) & fs
        if np.any(i_update):
            p[i_update,:] = x[i_update,:].copy()
            fp[i_update]  = fx[i_update]

        # Compare swarm's best position with global best position
        i_min = np.argmin(fp)
        
        # Stop if minimum found
        if fp[i_min] < minobj:
            g  = p[i_min,:].copy()
            fg = fp[i_min]
            break
        
        # Set new swarm best
        if fp[i_min] < fg:
            if verbose>=1: print('New best for swarm at iteration {:}: {:} {:}'.format(it, p[i_min, :], fp[i_min]))
            g  = p[i_min,:].copy()
            fg = fp[i_min]
            # stepsize = np.sqrt(np.sum((g - p_min)**2))
            # if stepsize <= minstep: break
            # if np.abs(fg - fp_min) <= minobj: break

        if verbose==2: print('Best after iteration {:}: {:} {:}'.format(it, g, fg))
        it += 1

    if (it == (maxiter+1)) and (verbose>=1): print('Maximum iterations reached --> {:}.'.format(maxiter))

    if not is_feasible(g): print("PSO could not find any feasible point in the search space.")

    if particle_output:
        return [g, fg, p, fp]
    else:
        return [g, fg]


if __name__ == '__main__':
    # import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    from jams.functions import ackley, griewank, goldstein_price, rastrigin, rosenbrock, six_hump_camelback
    '''
    This is the Ackley Function
    Global Optimum (n>=2): 0.0 at origin
    '''
    npara = 10
    bl = -10*np.ones(npara)
    bu = 10*np.ones(npara)
    bestx, bestf = pso(ackley, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Ackley ', bestx, bestf)
    '''
        This is the Griewank Function (2-D or 10-D)
        Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
           Global Optimum: 0, at origin
    '''
    npara = 10
    bl = -600*np.ones(npara)
    bu = 600*np.ones(npara)
    bestx, bestf = pso(griewank, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Griewank ', bestx, bestf)
    '''
    This is the Goldstein-Price Function
    Bound X1=[-2,2], X2=[-2,2]
    Global Optimum: 3.0,(0.0,-1.0)
    '''
    npara = 2
    bl = -2*np.ones(npara)
    bu = 2*np.ones(npara)
    bestx, bestf = pso(goldstein_price, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Goldstein ', bestx, bestf)
    '''
    This is the Rastrigin Function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, (0,0)
    '''
    npara = 2
    bl = -1*np.ones(npara)
    bu = 1*np.ones(npara)
    bestx, bestf = pso(rastrigin, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Rastrigin ', bestx, bestf)
    '''
    This is the Rosenbrock Function
    Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
           bl=[-5 -5]; bu=[5 5]; x0=[1 1];
    '''
    npara = 2
    bl = -2*np.ones(npara)
    bu = 5*np.ones(npara)
    bestx, bestf = pso(rosenbrock, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Rosenbrock ', bestx, bestf)
    '''
    This is the Six-hump Camelback Function.
    Bound: X1=[-5,5], X2=[-5,5]
    True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    npara = 2
    bl = -5*np.ones(npara)
    bu = 5*np.ones(npara)
    bestx, bestf = pso(six_hump_camelback, bl, bu, processes=4, init='lhs', psotype='fips', verbose=0, maxiter=250)
    print('Six_hump_camelback ', bestx, bestf)

