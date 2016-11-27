#!/usr/bin/env python
from functools import partial
import numpy as np
# ToDo:
#   maximum velocity
#   LHS / Sobol
#   local best / topologies
#   crossover with quadratic function
#   external function
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
        swarmsize=40, omega=0.5, phip=2., phig=2., maxiter=250,
        minstep=1e-8, minfunc=1e-8, debug=False, processes=1,
        particle_output=False):
    """
        Perform a particle swarm optimization (PSO)
       
        Parameters
        ==========
        func : function
            The function to be minimized
        lb : array
            The lower bounds of the design variable(s)
        ub : array
            The upper bounds of the design variable(s)
       
        Optional
        ========
        ieqcons : list
            A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in 
            a successfully optimized problem (Default: [])
        f_ieqcons : function
            Returns a 1-D array in which each element must be greater or equal 
            to 0.0 in a successfully optimized problem. If f_ieqcons is specified, 
            ieqcons is ignored (Default: None)
        args : tuple
            Additional arguments passed to objective and constraint functions
            (Default: empty tuple)
        kwargs : dict
            Additional keyword arguments passed to objective and constraint 
            functions (Default: empty dict)
        swarmsize : int
            The number of particles in the swarm (Default: 40)
        omega : scalar
            Particle velocity scaling factor (Default: 0.5)
        phip : scalar
            Scaling factor to search away from the particle's best known position
            (Default: 2.0)
        phig : scalar
            Scaling factor to search away from the swarm's best known position
            (Default: 2.0)
        maxiter : int
            The maximum number of iterations for the swarm to search (Default: 250)
        minstep : scalar
            The minimum stepsize of swarm's best position before the search
            terminates (Default: 1e-8)
        minfunc : scalar
            The minimum change of swarm's best objective value before the search
            terminates (Default: 1e-8)
        debug : boolean
            If True, progress statements will be displayed every iteration
            (Default: False)
        processes : int
            The number of processes to use to evaluate objective function and 
            constraints (default: 1)
        particle_output : boolean
            Whether to include the best per-particle position and the objective
            values at those.
       
        Returns
        =======
        g : array
            The swarm's best known position (optimal design)
        f : scalar
            The objective value at ``g``
        p : array
            The best known position per particle
        pf: arrray
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
    """
   
    assert len(lb)==len(ub), 'Lower- and upper-bounds must be the same length'
    assert hasattr(func, '__call__'), 'Invalid function handle'
    lb = np.array(lb)
    ub = np.array(ub)
    assert np.all(ub>lb), 'All upper-bound values must be greater than lower-bound values'
   
    vhigh = np.abs(ub - lb)
    vlow = -vhigh

    # Initialize objective function
    obj = partial(_obj_wrapper, func, args, kwargs)
    
    # Check for constraint function(s) #########################################
    if f_ieqcons is None:
        if not len(ieqcons):
            if debug:
                print('No constraints given.')
            cons = _cons_none_wrapper
        else:
            if debug:
                print('Converting ieqcons to a single constraint function')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        if debug:
            print('Single constraint function given in f_ieqcons')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the multiprocessing module if necessary
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)
        
    # Initialize the particle swarm ############################################
    S = swarmsize
    D = len(lb)  # the number of dimensions each particle has
    x = np.random.rand(S, D)  # particle positions
    v = np.zeros_like(x)  # particle velocities
    p = np.zeros_like(x)  # best particle positions
    fx = np.zeros(S)  # current particle function values
    fs = np.zeros(S, dtype=bool)  # feasibility of each particle
    fp = np.ones(S)*np.inf  # best particle function values
    g = []  # best swarm position
    fg = np.inf  # best swarm position starting value
    
    # Initialize the particle's position
    x = lb + x*(ub - lb)

    # Calculate objective and constraints for each particle
    if processes > 1:
        fx = np.array(mp_pool.map(obj, x))
        fs = np.array(mp_pool.map(is_feasible, x))
    else:
        for i in range(S):
            fx[i] = obj(x[i, :])
            fs[i] = is_feasible(x[i, :])
       
    # Store particle's best position (if constraints are satisfied)
    i_update = np.logical_and((fx < fp), fs)
    p[i_update, :] = x[i_update, :].copy()
    fp[i_update] = fx[i_update]

    # Update swarm's best position
    i_min = np.argmin(fp)
    if fp[i_min] < fg:
        fg = fp[i_min]
        g = p[i_min, :].copy()
    else:
        # At the start, there may not be any feasible starting point, so just
        # give it a temporary "best" point since it's likely to change
        g = x[0, :].copy()
       
    # Initialize the particle's velocity
    v = vlow + np.random.rand(S, D)*(vhigh - vlow)
       
    # Iterate until termination criterion met ##################################
    it = 1
    while it <= maxiter:
        rp = np.random.uniform(size=(S, D))
        rg = np.random.uniform(size=(S, D))
        '''
                        /*-------PSO variant with constriction coefficients------------------------------------*/
			/*  ''Clerc's analysis of the iterative system led him to propose a strategy for the
			 *  placement of "constriction coefficients" on the terms of the formulas; these
			 *  coefficients controlled the convergence of the particle and allowed an elegant and
			 *  well-explained method for preventing explosion, ensuring convergence, and
			 *  eliminating the arbitrary Vmax parameter. The analysis also takes the guesswork
			 *  out of setting the values of phi_1 and phi_2.''
			 *  ''this is the canonical particle swarm algorithm of today.''
			 *  [Poli et al., 2007] http://dx.doi.org/10.1007/s11721-007-0002-0
			 *  [Clerc and Kennedy, 2002] http://dx.doi.org/10.1109/4235.985692
			 *  
			 *  This being the canonical PSO of today, this variant is set as the default in PaGMO.
			 *-------------------------------------------------------------------------------------*/
			else if( m_variant == 5 ){
				for( d = 0; d < Dc; d++ ){
					r1 = m_drng();
					r2 = m_drng();
					V[p][d] = m_omega * ( V[p][d] + m_eta1 * r1 * (lbX[p][d] - X[p][d]) + m_eta2 * r2 * (best_neighb[d] - X[p][d]) );
				}
        '''
                                
        # Update the particles velocities
        v = omega*v + phip*rp*(p - x) + phig*rg*(g - x)
        # Update the particles' positions
        x = x + v
        # Correct for bound violations
        maskl = x < lb
        masku = x > ub
        x = x*(~np.logical_or(maskl, masku)) + lb*maskl + ub*masku

        # Update objectives and constraints
        if processes > 1:
            fx = np.array(mp_pool.map(obj, x))
            fs = np.array(mp_pool.map(is_feasible, x))
        else:
            for i in range(S):
                fx[i] = obj(x[i, :])
                fs[i] = is_feasible(x[i, :])

        # Store particle's best position (if constraints are satisfied)
        i_update = np.logical_and((fx < fp), fs)
        p[i_update, :] = x[i_update, :].copy()
        fp[i_update] = fx[i_update]

        # Compare swarm's best position with global best position
        i_min = np.argmin(fp)
        if fp[i_min] < fg:
            if debug:
                print('New best for swarm at iteration {:}: {:} {:}'\
                    .format(it, p[i_min, :], fp[i_min]))

            p_min = p[i_min, :].copy()
            stepsize = np.sqrt(np.sum((g - p_min)**2))

            if np.abs(fg - fp[i_min]) <= minfunc:
                print('Stopping search: Swarm best objective change less than {:}'\
                    .format(minfunc))
                if particle_output:
                    return p_min, fp[i_min], p, fp
                else:
                    return p_min, fp[i_min]
            elif stepsize <= minstep:
                print('Stopping search: Swarm best position change less than {:}'\
                    .format(minstep))
                if particle_output:
                    return p_min, fp[i_min], p, fp
                else:
                    return p_min, fp[i_min]
            else:
                g = p_min.copy()
                fg = fp[i_min]

        if debug:
            print('Best after iteration {:}: {:} {:}'.format(it, g, fg))
        it += 1

    print('Stopping search: maximum iterations reached --> {:}'.format(maxiter))
    
    if not is_feasible(g):
        print("However, the optimization couldn't find a feasible design. Sorry")
    if particle_output:
        return g, fg, p, fp
    else:
        return g, fg


if __name__ == '__main__':
    # import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    from functions import goldstein_price, griewank, rastrigin, rosenbrock, six_hump_camelback
    '''
        This is the Griewank Function (2-D or 10-D)
        Bound: X(i)=[-600,600], for i=1,2,...,10  !for visualization only 2!
           Global Optimum: 0, at origin
    '''
    npara = 10 # dimension of griewank function
    bl = -600*np.ones(npara)
    bu = 600*np.ones(npara)
    bestx, bestf = pso(griewank, bl, bu, processes=4)
    print(bestx, bestf)
