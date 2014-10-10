#!/usr/bin/env python
"""
    Defines common functions that are used in curve_fit or fmin parameter estimations.

    There are also some common test functions for parameter estimations such as Rosenbrock and Griewank.

    Defines the functions (execpt test functions) in two forms (ex. of 3 params):
        1. func(x, p1, p2, p3)
        2. func_p(x, p)  with p[0:3]
    These cost functions can be used for example with curve_fit
        p, cov  = opt.curve_fit(ufz.functions.f1x, x, y, p0=[p0,p1])

    Defines also two cost functions, one with absolute sum, one with squared sum of deviations.
        3. cost_func    sum(abs(obs-func))
        4. cost2_func   sum((obs-func)**2)
    These cost functions can be used for example with fmin
        p = opt.fmin(ufz.functions.cost_f1x, np.array([p1,p2]), args=(x,y), disp=False)
    or
        p, nfeval, rc = opt.fmin_tnc(ufz.functions.cost_f1x, [p1,p2], bounds=[[None,None],[None,None]],
                                     args=(x,y), approx_grad=True, disp=False)

    Note the different argument orders:
    curvefit wants f(x,*args) with the independent variable as the first argument
             and the parameters to fit as separate remaining arguments.
    fmin is a general minimiser with respect to the first argument, i.e. func(p,*args).

    There are also two common cost functions (absolute and squared deviations) where any function
    in the form func(x, p) can be used as second argument:
        5. cost_abs(p, func, x, y)
        6. cost_square(p, func, x, y)
    Used for example as
        p = opt.fmin(ufz.functions.cost_abs, np.array([p1,p2]), args=(ufz.functions.f1x_p,x,y), disp=False)
    or
        p, nfeval, rc = opt.fmin_tnc(ufz.functions.cost_square, [p1,p2], bounds=[[None,None],[None,None]],
                                     args=(ufz.functions.f1x_p,x,y), approx_grad=True, disp=False)


    Definition
    ----------
    Current functions are (there is always the second form with the name appended by _p;
    these are used in the cost functions).

        arrhenius         1 param:  Arrhenius temperature dependence of biochemical rates:
                                    exp((T-TC25)*E/(T25*R*(T+T0))), parameter: E
        f1x               2 params: General 1/x function: a + b/x
        fexp              3 params: General exponential function: a + b * exp(c*x)
        gauss             2 params: Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2)), parameter: mu, sig
        lasslop           6 params: Lasslop et al. (2010) a rectangular, hyperbolic light-response GPP
                                    with Lloyd & Taylor (1994) respiration and the maximum canopy uptake
                                    rate at light saturation decreases exponentially with VPD as in Koerner (1995)
        line0             1 params: Straight line: a*x
        line              2 params: Straight line: a + b*x
        lloyd_fix         2 params: Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        lloyd_only_rref   1 param:  Lloyd & Taylor (1994) Arrhenius type with fixed exponential term
        poly              n params: General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        sabx              2 params: sqrt(f1x), i.e. general sqrt(1/x) function: sqrt(a + b/x)

    Current test functions are
    ackley                >=2 params:     Ackley function, global optimum: 0.0 at origin
    goldstein_price       2 params:       Goldstein-Price function, global optimum: 3.0 (0.0,-1.0)
    griewank              2 or 10 params: Griewank function, global optimum: 0 at origin
    rastrigin             2 params:       Rastrigin function, global optimum: -2 (0,0)
    rosenbrock            2 params:       Rosenbrock function, global optimum: 0 (1,1)
    six_hump_camelback    2 params:       Six-hump Camelback function
                                          True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)


    Input / Output
    --------------
    See the help of the individual functions for explanations of in/out, etc.


    Examples
    --------
    p, cov  = opt.curve_fit(ufz.functions.f1x, x, y, p0=[p0,p1])
    p = opt.fmin(ufz.functions.cost_f1x, np.array([p1,p2]), args=(x,y), disp=False)


    License
    -------
    This file is part of the UFZ Python package.

    The UFZ Python package is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The UFZ Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2014 Matthias Cuntz


    History
    -------
    Written,  MC, Oct 2014
"""
from .functions import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1796"
__date__     = 'Date: 05.10.2014'
