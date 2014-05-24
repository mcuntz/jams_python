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
    >>> Rref = 1.0
    >>> E0   = 126.
    >>> T    = 293.15
    >>> resp = 2.0
    >>> from autostring import astr
    >>> print(astr(lloyd_fix(T, Rref, E0),3,pp=True))
    1.406
    >>> print(astr(lloyd_fix_p(T, [Rref, E0]),3,pp=True))
    1.406
    >>> print(astr(cost_lloyd_fix([Rref, E0], T, resp),3,pp=True))
    0.594
    >>> print(astr(cost2_lloyd_fix([Rref, E0], T, resp),3,pp=True))
    0.353

    >>> print(astr(poly(T,2,1),3,pp=True))
    295.150
    >>> print(astr(poly_p(T,[2,1]),3,pp=True))
    295.150

    >>> print(astr(griewank([0,0]),3,pp=True))
    0.000
    >>> print(astr(goldstein_price([0,-1]),3,pp=True))
    3.000


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
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2012-2013 Matthias Cuntz


    History
    -------
    Written,  MC, Dec 2012
    Modified, MC, Feb 2013 - ported to Python 3
              MC, May 2013 - general cost function cost_abs, cost_square
              MC, Oct 2013 - test functions such as Rosenbrock, Griewank, etc.
              MC, Feb 2014 - line0
"""
from __future__ import print_function
import numpy as np
import const # from ufz

# -----------------------------------------------------------
# general cost functions
def cost_abs(p, func, x, y):
    """ General cost function for robust optimising func(p, x) vs y with sum of absolute deviations"""
    return np.sum(np.abs(y-func(x,p)))

def cost_square(p, func, x, y):
    """ General cost function for least square optimising func(p, x) vs y"""
    return np.sum((y-func(x,p))**2)


# -----------------------------------------------------------
# arrhenius
def arrhenius(T, E):
  '''Arrhenius temperature dependence of rates
       T in C
       E in J
  '''
  return np.exp((T-(const.T25-const.T0))*E/(const.T25*const.R*(T+const.T0)))

def arrhenius_p(T, p):
  '''Arrhenius temperature dependence of rates
       T    in C
       p[0] in J
  '''
  return np.exp((T-(const.T25-const.T0))*p[0]/(const.T25*const.R*(T+const.T0)))

def cost_arrhenius(p, T, rate):
    """ Cost function for arrhenius with sum of absolute deviations """
    return np.sum(np.abs(rate-arrhenius_p(T,p)))

def cost2_arrhenius(p, T, rate):
    """ Cost function for arrhenius with sum of squared deviations """
    return np.sum((rate-arrhenius_p(T,p))**2)


# -----------------------------------------------------------
# a+b/x
def f1x(x,a,b):
  ''' General 1/x function: a + b/x
        x    independent variable
        a    1. parameter
        b    2. parameter
  '''
  return a+b/x

def f1x_p(x,p):
  ''' General 1/x function: a + b/x
        x    independent variable
        p    array of size 2, parameters
  '''
  return p[0]+p[1]/x

def cost_f1x(p,x,y):
  ''' Sum of absolut errors between obs and general 1/x function: a + b/x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-f1x_p(x,p)))

def cost2_f1x(p,x,y):
  ''' Sum of squared errors between obs and general 1/x function: a + b/x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-f1x_p(x,p))**2)


# -----------------------------------------------------------
# a+b*exp(c*x)
def fexp(x,a,b,c):
  ''' General exponential function: a + b * exp(c*x)
        x    independent variable in exp
        a    1. parameter
        b    2. parameter
        c    3. parameter
  '''
  return a+b*np.exp(c*x)

def fexp_p(x,p):
  ''' General exponential function: a + b * exp(c*x)
        x    independent variable in exp
        p    array of size 3, parameters
  '''
  return p[0]+p[1]*np.exp(p[2]*x)

def cost_fexp(p,x,y):
  ''' Sum of absolut errors between obs and general exponential function: a + b * exp(c*x)
        p    array of size 3, parameters
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-fexp_p(x,p)))

def cost2_fexp(p,x,y):
  ''' Sum of squared errors between obs and general exponential function: a + b * exp(c*x)
        p    array of size 3, parameters
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum((y-fexp_p(x,p))**2)


# -----------------------------------------------------------
# Gauss: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
def gauss(x,mu,sig):
  ''' Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        x    independent variable
        mu   mean
        sig  width
  '''
  return np.exp(-(x-mu)**2/(2.*sig**2))/(sig*np.sqrt(2.*np.pi))

def gauss_p(x,p):
  ''' Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        x    independent variable in exp
        p    array of size 2, parameters mu and sig
  '''
  return np.exp(-(x-p[0])**2/(2.*p[1]**2))/(p[1]*np.sqrt(2.*np.pi))

def cost_gauss(p,x,y):
  ''' Sum of absolut errors between obs and Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        p    array of size 2, parameters mu and sig
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-gauss_p(x,p)))

def cost2_gauss(p,x,y):
  ''' Sum of squared errors between obs and Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        p    array of size 2, parameters mu and sig
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum((y-gauss_p(x,p))**2)

# -----------------------------------------------------------
# lasslop
def lasslop(Rg, et, VPD, alpha, beta0, k, Rref):
    """ Lasslop et al. (2010) is basically the rectangular, hyperbolic
        light-response of NEE as by Falge et al. (2001), where the
        respiration is calculated with Lloyd & Taylor (1994), and the
        maximum canopy uptake rate at light saturation decreases
        exponentially with VPD as in Koerner (1995).
        Rg      Global radiation [W m-2]
        et      Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
        VPD     Vapour Pressure Deficit [Pa]
        alpha   Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
        beta0   Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
        k       e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
        Rref    Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = Rref*et
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.maximum(np.minimum(-k*(VPD-VPD0), 600.), -600.)
    beta  = np.where(VPD > VPD0, beta0*np.exp(kk), beta0)
    return -alpha*beta*Rg/(alpha*Rg+beta) + gamma

def lasslop_p(Rg, p):
    """ Lasslop et al. (2010) is basically the rectangular, hyperbolic
        light-response of NEE as by Falge et al. (2001), where the
        respiration is calculated with Lloyd & Taylor (1994), and the
        maximum canopy uptake rate at light saturation decreases
        exponentially with VPD as in Koerner (1995).
        Rg     Global radiation [W m-2]
        p[0]   Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
        p[1]   Vapour Pressure Deficit [Pa]
        p[2]   Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
        p[3]   Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
        p[4]   e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
        p[5]   Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = p[5]*p[0]
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.maximum(np.minimum(-p[4]*(p[1]-VPD0), 600.), -600.)
    beta  = np.where(p[1] > VPD0, p[3]*np.exp(kk), p[3])
    return -p[2]*beta*Rg/(p[2]*Rg+beta) + gamma

def cost_lasslop(p, Rg, et, VPD, NEE):
    """ Cost function for Lasslop with sum of absolute deviations """
    return np.sum(np.abs(NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3])))

def cost2_lasslop(p, Rg, et, VPD, NEE):
    """ Cost function for Lasslop with sum of squared deviations """
    return np.sum((NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3]))**2)


# -----------------------------------------------------------
# a+b*x
def line(x,a,b):
  ''' Straight line: a + b*x
        x    independent variable
        a    1. parameter
        b    2. parameter
  '''
  return a+b*x

def line_p(x,p):
  ''' Straight line: a + b*x
        x    independent variable
        p    array of size 2, parameters
  '''
  return p[0]+p[1]*x

def cost_line(p,x,y):
  ''' Sum of absolut errors between obs and straight line: a + b*x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-line_p(x,p)))

def cost2_line(p,x,y):
  ''' Sum of squared errors between obs and straight line: a + b*x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-line_p(x,p))**2)


# -----------------------------------------------------------
# b*x
def line0(x,a):
  ''' Straight line through origin: a*x
        x    independent variable
        a    parameter
  '''
  return a*x

def line0_p(x,p):
  ''' Straight line through origin: a*x
        x    independent variable
        a    parameter
  '''
  return p*x

def cost_line0(p,x,y):
  ''' Sum of absolut errors between obs and straight line though origin: a*x
        p    parameter
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-line0_p(x,p)))

def cost2_line0(p,x,y):
  ''' Sum of squared errors between obs and straight line through origin: a*x
        p    parameter
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-line0_p(x,p))**2)


# -----------------------------------------------------------
# lloyd_fix
def lloyd_fix(T, Rref, E0):
    """ Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        T       Temperature [k]
        Rref    Respiration at Tref=10 degC [umol(C) m-2 s-1]
        E0      Activation energy [K]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return Rref*np.exp(E0*(1./(Tref-T0)-1./(T-T0)))

def lloyd_fix_p(T, p):
    """ Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        T       Temperature [k]
        p[0]    Respiration at Tref=10 degC [umol(C) m-2 s-1]
        p[1]    Activation energy [K]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return p[0]*np.exp(p[1]*(1./(Tref-T0)-1./(T-T0)))

def cost_lloyd_fix(p, T, resp):
    """ Cost function for Lloyd with sum of absolute deviations """
    return np.sum(np.abs(resp-lloyd_fix_p(T,p)))

def cost2_lloyd_fix(p, T, resp):
    """ Cost function for Lloyd with sum of squared deviations """
    return np.sum((resp-lloyd_fix_p(T,p))**2)


# -----------------------------------------------------------
# lloyd_only_rref
def lloyd_only_rref(et, Rref):
    """ If E0 is know in Lloyd & Taylor (1994) then one can calc
        the exponential term outside the routine and the fitting
        becomes linear.
        et      exp-term in Lloyd & Taylor
        Rref    Respiration at Tref=10 degC [umol(C) m-2 s-1]
    """
    return Rref*et

def lloyd_only_rref_p(et, p):
    """ If E0 is know in Lloyd & Taylor (1994) then one can calc
        the exponential term outside the routine and the fitting
        becomes linear.
        et      exp-term in Lloyd & Taylor
        p[0]    Respiration at Tref=10 degC [umol(C) m-2 s-1]
    """
    return p[0]*et

def cost_lloyd_only_rref(p, et, resp):
    """ Cost function for rref with sum of absolute deviations """
    return np.sum(np.abs(resp-lloyd_only_rref_p(et,p)))

def cost2_lloyd_only_rref(p, et, resp):
    """ Cost function for rref  with sum of squared deviations """
    return np.sum((resp-lloyd_only_rref_p(et,p))**2)


# -----------------------------------------------------------
# c0*x[0] + c1*x[1] + c2*x[2] + ... + cn*x[n]
def multiline_p(p,*args):
  ''' Multiple linear regression line: c0*x[0] + c1*x[1] + c2*x[2] + ... + cn*x[n]
        p      constants
        *args  independent variables
  '''
  return np.polynomial.polynomial.polyval(x, list(args))


# -----------------------------------------------------------
# sqrt(a + b/x) - theoretical form of Jackknife-after-bootstrap

def sabx(x, a, b):
  """ sqrt(a + b/x)
        a, b  parameters
        x     independent variable
  """
  return np.sqrt(a+b/x)

def sabx_p(x, p):
  """ sqrt(a + b/x)
        p    array of size 2, parameters
        x    independent variable
  """
  return np.sqrt(p[0]+p[1]/x)

def cost_sabx(p,x,y):
  ''' Cost function for sqrt of general 1/x function with sum of absolute deviations
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-sabx_p(x,p)))

def cost2_sabx(p,x,y):
  ''' Cost function for sqrt of general 1/x function with sum of squared deviations
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-sabx_p(x,p))**2)


# -----------------------------------------------------------
# c0 + c1*x + c2*x**2 + ... + cn*x**n
def poly(x,*args):
  ''' General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        x    independent variable
        a    1. parameter
        b    2. parameter
        ...
  '''
  return np.polynomial.polynomial.polyval(x, list(args))

def poly_p(x,p):
  ''' General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        x    independent variable
        p    array of size n, parameters
  '''
  return np.polynomial.polynomial.polyval(x, p)

def cost_poly(p,x,y):
  ''' Sum of absolut errors between obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        p    array of size n, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-poly_p(x,p)))

def cost2_poly(p,x,y):
  ''' Sum of squared errors between obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        p    array of size n, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-poly_p(x,p))**2)


# -----------------------------------------------------------
# Test functions
def ackley(x):
    '''
    This is the Ackley Function
    Global Optimum (n>=2): 0.0 at origin
    '''
    a = 20.0
    b = 0.2
    c = 2.0*np.pi

    n  = np.size(x)
    s1 = np.sum(x**2)
    s2 = np.sum(np.cos(c*x))
    f  = -a * np.exp(-b*np.sqrt(1.0/n*s1)) - np.exp(1.0/n*s2) + a + np.exp(1.0)

    return f


def goldstein_price(x):
    '''
    This is the Goldstein-Price Function
    Bound X1=[-2,2], X2=[-2,2]
    Global Optimum: 3.0,(0.0,-1.0)
    '''
    x1 = x[0]
    x2 = x[1]
    u1 = (x1 + x2 + 1.0)**2
    u2 = 19. - 14.*x1 + 3.*x1**2 - 14.*x2 + 6.*x1*x2 +3.*x2**2
    u3 = (2.*x1 - 3.*x2)**2
    u4 = 18. - 32.*x1 + 12.*x1**2 + 48.*x2 -36.*x1*x2 + 27.*x2**2
    u5 = u1 * u2
    u6 = u3 * u4
    f = (1. + u5) * (30. + u6)
    return f


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


def rastrigin(x):
    '''
    This is the Rastrigin Function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, (0,0)
    '''
    x1 = x[0]
    x2 = x[1]
    f = x1**2 + x2**2 - np.cos(18.0*x1) - np.cos(18.0*x2)
    return f


def rosenbrock(x):
    '''
    This is the Rosenbrock Function
    Bound: X1=[-5,5], X2=[-2,8]; Global Optimum: 0,(1,1)
           bl=[-5 -5]; bu=[5 5]; x0=[1 1];
    '''

    x1 = x[0]
    x2 = x[1]
    a = 100.0
    f = a * (x2 - x1**2)**2 + (1 - x1)**2
    return f


def six_hump_camelback(x):
    '''
    This is the Six-hump Camelback Function.
    Bound: X1=[-5,5], X2=[-5,5]
    True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    x1 = x[0]
    x2 = x[1]
    f = (4 - 2.1*x1**2 + x1**4/3)*x1**2 + x1*x2 + (-4 + 4*x2**2)*x2**2
    return f


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # Rref = 1.0
    # E0   = 126.
    # T    = 293.15
    # resp = 2.0
    # print(lloyd_fix(T, Rref, E0))
    # #1.40590910521
    # print(lloyd_fix_p(T, [Rref, E0]))
    # #1.40590910521
    # print(cost_lloyd_fix([Rref, E0], T, resp))
    # #0.59409089479
    # print(cost2_lloyd_fix([Rref, E0], T, resp))
    # #0.352943991272

    # print(poly(T,2,1))
    # #295.15
    # print(poly_p(T,[2,1]))
    # #295.15
    # print(griewank([0,0]))
    # #0.0
    # print(goldstein_price([0,-1]))
    # #3.0
