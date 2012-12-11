#!/usr/bin/env python
import numpy as np
import const # from ufz

"""
    Defines common functions that are used in curve_fit or fmin parameter estimations.

    Defines the functions in two forms (ex. of 3 params):
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


    Definition
    ----------
    Current functions are (there is always the second form with the name appended by _p;
    these are used in the cost functions).
        arrhenius         1 param:  Arrhenius temperature dependence of biochemical rates
        f1x               2 params: General 1/x function: a + b/x
        fexp              3 params: General exponential function: a + b * exp(c*x)
        lasslop           6 params: Lasslop et al. (2010) a rectangular, hyperbolic light-response GPP
                                    with Lloyd & Taylor (1994) respiration and the maximum canopy uptake
                                    rate at light saturation decreases exponentially with VPD as in Koerner (1995)
        lloyd_fix         2 params: Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        lloyd_only_rref   1 param:  Lloyd & Taylor (1994) Arrhenius type with fixed exponential term


    Input / Output
    --------------
    See the help of the individual functions for explanations of in/out, etc.


    Examples
    --------
    >>> Rref = 1.0
    >>> E0   = 126.
    >>> T    = 293.15
    >>> resp = 2.0
    >>> print lloyd_fix(T, Rref, E0)
    1.40590910521
    >>> print lloyd_fix_p(T, [Rref, E0])
    1.40590910521
    >>> print cost_lloyd_fix([Rref, E0], T, resp)
    0.59409089479
    >>> print cost2_lloyd_fix([Rref, E0], T, resp)
    0.352943991272


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
    along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2012 Matthias Cuntz


    History
    -------
    Written, MC, Dec 2012
"""

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
        a    1. const
        b    2. const
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
        a    1. const
        b    2. const
        c    3. const
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # Rref = 1.0
    # E0   = 126.
    # T    = 293.15
    # resp = 2.0
    # print lloyd_fix(T, Rref, E0)
    # #1.40590910521
    # print lloyd_fix_p(T, [Rref, E0])
    # #1.40590910521
    # print cost_lloyd_fix([Rref, E0], T, resp)
    # #0.59409089479
    # print cost2_lloyd_fix([Rref, E0], T, resp)
    # #0.352943991272
