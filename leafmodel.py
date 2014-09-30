#!/usr/bin/env python
import numpy as np
from const import R, T0
from esat import esat
from errormeasures import mae
from scipy.optimize import fmin
from sce import sce

'''
    leafmodel: contains routines to model photosynthesis and stomatal conductance
    of canopies or leaves.
    
    
    major functions contained:
    --------------------------
    opt_leafmodel:  optimization against observed transpiration and
                    gross primary production
    cost_leafmodel: cost function for optimization
    leafmodel:      model for photosynthesis and stomatal conductance
    twoleafmodel:   uses leafmodel in a sunlit and shaded scheme
    farquhar:       photosysnthesis model
    leuning:        stomatal conductance modeling
    ball_berry:     stomatal conductance modeling


    License
    -------
    This file is part of the UFZ Python library.

    It is NOT released under the GNU Lesser General Public License, yet.

    If you use this routine, please contact Arndt Piayda.

    Copyright 2014 Arndt Piayda, Matthias Cuntz


    History
    -------
    Written,  AP+MC, Jul 2014
'''
###############################################################################
# Constants

# Michaelis-Menten constant for CO2 at 25 degC
Kc25    = 460.e-6 # [mol/mol]
# Activation energy for CO2
E_Kc    = 59356.  # [J]
# O2 partial pressure
O       = 0.21    # [mol/mol]
# Activation energy for O
E_O     = 35948.  # [J]
# Michaelis-Menten constant for O2 at 25 degC
Ko25    = 0.33    # [mol/mol]
# Activation energy for O2
E_Ko    = 35948.  # [J]
# effiency of photon capture
alpha   = 0.28    # [-]
# Activation energy for dark respiration Rd
E_Rd    = 50967.  # [J]
# decrease rate of Vcmax above temperature optimum (Medlyn approach)
HdV     = 200000. # [J/mol]
# decrease rate of Jmax above temperature optimum (Medlyn approach)
HdJ     = 220000. # [J/mol]
# Activation energy for Vcmax
E_Vcmax = 58520.  # [J]
# Activation energy for Jcmax
E_Jmax  = 35870.  # [J]

###############################################################################
def opt_leafmodel(params, Eobs, GPPobs, ci_ini, Tl, PAR, ea, ca, ga, gb, P,
                  PARdiff=None, fsun=None, dooptim=False, loc_search=True,
                  silent=True, n=5, eta=0.9):
    '''
    Optimization routine for leafmodel or twoleafmodel to infer model parameters
    fitting gross primary productivity and transpiration to observed values.
    
    
    Definition
    ----------
    def opt_leafmodel(params, Eobs, GPPobs, ci_ini, Tl, PAR, ea, ca, ga, gb, P,
                  PARdiff=None, fsun=None, dooptim=False, loc_search=True,
                  n=5, eta=0.9):
    
    
    Input 
    -----
    params        dict, parameters for leafmodel. Parameters given here decide
                  which conductance model and temperature dependence model is
                  used and which parameters shall be optimized and which not.
                  params CAN contain following keys with its respective initial
                  values: 'm','b','d0','Vcmax25','delSV','delSJ','Topt','omega'.
                  All parameters given in this way will be optimized. You can
                  give parameters also with the ending 'fix' like 'mfix' or
                  Toptfix', then they will not be optimized but kept fix during
                  optimization.
                  params MUST contain 'm' and 'b' (or with the 'fix' ending).
                  If you give 'd0' or 'd0fix', the Leuning model will be used,
                  else the Ball&Berry model will be used.
                  params MUST contain 'Vcmax25' or 'Vcamx25fix'. If you give
                  'Topt' and 'omega' (or with the fix ending), the June model
                  is used. If you give 'delSV' and 'delSJ' (or with the fix
                  ending), the Medlyn model is used. If you only give 'delSJ' or
                  'delSJfix', the von Caemmerer model is used. For model
                  descriptions see leafmodel. 
    Eobs          np.array(N), observed transpiration [mol/m2s]
    GPPobs        np.array(N), observed gross primary productivity [mol/m2s]
    ci_ini, Tl, PAR, ea, ca, ga, gb, P: same as for leafmodel
    
    
    Optional Input
    --------------
    PARdiff, fsun if both set to None (default) leafmodel is used, else
                  twoleafmodel is used and input is like twoleafmodel
    dooptim       bool, if True optimization is conducted, is False (default)
                  no optimization is conducted and output is calculated using
                  the input params
    loc_search    bool, if True optimization uses Nelder-Mead local search
                  (scipy.optimize.fmin) (default), if False shuffled complex
                  evolution global (ufz.sce) search is used. fmin is faster, but
                  params must be close to the optimum parameters. sce is slower
                  but always finds the minimum, even if params is far away from
                  optimum. The parameter limits of sce are set to 0 and 2*params
                  for each parameter.
    silent        bool, if True (default), no output during optimization, if
                  False objective function value is printed to the console at
                  each cost function evaluation
    n, eta        same as for leafmodel
    
    
    Output
    ------
    opt_params    list, optimized parameters (does not contain the fixed values)
    ci_conv, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod:
                  same as leafmodel
    
    
    Examples
    --------
    >>> # Define input
    >>> ci_ini = np.array([3.361e-4, 3.354e-4, 3.350e-4, 3.345e-4, 3.337e-4])
    >>> Tl     = np.array([21.041, 21.939, 21.997, 21.026, 20.237])
    >>> PAR    = np.array([1.144e-3, 1.294e-3, 1.445e-3, 1.733e-3, 1.752e-3])
    >>> ea     = np.array([1636.967, 1609.575, 1569.418, 1561.334, 1460.839])
    >>> ca     = np.array([4.202e-4, 4.193e-4, 4.188e-4, 4.182e-4, 4.171e-4])
    >>> ga     = np.array([3.925, 5.998, 3.186, 7.924, 5.501])
    >>> gb     = np.array([1.403, 1.678, 1.732, 1.691, 1.782])
    >>> P      = np.array([99893.7, 99907.7, 99926.5, 99928.0, 99924.7])
    >>> PARdiff= np.array([3.746e-4, 4.059e-4, 4.428e-4, 4.663e-4, 4.849e-4])
    >>> fsun   = np.array([0.682, 0.703, 0.718, 0.729, 0.737])
    >>> Eobs   = np.array([3.552e-3, 3.854e-3, 3.828e-3, 4.455e-3, 5.322e-3])
    >>> GPPobs = np.array([1.293e-5, 1.351e-5, 1.339e-5, 1.318e-5, 1.139e-5])
    
    >>> # use default methods Leuning + June with fixed b, Topt and omega
    >>> params = {'m':15., 'bfix':0.0041, 'd0':5.e3,\
    'Vcmax25':60.e-6, 'Toptfix':25., 'omegafix':18.}
    >>> opt_params, ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod = \
    opt_leafmodel(params, Eobs, GPPobs, ci_ini, Tl, PAR, ea, ca, ga, gb, P,\
    PARdiff, fsun, dooptim=True)
    >>> # print GPP [mumol(m2s] and E [mmol/m2s]
    >>> print np.ma.round(GPPmod*1.e6, 3)
    [12.877 13.29 13.384 13.273 13.013]
    >>> print np.ma.round(Emod*1.e3, 3)
    [3.546 4.442 4.306 4.31 4.157]
    >>> # print optimized parameters
    >>> m, d0, Vcmax25 = opt_params
    >>> print np.round(m, 3)
    48.346
    >>> print np.round(d0/100., 3) #[hPa]
    4.343
    >>> print np.round(Vcmax25*1.e6, 3) #[mumol/m2s]
    65.763
    
    
    License
    -------
    This file is part of the UFZ Python library.

    It is NOT released under the GNU Lesser General Public License, yet.

    If you use this routine, please contact Arndt Piayda.

    Copyright 2014 Arndt Piayda, Matthias Cuntz


    History
    -------
    Written,  AP+MC, Jul 2014
    '''
    # select methods according to given paramters and sort them according to
    # fixed or non fixed
    gs_method = 'BB'
    fix       = {}
    
    if ('m' in params):
        guess = [params['m']]
    else:
        fix.update({'mfix':params['mfix']})
    if ('b' in params):
        guess += [params['b']]
    else:
        fix.update({'bfix':params['bfix']})
    
    if ('d0' in params):
        guess += [params['d0']]
        gs_method = 'Leuning'
    elif('d0fix' in params):
        fix.update({'d0fix':params['d0fix']})
        gs_method = 'Leuning'
    
    if ('Vcmax25' in params):
        guess += [params['Vcmax25']]
    else:
        fix.update({'Vcmax25fix':params['Vcmax25fix']})
                            
    temp_method = 'June'
    if ('Topt' in params):
        guess += [params['Topt']]
        assert ('omega' in params) or ('omegafix' in params)
    elif ('Toptfix' in params):
        fix.update({'Toptfix':params['Toptfix']})
        assert ('omega' in params) or ('omegafix' in params)
    if ('omega' in params):
        guess += [params['omega']]
    elif ('omegafix' in params):
        fix.update({'omegafix':params['omegafix']})
    else:
        if ('delSV' in params):
            guess += [params['delSV']]
            temp_method = 'Medlyn'
        elif ('delSVfix' in params):
            fix.update({'delSVfix':params['delSVfix']})
            temp_method = 'Medlyn'
        else:
            temp_method = 'Caemmerer'
    
        if ('delSJ' in params):
            guess += [params['delSJ']]
        elif ('delSJfix' in params):
            fix.update({'delSJfix':params['delSJfix']})
    
    # optimization
    if dooptim:
        # use local nelmin optimization
        if loc_search:
            opt_params = fmin(cost_leafmodel, guess,
                              args=(Eobs, GPPobs, ci_ini, Tl, PAR, ea, ca, ga,
                                    gb, P, PARdiff, fsun, n, eta, gs_method,
                                    temp_method, fix, silent), disp=0)
        # use global sce optimization
        else:
            def wrap_opt(opt_params):
                obj = cost_leafmodel(opt_params, Eobs, GPPobs, ci_ini, Tl, PAR,
                                     ea, ca, ga, gb, P, PARdiff,fsun, n, eta,
                                     gs_method, temp_method, fix, silent)
                return obj
            guess=np.array(guess)
            opt_params = sce(wrap_opt, guess, np.zeros_like(guess) , guess*2.,
                             maxn=10000, kstop=15)    
    else:
        opt_params = guess
        
    # distribute parameters
    par1, par2 = dist_params(opt_params, fix, gs_method, temp_method)
    
    # calculate output
    if (PARdiff==None) and (fsun==None):
        ci, gcmol_c, gcmol_h, gsmol_c, gsmol_h, no_conv, Jc, Je, Rd, GPPmod, Emod =\
            leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,
                      gs_method=gs_method, temp_method=temp_method, n=n, eta=eta)
    elif (PARdiff!=None) and (fsun!=None):
        ci, gcmol_c, gcmol_h, gsmol_c, gsmol_h, no_conv, Jc, Je, Rd, GPPmod, Emod =\
            twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun, par1, par2,
                         gs_method=gs_method, temp_method=temp_method, n=n, eta=eta)
    else:
        raise ValueError('opt_leafmodel: if two leaf model is desired, PARdiff and fsun must be given')
    
    return opt_params, ci, gcmol_c, gcmol_h, gsmol_c, gsmol_h, no_conv, Jc, Je, Rd, GPPmod, Emod

###############################################################################
def cost_leafmodel(params, Eobs, GPPobs, ci_ini, Tl, PAR, ea, ca, ga, gb, P,
                   PARdiff, fsun, n, eta, gs_method, temp_method, fix,
                   silent=True):
    '''
    opt_leafmodel: cost function for optimization
    (input see leafmodel and twoleafmodel)
    
    returns: obj objective function value
    '''
    # distribute parameters
    par1, par2 = dist_params(params, fix, gs_method, temp_method)
    
    if (PARdiff==None) and (fsun==None):
        # run the model with current parameters
        ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod =\
            leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,
                      gs_method=gs_method, temp_method=temp_method, n=n,
                      eta=eta)
    elif (PARdiff!=None) and (fsun!=None):
        ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod =\
            twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun,
                         par1, par2, gs_method=gs_method,
                         temp_method=temp_method, n=n, eta=eta)
    else:
        raise ValueError('cost_leafmodel: if two leaf model is desired, PARdiff and fsun must be given')
        
    # mask values where ci did not converged
    Emod.mask[no_conv], GPPmod.mask[no_conv] = True, True
    # calculate range of observations
    Eobs_r   = np.ma.ptp(Eobs)
    GPPobs_r = np.ma.ptp(GPPobs)
    # if for some time steps ci did not converge, tell fmin that's bad! 
    if no_conv.size>0:
        return 1000.
    # else calculate objective function
    else:
        # combined mean absolute error
        obj = (
               mae(Eobs/Eobs_r, Emod/Eobs_r)**6 +
               mae(GPPobs/GPPobs_r, GPPmod/GPPobs_r)**6
               )**(1./6.)
        if not silent:
            print obj
        
        return obj

###############################################################################
def dist_params(params, fix, gs_method, temp_method):
    '''
    opt_leafmodel+cost_leafmodel: distibutes parameters to variables according
    to selected methods
    '''
    params = np.insert(params, 0, fix['mfix']) if 'mfix' in fix else params
    params = np.insert(params, 1, fix['bfix']) if 'bfix' in fix else params
    
    if gs_method=='Leuning':
        params = np.insert(params, 2, fix['d0fix']) if 'd0fix' in fix else params
        params = np.insert(params, 3, fix['Vcmaxfix']) if 'Vcmaxfix' in fix else params        
        if temp_method=='June':
            params = np.insert(params, 4, fix['Toptfix']) if 'Toptfix' in fix else params
            params = np.insert(params, 5, fix['omegafix']) if 'omegafix' in fix else params
        elif temp_method=='Medlyn':
            params = np.insert(params, 4, fix['delSVfix']) if 'delSVfix' in fix else params
            params = np.insert(params, 5, fix['delSJfix']) if 'delSJfix' in fix else params
        elif temp_method=='Caemmerer':
            params = np.insert(params, 4, fix['delSJfix']) if 'delSJfix' in fix else params
    else:
        params = np.insert(params, 2, fix['Vcmaxfix']) if 'Vcmaxfix' in fix else params
        if temp_method=='June':
            params = np.insert(params, 3, fix['Toptfix']) if 'Toptfix' in fix else params
            params = np.insert(params, 4, fix['omegafix']) if 'omegafix' in fix else params
        elif temp_method=='Medlyn':
            params = np.insert(params, 3, fix['delSVfix']) if 'delSVfix' in fix else params
            params = np.insert(params, 4, fix['delSJfix']) if 'delSJfix' in fix else params
        elif temp_method=='Caemmerer':
            params = np.insert(params, 3, fix['delSJfix']) if 'delSJfix' in fix else params
    
    if gs_method=='Leuning':
        par1 = params[:3]
        par2 = params[3:]
    if gs_method=='BB':
        par1 = params[:2]
        par2 = params[2:]
    
    return par1, par2

###############################################################################
def leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,
              gs_method='Leuning', temp_method='June', n=5, eta=0.9):
    '''
    Model to calculate photosynthesis and stomatal conductance of canopies. The
    Farquhar model (Farquhar et al, 1980) is used for photosynthesis. Stomatal
    conductance can either be modeled with Ball & Berry (1987) or with
    Leuning (1995). Aerodynamic and leaf boundary layer conductances are used
    to calculate leaf surface concentrations and the model iterates until a
    conversion of leaf internal CO2 concentration. Canopy transpiration and
    gross primary production are among the results. The temperature dependency
    of the photosynthesis model can be modeled according to June et al. (2004),
    Medlyn et al. (2002) and von Caemmerer (2000).
     
    
    Definition
    ----------
    def leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,
              gs_method='Leuning', temp_method='June', n=5, eta=0.9):
    
    
    Input
    -----
    ci_ini      np.array(N), initial leaf internal co2 concentration [mol/mol]
                Can be approximated with 0.8*ca
    Tl          np.array(N), leaf temperature [degC]
    PAR         np.array(N), total (direct+diffuse) photon flux density
                [mol/m2s]
    ea          np.array(N), vapour pressure of the atmosphere [Pa]
    ca          np.array(N), atmospheric CO2 concentration [mol/mol]
    ga          np.array(N), aerodynamic conductivity [mol/m2s] (of heat and
                water, carbon is calulated in the model)
    gb          np.array(N), leaf boundary layer conductivity [mol/m2s] (of heat
                and water, carbon is calulated in the model)
    P           np.array(N), atmospheric pressure [Pa]
    par1        if gs_method='Leuning' (default): list(3)
                    [0] m    slope of the Leuning model [mol(H2O)/mol(air)]
                    [1] b    offset of the Leuning model [mol(H2O)/m2leaf s]
                    [2] d0   sensitivity parameter to vpds [Pa]
                if gs_method='BB': list(2)
                    [0] m    slope of the BB model [mol(H2O)/mol(air)]
                    [1] b    offset of the BB model [mol(H2O)/m2leaf s]
    par2        if temp_method is June (default): list(3) 
                    [0] Vcmax25 maximum carboxylation rate at 25 degC [mol/m2s]
                    [1] Topt    optimum temperature of Jmax [degC]
                    [2] omega   the difference in temperature from Topt at
                                which Jmax falls to e-1 (0.37) of its value
                                at Topt. (In the paper its on average
                                18+-0.6 degC)
                if temp_method is Medlyn: list(3)
                    [0] Vcmax25 maximum carboxylation rate at 25 degC [mol/m2s]
                    [1] delSV   entropy factor for Vcmax25 [J/mol/K]
                    [2] delSJ   entropy factor for Jmax25 [J/mol/K]
                if temp_method is Caemmerer: list(2)
                    [0] Vcmax25 maximum carboxylation rate at 25 degC [mol/m2s]
                    [1] delSJ   entropy factor for Jmax25 [J/mol/K]

    
    Optional Input
    --------------
    gs_method   str, stomatal conductance method: 'Leuning' (default), 'BB'
                (Ball&Berry)
    temp_method str, temperature dependecy method: 'June', 'Medlyn', 'Caemmerer'
    n           int, number of iteration for leaf internal CO2 convergence
    eta         float, smoothing factor for transition between Jmax and Vcmax
                dominated assimilation rate [-]
                
    
    Output
    ------
    ci_conv     np.array(N), leaf internal CO2 concentration [mol/mol]
    gc_c        np.array(N), canopy conductivity for carbon [mol/m2s]
    gc_h        np.array(N), canopy conductivity for water and heat [mol/m2s]
    gs_c        np.array(N), stomatal conductivity for carbon [mol/m2s]
    gs_h        np.array(N), stomatal conductivity for water and heat [mol/m2s]
    no_conv     np.array(x), indices of values where no ci conversion was
                achieved with n iterations
    Jc          np.array(N), carboxylation limited assimilation rate [mol/m2s]
    Je          np.array(N), electron transport limited assimilation rate
                [mol/m2s]
    Rd          np.array(N), leaf dark respiration rate [mol/m2s]
    GPPmod      np.array(N), gross primary productivity [mol/m2s]
    Emod        np.array(N), transpiration [mol/m2s]

                         
    Examples
    --------
    >>> # Define input
    >>> ci_ini = np.array([3.361e-4, 3.354e-4, 3.350e-4, 3.345e-4, 3.337e-4])
    >>> Tl     = np.array([21.041, 21.939, 21.997, 21.026, 20.237])
    >>> PAR    = np.array([1.144e-3, 1.294e-3, 1.445e-3, 1.733e-3, 1.752e-3])
    >>> ea     = np.array([1636.967, 1609.575, 1569.418, 1561.334, 1460.839])
    >>> ca     = np.array([4.202e-4, 4.193e-4, 4.188e-4, 4.182e-4, 4.171e-4])
    >>> ga     = np.array([3.925, 5.998, 3.186, 7.924, 5.501])
    >>> gb     = np.array([1.403, 1.678, 1.732, 1.691, 1.782])
    >>> P      = np.array([99893.7, 99907.7, 99926.5, 99928.0, 99924.7])
    
    >>> # use default methods Leuning + June
    >>> par1   = [15., 0.0041, 5.e3]
    >>> par2   = [60.e-6, 25., 18.]
    >>> ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod = \
    leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2)
    >>> # print GPP [mumol(m2s] and E [mmol/m2s]
    >>> print np.ma.round(GPPmod*1.e6, 3)
    [12.333 12.685 12.696 12.462 12.157]
    >>> print np.ma.round(Emod*1.e3, 3)
    [3.339 4.338 4.287 4.02 3.845]
    
    >>> # use methods Ball&Berry + Medlyn
    >>> par1   = [15., 0.0041]
    >>> par2   = [60.e-6, 700., 700.]
    >>> ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod = \
    leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,\
    gs_method='BB', temp_method='Medlyn')
    >>> # print GPP [mumol(m2s] and E [mmol/m2s]
    >>> print np.ma.round(GPPmod*1.e6, 3)
    [14.286 14.496 14.499 14.526 14.287]
    >>> print np.ma.round(Emod*1.e3, 3)
    [2.718 3.326 3.22 3.167 2.994]
    
    
    License
    -------
    This file is part of the UFZ Python library.

    It is NOT released under the GNU Lesser General Public License, yet.

    If you use this routine, please contact Arndt Piayda.

    Copyright 2014 Arndt Piayda, Matthias Cuntz


    History
    -------
    Written,  AP+MC, Jul 2014
    '''
    
    # array for convergence of leaf internal carbon concentration ci
    ci_conv = np.ma.empty((n,ci_ini.size))
    # initial conditions for ci convergence
    ci_conv[0,:] = ci_ini # [mol(CO2)/mol(air)]
    Emod         = 0.     # [mol(H2O)/m2leaf s]
    gs_h         = 1.     # [mol(H2O)/m2leaf s]
    
    # atmospheric (wa) and leaf internal (wi) water concentration
    wa = ea/P    # [mol(H2O)/mol(air)]
    ei = esat(Tl+T0)
    wi = ei/P # [mol(H2O)/mol(air)]
    
    # loop for ci convergence
    for i in range(n-1):
        # assimilation A, carboxylation limited assimilation Jc,
        # electron transport limited assimilation Je and leaf dark respiration
        # rate Rd, all in [mol/m2leaf s]
        A, Jc, Je, Rd = farquhar(Tl, ci_conv[i,:], PAR, par2,
                                 temp_method=temp_method, eta=eta)
        
        # aerodynamic + leaf boundary layer conductances
        gab_c         = 1./(1.37/gb + 1./ga) # [mol(CO2)/m2leaf s] 
        gab_h         = 1./(1./gb + 1./ga) # [mol(H2O)/m2leaf s]
        
        # carbon concentration at the leaf surface 
        cs            = ca - A/gab_c #[mol(CO2)/mol(air)]
        # water concentration at the leaf surface
        ws            = wa - Emod/gab_h # [mol(H2O)/mol(air)]
        # vapour pressure at the leaf surface 
        es            = ws*P         # [Pa]
        # relative humidity at the leaf surface 
        rHs           = es/ei     # [-]
        
        # stomatal conductance for carbon [mol(CO2)/m2leaf s] 
        if gs_method=='Leuning': 
            gs_c = leuning(A, ei-es, cs, Tl, par1)
            #gs_c = leuning(A, wi-ws, cs, Tl, par1)
        if gs_method=='BB':
            gs_c = ball_berry(A, rHs, cs, par1)
        # stomatal conductance for water [mol(H2O)/m2leaf s]
        gs_h = gs_c * 1.6
        
        # canopy conductances 
        gc_c           = 1./(1./gs_c + 1./gab_c) # [mol(CO2)/m2leaf s] 
        gc_h           = 1./(1./gs_h + 1./gab_h) # [mol(H2O)/m2leaf s]
        # new ci value
        ci_conv[i+1,:] = ca - A/(gc_c) # [mol(CO2)/mol(air)]
        # new modelled transpiration
        Emod           = gc_h * (wi-wa) # [mol(H2O)/m2leaf s]
    
    # number of non-coverging ci values
    no_conv = np.ma.where(np.ma.abs(ci_conv[-1,:]-ci_conv[-2,:])>10.e-6)[0]
    
    # modeled transpiration
    Emod   = gc_h * ((ei-ea)/P) # [mol(H2O)/m2leaf s]
    # modeled gross primary production
    GPPmod = gc_c * (ca - ci_conv[-1,:]) + Rd       # [mol(CO2)/m2leaf s]
    
    return ci_conv[-1,:], gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod

###############################################################################
def twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun, par1, par2,
                 gs_method='Leuning', temp_method='June', n=5, eta=0.9):
    '''
    Model to calculate photosynthesis and stomatal conductance of canopies in a
    two leaf scheme. It calculates leafmodel for a sunlit and shaded leaf
    fraction of the canopy. Inputs are the same as for leafmodel plus PARdiff
    and fsun.
    
    Definition
    ----------
    def twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun, par1, par2,
                 gs_method='Leuning', temp_method='June', n=5, eta=0.9):
    
    
    Input 
    -----
    same as for leafmodel except:
    PARdiff     np.array(N), diffuse photon flux density [mol/m2s]
    fsun        np.array(N), fraction of sunlit leaf area [-]
    
    
    Optional Input
    --------------
    same as for leafmodel
    
    
    Output
    ------
    same as for leafmodel
    
    
    Examples
    --------
    >>> # Define input
    >>> ci_ini = np.array([3.361e-4, 3.354e-4, 3.350e-4, 3.345e-4, 3.337e-4])
    >>> Tl     = np.array([21.041, 21.939, 21.997, 21.026, 20.237])
    >>> PAR    = np.array([1.144e-3, 1.294e-3, 1.445e-3, 1.733e-3, 1.752e-3])
    >>> ea     = np.array([1636.967, 1609.575, 1569.418, 1561.334, 1460.839])
    >>> ca     = np.array([4.202e-4, 4.193e-4, 4.188e-4, 4.182e-4, 4.171e-4])
    >>> ga     = np.array([3.925, 5.998, 3.186, 7.924, 5.501])
    >>> gb     = np.array([1.403, 1.678, 1.732, 1.691, 1.782])
    >>> P      = np.array([99893.7, 99907.7, 99926.5, 99928.0, 99924.7])
    >>> PARdiff= np.array([3.746e-4, 4.059e-4, 4.428e-4, 4.663e-4, 4.849e-4])
    >>> fsun   = np.array([0.682, 0.703, 0.718, 0.729, 0.737])
    
    >>> # use default methods Leuning + June
    >>> par1   = [15., 0.0041, 5.e3]
    >>> par2   = [60.e-6, 25., 18.]
    >>> ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod = \
    twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun, par1, par2)
    >>> # print GPP [mumol(m2s] and E [mmol/m2s]
    >>> print np.ma.round(GPPmod*1.e6, 3)
    [11.861 12.255 12.344 12.176 11.924]
    >>> print np.ma.round(Emod*1.e3, 3)
    [3.253 4.232 4.205 3.952 3.792]
    
    >>> # use methods Ball&Berry + Medlyn
    >>> par1   = [15., 0.0041]
    >>> par2   = [60.e-6, 700., 700.]
    >>> ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod = \
    twoleafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, PARdiff, fsun, par1, par2,\
    gs_method='BB', temp_method='Medlyn')
    >>> # print GPP [mumol(m2s] and E [mmol/m2s]
    >>> print np.ma.round(GPPmod*1.e6, 3)
    [13.482 13.785 13.904 14.002 13.839]
    >>> print np.ma.round(Emod*1.e3, 3)
    [2.612 3.206 3.127 3.083 2.926]
    
    
    License
    -------
    This file is part of the UFZ Python library.

    It is NOT released under the GNU Lesser General Public License, yet.

    If you use this routine, please contact Arndt Piayda.

    Copyright 2014 Arndt Piayda, Matthias Cuntz


    History
    -------
    Written,  AP+MC, Jul 2014
    '''
    # sunlit leaves
    ci_l, gc_c_l, gc_h_l, gs_c_l, gs_h_l, no_conv_l, Jc_l, Je_l, Rd_l,\
    GPP_l, E_l = leafmodel(ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2,
                           gs_method=gs_method, temp_method=temp_method,
                           n=n, eta=eta)
        
    # shaded leaves
    ci_s, gc_c_s, gc_h_s, gs_c_s, gs_h_s, no_conv_s, Jc_s, Je_s, Rd_s,\
    GPP_s, E_s = leafmodel(ci_ini, Tl, PARdiff, ea, ca, ga, gb, P, par1, par2,
                           gs_method=gs_method, temp_method=temp_method,
                           n=n, eta=eta)
    
    # combine
    ci      = fsun*ci_l   + (1.-fsun)*ci_s
    gc_c    = fsun*gc_c_l + (1.-fsun)*gc_c_s
    gc_h    = fsun*gc_h_l + (1.-fsun)*gc_h_s
    gs_c    = fsun*gs_c_l + (1.-fsun)*gs_c_s
    gs_h    = fsun*gs_h_l + (1.-fsun)*gs_h_s
    no_conv = np.concatenate((no_conv_l,no_conv_s))
    Jc      = fsun*Jc_l   + (1.-fsun)*Jc_s
    Je      = fsun*Je_l   + (1.-fsun)*Je_s
    Rd      = fsun*Rd_l   + (1.-fsun)*Rd_s
    GPPmod  = fsun*GPP_l  + (1.-fsun)*GPP_s
    Emod    = fsun*E_l    + (1.-fsun)*E_s
    
    return ci, gc_c, gc_h, gs_c, gs_h, no_conv, Jc, Je, Rd, GPPmod, Emod

###############################################################################
def farquhar(Tl, ci, PAR, par, temp_method='June', eta=0.9):
    '''
    leafmodel: photosynthesis model
    Tl : leaf temperature [degC]
    ci : leaf intercellular co2 concentration [mol/mol]
    PAR: photosynthetic active photon flux density [mol/m2leaf s]
    par: if temp_method is June: 
         [0] VCmax25: maximum carboxylation rate at 25 degC [mol/m2leaf s]
         [1] Topt:    optimum temperature of Jmax [degC]
         [2] omega:   the difference in temperature from Topt at which Jmax
                      falls to e-1 (0.37) of its value at Topt. (In the paper
                      its on average 18+-0.6 degC)
         if temp_method is Medlyn:
         [0] VCmax25: maximum carboxylation rate at 25 degC [mol/m2leaf s]
         [1] delSV:   entropy factor for VCmax [J/mol/K]
         [2] delSJ:   entropy factor for Jmax [J/mol/K]
         if temp_method is Caemmerer:
         [0] VCmax25: maximum carboxylation rate at 25 degC [mol/m2leaf s]
         [1] delSJ:   entropy factor for Jmax [J/mol/K]
    temp_method: temperature dependecy methods: 'June', 'Medlyn', 'Caemmerer'
    eta: smoothiing factor for transition between Jmax and Vcmax dominated
         assimilation rate [-]
         
    returns:
    A : assimilation rate [mol/m2leaf s]
    Jc: carboxylation limited assimilation flux [mol/m2leaf s]
    Je: electron transport limited assimilation flux [mol/m2leaf s]
    Rd: leaf dark respiration [mol/m2leaf s]
    '''
    if temp_method=='June':
        assert len(par)==3, 'farquhar: wrong parameters for method June'
        Vcmax25 = par[0]
        Topt    = par[1]
        omega   = par[2]
    elif temp_method=='Medlyn':
        assert len(par)==3, 'farquhar: wrong parameters for method Medlyn'
        Vcmax25 = par[0]
        delSV   = par[1]
        delSJ   = par[2]
    elif temp_method=='Caemmerer':
        assert len(par)==2, 'farquhar: wrong parameters for method Caemmerer'
        Vcmax25 = par[0]
        delSJ   = par[1]
    else:
        raise ValueError('farquhar: unknown temperature method')
        
    # CO2 compensation point
    gamma_star = compensationpoint(Tl) # [mol/mol]
    # maximum electron transport rate at 25 degC
    Jm25 = 1.67*Vcmax25 # [mol/m-2s-1] Medlyn et al. 2002
    # leaf or dark respiration
    Rd25 = 0.011*Vcmax25 # [mol/m-2s-1]

    # calculate temperature dependencies
    Kc    = ArrheniusTemp(Kc25, E_Kc, Tl)
    Ko    = ArrheniusTemp(Ko25, E_Ko, Tl)    
    Rd    = ArrheniusTemp(Rd25, E_Rd, Tl)
    if temp_method=='June':
        Vcmax = ArrheniusTemp(Vcmax25, E_Vcmax, Tl)
        Jm    = JuneTemp(Jm25, Topt, Tl, omega)
    elif temp_method=='Medlyn':
        Vcmax = MedlynTemp(Vcmax25, E_Vcmax, delSV, Tl, HdV)
        Jm    = MedlynTemp(Jm25, E_Jmax, delSJ, Tl, HdJ) 
    elif temp_method=='Caemmerer':
        Vcmax = ArrheniusTemp(Vcmax25, E_Vcmax, Tl)
        Jm    = MedlynTemp(Jm25, E_Jmax, delSJ, Tl, HdJ) 

    # carboxylation limited assimilation flux [mol/m2leaf s]
    Jc = Vcmax * (ci - gamma_star) / (ci + Kc*(1. + O/Ko))
    # electron transport limited assimilation flux [mol/m2leaf s]
    Je = (alpha * PAR * Jm / np.ma.sqrt(alpha**2 * PAR**2 + Jm**2)) * \
         ((ci - gamma_star) / (4.*(ci + 2. * gamma_star)))
    # assimilation [mol/m2leaf s]
    A = smoothmin(Jc, Je, eta) - Rd
    
    return A , Jc, Je, Rd# [mol/m2leaf s]

###############################################################################
def compensationpoint(T):
    '''
    leafmodel: CO2 compensation point
    T: temperature [degC]
    
    returns:
    gamma_star: CO2 compensation point [mol/mol] 
    '''
    gamma_star = 1.7e-6 * T
    return gamma_star

###############################################################################
def ball_berry(A, rH, cs, par):
    '''
    leafmodel: stomatal conductance
    A: carbon assimilation rate [mol/m2leaf s]
    rH: relative humidity [-]
    cs: carbon mixing ratio at the leaf surface [mol/mol]
    par: [0] m: slope [mol(H2O)/mol(air)]
         [1] b: offset [mol(H2O)/m2leaf s]
    
    returns:
    gs: stomatal conductance [mol/m2leaf s] 
    '''
    m = par[0]
    b = par[1]
    return m * A * rH/cs + b

###############################################################################
def leuning(A, vpds, cs, Tl, par):
    '''
    leafmodel: stomatal conductance
    A: carbon assimilation rate [mol/m2leaf s]
    vpds: vapour pressure deficit at the leaf surface [Pa]
    cs: carbon mixing ratio at the leaf surface [mol/mol]
    Tl: temperature at the leaf surface [degC]
    par: [0] m:  slope [mol(H2O)/mol(air)]
         [1] b:  offset [mol(H2O)/m2leaf s]
         [2] d0: sensitivity parameter to vpds [mol(H2O)/mol(air) or Pa]
    
    returns:
    gs: stomatal conductance [mol/m2leaf s] 
    '''
    m  = par[0]
    b  = par[1]
    d0 = par[2]
    # CO2 compensation point
    gamma_star = compensationpoint(Tl) # [mol/mol]
    return m * A/(cs-gamma_star) * 1./(1.+vpds/d0) + b

###############################################################################
def smoothmin(x,y,eta):
    '''
    leafmodel: smooth minimum
    Smooth minumum function for the transition from Jmax to Vcmax dominated
    assimilation rate
    '''
    z = (x+y)**2 - 4.*eta*x*y
    z = np.ma.maximum(z,1.e-18)
    return (x+y-np.ma.sqrt(z)) / (2.*eta)

###############################################################################
def ArrheniusTemp(r25, E, T):
    '''
    leafmodel: temperature dependency
    Arrhenius style temperature dependency
    r25: activity rate at 25 degC [unit depends on input]
    E: activation energy [J]
    T: temperature [degC]
    
    returns: r activity rate at T [unit depends on input]
    '''
    return r25 * np.ma.exp(((T-25.)*E)/(298.*R*(T+273.)))

###############################################################################
def JuneTemp(r25, Topt, T, omega=18.):
    '''
    leafmodel: temperature dependency
    Temperature dependency after June et al. 2004
    r25 : activity rate at 25 degC [unit depends on input]
    Topt: optimum temperature [degC]
    T: temperature [degC]
    omega: the difference in temperature from Topt at which activity rate falls
           to e-1 (0.37) of its value at Topt. (In the paper its on average
           18+-0.6 degC for Jmax)
           
    returns:
    r: activity rate at temperature T [unit depends on input]
    '''
    return r25 * np.ma.exp((-(T-Topt)**2 + (25.-Topt)**2)/(omega**2))

###############################################################################
def MedlynTemp(r25, E, delS, T, Hd):
    '''
    leafmodel: temperature dependency
    Temperature dependency after Medlyn et al. 2002
    r25 : activity rate at 25 degC [unit depends on input]
    E: activation energy [J]
    delS: entropy factor [J/mol/K]
    T: temperature [degC]
    Hd : decrease rate above optimum [J/mol]
    
    returns:
    r: activity rate at temperature T [unit depends on input]
    '''
    return r25 * np.ma.exp(((T-25.)*E)/(298.*R*(T+273.))) * \
                      (1. + np.ma.exp((298.*delS-Hd)/(298.*R))) / \
                      (1. + np.ma.exp(((T+273.)*delS-Hd)/((T+273.)*R)))
                      
###############################################################################
def Topt(Hd, E, delS):
    '''
    leafmodel: temperature dependency
    Optimum temperature from entropy value after Medlyn et al. 2002
    Hd : decrease rate above optimum [J/mol]
    E: activation energy [J]
    delS: entropy factor [J/mol/K]
    
    returns:
    Topt: optimal temperature [degC]
    '''
    return Hd / (delS - R*np.ma.log(E/(Hd-E))) - 273.

###############################################################################
def sunfrac(G_of_theta, theta, lai):
    '''
    twoleafmodel: sunlit fraction of the canopy
    G_of_theta: leaf projection function at zenith angle theta [-]
    theta: zenith angle, 90. at horizon, 0. at noon [deg]
    lai: leaf area index [-]
    
    returns: sun lit fraction of lai [-]
    
    '''
    K = G_of_theta/np.ma.cos(np.deg2rad(theta))
    f_sun = (1.-np.ma.exp(-K*lai))/(K*lai) 
    return f_sun

###############################################################################
def gb(u, d):
    '''
    leafmodel: leaf boundayr layer conductance
    Leaf boundary layer conductance after Bonan et al. 2002
    u: wind speed [m/s]
    d: leaf size [m]
    
    returns: leaf boundary layer conductance [m/s]
    '''
    return 1./(200.*np.ma.sqrt(d/np.ma.abs(u)))

###############################################################################
def ga(u, ustar):
    '''
    leafmodel: aerodynamic conductance
    Aerodynamic conductance
    u: wind speed [m/s]
    ustar: friction velocity [m/s]
    
    returns: aerodynamic conductance [m/s]
    '''
    return ustar**2/u

###############################################################################
def g2mol(g, Ta, P):
    '''
    leafmodel: conductance unit conversion
    Converts conductivity, e.g. aerodynamic or leaf boundary layer conductivity,
    from [m/s] to [mol/m2s]
    g: conductivity [m/s]
    Ta: air temperature [degC]
    P: air pressure [Pa]
    
    returns: conductivity [mol/m2s]
    '''
    return g / (R * (Ta+T0)/P)

###############################################################################
def leaftemp(Ta, H, ga, rho, cp):
    '''
    leafmodel: leaf temperature estimation
    Leaf temperature
    Ta: air temperature [degC]
    H: sensible heat flux [W/m2]
    ga: aerodynamic conductance [m/s]
    rho: air density [kg/m3]
    cp: heat capacity of dry air [J/kg/K]
    
    returns: leaf temperature [degC]
    '''
    return Ta+(H/(ga*rho*cp))


if __name__ == '__main__':
    import doctest
    doctest.testmod()