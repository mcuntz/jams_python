#!/usr/bin/env python
"""
    Test functions for parameter sensitivity analysis from
        Ishigami and Homma (1990) An importance qualification technique in uncertainty analysis for computer models,
                                  Proceedings of the isuma '90, First International Symposium on Uncertainty
                                  Modelling and Analysis, University of Maryland, Dec. 03 - Dec 05 1990, 398-403
        Oakley and O'Hagan (2004) Probabilistic sensitivity analysis of complex models: a Bayesian approach
                                  J. R. Statist. Soc. B 66, Part 3, 751-769.
        Morris (1991)             Factorial sampling plans for preliminary computational experiments,
                                  Technometrics 33, 161-174.
        Saltelli et al. (2008)    Global Sensitivity Analysis. The Primer, John Wiley & Sons, pp. 292
        Saltelli et al. (2010)    Variance based sensitivity analysis of model output, Design and estimator
                                  for the total sensitivity index, Comp. Phys. Comm. 181, 259-270.
        Sobol' (1990),            Sensitivity estimates for nonlinear mathematical models,
                                  Matematicheskoe Modelirovanie 2, 112-118 (in Russian),
                                  translated in English in Sobol' (1993).
        Sobol' (1993)             Sensitivity analysis for non-linear mathematical models,
                                  Mathematical Modelling and Computational Experiment 1, 407-414,
                                  English translation of Russian original paper Sobol' (1990).


    Current functions are:
    B                     B of Saltelli et al. (2010)
    G / g                 G-function attributed to Sobol' (1990, 1993), given by Saltelli et al. (2008, 2010)
    Gstar                 G* of Saltelli et al. (2010)
    ishigami_homma        Ishigami and Homma (1990), given by Saltelli et al. (2008, page 179)
    K                     K  of Saltelli et al. (2010)
    morris                After Morris (1991)
    oakley_ohagan         Oakley and O'Hagan (2004), parameters given in Saltelli et al. (2008)
                          or on http://www.jeremy-oakley.staff.shef.ac.uk/psa_example.txt


    Input / Output
    --------------
    See the help of the individual functions for explanations of in/out, etc.


    Examples
    --------
    >>> from ufz.autostring import astr

    !!! ToDo !!!
    print(astr(griewank([0,0]),3,pp=True))
    0.000
    print(astr(goldstein_price([0,-1]),3,pp=True))
    3.000


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

    Copyright 2012-2015 Matthias Cuntz


    History
    -------
    Written,  MC, Mar 2015
"""
from __future__ import print_function

__all__ = ['B', 'g', 'G', 'Gstar', 'K', 'morris', 'oakley_ohagan', 'ishigami_homma']

import numpy as np

# -----------------------------------------------------------

def B(X):
    '''
        B function, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

        
        Input
        -----
        X        (nX,) or (nX,npoints) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    '''
    # Parameter sets are assumed to be in following ordering: (x_1, x_2, ..., X_m, w_1, w_2, ..., w_m)
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] % 2 == 0, 'X.shape[0] must be even.'

    m  = iX.shape[0]//2
    y  = np.sum(iX[:m,:]*iX[m:,:], axis=0)

    if isone:
        return y[0]
    else:
        return y

# -----------------------------------------------------------

def g(X, a):
    """
        G-function: Sobol' (1990) Matematicheskoe Modelirovanie 2, 112-118 (in Russian)
                    Sobol' (1993) Mathematical Modelling and Computational Experiment 1, 407-414 (English translation)

        
        Input
        -----
        X        (nX,) or (nX,npoints) array of floats
        a        (nX,) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    return Gstar(X, np.ones(len(a)), np.zeros(len(a)), a)

# -----------------------------------------------------------

def G(X, a):
    """
        G-function: Sobol' (1990) Matematicheskoe Modelirovanie 2, 112-118 (in Russian)
                    Sobol' (1993) Mathematical Modelling and Computational Experiment 1, 407-414 (English translation)

        
        Input
        -----
        X        (nX,) or (nX,npoints) array of floats
        a        (nX,) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    return Gstar(X, np.ones(len(a)), np.zeros(len(a)), a)

# -----------------------------------------------------------

def Gstar(X, alpha, delta, a):
    """
        G* example, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

        
        Input
        -----
        X        (nX,) or (nX,npoints) array of floats
        alpha    (nX,) array of floats
        delta    (nX,) array of floats
        a        (nX,) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 10
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    if isinstance(alpha, (list, tuple)):
        nalpha = len(alpha)
    else:
        nalpha = alpha.size

    assert iX.shape[0] == nalpha, 'X.shape[0] must len(alpha).'

    ialpha = np.array(alpha).reshape((nalpha,1))
    idelta = np.array(delta).reshape((nalpha,1))
    ia     = np.array(a).reshape((nalpha,1))
    yi = ((1.+ialpha)*np.abs(2.*(iX+idelta-np.trunc(iX+idelta))-1.)**ialpha+ia)/(1.+ia)
    y  = np.prod(yi, axis=0)

    if isone:
        return y[0]
    else:
        return y

# -----------------------------------------------------------

def ishigami_homma(X, a, b):
    """
        Ishigami and Homma (1990), given by Saltelli et al. (2008, page 179)

        
        Input
        -----
        X    (nX,) or (nX,npoints) array of floats
        a    (nX,) array of floats
        b    (nX,) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 3
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 3, 'X.shape[0] must 3.'

    y = np.sin(iX[0,:])  + a*(np.sin(iX[1,:]))**2  + b*iX[2,:]**4  * np.sin(iX[0,:])
    
    if isone:
        return y[0]
    else:
        return y

# -----------------------------------------------------------

def K(X):
    """
        K example, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

        
        Input
        -----
        X        (nX,) or (nX,npoints) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 10
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] > 1, 'X.shape[0] must be > 1.'

    nX = iX.shape[0]
    for ii in range(1,nX):
        iX[ii,:] = iX[ii-1,:] * iX[ii,:]
    for ii in range(nX):
        iX[ii,:] = (-1)**(ii+1) * iX[ii,:]
    y = iX.sum(axis=0)
    
    if isone:
        return y[0]
    else:
        return y

# -----------------------------------------------------------

def morris(X, beta0, beta1, beta2, beta3, beta4):
    """
        Morris-function: Morris (1991) Technometrics 33, 161-174

        
        Input
        -----
        X        (20,) or (20,npoints) array of floats
        beta0    float
        beta1    (20,) array of floats
        beta2    (20,20) array of floats
        beta3    (20,20,20) array of floats
        beta4    (20,20,20,20) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    om = 2.*(iX - 0.5)
    ii = np.array((2, 4, 6))
    om[ii,:] = 2.*(1.1*iX[ii,:]/(iX[ii,:]+0.1) - 0.5)

    nn  = 20

    assert iX.shape[0] == nn, 'X.shape[0] must '+str(nn)+'.'
    assert np.all(np.array(beta1.shape) == nn), 'beta1.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta2.shape) == nn), 'beta2.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta3.shape) == nn), 'beta3.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta4.shape) == nn), 'beta4.shape must be all '+str(nn)+'.'

    y    = np.empty(iX.shape[1])
    y[:] = beta0
    for i in range(0,nn):
        y[:] += beta1[i]*om[i,:]
        for j in range(i+1,nn):
            y[:] += beta2[i,j]*om[i,:]*om[j,:]
            for l in range(j+1,nn):
                y[:] += beta3[i,j,l]*om[i,:]*om[j,:]*om[l,:]
                for s in range(l+1,nn):
                    y[:] += beta4[i,j,l,s]*om[i,:]*om[j,:]*om[l,:]*om[s,:]

    if isone:
        return y[0]
    else:
        return y

# -----------------------------------------------------------
    
def oakley_ohagan(X):
    """
        Oakley and O'Hagan (2004) J. R. Statist. Soc. B 66, Part 3, 751-769

        
        Input
        -----
        X        (15,) or (15,npoints) array of floats


        Output
        ------
        float or (npoints,) array of floats

        
        History
        -------
        Written,  MC & JM, Mar 2015
    """
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    nn = 15

    assert iX.shape[0] == nn, 'X.shape[0] must '+str(nn)+'.'

    a1 = np.array([0.01, 0.05, 0.23, 0.04, 0.12, 0.39, 0.39, 0.61, 0.62, 0.40, 1.07, 1.15, 0.79, 1.12, 1.20])
    a2 = np.array([0.43, 0.09, 0.05, 0.32, 0.15, 1.04, 0.99, 0.97, 0.90, 0.81, 1.84, 2.47, 2.39, 2.00, 2.26])
    a3 = np.array([0.10, 0.21, 0.08, 0.27, 0.13, 0.75, 0.86, 1.03, 0.84, 0.80, 2.21, 2.04, 2.40, 2.05, 1.98])
    M  = np.array([
        [-0.02, -0.19,   0.13,  0.37,  0.17,  0.14, -0.44,  -0.08,  0.71, -0.44,  0.5,  -0.02, -0.05,  0.22,  0.06],
        [ 0.26,  0.05,   0.26,  0.24, -0.59, -0.08, -0.29,   0.42,  0.5,   0.08, -0.11,  0.03, -0.14, -0.03, -0.22],
        [-0.06,  0.2,    0.1,  -0.29, -0.14,  0.22,  0.15,   0.29,  0.23, -0.32, -0.29, -0.21,  0.43,  0.02,  0.04],
        [ 0.66,  0.43,   0.3,  -0.16, -0.31, -0.39,  0.18,   0.06,  0.17,  0.13, -0.35,  0.25, -0.02,  0.36, -0.33],
        [-0.12,  0.12,   0.11,  0.05, -0.22,  0.19, -0.07,   0.02, -0.1,   0.19,  0.33,  0.31, -0.08, -0.25,  0.37],
        [-0.28, -0.33,  -0.1,  -0.22, -0.14, -0.14, -0.12,   0.22, -0.03, -0.52,  0.02,  0.04,  0.36,  0.31,  0.05],
        [-0.08,  0.004,  0.89, -0.27, -0.08, -0.04, -0.19,  -0.36, -0.17,  0.09,  0.4,  -0.06,  0.14,  0.21, -0.01],
        [-0.09,  0.59,   0.03, -0.03, -0.24, -0.1,   0.03,   0.1,  -0.34,  0.01, -0.61,  0.08,  0.89,  0.14,  0.15],
        [-0.13,  0.53,   0.13,  0.05,  0.58,  0.37,  0.11,  -0.29, -0.57,  0.46, -0.09,  0.14, -0.39, -0.45, -0.15],
        [ 0.06, -0.32,   0.09,  0.07, -0.57,  0.53,  0.24,  -0.01,  0.07,  0.08, -0.13,  0.23,  0.14, -0.45, -0.56],
        [ 0.66,  0.35,   0.14,  0.52, -0.28, -0.16, -0.07,  -0.2,   0.07,  0.23, -0.04, -0.16,  0.22,  0,    -0.09],
        [ 0.32, -0.03,   0.13,  0.13,  0.05, -0.17,  0.18,   0.06, -0.18, -0.31, -0.25,  0.03, -0.43, -0.62, -0.03],
        [-0.29,  0.03,   0.03, -0.12,  0.03, -0.34, -0.41,   0.05, -0.27, -0.03,  0.41,  0.27,  0.16, -0.19,  0.02],
        [-0.24, -0.44,   0.01,  0.25,  0.07,  0.25,  0.17,   0.01,  0.25, -0.15, -0.08,  0.37, -0.3,   0.11, -0.76],
        [ 0.04, -0.26,   0.46, -0.36, -0.95, -0.17,  0.003,  0.05,  0.23,  0.38,  0.46, -0.19,  0.01,  0.17,  0.16] ])

    y = np.dot(a1,iX) + np.dot(a2,np.sin(iX)) + np.dot(a3,np.cos(iX))
    for i in range(iX.shape[1]):
        y[i] += np.dot(iX[:,i].T,np.dot(M,iX[:,i]))
    
    if isone:
        return y[0,0]
    else:
        return y.squeeze()

# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # print(griewank([0,0]))
    # #0.0
    # print(goldstein_price([0,-1]))
    # #3.0
