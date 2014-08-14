#!/usr/bin/env python
from __future__ import print_function
import numpy as np

"""
    Defines common error measures.

    
    Definition
    ----------
    def bias(y_obs,y_mod):    bias
    def mae(y_obs,y_mod):     mean absolute error
    def mse(y_obs,y_mod):     mean squared error
    def rmse(y_obs,y_mod):    root mean squared error
    def nse(y_obs,y_mod):     Nash-Sutcliffe-Efficiency
    def pear2(y_obs,y_mod):   Squared Pearson correlation coefficient


    Input
    -----
    y_obs        np.array(N) or np.ma.array(N)
    y_mod        np.array(N) or np.ma.array(N)


    Output
    ------
    measure      float: error measure of respective error function
                        (see definitions for details)

    Restrictions
    ------------
    Deals with masked and unmasked arrays. When nan is found in an unmasked
    array, it will be masked. All measures are applied only on values where
    both, y_obs and y_mod, have valid entries (not masked and not nan)
    
    
    Examples
    --------
    -> see respective function


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

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Jul 2014
"""

def bias(y_obs,y_mod):
    """
    calculates bias = mean(y_obs) - mean(y_mod)
    
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate bias
    >>> print(np.round(bias(y_obs, y_mod),2))
    -0.43

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('bias: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('bias: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)
    return np.ma.mean(y_obsr) - np.ma.mean(y_modr)

def mae(y_obs,y_mod):
    """
    calculates mean absolute error = mean(abs(y_obs - y_mod))
        
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate mean absolute error
    >>> print(np.round(mae(y_obs, y_mod),2))
    0.55

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('mae: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('mae: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)
        return np.ma.mean(np.ma.abs(y_obsr-y_modr))

def mse(y_obs,y_mod):
    """
    calculates mean squared error = mean((y_obs - y_mod)**2)
        
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate mean squared error
    >>> print(np.round(mse(y_obs, y_mod),2))
    0.54

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('mse: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('mse: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)
        return np.ma.mean((y_obsr-y_modr)**2)

def rmse(y_obs,y_mod):
    """
    calculates root mean squared error = sqrt(mean((y_obs - y_mod)**2))
        
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate root mean squared error
    >>> print(np.round(rmse(y_obs, y_mod),2))
    0.73

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('rmse: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('rmse: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)
        return np.ma.sqrt(np.ma.mean((y_obsr-y_modr)**2))

def nse(y_obs,y_mod):
    """
    calculates Nash-Sutcliffe-Efficiency = 1 - (sum((y_obs - y_mod)**2) / sum((y_obs - mean(y_obs))**2))
       
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate Nash-Sutcliffe-Efficiency
    >>> print(np.round(nse(y_obs, y_mod),2))
    0.71

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('r2: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('r2: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)
        a = np.ma.sum((y_obsr - y_modr)**2)
        b = np.ma.sum((y_obsr - np.ma.mean(y_obsr))**2)
    return 1. - (a / b)

def pear2(y_obs,y_mod):
    """
    calculates squared Pearson correlation coeffcient
     
    Examples
    --------
    >>> # Create some data
    >>> y_obs = np.array([12.7867, 13.465, 14.1433, 15.3733, 16.6033])
    >>> y_mod = np.array([12.8087, 13.151, 14.3741, 16.2302, 17.9433])
    >>> # calculate Squared Pearson correlation coefficient
    >>> print(np.round(pear2(y_obs, y_mod),2))
    0.99

    """
    # check
    if (y_obs.ndim!=1) or (y_mod.ndim!=1):
        raise ValueError('pear2: input must be 1D')
    elif y_obs.size!=y_mod.size:
        raise ValueError('pear2: input must be of same size')
    # calc
    else:
        # check if masked or not
        try:
            temp = y_obs.mask
            temp = y_mod.mask
        except AttributeError:
            y_obs=np.ma.array(y_obs, mask=np.isnan(y_obs))
            y_mod=np.ma.array(y_mod, mask=np.isnan(y_mod))
        y_modr = np.ma.array(y_mod, mask=y_mod.mask | y_obs.mask)    
        y_obsr = np.ma.array(y_obs, mask=y_mod.mask | y_obs.mask)    
        return np.corrcoef(y_obsr.compressed(), y_modr.compressed())[0,1]**2

if __name__ == '__main__':
    import doctest
    doctest.testmod()