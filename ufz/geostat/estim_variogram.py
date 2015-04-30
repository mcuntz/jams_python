# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 19:41:36 2015

@author: fhesse
"""

import numpy as np
import random
import matplotlib.pyplot as plt

def estim_variogram( field_data ):
    
    sample_num = 10
    bin_num = 50 + 1
    
    lag_mean = np.zeros(( bin_num ))
    gamma_mean = np.zeros(( bin_num ))
    
    for sample_i in range(0, sample_num):
        
        gamma_i, lag_i = estim_single_variogram( field_data, bin_num )  
        lag_mean = lag_mean + lag_i/sample_num
        gamma_mean = gamma_mean + gamma_i/sample_num
    
#    plt.plot(lag_mean, gamma_mean)
#    plt.show( )
    
    return gamma_mean, lag_mean

def estim_single_variogram( field_data, bin_num ):

    x_vec = field_data[0, :]
    y_vec = field_data[1, :]
    u_vec = field_data[3, :]
    
    # sampling the data vector
    sample_num = 100
    pos_num = len( x_vec )
    sample_ind = [0]*sample_num
    for sample_i in range(0, sample_num):
        sample_ind[sample_i] = int( random.randint(0, pos_num) - 1 )

    x_sample = x_vec[sample_ind]
    y_sample = y_vec[sample_ind]
    u_sample = u_vec[sample_ind]
    
    #-- calculate empirical variogram cloud -----------------------------------
    
    index = 1
    lag = np.zeros(( sample_num*sample_num - sample_num + 1 ))
    gamma = np.zeros(( sample_num*sample_num - sample_num + 1 ))
    for i in range(0, sample_num - 1):
        for j in range(0, sample_num):
            
            delta_x = x_sample[j] - x_sample[i]
            delta_y = y_sample[j] - y_sample[i]
            lag[index] = np.sqrt( delta_x**2 + delta_y**2)
            gamma[index] = 0.5*abs( u_sample[j] - u_sample[i] )**2
            index = index + 1
    
    #-- calculate emperical variogram distribution ----------------------------

    bin_max = np.max( x_vec )
    bin_array = np.linspace(0, bin_max, bin_num)
    
    lag_sort = np.sort(lag)
    sort_index = np.argsort(lag)
    gamma_sort = gamma[sort_index]
    
    lag_index = np.digitize(lag_sort, bin_array)
    lag_mean = np.zeros(( bin_num ))
    gamma_mean = np.zeros(( bin_num ))
    
    for bin_i in range(1, bin_num+1):
        
        index_i = np.where( lag_index == bin_i )
        lag_mean[ bin_i-1 ] = np.mean(lag_sort[ index_i ])
        gamma_mean[ bin_i-1 ] = np.mean(gamma_sort[ index_i ])

#    print(np.where( lag_index == 34 ))
    
#    print( lag )
#    sample_pos = floor(rand(sample_num, 1)*pos_num) + 1
    
#    plt.plot(lag_mean, gamma_mean)
#    plt.show( )
    
#    print( container )
    
    return gamma_mean, lag_mean