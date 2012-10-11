#!/usr/bin/env python
import numpy as np

def jab(arr_stat=None, arr_ind=None, nB=None):
    """
        Computes the Jackknife-after-Bootstrap statistics which is the standard error of the
        Bootstrap standard error, depending on the number of bootstrap (B)
        The statistic should converge to a very small number (error) with B --> inf. 

        Definition
        ----------
        def jab():


        Optional Input
        --------------
        arr_stat            array that holds in the columns the statistic of interest per Bootstrap and in the rows the number of parameters 
        
        arr_ind             array that holds the bootstrap indices, the number of rows corresponds to the number of bootstraps, the number of rows to the number of data points from which the statistic is computed
       
        nB                  number of Bootstraps 


        Output
        ------
        standard error of bootstrap standard error 


        Restrictions
        ------------
        ATTENTION: be careful with the dimensions of the arrays


        References
        ----------
        Jackknife-after-Bootstrap standard errors and influence functions, Efron B, 1990
        Weighted Jackknife-after-bootstrap: A heuristic approach, Wang J, 1997

        Example
        --------
        *computation of JAB error for the means for two parameters(Bstat columns)
        for five bootstraps (Bstat rows)
        *indices per bootstrap (indices colums) for five bootstraps (indices rows)
        *eg. first bootstrap with indices [0, 0, 4, 5, 5, 3, 2, 6] correponds to means
        [ 3.625,  5.375] for parameter 1 and 2 respectively
        *[ 0.36657688  0.59590163] are the corresponding errors (similar to standard deviation) of the
        errors (standard deviation) of the bootstrap samples (over all five boots)
         >>> import numpy as np
         >>> Bstat   = np.array([[ 3.625,  5.375],\
                               [ 3.625,  4.125],\
                               [ 4.625,  3.75 ],\
                               [ 4.   ,  4.5  ],\
                               [ 4.25 ,  3.625]])
         >>> indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],\
                                 [1, 2, 3, 0, 5, 5, 4, 1],\
                                 [3, 0, 3, 5, 4, 4, 1, 7],\
                                 [0, 3, 5, 2, 2, 1, 4, 4],\
                                 [4, 4, 5, 5, 3, 3, 2, 1]])
         >>> nB      = 5                
         >>> tmp     = jab(Bstat,indices,nB) 
         >>> print tmp
         [ 0.36657688  0.59590163]

        
        History
        -------
        Written, MG, August 2012
        
    """

    # Check input
    if (arr_ind==None):
        raise ValueError('No array of bootstrap indices given.')
    if (arr_stat==None):
        raise ValueError('No array of statistic given.')
    if (nB==None):
        raise ValueError("Number of B not given!")

    # errors in bootstrap samples
    se_boot     = np.empty([arr_ind.shape[1],arr_stat.shape[1]])

    # convert indices to integer for mask
    arr_ind = arr_ind.astype(int)
    #how much samples do we have(number of lines from indices matrix)
    n       = arr_ind.shape[1]
    n_float = float(n)


    # boolean array to mark the indices in the boot samples
    # mark in the rows, which index was chosen (true)
    #mask      = np.empty([n,np.shape(arr_ind)[1]],dtype=bool)
    mask      = np.empty([np.shape(arr_ind)[1],np.shape(arr_ind)[0]],dtype=bool)
    mask[:,:] = False

    for i in xrange(arr_ind.shape[0]):
        mask[arr_ind[i,:],i]=True

    seboot_i = np.empty([arr_ind.shape[1],arr_stat.shape[1]])
    for k in xrange(n):
        # indexes of bootstraps of kth line that doesn't exist
        # use that samples to estimate error se_Bi
        ind = np.where(mask[k,:]==False)[0]
        # stats of bootstraps where kth line didn't appear
        tmp = np.empty([ind.shape[0],arr_stat.shape[1]])
        tmp[:,:] = arr_stat[ind,:]        
        # compute std of these missing samples
        seboot_i[k,:] = np.std(tmp[:,:],axis=0)

    # initialize array for jab (len of data stat vector)
    se_jab_boot = np.empty([arr_stat.shape[1]])
    for i in xrange(arr_stat.shape[1]):
        se_jab_boot[i] = np.sqrt((n_float-1.)/n_float * \
                         np.sum((seboot_i[:,i] - np.mean(seboot_i[:,i]))**2) )

    return se_jab_boot

    

"""
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    Bstat = np.array([[ 3.625,  5.375],
                      [ 3.625,  4.125],
                      [ 4.625,  3.75 ],
                      [ 4.   ,  4.5  ],
                      [ 4.25 ,  3.625]])
    indices = np.array([[0, 0, 4, 5, 5, 3, 2, 6],
                        [1, 2, 3, 0, 5, 5, 4, 1],
                        [3, 0, 3, 5, 4, 4, 1, 7],
                        [0, 3, 5, 2, 2, 1, 4, 4],
                        [4, 4, 5, 5, 3, 3, 2, 1]])
    nB      = 5                
    tmp     = jab(Bstat,indices,nB)
    print tmp
"""

