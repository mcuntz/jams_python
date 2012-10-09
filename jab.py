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
        arr_stat array that holds the statistic of interest per Bootstrap
        
        arr_ind  array that holds the bootstrap indices, dim.: nB * indices
       
        nB       number of Bootstraps 


        Output
        ------
        standard error of bootstrap error 


        Restrictions
        ------------



        References
        ----------
        Jackknife-after-Bootstrap standard errors and influence functions, Efron B, 1990
        Weighted Jackknife-after-bootstrap: A heuristic approach, Wang J, 1997

        Examples
        --------
         >>> import numpy as np
         >>> arr_stat = np.array([[ 0.24132883,  0.09431899,  0.09271462,  0.1022947 ,  0.12487993],\
                                  [ 0.06426133,  0.04065362,  0.07352835,  0.05208864,  0.07569657], \
                                  [ 0.7796398 ,  0.86108454,  0.94973524,  0.79139866,  0.82595317],\
                                  [ 0.16941597,  0.12972866,  0.01450686,  0.15166605,  0.17000265],\
                                  [ 0.3450626 ,  0.24887742,  0.04347608,  0.45702632,  0.30438783]])
         >>> arr_ind  = np.array([[ 4.,  0.,  1.,  0.,  4.],\
                                  [ 0.,  0.,  1.,  3.,  1.],\
                                  [ 3.,  4.,  2.,  0.,  0.],\
                                  [ 0.,  2.,  2.,  4.,  3.],\
                                  [ 4.,  3.,  4.,  3.,  0.]])
         >>> nB       = 5                
         >>> tmp      = jab(arr_stat,arr_ind,nB)
         >>> print tmp
         [ 0.06304995  0.00944466  0.02912111  0.01306797  0.07442766]
        
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

    # initialize array for jab  ([number of paras,number of bootstraps])
    se_jab_boot = np.empty([arr_stat.shape[0]])

    ind         = np.ones([arr_ind.shape[1]])
    se_boot     = np.empty([arr_ind.shape[0],arr_stat.shape[0]])

    arr_ind = arr_ind.astype(int)
    # number of lines from indices matrix
    n       = arr_ind.shape[0]
    n_float = float(n)


    # boolean array to mark the indices in the boot samples
    # mark in the rows, which index was chosen (true)
    a      = np.empty([n,np.shape(arr_ind)[1]],dtype=bool)
    a[:,:] = False

    for i in xrange(arr_ind.shape[1]):
        a[arr_ind[:,i],i]=True

    seboot_i = np.zeros([n,arr_stat.shape[0]])
    for k in xrange(n):
        # indexes of bootstraps of kth line that doesn't exist
        # use that samples to estimate error se_Bi
        ind = np.where(a[k,:]==False)[0]
        # stats of bootstraps where kth line didn't appear

        tmp           = np.empty([arr_stat.shape[0],ind.shape[0]])
        tmp           = arr_stat[:,ind]
        # compute std of these missing samples
        seboot_i[k,:] = (np.std(tmp[:,:],axis=1))
        
    for kk in xrange(arr_stat.shape[0]):
        se_jab_boot[kk] = np.sqrt((n_float-1.)/n_float * \
                         np.sum((seboot_i[:,kk] - np.mean(seboot_i[:,kk]))**2) )




    return se_jab_boot 


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    #arr_stat = np.array([[ 0.24132883,  0.09431899,  0.09271462,  0.1022947 ,  0.12487993],\
    #                     [ 0.06426133,  0.04065362,  0.07352835,  0.05208864,  0.07569657], \
    #                     [ 0.7796398 ,  0.86108454,  0.94973524,  0.79139866,  0.82595317],\
    #                     [ 0.16941597,  0.12972866,  0.01450686,  0.15166605,  0.17000265],\
    #                     [ 0.3450626 ,  0.24887742,  0.04347608,  0.45702632,  0.30438783]])
    #arr_ind  = np.array([[ 4.,  0.,  1.,  0.,  4.],\
    #                     [ 0.,  0.,  1.,  3.,  1.],\
    #                     [ 3.,  4.,  2.,  0.,  0.],\
    #                     [ 0.,  2.,  2.,  4.,  3.],\
    #                     [ 4.,  3.,  4.,  3.,  0.]])
    #arr_B    = 5                
    #tmp      = jab(arr_stat,arr_ind,nB)
    #print tmp


