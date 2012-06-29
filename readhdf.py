#!/usr/bin/env python
import numpy as np
import os
from pyhdf.SD import SD
import h5py as hdf5
from readhdf4 import *
from readhdf5 import *

def readhdf(fName, var='', reform=False, squeeze=False, variables=False,
            attributes=False, fileattributes=False, sort=False):
    """
        Get variables or prints information of hdf5 file.


        Definition
        ----------
        def readhdf5(fName, var='', reform=False, squeeze=False, variables=False,
                     attributes=False, fileattributes=False,  sort=False):

        Input
        -----
        fName            hdf5 file name


        Optional Input Parameters
        -------------------------
        var              name of variable in hdf5 file


        Options
        -------
        reform           if output is array then squeeze(array)
        squeeze          same as reform
        variables        get list of variables in hdf5 file
        attributes       get dictionary of all attributes of specific variable
        fileattributes   get dictionary of all attributes of the file
        sort             sort variable names.
        quiet            quietly return None if error occurs


        Output
        ------
        Either variable array or information such as list of all variables in hdf5 file.


        Examples
        --------
        >>> a = readhdf5('test_readhdf5.hdf5', fileattributes=True)
        >>> print a['NB_PARAMETERS']
        9

        >>> print readhdf5('test_readhdf5.hdf5', variables=True)
        [u'chs']

        >>> a = readhdf5('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> for i in a: print i, a[i]
        Double [ 1.1]
        Inttest [99]
        What the hell ['']
        LLLLLL [528040]

        >>> print readhdf5('test_readhdf5.hdf5', var='chs')
        [[   1.            2.            3.            3.            2.        ]
         [   1.           23.          254.            5.            4.65399981]
         [ 654.6539917    54.54000092  546.53997803  564.54602051    5.5       ]
         [   1.10000002    2.20000005    0.            3.29999995    4.4000001 ]
         [   0.            1.            2.            3.            0.        ]]


        History
        -------
        Written, MZ, Jun 2012
        """
    # Open hdf5 file
    if not os.path.isfile(fName):
        raise IOError('Fiel does not exist: '+fName)
    try:
        f = hdf5.File(fName, 'r')
        f.close()
        return readhdf5(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                        attributes=attributes, fileattributes=fileattributes, sort=sort)
    except IOError:
        try:
            f = SD(fName)
            f.close()
            return readhdf4(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                            attributes=attributes, fileattributes=fileattributes, sort=sort)
        except HDF4Error:
            raise IOError('File is not in hdf4 or hdf5 format: '+fName)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # a = readhdf('test_readhdf5.hdf5', fileattributes=True)
    # print a['NB_PARAMETERS']
    # #    9

    # print readhdf('test_readhdf5.hdf5', variables=True)
    # #    [u'chs']

    # a = readhdf('test_readhdf5.hdf5', var='chs', attributes=True)
    # for i in a: print i, a[i]
    #     # Double [ 1.1]
    #     # Inttest [99]
    #     # What the hell ['']
    #     # LLLLLL [528040]

    # print readhdf('test_readhdf5.hdf5', var='chs')
    #     # [[   1.            2.            3.            3.            2.        ]
    #     #  [   1.           23.          254.            5.            4.65399981]
    #     #  [ 654.6539917    54.54000092  546.53997803  564.54602051    5.5       ]
    #     #  [   1.10000002    2.20000005    0.            3.29999995    4.4000001 ]
    #     #  [   0.            1.            2.            3.            0.        ]]
