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
        Get variables or print information of hdf4 or hdf5 file.


        Definition
        ----------
        def readhdf(fName, var='', reform=False, squeeze=False, variables=False,
                    attributes=False, fileattributes=False, sort=False):

        Input
        -----
        fName            hdf4/5 file name


        Optional Input Parameters
        -------------------------
        var              name of variable in hdf file


        Options
        -------
        reform           if output is array then squeeze(array)
        squeeze          same as reform
        variables        get list of variables in hdf file
        attributes       get dictionary of all attributes of specific variable
        fileattributes   get dictionary of all attributes of the file
        sort             sort variable names.


        Output
        ------
        Either variable array or information of file aor variable
        such as list of all variables or dictionary of variable attributes.


        Examples
        --------
        >>> # HDF4
        >>> var = readhdf('test_readhdf4.hdf4', fileattributes=True)
        >>> print var.keys()
        ['OldCoreMetadata.0', 'HDFEOSVersion', 'OldArchiveMetadata.0', 'OldStructMetadata.0', 'StructMetadata.0']
        >>> print var['HDFEOSVersion']
        ('HDFEOS_V2.14', 0, 4, 12)

        >>> var = readhdf('test_readhdf4.hdf4', variables=True)
        >>> print var
        ['QC_250m_1', 'sur_refl_b02_1', 'sur_refl_b01_1', 'num_observations']

        >>> var = readhdf('test_readhdf4.hdf4', variables=True, sort=True)
        >>> print var
        ['QC_250m_1', 'num_observations', 'sur_refl_b01_1', 'sur_refl_b02_1']
    
        >>> var = readhdf('test_readhdf4.hdf4', var='sur_refl_b01_1')
        >>> print var
        [[7492 7327 7327 7131 7187]
         [6604 6604 7423 7131 7131]
         [7441 7441 7423 7423 7507]]

        >>> var = readhdf('test_readhdf4.hdf4', var='sur_refl_b01_1', attributes=True)
        >>> print var.keys()
        ['_FillValue', 'Nadir Data Resolution', 'scale_factor', 'valid_range', 'add_offset', 'long_name', 'calibrated_nt', 'units', 'scale_factor_err', 'add_offset_err', 'HorizontalDatumName']
        >>> print var['_FillValue']
        (-28672, 3, 22, 1)


        >>> # HDF5
        >>> a = readhdf('test_readhdf5.hdf5', fileattributes=True)
        >>> print a['NB_PARAMETERS']
        9

        >>> print readhdf('test_readhdf5.hdf5', variables=True)
        [u'chs']

        >>> a = readhdf('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> print a.keys()
        [u'Double', u'Inttest', u'What the hell', u'LLLLLL']
        >>> print a['Double']
        [ 1.1]

        >>> print readhdf('test_readhdf5.hdf5', var='chs')
        [[   1.            2.            3.            3.            2.        ]
         [   1.           23.          254.            5.            4.65399981]
         [ 654.6539917    54.54000092  546.53997803  564.54602051    5.5       ]
         [   1.10000002    2.20000005    0.            3.29999995    4.4000001 ]
         [   0.            1.            2.            3.            0.        ]]


        History
        -------
        Written, MC, Jun 2012
        """
    # Open hdf5 file
    if not os.path.isfile(fName):
        raise IOError('File does not exist: '+fName)
    try:
        f = hdf5.File(fName, 'r')
        f.close()
        return readhdf5(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                        attributes=attributes, fileattributes=fileattributes, sort=sort)
    except IOError:
        try:
            f = SD(fName)
            f.end()
            return readhdf4(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                            attributes=attributes, fileattributes=fileattributes, sort=sort)
        except HDF4Error:
            raise IOError('File is not in hdf4 or hdf5 format: '+fName)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # print readhdf('test_readhdf4.hdf4', variables=True)
    # print readhdf('test_readhdf5.hdf5', variables=True)
