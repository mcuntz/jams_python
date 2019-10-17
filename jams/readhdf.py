#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import os

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
        # Do only HDF5 tests. HDF4 tests are in wrapper function hdfread
        # so that al least the HDF5 doctests can be run for Python 3.
        # HDF5
        >>> a = readhdf('test_readhdf5.hdf5', fileattributes=True)
        >>> print(a['NB_PARAMETERS'])
        9

        >>> print([ str(i) for i in readhdf('test_readhdf5.hdf5', variables=True) ])
        ['chs']

        >>> a = readhdf('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> print([ str(i) for i in sorted(list(a.keys())) ])
        ['Double', 'Inttest', 'LLLLLL', 'What the hell']
        >>> print(a['Double'])
        [1.1]

        >>> from autostring import astr
        >>> print(astr(readhdf('test_readhdf5.hdf5', var='chs'),3,pp=True))
        [['  1.000' '  2.000' '  3.000' '  3.000' '  2.000']
         ['  1.000' ' 23.000' '254.000' '  5.000' '  4.654']
         ['654.654' ' 54.540' '546.540' '564.546' '  5.500']
         ['  1.100' '  2.200' '  0.000' '  3.300' '  4.400']
         ['  0.000' '  1.000' '  2.000' '  3.000' '  0.000']]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2013 Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written,  MC, Jun 2012
        Modified, MC, Feb 2013 - ported to Python 3
    """
    # Open hdf5 file
    if not os.path.isfile(fName):
        raise IOError('File does not exist: '+fName)

    try:
        from pyhdf.SD import SD
        haveh4 = True
        from readhdf4 import readhdf4
    except ImportError:
        haveh4 = False

    try:
        import h5py as hdf5
        haveh5 = True
        from readhdf5 import readhdf5
    except ImportError:
        haveh5 = False

    if haveh5:
        try:
            f = hdf5.File(fName, 'r')
            f.close()
            return readhdf5(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                            attributes=attributes, fileattributes=fileattributes, sort=sort)
        except:
            if haveh4:
                try:
                    f = SD(fName)
                    f.end()
                    return readhdf4(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                                attributes=attributes, fileattributes=fileattributes, sort=sort)
                except:
                    raise IOError('File is not in hdf4 nor in hdf5 format: '+fName)
            else:
                raise IOError('File is not in hdf5 format: '+fName)
    else:
        if haveh4:
            try:
                f = SD(fName)
                f.end()
                return readhdf4(fName, var=var, reform=reform, squeeze=squeeze, variables=variables,
                                attributes=attributes, fileattributes=fileattributes, sort=sort)
            except:
                raise IOError('File is not in hdf4 format: '+fName)



def hdfread(*args, **kwargs):
    """
        Wrapper for readhdf.
        def readhdf(fName, var='', reform=False, squeeze=False, variables=False,
                    attributes=False, fileattributes=False, sort=False):


        Examples
        --------
        # HDF4
        >>> var = hdfread('test_readhdf4.hdf4', fileattributes=True)
        >>> print(list(var.keys()))
        ['OldCoreMetadata.0', 'HDFEOSVersion', 'OldArchiveMetadata.0', 'OldStructMetadata.0', 'StructMetadata.0']
        >>> print(var['HDFEOSVersion'])
        ('HDFEOS_V2.14', 0, 4, 12)

        >>> var = hdfread('test_readhdf4.hdf4', variables=True)
        >>> print(var)
        ['QC_250m_1', 'sur_refl_b02_1', 'sur_refl_b01_1', 'num_observations']

        >>> var = hdfread('test_readhdf4.hdf4', variables=True, sort=True)
        >>> print(var)
        ['QC_250m_1', 'num_observations', 'sur_refl_b01_1', 'sur_refl_b02_1']

        >>> var = hdfread('test_readhdf4.hdf4', var='sur_refl_b01_1')
        >>> print(var)
        [[7492 7327 7327 7131 7187]
         [6604 6604 7423 7131 7131]
         [7441 7441 7423 7423 7507]]

        >>> var = hdfread('test_readhdf4.hdf4', var='sur_refl_b01_1', attributes=True)
        >>> print(list(var.keys()))
        ['_FillValue', 'Nadir Data Resolution', 'scale_factor', 'valid_range', 'add_offset', 'long_name', 'calibrated_nt', 'units', 'scale_factor_err', 'add_offset_err', 'HorizontalDatumName']
        >>> print(var['_FillValue'])
        (-28672, 3, 22, 1)
    """
    return readhdf(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    try:
        from pyhdf.SD import SD
        doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    except ImportError:
        doctest.run_docstring_examples(readhdf, None, optionflags=doctest.NORMALIZE_WHITESPACE)

    # print readhdf('test_readhdf4.hdf4', variables=True)
    # print readhdf('test_readhdf5.hdf5', variables=True)

