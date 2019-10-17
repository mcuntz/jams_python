#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def readhdf5(fName, var='', reform=False, squeeze=False, variables=False,
             attributes=False, fileattributes=False, sort=False):
    """
        Get variables or print information of hdf5 file.


        Definition
        ----------
        def readhdf5(fName, var='', reform=False, squeeze=False, variables=False,
                     attributes=False, fileattributes=False, sort=False):

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
        >>> print(a['NB_PARAMETERS'])
        9

        >>> print([str(i) for i in readhdf5('test_readhdf5.hdf5', variables=True)])
        ['chs']

        >>> a = readhdf5('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> print([ str(i) for i in sorted(a)])
        ['Double', 'Inttest', 'LLLLLL', 'What the hell']
        >>> print(a['Double'])
        [1.1]

        >>> from autostring import astr
        >>> print(astr(readhdf5('test_readhdf5.hdf5', var='chs'),3,pp=True))
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

        Copyright (c) 2012-2013 Matthias Zink, Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MZ, Jun 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Oct 2013 - hdf5read
    """
    try:
        import h5py as hdf5
    except:
        raise Error('No HDF5 support available, i.e. h5py')
    # Open hdf5 file
    try:
        f = hdf5.File(fName, 'r')
    except IOError:
        raise IOError('Cannot open file: '+fName)
    # Get attributes of the file
    if fileattributes:
        attrs = dict()
        attr = list(f.attrs.keys())
        for a in attr:
            attrs[a] = f.attrs.get(a)
        f.close()
        return attrs
    # Variables
    vars = list(f.keys())
    nvars = len(vars)
    # Sort and get sort indices
    if sort:
      svars = sorted(vars)
      ivars = list()
      for v in svars:
        ivars.append(vars.index(v))
    if variables:
      f.close()
      if sort:
        return svars
      else:
        return vars
    # Get attributes of variables
    if attributes:
      if var not in vars:
          f.close()
          raise ValueError('variable '+var+' not in file '+fname)
      attrs = dict()
      attr = list(f[var].attrs.keys())
      for a in attr:
        attrs[a] = f[var].attrs.get(a)
      f.close()
      return attrs
    # Get variable
    if var == '':
        f.close()
        raise ValueError('Variable name has to be given.')
    if var != '':
      if var not in vars:
          f.close()
          raise ValueError('variable '+var+' not in file '+fname)
      try:
          arr = f[var][:]
      except IOError:
          f.close()
          raise IOError('Cannot read variable '+var+' in file '+fname)
      if reform or squeeze:
        f.close()
        return arr.squeeze()
      else:
        f.close()
        return arr


def hdf5read(*args, **kwargs):
    """
        Wrapper for readhdf5.
        def readhdf5(fName, var='', reform=False, squeeze=False, variables=False,
                     attributes=False, fileattributes=False, sort=False):


        Examples
        --------
        >>> a = hdf5read('test_readhdf5.hdf5', fileattributes=True)
        >>> print(a['NB_PARAMETERS'])
        9

        >>> print([str(i) for i in hdf5read('test_readhdf5.hdf5', variables=True)])
        ['chs']

        >>> a = hdf5read('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> print([ str(i) for i in sorted(a)])
        ['Double', 'Inttest', 'LLLLLL', 'What the hell']
        >>> print(a['Double'])
        [1.1]

        >>> from autostring import astr
        >>> print(astr(hdf5read('test_readhdf5.hdf5', var='chs'),3,pp=True))
        [['  1.000' '  2.000' '  3.000' '  3.000' '  2.000']
         ['  1.000' ' 23.000' '254.000' '  5.000' '  4.654']
         ['654.654' ' 54.540' '546.540' '564.546' '  5.500']
         ['  1.100' '  2.200' '  0.000' '  3.300' '  4.400']
         ['  0.000' '  1.000' '  2.000' '  3.000' '  0.000']]
    """
    return readhdf5(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    try:
        import h5py as hdf5
    except ImportError:
        raise ImportError('No HDF5 support available, i.e. h5py')
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # var = readhdf5('test_readhdf5.hdf5', variables=True)
    # print var

    # varAtt = readhdf5('test_readhdf5.hdf5', var='EM_BB', attributes=True)
    # print varAtt
    # #{u'CAL_SLOPE': 999.0, u'PRODUCT': 'EM_BB', u'PRODUCT_ID': 999, u'CAL_OFFSET': 999.0, u'N_LINES': 651, u'N_COLS': 1701, u'SCALING_FACTOR': 10000.0, u'NB_BYTES': 2, u'OFFSET': 0.0, u'UNITS': 'Dimensionless', u'MISSING_VALUE': -8000, u'CLASS': 'Data'}
    # varAtt['SCALING_FACTOR']
    # #1000.0
    # data = readhdf5('test_readhdf5.hdf5', var='EM_BB')
    # print data
    # #array([[-8000, -8000, -8000, ..., -8000, -8000, -8000],
    # #       [-8000, -8000, -8000, ..., -8000, -8000, -8000],
    # #       [-8000, -8000, -8000, ..., -8000, -8000, -8000],
    # #       ...,
    # #        [-8000, -8000, -8000, ..., -8000, -8000, -8000],
    # #       [-8000, -8000, -8000, ..., -8000, -8000, -8000],
    # #       [-8000, -8000, -8000, ..., -8000, -8000, -8000]], dtype=int16)

