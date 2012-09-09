#!/usr/bin/env python
import numpy as np
from pyhdf.SD import SD

def readhdf4(fName, var='', reform=False, squeeze=False, variables=False,
            attributes=False, fileattributes=False, sort=False):
    """
        Get variables or print information of hdf4 file.


        Definition
        ----------
        def readhdf4(fName, var='', reform=False, squeeze=False, variables=False,
                     attributes=False, fileattributes=False, sort=False):

        Input
        -----
        fName            hdf4 file name


        Optional Input Parameters
        -------------------------
        var              name of variable in hdf4 file


        Options
        -------
        reform           if output is array then squeeze(array)
        squeeze          same as reform
        variables        get list of variables in hdf4 file
        attributes       get dictionary of all attributes of specific variable
        fileattributes   get dictionary of all attributes of the file
        sort             sort variable names


        Output
        ------
        Either variable array or information of file or variable
        such as list of all variables or attributes of a variable.


        Examples
        --------
        >>> from pyhdf.SD import SD
        >>> var = readhdf4('test_readhdf4.hdf4', fileattributes=True)
        >>> print var.keys()
        ['OldCoreMetadata.0', 'HDFEOSVersion', 'OldArchiveMetadata.0', 'OldStructMetadata.0', 'StructMetadata.0']
        >>> print var['HDFEOSVersion']
        ('HDFEOS_V2.14', 0, 4, 12)

        >>> var = readhdf4('test_readhdf4.hdf4', variables=True)
        >>> print var
        ['QC_250m_1', 'sur_refl_b02_1', 'sur_refl_b01_1', 'num_observations']

        >>> var = readhdf4('test_readhdf4.hdf4', variables=True, sort=True)
        >>> print var
        ['QC_250m_1', 'num_observations', 'sur_refl_b01_1', 'sur_refl_b02_1']
    
        >>> var = readhdf4('test_readhdf4.hdf4', var='sur_refl_b01_1')
        >>> print var
        [[7492 7327 7327 7131 7187]
         [6604 6604 7423 7131 7131]
         [7441 7441 7423 7423 7507]]

        >>> var = readhdf4('test_readhdf4.hdf4', var='sur_refl_b01_1', attributes=True)
        >>> print var.keys()
        ['_FillValue', 'Nadir Data Resolution', 'scale_factor', 'valid_range', 'add_offset', 'long_name', 'calibrated_nt', 'units', 'scale_factor_err', 'add_offset_err', 'HorizontalDatumName']
        >>> print var['_FillValue']
        (-28672, 3, 22, 1)


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
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012 Matthias Cuntz


        History
        -------
        Written, MC, Jun 2012
        """
    # Open hdf4 file
    try:
        f = SD(fName)
    except HDF4Error:
        raise IOError('Cannot open file: '+fName)
    # Get attributes of the file
    if fileattributes:
        attr = f.attributes(full=1)
        f.end()
        return attr
    # Variables
    svars = f.datasets().keys()
    # Sort and get sort indices
    if variables:
      f.end()
      if sort:
        svars.sort()
        return svars
      else:
        return svars
    # Get attributes of variables
    if attributes:
      if var not in svars:
          f.end()
          raise ValueError('Variable '+var+' not in file '+fname)
      attrs = f.select(var).attributes(full=1)
      f.end()
      return attrs
    # Get variable
    if var == '':
        f.end()
        raise ValueError('Variable name has to be given.')
    if var != '':
      if var not in svars:
          f.end()
          raise ValueError('Variable '+var+' not in file '+fname)
      try:
          arr = np.array(f.select(var).get())
      except ValueError:
          f.end()
          raise IOError('Cannot read variable '+var+' in file '+fName)
      f.end()
      if reform or squeeze:
        return arr.squeeze()
      else:
        return arr


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # var = readhdf4('test_readhdf4.hdf4', fileattributes=True)
    # print var.keys()
    # # ['OldCoreMetadata.0',
    # #  'HDFEOSVersion',
    # #  'OldArchiveMetadata.0',
    # #  'OldStructMetadata.0',
    # #  'StructMetadata.0']
    # print var['HDFEOSVersion']
    # # ('HDFEOS_V2.14', 0, 4, 12)

    # var = readhdf4('test_readhdf4.hdf4', variables=True)
    # print var
    # # ['QC_250m_1', 'sur_refl_b02_1', 'sur_refl_b01_1', 'num_observations']

    # var = readhdf4('test_readhdf4.hdf4', variables=True, sort=True)
    # print var
    # # ['QC_250m_1', 'num_observations', 'sur_refl_b01_1', 'sur_refl_b02_1']
    
    # var = readhdf4('test_readhdf4.hdf4', var='sur_refl_b01_1')
    # print var
    # # [[7492 7327 7327 7131 7187]
    # #  [6604 6604 7423 7131 7131]
    # #  [7441 7441 7423 7423 7507]]

    # var = readhdf4('test_readhdf4.hdf4', var='sur_refl_b01_1', attributes=True)
    # print var.keys()
    # # ['_FillValue',
    # #  'Nadir Data Resolution',
    # #  'scale_factor',
    # #  'valid_range',
    # #  'add_offset',
    # #  'long_name',
    # #  'calibrated_nt',
    # #  'units',
    # #  'scale_factor_err',
    # #  'add_offset_err',
    # #  'HorizontalDatumName']
    # print var['_FillValue']
    # # (-28672, 3, 22, 1)