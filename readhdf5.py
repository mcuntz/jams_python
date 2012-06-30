#!/usr/bin/env python
import numpy as np
import h5py as hdf

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
        >>> print a['NB_PARAMETERS']
        9

        >>> print readhdf5('test_readhdf5.hdf5', variables=True)
        [u'chs']

        >>> a = readhdf5('test_readhdf5.hdf5', var='chs', attributes=True)
        >>> print a.keys()
        [u'Double', u'Inttest', u'What the hell', u'LLLLLL']
        >>> print a['Double']
        [ 1.1]

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
    try:
        f = hdf.File(fName, 'r')
    except IOError:
        raise IOError('Cannot open file: '+fName)
    # Get attributes of the file
    if fileattributes:
        attrs = dict()
        attr = f.attrs.keys()
        for a in attr:
            attrs[a] = f.attrs.get(a)
        f.close()
        return attrs
    # Variables
    vars = f.keys()
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
      attr = f[var].attrs.keys()
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
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
