#!/usr/bin/env python
from __future__ import print_function
import numpy as np                       # array manipulation
import netCDF4 as nc
from jams.readnetcdf import readnetcdf

def writenetcdf(fhandle, vhandle=None, var=None, time=None, isdim=False, name=None, dims=None,
                attributes=None, fileattributes=None, comp=False, vartype=None, create_var=True ):
    """
        Writes dimensions, variables, dependencies and attributes to NetCDF file
        for 1D to 6D data.

        All attribues must be lists; except the data itself.


        Definition
        ----------
        def writenetcdf(fhandle, vhandle=None, var=None, time=None, isdim=False, name=None, dims=None,
                        attributes=None, fileattributes=None, comp=False):


        Input           Format                  Description
        -----           -----                   -----------
        fhandle         string                  file handle of nc.Dataset(FileName, 'w')
        vhandle         string                  variable handle of the particular variable
        var             array like              data (assumed to be netcdf float4)
        time            integer or array like   particular time step of the data
        isdim           boolean                 defines if current var is a dimension
        name            string                  name of the variable
        dims            1D list                 variable dependencies (e.g. [time, x, y])
        attributes      2D list or dictionary   variable attributes
        fileattributes  2D list or dictionary   global attributes of NetCDF file (history, description, ...)
        comp            boolean                 compress data on the fly using zlib
        vartype         string                  netcdf variable type
                                                default: 'f4' for normal variables
                                                         'f8' for variable with isdim=True and dims=None (=unlimited)
        create_var      boolean                 create variable for dimension although var is None

        Description
        -----------
        Open the NetCDF:
          fhandle =  nc.Dataset("Filename", 'w', format='NETCDF4')

        Then call writenetcdf with specifications of the dimensions:
          define
          - the dimension (isdim=True)
          - attributes of the dimensions (e.g. attributes=['units', 'm'])
          - the data of the dimensions as var=numpy.array or None (=unlimited dimension),

        Now write the data as variables:
        - when creating a variable, save the variable handle (vhandle)
        - specify the dimensions of data with dims
        - specify the timestep or timesteps with time


        Examples
        --------
        >>> import numpy as np
        >>> import netCDF4 as nc
        >>> fhandle =  nc.Dataset('writenetcdf_test.nc', 'w', format='NETCDF4')
        >>> dat    = np.array([[-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. , -9.,-9.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5., 5.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. ,  5.,-9., -9.],
        ...                     [-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.] ])
        >>> FiAtt   = ([['description', 'test writing with writenetcdf.py'],
        ...            ['history'    , 'Created by Matthias Zink']])
        >>> handle  = writenetcdf(fhandle, fileattributes=FiAtt)

        # dimensions
        >>> varName = 'time'
        >>> dims    = None
        >>> varAtt  = ([['units',    'hours since 2011-01-01 00:00:00'],
        ...             ['calendar', 'gregorian']])
        >>> thand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, isdim=True)

        >>> varName = 'lon'
        >>> varAtt  = ([['units'         , 'degrees_east'],
        ...             ['standard_name' , 'longitude'   ],
        ...             ['missing_value' , -9.]])
        >>> dims    = np.shape(dat)[0]
        >>> var     = np.arange(np.shape(dat)[0])+1
        >>> handle  = writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)

        >>> varName = 'lat'
        >>> varAtt  = ([['units'         , 'degrees_north'],
        ...             ['standard_name' , 'latitude'     ],
        ...             ['missing_value' , -9.]])
        >>> dims    = np.shape(dat)[1]
        >>> var     = np.arange(np.shape(dat)[1])+1
        >>> handle  = writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)

        # times & variable
        >>> times   = np.array([0.5, 1., 1.5, 2.])
        >>> handle  = writenetcdf(fhandle, thand, time=list(range(np.size(times))), var=times)

        >>> varName = 'TESTING'
        >>> varAtt  = {'units': 'm',
        ...            'long_name': 'Does this writing routine work?',
        ...            'missing_value': -9.}
        >>> dims    = ['time','lon','lat']
        >>> vhand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, comp=True)
        >>> for i in range(2):
        ...     handle  = writenetcdf(fhandle, vhand, time=i, var=dat*(i+1))
        >>> handle  = writenetcdf(fhandle, vhand, time=[2,3], var=np.array([dat*2,dat*3]))

        # type other than float
        >>> varName = 'TESTING2'
        >>> varAtt  = ([['units'        , 'm'                              ],
        ...             ['long_name'     ,'Does this writing routine work?'],
        ...             ['missing_value' , -9]])
        >>> dims    = ['time','lon','lat']
        >>> typ     = 'i4'
        >>> vhand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, comp=True, vartype=typ)
        >>> for i in range(2):
        ...     handle  = writenetcdf(fhandle, vhand, time=i, var=np.array(dat,dtype=np.int)*(i+1))
        >>> fhandle.close()


        # check file
        >>> from readnetcdf import readnetcdf
        >>> print([str(i) for i in readnetcdf('writenetcdf_test.nc', variables=True)])
        ['time', 'lon', 'lat', 'TESTING', 'TESTING2']
        >>> readdata = readnetcdf('writenetcdf_test.nc', var='TESTING')
        >>> print(np.any((readdata[0,:,:] - dat) != 0.))
        False
        >>> readdata2 = readnetcdf('writenetcdf_test.nc', var='TESTING2')
        >>> print(np.any((readdata2[0,:,:] - np.array(dat,dtype=np.int)) != 0.))
        False
        >>> if readdata.dtype == np.dtype('float32'): print('Toll')
        Toll
        >>> if readdata2.dtype == np.dtype('int32'): print('Toll')
        Toll

        >>> import os
        >>> os.remove('writenetcdf_test.nc')


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012-2013 Matthias Cuntz, Matthias Zink, Stephan Thober


        History
        -------
        Written,  MZ & MC, Feb 2012
        Modified, ST,      May 2012 - type 'f8' for time dimension
                  MC,      Jun 2012 - vartype
                  MC,      Feb 2013 - ported to Python 3
                  MC,      Apr 2014 - attributes can be given as dictionary e.g. from readnetcdf with attributes=True
                  ST,      May 2015 - added create_var flag that allows to disable automatic creation of variables for dimensions
    """
    # create File attributes
    if fileattributes is not None:
        if type(fileattributes) is dict:
            for k in fileattributes:
                fhandle.setncattr(k, fileattributes[k])
        else:
            for i in range(len(fileattributes)):
                fhandle.setncattr(fileattributes[i][0], fileattributes[i][1])
        return None

    # create dimensions
    if vhandle is not None:
        hand = vhandle
    else:
        if vartype is None:
            typ = 'f4'
            if isdim:
                if dims is None:
                    typ = 'f8'
        else:
            typ = vartype
        if isdim:
            idim = fhandle.createDimension(name, dims)
            # create variable for the dimension
            if not var is None or create_var:
                hand = fhandle.createVariable(name, typ, (name,))
        else:
            keys   = list(fhandle.dimensions.keys())
            for i in range(len(dims)):
                if dims[i] not in keys:
                    raise ValueError('Dimension '+str(dims[i])+' not in file dimensions: '+''.join([i+' ' for i in keys]))
            hand = fhandle.createVariable(name, typ, tuple(dims), zlib=comp)

    if attributes is not None:
        if type(attributes) is dict:
            for k in attributes:
                hand.setncattr(k, attributes[k])
        else:
            for i in range(len(attributes)):
                hand.setncattr(attributes[i][0], attributes[i][1])

    if var is not None:
        shand = hand.shape
        if time is not None:
            svar = np.shape(var)
            if np.size(np.shape(time)) == 0:
                if np.size(svar) != (np.size(shand)-1):
                    raise ValueError('Variable and handle dimensions do not agree for variable time: '+str(svar)+' and '+str(shand))
            elif np.size(np.shape(time)) == 1:
                if np.size(svar) != np.size(shand):
                    raise ValueError('Variable and handle dimensions do not agree for variable time vector: '+str(svar)+' and '+str(shand))
            else:
                raise ValueError('Time must be scalar or index vector.')
            if np.size(shand) == 1:
                hand[time] = var
            elif np.size(shand) == 2:
                hand[time,:] = var
            elif np.size(shand) == 3:
                hand[time,:,:] = var
            elif np.size(shand) == 4:
                hand[time,:,:,:] = var
            elif np.size(shand) == 5:
                hand[time,:,:,:,:] = var
            elif np.size(shand) == 6:
                hand[time,:,:,:,:,:] = var
            else:
                raise ValueError('Number of dimensions not supported (>6): '+str(np.size(shand)))
        else:
            if np.size(var) != np.size(hand):
                raise ValueError('Variable and handle elements do not agree: '+str(np.size(var))+' and '+str(np.size(hand)))
            if np.size(shand) == 1:
                hand[:] = var
            elif np.size(shand) == 2:
                hand[:,:] = var
            elif np.size(shand) == 3:
                hand[:,:,:] = var
            elif np.size(shand) == 4:
                hand[:,:,:,:] = var
            elif np.size(shand) == 5:
                hand[:,:,:,:,:] = var
            elif np.size(shand) == 6:
                hand[:,:,:,:,:,:] = var
            else:
                raise ValueError('Number of dimensions not supported (>6): '+str(np.size(shand)))
    if not var is None or create_var:
        return hand

# write to file
def dumpnetcdf( fname, dims=None, fileattributes=None, vnames=None, create=True, **variables ):
    """
        Writes variables with the same dimension to a netcdf file (1D-5D).
        It is a wrapper around writenetcdf.

        Variables must be given as keyword arguments. The key becomes the name of the variable
        in the netcdf file. The value becomes the array.


        Definition
        ----------
        def dumpnetcdf( fname, dims=None, fileattributes=None, vnames=None, create=True, **variables )


        Input           Format                  Description
        -----           -----                   -----------
        fname           string                  filename of netcdf file
        dims            1D list of strings      dimension names (e.g., ['time', 'x', 'y']),
                                                default is [ 'x', 'y', 'z', 'u', 'v' ]
        fileattributes  2D list or dictionary   global attributes of NetCDF file
                                                (e.g., history, description, ...)
        create          bool                    flag for creating file or appending variables
        vnames          dictionary              dictionary that specifies the variable name
                                                for each given variable
        variables       keyword arguments       keys become variable names, values are either
                                                numpy arrays to write OR list/tuple of
                                                [ numpy array, attributes ] where attributes
                                                must be a dictionary

        Description
        -----------
        write netcdf file in one single line
        dump_netcdf( 'test.nc', dims = [ 'time', 'northing', 'easting' ],
                     fileattributes
                     sm = sm_array, smi = [ smi_array, {'long_name': 'soil moisture index'} ] )

        Restrictions
        ------------
        all arrays must have the same dimensions and must be numpy arrays

        Examples
        --------
        >>> import numpy as np
        >>> dat    = np.array([[-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. , -9.,-9.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5., 5.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. ,  5.,-9., -9.],
        ...                     [-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.] ])
        >>> dims    = [ 'xx', 'yy' ]
        >>> FiAtt   = ([['description', 'test dump_netcdf'],
        ...             ['history'    , 'Created by Stephan Thober']])
        >>> dumpnetcdf( 'test_dump.nc', dims = dims, fileattributes = FiAtt,
        ...            data=( dat, {'long_name':'CHS'} ) )

        # check file
        >>> from readnetcdf import readnetcdf
        >>> print([str(i) for i in readnetcdf('test_dump.nc', variables=True)])
        ['data']
        >>> readdata = readnetcdf('test_dump.nc', var='data')
        >>> print(np.any((readdata[:,:] - dat) != 0.))
        False
        >>> if readdata.dtype == np.dtype('float64'): print('Toll')
        Toll

        >>> import os
        >>> os.remove('test_dump.nc')


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012-2013 Stephan Thober


        History
        -------
        Written,  ST, Sep 2014
        Modified  ST, Jan 2015 - bug fix, only parse ndim dimensions to writenetcdf
                  ST & MZ, Mar 2015 - added flag to append variables
                  ST, Jun 2016 - included dump for single numbers without attributes
                  ST, Aug 2016 - included separate variable names
                  MC, Nov 2016 - ported to Python 3, mostly dictionary behaviour
    """
    # check that variables are given
    if variables == dict():
        raise ValueError('Variables must be given as keyword arguments.')
    # check if dims are given
    if dims is None:
        dims = [ 'x', 'y', 'z', 'u', 'v' ]
    # check if variable names are given
    if vnames is not None:
        if type(vnames) != dict:
            raise ValueError('Variable names must be dictionary with the same keys as given variables')
        if len(vnames) != len(variables):
            raise ValueError('Number of given variable names does not match number of variables')
    # get dimensions from first variable
    vv1 = variables[list(variables.keys())[0]]
    if isinstance(vv1,(tuple,list)):
        arr_shape = vv1[0].shape
    else:
        arr_shape = vv1.shape
    # open netcdf file
    if create:
        fh = nc.Dataset( fname, 'w', 'NETCDF4' )
    else:
        fh = nc.Dataset( fname, 'a', 'NETCDF4' )
    # write file attributes
    if fileattributes is not None:
        writenetcdf(fh, fileattributes=fileattributes)
    # create dimensions according to first variable
    if create:
        for dd in range(len(arr_shape)):
            writenetcdf( fh, name = dims[dd], dims = arr_shape[dd],
                         var = None, isdim = True, create_var = False )
    else:
        # read dimensions from file
        file_dims = get_dims( fname )
        for dd in range(len(arr_shape)):
            # write dimensions if they do not exist
            if not dims[dd] in file_dims:
                writenetcdf( fh, name = dims[dd], dims = arr_shape[dd],
                            var = None, isdim = True, create_var = False )
    # loop over variables
    cc = 0
    for key in variables:
        value = variables[key]
        # update vnames if required
        if vnames is not None:
            key = vnames[key]
            cc += 1
        if len(value) == 1:
            # WRITE SINGLE NUMBER WITHOUT ATTRIBUTE
            if  value.ndim > len( dims ):
                raise ValueError( '***size mismatch: variable '+key )
            writenetcdf( fh, name = key, dims = dims[: value.ndim], var = value,
                         vartype=value.dtype, comp = True )
        else:
            if type( value[-1] ) == dict:
                # WRITE WITH ATTRIBUTE
                if  value[0].ndim > len( dims ):
                    raise ValueError( '***size mismatch: variable '+key )
                writenetcdf( fh, name = key, dims = dims[: value[0].ndim ], var = value[0],
                             vartype = value[0].dtype, attributes = value[-1], comp = True )
            else:
                # WRITE ARRAY WITHOUT ATTRIBUTE
                if  value.ndim > len( dims ):
                    raise ValueError( '***size mismatch: variable '+key )
                writenetcdf( fh, name = key, dims = dims[: value.ndim], var = value,
                            vartype=value.dtype, comp = True )
    fh.close()

# returns list of all dimensions given a filename
def get_dims( fname ):
    """
        This functions returns all dimension names contained in a netcdf file as list

        Definition
        ----------
        def get_dims( fname )


        Input           Format                  Description
        -----           -----                   -----------
        fname           string                  filename of netcdf file

        Description
        -----------
        reads all dimension names from netcdf file
        dnames = get_dims( 'test.nc' )

        Restrictions
        ------------
        None

        Examples
        --------
        >>> import numpy as np
        >>> dat    = np.array([[-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. , -9.,-9.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5., 5.,5., -9. ,  5., 5.,  5.],
        ...                     [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. ,  5.,-9., -9.],
        ...                     [-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.] ])
        >>> dims    = [ 'xx', 'yy' ]
        >>> FiAtt   = ([['description', 'test dump_netcdf'],
        ...            ['history'    , 'Created by Stephan Thober']])
        >>> dumpnetcdf( 'test_dump.nc', dims = dims, fileattributes = FiAtt,
        ...            data = ( dat, {'long_name':'CHS'} ) )

        # check file
        >>> import sys
        >>> pyver = sys.version_info
        >>> from readnetcdf import readnetcdf
        >>> out = get_dims( 'test_dump.nc' )
        >>> print(out) if pyver > (3,0) else print([ i.encode('UTF-8') for i in out ])
        ['xx', 'yy']

        >>> import os
        >>> os.remove('test_dump.nc')


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2015-2015 Stephan Thober


        History
        -------
        Written,  ST, Jun 2015
                  MC, Nov 2016 - adapted docstring to Python 2 and 3
    """
    file_vars = readnetcdf( fname, var = '', variables = True )
    file_dims = []
    for var in file_vars:
        tmp_dims = readnetcdf( fname, var = var, dims = True )
        for dims in tmp_dims:
            if not dims in file_dims:
                file_dims.append( dims )
    return file_dims


if __name__ == '__main__':
    import doctest
    try:
        import netCDF4 as nc
        doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    except:
        raise IOError('No NetCDF4 support available.')
