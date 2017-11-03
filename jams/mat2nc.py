#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import netCDF4 as nc
import scipy.io as sio
import datetime
import os
import warnings
from collections import OrderedDict

from jams.writenetcdf import writenetcdf
from jams.readnetcdf  import readnetcdf

def mat2nc(fname, overwrite=False, fname_out=None, verbose=True, squeeze=True,
           varname=None,
           lat=None, lon=None, lat_file=None, lon_file=None, transposed=False, regular_grid=True,
           anchor_time=None, data_time=None):
    """
        This functions writes the content of a Matlab file into a NetCDF file.
        It allows for specifying longitudes and latitudes as well as a validity time.

        Definition
        ----------
        def mat2nc( fname )


        Input           Format                  Description
        -----           -----                   -----------
        fname           (array-like of) string  Filename of matlab file.
                                                If a list of file names is given, multiple files will be read and written in one netcdf file.
                                                This feature is implemented such that multiple time steps can be written to one file.
                                                Make sure that all files specified have the same structure (e.g. same grid) and that only
                                                one variable is stored per file. The "varname", "fname_out", "anchor_time" and an array of
                                                "data_time" should be specified in this case.
        overwrite       bool                    Optional argument.
                                                True:  existing nc-file will be overwritten,
                                                False: if nc-file exists, script will be stopped
                                                Default: False
        fname_out       string                  Optional argument.
                                                Name of output NetCDF file
                                                Default: abc.mat --> abc.nc
        verbose         bool                    Optional argument.
                                                True:  variable names found etc. are print on screen
                                                False: silent
                                                Default: True
        squeeze         bool                    Optional argument.
                                                True:  all axis with shape 1 are removed
                                                False: axis with shape 1 are written (arrays will have same shape like in Matlab)
                                                Default: True
        varname         (array-like of) string  Optional argument.
                                                Use this name for variable in NetCDF file.
                                                If there are multiple variables in one Matlab file a list of
                                                varnames should be specified.
                                                Default: Name of variable in Matlab file
        lat             array-like              Optional argument for input specifying latitude
                                                (needs to match dimensions of all variables read from mat-file)
        lon             array-like              Optional argument for input specifying longitude
                                                (needs to match dimensions of all variables read from mat-file)
        lat_file        string                  Optional argument for input file (*.mat) specifying latitude
        lon_file        string                  Optional argument for input file (*.mat) specifying longitude
        transposed      boolean                 Optional argument. If true all gridded data will be transposed.
                                                This might be necessary if data are stored with swapped lon and lat,
                                                i.e. longitude is y-axis and latitude is x-axis.
        anchor_time     string                  Time where data in Matlab file are valid. Has to follow the pattern
                                                    'hours since YYYY-MM-DD HH:MM:SS'
        data_time       array-like              Time of each data set written, i.e. if there are 4 precipitations fields and
                                                they are valid 6, 12, 18, 24 hrs after anchor time, data_time=[6,12,18,24].


        Description
        -----------
        None

        Restrictions
        ------------
        None

        Examples
        --------
        >>> import scipy.io as sio
        >>> ss = np.array([[10.0]])
        >>> tt = np.array([[10.0, 1.0, 1.0],[20.0, 2.0, 2.0]])
        >>> uu = np.array([[[10.0], [1.0], [1.0]],[[20.0], [2.0], [2.0]]])
        >>> vv = np.array([10.0,1.0, 1.0, 20.0, 2.0, 2.0])
        >>> sio.savemat('test.mat', {'ss':ss, 'tt':tt, 'uu':uu, 'vv':vv})

        # --------------------------------------------
        # Write squeezed arrays (all axis with shape 1 are removed)
        # --------------------------------------------
        >>> if os.path.isfile('test.nc'): os.remove('test.nc')
        >>> mat2nc('test.mat')
        writes  'ss'                  shape =  (1,)                  to *.nc
        writes  'tt'                  shape =  (2, 3)                to *.nc
        writes  'uu'                  shape =  (2, 3)                to *.nc
        writes  'vv'                  shape =  (6,)                  to *.nc
        wrote file  test.nc

        >>> from jams.readnetcdf import readnetcdf
        >>> readnetcdf('test.nc',var='ss')
        array([ 10.])
        >>> readnetcdf('test.nc',var='tt')
        array([[ 10.,   1.,   1.],
               [ 20.,   2.,   2.]])
        >>> readnetcdf('test.nc',var='uu')
        array([[ 10.,   1.,   1.],
               [ 20.,   2.,   2.]])
        >>> readnetcdf('test.nc',var='vv')
        array([ 10.,   1.,   1.,  20.,   2.,   2.])

        # --------------------------------------------
        # Write un-squeezed, original arrays (same shape as in Matlab)
        # --------------------------------------------
        >>> if os.path.isfile('test.nc'): os.remove('test.nc')
        >>> mat2nc('test.mat',squeeze=False,verbose=False)

        >>> readnetcdf('test.nc',var='ss')
        array([[ 10.]])
        >>> readnetcdf('test.nc',var='tt')
        array([[ 10.,   1.,   1.],
               [ 20.,   2.,   2.]])
        >>> readnetcdf('test.nc',var='uu')
        array([[[ 10.],
                [  1.],
                [  1.]],
               [[ 20.],
                [  2.],
                [  2.]]])
        >>> readnetcdf('test.nc',var='vv')
        array([[ 10.,   1.,   1.,  20.,   2.,   2.]])

        >>> if os.path.isfile('test.nc'): os.remove('test.nc')
        >>> if os.path.isfile('test.mat'): os.remove('test.mat')


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

        Copyright 2016-2016 Juliane Mai


        History
        -------
        Written,  JM, Sep 2016
        Modified, JM, Oct 2016 - add feature xy
                  MC, Nov 2016 - single file not working in Python 3, Order .mat content by name (for docstrings).
    """

    # make sure that fname is an array
    import sys
    if sys.version_info > (3,0): basestring = str
    if hasattr(fname, "__iter__") and not isinstance(fname, basestring):
        filename_mat = np.array(fname)
    else:
        filename_mat = np.array([fname])

    # make sure that varname is an array
    if varname is not None:
        if hasattr(varname, "__iter__"):
            varname = np.array(varname)
        else:
            varname = np.array([varname])

    # make sure that data_time is an array
    if data_time is not None:
        if hasattr(data_time, "__iter__"):
            data_time = np.array(data_time)
        else:
            data_time = np.array([data_time])

    # make sure that output filename is a valid one
    if fname_out is None:
        filename_nc  = '.'.join(filename_mat[0].split('.')[:-1])+'.nc'
        if filename_nc == '.nc': # in case input filename is missing file ending ".mat"
            filename_nc = 'new.nc'
    else:
        filename_nc  = fname_out

    # make sure that output file does not exist
    if os.path.isfile(filename_nc) :
        if overwrite:
            # remove existing nc-file
            os.remove(filename_nc)
        else:
            # stop script
            raise ValueError('mat2nc: NetCDF file already exists')

    # assure that time points are sorted
    if ( not (data_time is None) ):
        idx = np.argsort(data_time)
        data_time = data_time[idx]
        filename_mat = filename_mat[idx]

    # checking of input arguments
    if (not(lon is None) and lat is None) or (lon is None and not(lat is None)):
        raise ValueError('mat2nc: either both lat and lon have to be specified or none of it')

    if (not(lon_file is None) and lat_file is None) or (lon_file is None and not(lat_file is None)):
        raise ValueError('mat2nc: either both lat_file and lon_file have to be specified or none of it')

    if (not(lat_file is None) and not(lat is None)):
        raise ValueError('mat2nc: grid information (lat/lon) can either be given as file or as array but not both')

    # initialize counters and variables
    key_content  = None
    grid_avail   = (not (lon is None)) or (not (lon_file is None))
    time_avail   = not (anchor_time is None)
    dim_name     = np.array([ 'dim_'+str(ii) for ii in range(100) ])  # ['dim_0', 'dim_1', 'dim_2', ...., 'dim_99']
    dims         = np.zeros(np.shape(dim_name),dtype=int)
    idim         = 0 # index where to write next dim

    # Read longitude info from file
    cc = 0
    if (not lon_file is None):
        if not os.path.isfile(lon_file):
            raise ValueError('mat2nc: Longitude file not found')
        lon_contents = sio.loadmat(lon_file)
        for ikey,key in enumerate(lon_contents):
            if key[0:2] != '__':
                cc += 1
                lon = np.array(lon_contents[key])

                if transposed:
                    lon = np.transpose(lon)

        if regular_grid:
            # look if longitudes are constant per column
            lon_tmp   = np.round(lon,3)
            lon_const = np.all([ list(lon_tmp[:,ii]).count(lon_tmp[0,ii]) == len(lon_tmp[:,ii]) for ii in range(np.shape(lon_tmp)[1])])

            if not lon_const:
                print('lon[:,0]: ',lon_tmp[:,0])
                raise ValueError('mat2nc: Longitudes are not constant per column. Use either irregular=True or try transpose=True')

    if (cc > 1):
        raise ValueError('mat2nc: More than one variable found in longitude file')


    # Read latitude info from file
    cc = 0
    if (not lat_file is None):
        if not os.path.isfile(lat_file):
            raise ValueError('mat2nc: Latitude file not found')
        lat_contents = sio.loadmat(lat_file)
        for ikey,key in enumerate(lat_contents):
            if key[0:2] != '__':
                cc += 1
                lat = np.array(lat_contents[key])

                if transposed:
                    lat = np.transpose(lat)

        if regular_grid:
            # look if longitudes are constant per column
            lat_tmp   = np.round(lat,3)
            lat_const = np.all([ list(lat_tmp[ii,:]).count(lat_tmp[ii,0]) == len(lat_tmp[ii,:]) for ii in range(np.shape(lat_tmp)[0])])

            if not lat_const:
                print('lat[0,:]: ',lat_tmp[0,:])
                raise ValueError('mat2nc: Latitudes are not constant per row. Use either irregular=True or try transpose=True')

    if (cc > 1):
        raise ValueError('mat2nc: More than one variable found in latitude file')

    # open netcdf file and add some general information
    fh = nc.Dataset( filename_nc, 'w', 'NETCDF4' )
    # File attributes
    FiAtt   = ([['description', 'Dump content of '+", ".join(filename_mat)+' to NetCDF'],
                ['history'    , 'Created by Juliane Mai (mat2nc.py)']])
    handle  = writenetcdf(fh, fileattributes=FiAtt)

    dd = 0
    # if time available, write this dimensions
    if (time_avail):
        varName     = 'ntime'
        varAtt  = ([['units',    anchor_time],
                    ['calendar', 'gregorian']])
        ntime        = np.shape(np.array(data_time))[0]
        dims[dd]     = 1 # means that only one time step is written
        dim_name[dd] = varName
        # create dimension
        th = writenetcdf(fh, name=dim_name[dd], dims=None, attributes=varAtt, isdim=True)
        # write values to dimension
        times   = np.array(data_time)
        writenetcdf(fh, th, time=range(ntime), var=times)

        dd += 1

    # if grid available, write this dimensions and data first
    if (grid_avail):

        # (1) check if variable shape matches given grid shape
        arr_shape_1 = np.shape(lat)
        arr_shape_2 = np.shape(lon)
        if (arr_shape_1 != arr_shape_2):
            raise ValueError('mat2nc: shapes of given lat and lon grids are not matching')

        # TODO: If dim_lat = dim_lon the identification of lon and lat axis in variables might not work properly
        if (np.shape(lon)[0] == np.shape(lat)[1]):
            print('WARNING: ')
            raise ValueError('mat2nc: dimension of lat and lon are equal. This is not supported yet.')

        # (2) write dimensions and lat/lon data
        #        lon is x-axis ==> number columns
        #        lat is y-axis ==> number of rows
        varName     = 'nlon'
        varAtt      = ([['units'         , 'degrees_east'],
                        ['standard_name' , 'longitude'   ],
                        ['missing_value' , -9.]])
        nlon         = np.shape(lon)[1]
        dims[dd]     = nlon
        dim_name[dd] = varName
        if regular_grid:
            # create dimension and write variable together
            var          = lon[0] #lon[:,0]
            writenetcdf(fh, name=dim_name[dd], dims=dims[dd], var=var, attributes=varAtt, isdim=True)
        else:
            # just create dimension but don't write variable yet
            writenetcdf(fh, name=dim_name[dd], dims=dims[dd],
                        var=None, isdim=True, create_var=False )
        dd += 1
        varName     = 'nlat'
        varAtt      = ([['units'         , 'degrees_north'],
                        ['standard_name' , 'latitude'     ],
                        ['missing_value' , -9.]])
        nlat         = np.shape(lat)[0]
        dims[dd]     = nlat
        dim_name[dd] = varName
        if regular_grid:
            var          = lat[:,0] #lat[0]
            writenetcdf(fh, name=dim_name[dd], dims=dims[dd], var=var, attributes=varAtt, isdim=True)
        else:
            # just create dimension but don't write variable yet
            writenetcdf(fh, name=dim_name[dd], dims=dims[dd],
                        var=None, isdim=True, create_var=False )
        dd += 1

        # now both dimensions exist and grid can be written as extra variable
        if not regular_grid:
            idx=[ np.where(dims==dd)[0][0] for dd in np.shape(lon) ]

            attributes = {'standard_name': 'longitude', 'long_name': 'longitudes of irregular grid', 'units': 'degrees_??', 'missing_value': -9.}
            writenetcdf(fh, name='lon', dims=dim_name[idx], var=lon,
                        comp=True, create_var=True,
                        attributes=attributes,
                        vartype=lon.dtype )

            idx=[ np.where(dims==dd)[0][0] for dd in np.shape(lat) ]
            attributes = {'standard_name': 'latitude', 'long_name': 'latitudes of irregular grid', 'units': 'degrees_??', 'missing_value': -9.}
            writenetcdf(fh, name='lat', dims=dim_name[idx], var=lat,
                        comp=True, create_var=True,
                        attributes=attributes,
                        vartype=lat.dtype )


    tt = 0
    for ff in filename_mat:
        vv = 0

        # read content from Matlab file
        mat_contents = sio.loadmat(ff)
        mat_contents = OrderedDict(sorted(mat_contents.items(), key=lambda t: t[0]))

        for ikey, key in enumerate(mat_contents):
            if key[0:2] != '__':

                key_content = np.array(mat_contents[key])
                tmp = key_content

                if transposed:
                    tmp = np.transpose(tmp)

                arr_shape = np.shape(tmp)

                # if there is exactly one entry for time_data given, i.e. the data read is the data for exactly one time point
                # but it might be that this dimension is missing, i.e. data_shape=[lat, lon] but should be data_shape=[1, lat, lon]
                if (time_avail):
                    if (not (1 in arr_shape)):
                        tmp = np.expand_dims(tmp,axis=0)
                        arr_shape = np.shape(tmp)

                if (not grid_avail and not time_avail):
                    # remove axis with shape 1, i.e.
                    #        shape=[1,23,30,1,4] --> shape=[23,30,4]
                    # but    shape=[1,1,1]       --> shape=[1] (otherwise cannot be written to NetCDF)
                    if (squeeze):
                        if  (hasattr(tmp, "__len__")):
                            tmp = np.squeeze(tmp)

                            # in case tmp contains only one scalar, make again proper array of it
                            if (tmp.shape == ()):
                                tmp = np.array([tmp])

                if (not (varname is None)):
                    key = varname[vv]

                if verbose:
                    print('writes ',repr(key).ljust(35),' shape = ',repr(np.shape(tmp)).ljust(20),' to *.nc')

                # write (missing) dimensions
                arr_shape = np.shape(tmp)
                for dd in range(len(arr_shape)):
                    if (not arr_shape[dd] in dims) and (not dim_name[idim] in ['nlat','nlon','ntime']):

                        dims[idim] = arr_shape[dd]
                        writenetcdf(fh, name=dim_name[idim], dims=dims[idim],
                                    var=None, isdim=True, create_var=False )
                        idim +=1
                        if (idim > len( dim_name )):
                            raise ValueError( 'mat2nc: Too many different dimensions found. You may have to increase size of dim_name.' )

                # find matching dim names for current shape
                idx=[ np.where(dims==arr_shape[dd])[0][0] for dd in range(len(arr_shape)) ]

                if (grid_avail):
                    # check if dimension matching lon and lat are found for this variable
                    if (not 'nlat' in dim_name[idx]) or (not 'nlon' in dim_name[idx]):
                        if verbose:
                            print('    WARNING: lat and/or lon dimension not found for variable')

                attributes = {'long_name': key, 'units': 'not defined'}
                if (not (varname is None)):
                    if (varname[vv].lower() == 'av'):
                        attributes = {'long_name': 'Accumulated surface flux of latent heat (accumulated FV)', 'units': 'J/m^2'}
                    elif (varname[vv].lower() == 'fb'):
                        attributes = {'long_name': 'Downward solar flux', 'units': 'W/m^2'}
                    elif (varname[vv].lower() == 'fi'):
                        attributes = {'long_name': 'Surface incoming infrared flux', 'units': 'W/m^2'}
                    elif (varname[vv].lower() == 'hu'):
                        attributes = {'long_name': 'Specific humidity', 'units': 'kg/kg'}
                    elif (varname[vv].lower() == 'p0'):
                        attributes = {'long_name': 'Surface pressure', 'units': 'mbar'}
                    elif (varname[vv].lower() == 'pr'):
                        attributes = {'long_name': 'Precipitation', 'units': 'mm'}
                    elif (varname[vv].lower() == 'tt'):
                        attributes = {'long_name': 'Air temperature', 'units': 'degree C'}
                    elif (varname[vv].lower() == 'uu'):
                        attributes = {'long_name': 'U-component of the wind (along the grid X axis)', 'units': 'kts'}
                    elif (varname[vv].lower() == 'vv'):
                        attributes = {'long_name': 'V-component of the wind (along the grid Y axis)', 'units': 'kts'}

                if (time_avail):
                    if (tt==0):
                        vh   = writenetcdf(fh, name=key, dims=dim_name[idx],
                                           comp=True,
                                           attributes=attributes,
                                           vartype=key_content.dtype)
                    writenetcdf(fh, vh, var = tmp[0], time = tt)
                    tt += 1
                else:
                    writenetcdf(fh, name=key, dims=dim_name[idx], var=tmp,
                                comp=True,
                                attributes=attributes,
                                vartype=key_content.dtype )
                vv += 1

    fh.close()

    if verbose:
        print('wrote file ',filename_nc)

if __name__ == '__main__':
    import doctest
    try:
        import netCDF4 as nc
        doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    except:
        raise IOError('No NetCDF4 support available.')
