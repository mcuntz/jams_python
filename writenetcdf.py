#!/usr/bin/env python
import numpy as np                       # array manipulation
import netCDF4 as nc                     # netCDF interphase

def writenetcdf(fhandle, vhandle=None, var=None, time=None, isdim=False, name=None, dims=None,
                attributes=None, fileattributes=None, comp=False):
    """ 
        Writes dimensions, variables, dependencies and attributes to NetCDF file
        for 1D or 9D data

        except the data, all variables given to the functions have to be python lists
        
        Definition
        ----------
        def writenetcdf(fhandle, vhandle=None, var=None, time=None, isdim=False, name=None, dims=None,
                        attributes=None, fileattributes=None, comp=False)

        Input          Format                 Description
        -----          -----                  -----------
        fhandle        string                 file handle of nc.Dataset(FileName, 'w')
        vHandle        string                 varaible handle of the particular variable
        var            array like             data (assumed to be float4)
        time           integer or array like  particular time step of the data
        isdim          boolean                defines if current var is a dimension
        name           string                 name of the variable
        dims           python list(1D)        variable dependencies (e.g. [time, x, y])
        attributes     python list(2D)        variable attributes
        fileattributes python list(2D)        global attributes of NetCDF file (history, description, ...)
        comp           boolean                compress data on the fly using zlib

        Description
        ------
        Please at first open the NetCDF with (save fhandle to write later on into the file)
        fhandle =  nc.Dataset("Filename", 'w', format='NETCDF4')

        Afterwards call writenetcdf with the specification of the dimensions:
        - define that the current objects is a deimension (<isdim>=True)
        - the data of the dimensions as (<var> = ARRAY(1D) or None(unlimited dimension)), 
        - if you like some attributes of the dimesnions (<attributes>=['units', 'm'])

        Subsequently you can put the data as variables into the file:
        - when creating a variable please save the variable handle (vhand) 
        - be aware to specify the dimensions of the data with <dims>
        - specify current timestep or timesteps with <time>

        EXAMPLES:
        ------
        >>> import numpy as np
        >>> import netCDF4 as nc
        >>> data    = np.array([[-9., 5., 5., -9. ,5.,-9.,5.,     -9. ,  5., 5.,  5.], \
                                [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. , -9.,-9.,  5.], \
                                [ 5.,-9.,-9., -9. ,5., 5.,5., -9. ,  5., 5.,  5.], \
                                [ 5.,-9.,-9., -9. ,5.,-9.,5., -9. ,  5.,-9., -9.], \
                                [-9., 5., 5., -9. ,5.,-9.,5., -9. ,  5., 5.,  5.] ])
        >>> fhandle =  nc.Dataset('writenetcdf_test.nc', 'w', format='NETCDF4')
        >>> FiAtt=([['description', 'test writing with writenetcdf.py'], \
                    ['history'    , 'Created by Matthias Zink']         ])
        >>> handle  = writenetcdf(fhandle, fileattributes=FiAtt)

        # Dimesnions
        >>> varName = 'time'
        >>> dims    = None
        >>> varAtt  =([['units',    'hours since 2011-01-01 00:00:00'], \
                       ['calendar', 'gregorian']                       ])
        >>> thand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, isdim=True)                

        >>> varName = 'lon'
        >>> varAtt  = ([['units'         , 'degrees_east'], \
                        ['standard_name' , 'longitude'   ], \
                        ['missing_value' , -9.]           ])
        >>> dims    = np.shape(data)[0]
        >>> var     = np.arange(np.shape(data)[0])+1
        >>> handle  = writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)

        >>> varName = 'lat'
        >>> varAtt  = ([['units'         , 'degrees_north'], \
                        ['standard_name' , 'latitude'     ], \
                        ['missing_value' , -9.]            ])
        >>> dims    = np.shape(data)[1]
        >>> var     = np.arange(np.shape(data)[1])+1
        >>> handle  = writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)

        # variable & time
        >>> varName = 'TESTING'
        >>> varAtt  = ([['units'        , 'm'                              ], \
                        ['long_name'     ,'Does this writing routine work?'], \
                        ['missing_value' , -9                              ]])
        >>> dims    = ['time','lon','lat']
        >>> vhand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, comp=True)
        >>> for i in xrange(2): \
                handle  = writenetcdf(fhandle, vhand, time=i, var=data*(i+1))
        >>> handle  = writenetcdf(fhandle, vhand, time=[2,3], var=np.array([data*2,data*3]))
        >>> times       = [0.5, 1., 1.5, 2.]
        >>> handle  = writenetcdf(fhandle, thand, time=range(4), var=times)

        >>> from readnetcdf import *
        >>> readnetcdf('writenetcdf_test.nc', variables=True)
        [u'time', u'lon', u'lat', u'TESTING']
        >>> readdata= readnetcdf('writenetcdf_test.nc', var='TESTING')
        >>> print np.any((np.ma.getdata(readdata[0,:,:]) - data) != 0.)
        False

        >>> import os
        >>> os.remove('writenetcdf_test.nc')

        History
        -------
        Written,  Matthias Zink,  Feb. 2012
        Modified, Stephan Thober, May  2012 - added double precision for time handle
    """

    # create File attributes
    if fileattributes != None:
        for i in xrange(len(fileattributes)):
            fhandle.setncattr(fileattributes[i][0], fileattributes[i][1])
        return None

    # create dimensions
    if vhandle != None:
        hand = vhandle
    else:
        if isdim:
            idim = fhandle.createDimension(name, dims)
            # create variable for the dimension
            if dims == None:
                hand = fhandle.createVariable(name, 'f8', (name,))
            else:
                hand = fhandle.createVariable(name, 'f4', (name,))
        else:
            keys   = fhandle.dimensions.keys()
            for i in xrange(len(dims)):
                if dims[i] not in keys:
                    raise ValueError('writenetcdf error: Dimension '+str(dims[i])+' not in file dimensions: '+''.join([i+' ' for i in keys]))
            hand = fhandle.createVariable(name, 'f4', (dims), zlib=comp)

    if attributes != None:
        for i in xrange(len(attributes)):
            hand.setncattr(attributes[i][0], attributes[i][1])

    if var != None:
        shand = hand.shape
        if time != None:
            svar = np.shape(var)
            if np.size(np.shape(time)) == 0:
                if np.size(svar) != (np.size(shand)-1):
                    raise ValueError('writenetcdf error: Variable and handle dimensions do not agree for variable time: '+str(svar)+' and '+str(shand))
            elif np.size(np.shape(time)) == 1:
                if np.size(svar) != np.size(shand):
                    raise ValueError('writenetcdf error: Variable and handle dimensions do not agree for variable time vector: '+str(svar)+' and '+str(shand))
            else:
                raise ValueError('writenetcdf error: Time must be scalar or index vector.')
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
                raise ValueError('writenetcdf error: Number of dimensions not supported (>6): '+str(np.size(shand)))
        else:
            if np.size(var) != np.size(hand):
                raise ValueError('writenetcdf error: Variable and handle elements do not agree: '+str(np.size(var))+' and '+str(np.size(hand)))
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
                raise ValueError('writenetcdf error: Number of dimensions not supported (>6): '+str(np.size(shand)))
    return hand

if __name__ == '__main__':
    import doctest
    doctest.testmod()


