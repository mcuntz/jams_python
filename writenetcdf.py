#!/usr/bin/env python
import numpy as np                       # array manipulation
import netCDF4 as nc                     # netCDF interphase

def writenetcdf(fhandle, vhandle=None, var=None, time=None, isdim=False, name=None, dims=None,
                attributes=None, fileattributes=None, comp=False):
    """ 
        Writes dimensions, variables, dependencies and attributes to NetCDF file
        only for 2D or 3D data

        except the data, all variables given to the functions have to be python lists
        
        Definition
        ----------
        def writenetcdf(fhandle, dim='', var=None, name='DATA', dims='',
                        attributes='', fileattributes=''):


        Input          Format         Description
        -----          -----          -----------
        fhandle         string        File handle of nc.Dataset(FileName, 'w')
        dim      python list   Dimensions of the NetCDF file
        var       array like    Data
        name         string        Name of the variable
        dims          python list   Variable Dependencies, var are assumed to be float4
        attributes   python list   Variable Attributes
        fileattributes  python list   Attributes of NetCDF file (history, description, ...)


        Description
        ------
        Please at first call the writenetcdf with the specification of the <dim>, the data 
        of the dimensions as <var> and the attributes of the variables. Please make sure to
        set <open> to 'True'.

        In the following you can put the variables into the file. 2D is always assigned ti time step 0.
        3D arrays should contain the field data in the first and second dimension, while the third dimension
        represents the time.

        EXAMPLES:
        ------
        varAtt=[['units', 'hours since 2000-01-01 00:00:00'],['calendar','gregorian']]
        writenetcdf('test.nc', dim=['time', 'None'], var =[0], attributes=varAtt, open=True)
        
        varAtt=[['units', 'm'],['long_name','easting'], ['missing_value', -9999.]]
        FiAtt=[['description', 'radiation data determined from HDF5 files'],['history','Created by me'],['source','LSA SAF http://landsaf.meteo.pt/']]
        writenetcdf('test.nc', dim=['xc', '4'], var = [1.,2.,3.,4.], attributes=varAtt, open=False)
        
        varAtt=[['units', 'm'],['long_name','northing'], ['missing_value', float(-9999.0)]]
        writenetcdf('test.nc', dim=['yc', '4'], var = [5.,6.,7.,8.], attributes=varAtt, open=False)
        
        varAtt=[['units', 'm3s-1']               ,\
        ['long_name','river discharge']  ,\
        ['missing_value', -9999.]]
        datas          = np.array([[1.,2.,3.,4.], \
        [1.,2.,3.,4.], \
        [1.,2.,3.,4.], \
        [1.,2.,3.,-9999.]]) #2D
        #datas          = np.array([  [ [1.,2.,3.,4.],[1.,2.,3.,4.],[1.,5.,6.,7.],[1.,8.,9.,10.]], [[1.,2.,3.,4.],[1.,2.,3.,4.],[1.,5.,6.,7.],[1.,8.,9.,10.]]]) #3d
        writenetcdf('test.nc', name = 'TEST', dims  = ('time', 'yc', 'xc'), var = datas, attributes=varAtt, open=False)
        
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
            print shand, np.size(shand)
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
                print 'hand: ',repr(hand[time]),var
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
