#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

__all__ = ['nc2nc']

def nc2nc(ifile, ofile, dvar=[], rvar={}, rname={}, ratt={}, hist=None):
    '''
        Copies netcdf file deleting variables, replacing values of other variables.


        Definition
        ----------
        def nc2nc(ifile, ofile, dvar=[], rvar={}, rname={}, ratt={}):


        Input
        -----
        ifile   string
                input netcdf file name
        ofile   string
                output netcdf file name


        Optional Input
        --------------
        dvar    list
                variable names (case sensitive) that will not be copied to output file
        rvar    dictionary
                dictionary keys are the variable names (case sensitive),
                which values will be set given by the dictionary values
                instead of the values from the input file
        rname   dictionary
                dictionary keys are the variable names (case sensitive),
                which will be renamed with dictionary values in the output file
        ratt    dictionary of dictionaries
                Replace/set attributes of variables in dictionary keys (case sensitive).
                Dictionary values are also dictionaries with {'attribute_name': attribute_value}


        Output
        ------
        netcdf file


        Notes
        -----
        Make NETCDF4 files, using compression.


        Examples
        --------
        nc2nc('in.nc', 'out.nc',
              dvar=['var1','var2'],
              rvar={'var3': np.ones(100)},
              rname={'var3':'variable3'}
              ratt={'var3':{'long_name':'New variable 3', units='kg/m^2/s'}})


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2019 Matthias Cuntz


        History
        -------
        Written,  Matthias Cuntz, Mar 2019 - from a code found by Vanessa. Source unknown.
        Modified, Matthias Cuntz, Apr 2019 - allow variable removal, renaming and attribute editing
                                           - use NETCDF4 and zlib variables.
    '''
    import numpy as np
    import netCDF4 as nc

    if not isinstance(dvar, (list, tuple, np.ndarray)): dvar = [dvar]
    rvlist = list(rvar.keys())
    rnlist = list(rname.keys())
    ralist = list(ratt.keys())
    
    with nc.Dataset(ifile) as src, nc.Dataset(ofile, 'w', format='NETCDF4') as dst:
        # copy file attributes
        for name in src.ncattrs():
            dst.setncattr(name, src.getncattr(name))
        # extend history
        if hist is not None:
            if 'history' in src.ncattrs():
                histo = dst.getncattr('history')+'\n'+hist
                dst.setncattr('history', histo)
            else:
                dst.setncattr('history', hist)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
        # copy all file data, deleting variables, replacing other variable data
        for name, variable in src.variables.items():
            if (name in dvar):
                # delete variable
                pass
            else:
                # create output variable
                if name in rnlist:
                    oname = rname[name]
                else:
                    oname = name
                x = dst.createVariable(oname, variable.datatype, variable.dimensions, zlib=True)
                x.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
                # Overwrite or set new attributes
                if name in ralist:
                    iatt = ratt[name]
                    for rn in iatt:
                        x.setncattr(rn, iatt[rn])
                if (name in rvlist):
                    # replace variable data
                    dst.variables[oname][:] = rvar[name]
                else:
                    # copy variable
                    dst.variables[oname][:] = src.variables[name][:]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
