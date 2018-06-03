#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from netcdf4 import NcDataset


def concatDimension(fnames, dim, dimvar=None, outfname=None, **kwargs):
    """
    Arguments:
    ----------
    datasets (List[str]):          datasets to concatenate
    dim (String):                  name of the dimension to concatenate
    dimvar (Optional[String]):     name of the variable defining values
                                   for the concatenation dimension, 
                                   default: dimvar = dim
    outfname (Optional[str]):      file name of an output dataset
    kwargs (Optional[Any]):        paramaters to pass to the createVariable
                                   Method of NcDataset
                                

    Return
    ------
    Optional[netcdf4.NcDataset]

    Purpose
    -------
    Concatenate the given datasets into a new one along the given dimension
    'dim'. Values of the dimensions need to be given in 'dimvar' (default
    variable name is the same as the concatenation dimension).
    If 'outfname' is given the concatenation is file based, otherwise an in-memory
    dataset will be returned.
    """

    def getDimVarValues(fnames, vardim):
        vals = []
        types = []
        for fname in fnames:
            with NcDataset(fname, "r") as nc:
                # TODO:
                # check if values are given, otherwise make up
                # some fake data simpling starting at 0
                data = nc.variables[vardim][:]
                types.append(data.dtype)
                vals.extend(data.tolist())
        return np.unique(vals).astype(dtype=np.find_common_type(types, []))
    
    def getCommonVarTypes(fnames):
        vars = {}
        for fname in fnames:
            with NcDataset(fname, "r") as nc:
                for vname, var in nc.variables.items():
                    try:
                        vars[vname].append(var.dtype)
                    except KeyError:
                        vars[vname] = [var.dtype]
        return {
            vname: np.find_common_type(vtypes, [])
            for vname, vtypes in vars.items()}
            
    if not dimvar:
        dimvar = dim

    dimvals = getDimVarValues(fnames, dimvar)
    dtypes = getCommonVarTypes(fnames)

    # the new concat dimension
    out = NcDataset(outfname, "w")
    out.createDimension(dim, None, fail=False)
    var = out.createVariable(dimvar, dimvals.dtype, (dim,), fail=False, **kwargs)
    var[:] = dimvals

    for i, fname in enumerate(fnames):

        with NcDataset(fname, "r") as nc:
        
            out.copyDimensions(nc.dimensions, skip=out.dimensions, fix=True)
            out.copyAttributes(nc.attributes)

            for vname, var in nc.variables.items():

                outvar = out.copyVariable(
                    var, dtype=dtypes[vname], data=False, fail=False, **kwargs)

                # build up the indices
                if dim in var.dimensions:
                    # idx = [dimvals.index(v) for v in nc.variables[dimvar]]
                    idx = [np.where(dimvals == v)[0][0] for v in nc.variables[dimvar]]
                    slices = [slice(None,)] * len(var.dimensions)
                    slices[var.dimensions.index(dim)] = idx
                else:
                    # writing them again is stupid...
                    slices = slice(None)

                data = var[:]

                # this seems to be necessary, because netCDF4 sometimes comes
                # up with a masked without a explicit user definition
                # (i.e. no fill_value, valid_range, whatever is set)
                if isinstance(data, np.ma.MaskedArray):
                    data = data.data
                outvar[slices] = data

    if outfname is None:
        return out
    out.close()
