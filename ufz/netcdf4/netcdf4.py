#! /usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset, Variable, Group, date2num, num2date
from collections import OrderedDict

def _tupelize(arg):
    if isinstance(arg,str):
        return (arg,)
    try:
        return tuple(arg)
    except TypeError:
        return (arg,)

def copyGroups(ncin, groups, skip=None):
    for k in groups:
        if k not in _tupelize(skip):
            group = copyGroup(ncin,ncout,k)
            copyGroups(ncin.groups[k],group)

def copyGroup(ncin, group):
    out = ncout.createGroup(group.name)
    group.copyDimensions(out)
    group.copyVariables(out)
    group.copyAttributes(out)
    return out
 
def copyDimension(ncin, dim):
    return ncin.createDimension(dim.name, None if dim.isunlimited() else len(dim))

def copyDimensions(ncin, dimensions, skip=None):
    for k,v in dimensions.items():
        if k not in _tupelize(skip):
            copyDimension(ncin, v)

def copyAttributes(ncin, attributes, skip=None):
    for k,v in attributes.items():
        if k not in _tupelize(skip):
            if k == "missing_value":
                v = ncin.dtype.type(v)
            ncin.createAttribute(k, v)
            
def copyGroup(ncout, ncin, skipdims=None, skipgroups=None, skipvars=None, skipattrs=None):
    ncout.set_fill_off()
    ncout.copyDimensions(ncin.dimensions, skipdims)
    ncout.copyVariables(ncin.variables, skipvars)
    ncout.copyAttributes(ncin.attributes, skipattrs)
    ncout.copyGroups(ncin.groups, skipgroups)
            
def copyVariable(ncin, var, data=True):
    """
    Arguments
    ---------
    ncin : NcDataset in write/append mode
    var  : NcVariable
    data : boolean

    Purpose
    -------
    Creates variables in Datset and also write data if data=True
    """
    invardef = var.definition
    invar = ncin.createVariable(
        invardef.pop("name"),invardef.pop("dtype"),**invardef
    )
    copyAttributes(invar,var.attributes)
    if data:
        invar[:] = var[:] 
    return invar

def copyVariables(ncin, variables, skip=None, data=True):
    for k,v in variables.items():
        if k not in _tupelize(skip):
            copyVariable(ncin, v, data)

def getVariables(ncin):
    out = OrderedDict()
    for v in getattr(ncin, "variables").values():
        out[v.name] = NcVariable(ncin,v.name,v.dtype,v.dimensions,id=v._varid)
    return out

def getAttributes(ncin):
    out = OrderedDict()
    for k in ncin.ncattrs():
        if not k.startswith("_"):
            out[k] = ncin.getncattr(k)
    return out

def getVariableDefinition(ncvar):
    out = ncvar.filters() if ncvar.filters() else {}
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : ncvar.dimensions,
        "chunksizes" : ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None,
        "fill_value" : ncvar._FillValue if "_FillValue" in dir(ncvar) else None,
    })
    return out

def getDates(nc, timesteps=None, timevar="time", units=None, calendar=None):
    """
    timesteps:  time_steps to return datetimeobjects for. If not given datetime objects
                for the entire file will be returned.
    units:      time units ('time units since reference time'). Must be given as argument
                if not given in time variable attributes. AnAttributeError is risen 
                otherwise.
    calendar:   calendar string. If not given it is tried to read the calendar from file.
                If not set in the time variable attributes a standard calendar is assumed        
                Valid calendars are:
                  'standard', 'gregorian', 'proleptic_gregorian' 'noleap',
                  '365_day', '360_day', 'julian', 'all_leap', '366_day'
    """
    var = nc.variables[timevar]        
    if not units:
        try:
            units = var.units
        except AttributeError:
            raise AttributeError(
                "Time variable does not specify an units attribute! Pass as argument."
            )

    if not calendar:
        try:
            calendar = var.calendar
        except AttributeError:
            calendar = "standard"

    if not timesteps:
        timesteps = var[:]

    dates = num2date(timesteps,units,calendar)
    
    try:
        return [d.date() for d in dates]
    except AttributeError:
        return dates

def getFillValue(ncin):
    try:
        return ncin.getncattr("_FillValue")
    except AttributeError:
        return None
        
def setAttribute(ncin, name, value):
    ncin.setncattr(name, value)

def setAttributes(ncin, attdict):
    ncin.setncatts(attdict)
 
def filterVariables(ncin, dims=None, ndim=None):
    """
    Arguments
    ---------
    dims (optional) : tuple/list of dimension strings
    ndim (optional) : int number of dimensions

    Return
    ------
    tuple of strings
    
    Purpose
    -------
    Return all Variables that are based on the dimension(s) given in dims
    """
    out = OrderedDict()
    dims = set(dims or {})
    
    for v in ncin.variables.values():
        if dims.issubset(set(v.dimensions)):
            if ndim:
                if ndim == len(v.dimensions):
                    out[v.name] = v
            else:
                out[v.name] = v
                
    return out

def filterDimensions(ncin, lengths):
    out = OrderedDict()
    dims = {len(d) : d  for d in ncin.dimensions.values()}
    for l in lengths:
        try:
            d = dims[l]
            out[d.name] = d
        except KeyError:
            pass
    return out
    
class NcVariable(Variable):
    def __init__(self,*args,**kwargs):
        super(NcVariable,self).__init__(*args,**kwargs)

    copyAttributes   = copyAttributes
    createAttribute  = setAttribute
    createAttributes = setAttributes
    attributes       = property(fget=getAttributes)
    definition       = property(fget=getVariableDefinition)
    fill_value       = property(fget=getFillValue)
    
class NcGroup(Group):
    def __init__(self,*args,**kwargs):
        super(NcGroup,self).__init__(*args,**kwargs)
        for k,v in zip( self.variables, getVariables(self).values() ):
            self.variables[k] = v 
      
    copyAttributes   = copyAttributes
    copyGroup        = copyGroup
    copyGroups       = copyGroups
    createAttribute  = setAttribute
    createAttributes = setAttributes
    filterVariables  = filterVariables
    filterDimensions = filterDimensions
    attributes       = property(fget=getAttributes)
    
class NcDataset(Dataset):
    def __init__(self,*args,**kwargs):
        super(NcDataset, self).__init__(*args,**kwargs)
        for k,v in zip( self.variables, getVariables(self).values() ):
            self.variables[k] = v 
        
    def createVariable(self,*args,**kwargs):
        var = NcVariable(self,*args,**kwargs)
        self.variables[var.name] = var
        return var
    
    def createGroup(self,name):
        grp = NcGroup(self,name)
        self.groups[name] = grp
        return grp

    def __enter__(self):
        return self
    
    def __exit__(self,*args,**kwargs):
        self.close()
            
    copyDataset      = copyGroup
    copyDimension    = copyDimension
    copyDimensions   = copyDimensions
    copyAttributes   = copyAttributes
    copyFile         = copyGroup
    copyVariable     = copyVariable
    copyVariables    = copyVariables
    copyGroup        = copyGroup
    copyGroups       = copyGroups
    createAttribute  = setAttribute
    createAttributes = setAttributes
    filterVariables  = filterVariables
    filterDimensions = filterDimensions
    getDates         = getDates
    attributes       = property(fget=getAttributes)

        
if __name__ == "__main__":
    pass
