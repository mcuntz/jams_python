#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author
------
David Schaefer

Purpose
-------
A sanitizing layer for the netCDF4 library. Adds a number of convenince methods
and aims for a cleaner user interface. All classes avaliable are children of their
associated netCDF4 counterparts.
"""

from netCDF4 import Dataset, Group, Dimension, Variable, date2num, num2date
from collections import OrderedDict

def _tupelize(arg):
    if isinstance(arg,str):
        return (arg,)
    try:
        return tuple(arg)
    except TypeError:
        return (arg,)

def copyGroup(ncin, group, skipdims=None, skipgroups=None, skipvars=None, skipattrs=None):
    """
    Arguments
    ---------
    ncin                  : Instance of an object with a createGroup method
                            (i.e. NcDataset, NcGroup)
    group                 : Instance of an object with dimensions/variables/attributes/groups attributes
                            (i.e. NcDataset, NcGroup)
    skipdims (optional)   : string or list/tuple of strings
                            Name(s) of dimension(s) to skip
    skipgroups (optional) : string or list/tuple of strings
                            Name(s) of group(s) to skip
    skipvars (optinal)    : string or list/tuple of strings
                            Name(s) of variable(s) to skip
    skipattrs (optinal)   : string or list/tuple of strings
                            Name(s) of attribute(s) to skip

    Return
    ------
    NcGroup

    Purpose
    -------
    Copy the given group to ncin
    """
    out = ncout.createGroup(group.name)
    out.set_fill_off()
    out.copyDimensions(group.dimensions, skipdims)
    out.copyVariables(group.variables, skipvars)
    out.copyAttributes(group.attributes, skipattrs)
    out.copyGroups(group.groups, skipgroups)
    return out

def copyDataset(ncin, group, skipdims=None, skipgroups=None, skipvars=None, skipattrs=None):
    """
    Arguments
    ---------
    ncin                  : Instance of an object with a createGroup method
                            (i.e. NcDataset, NcGroup)
    group                 : Instance of an object with dimensions/variables/attributes/groups attributes
                            (i.e. NcDataset, NcGroup)
    skipdims (optional)   : string or list/tuple of strings
                            Name(s) of dimension(s) to skip
    skipgroups (optional) : string or list/tuple of strings
                            Name(s) of group(s) to skip
    skipvars (optinal)    : string or list/tuple of strings
                            Name(s) of variable(s) to skip
    skipattrs (optinal)   : string or list/tuple of strings
                            Name(s) of attribute(s) to skip

    Return
    ------
    NcDataset/NcGroup

    Purpose
    -------
    Copy the content of given group to ncin
    """
    ncin.set_fill_off()
    ncin.copyDimensions(group.dimensions, skipdims)
    ncin.copyVariables(group.variables, skipvars)
    ncin.copyAttributes(group.attributes, skipattrs)
    ncin.copyGroups(group.groups, skipgroups)
    return ncin
 
def copyGroups(ncin, groups, skip=None):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createGroup method
                      (i.e. NcDataset, NcGroup)
    groups          : Dictionary
                      key   : group name (string)
                      value : instance of an object with dimensions/variables/attributes/groups attributes
    skip (optional) : string or list/tuple of strings
                      Name(s) of group(s) to skip
    
    Return
    ------
    None

    Purpose
    -------
    Copy the given groups to ncin
    """
    for g in groups.values():
        if g.name not in _tupelize(skip):
            ncin.copyGroup(g)

def copyDimension(ncin, dim):
    """
    Arguments
    ---------
    ncin  : Instance of an object with a createDimension method
            (i.e. NcDataset, NcGroup)
    group : Instance of NcDimension 
    
    Return
    ------
    netCDF4.Dimension

    Purpose
    -------
    Copy the given dimension to ncin
    """
    return ncin.createDimension(dim.name, None if dim.isunlimited() else len(dim))

def copyDimensions(ncin, dimensions, skip=None):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createDimension method
                      (i.e. NcDataset, NcGroup)
    dimension       : Dictionary
                      key   : dimension name (string)
                      value : instance of NcDimension
    skip (optional) : string or list/tuple of strings
                      Name(s) of dimension(s) to skip
    
    Return
    ------
    None

    Purpose
    -------
    Copy the given dimensions to ncin
    """
    for d in dimensions.values():
        if d.name not in _tupelize(skip):
            ncin.copyDimension(d)

def copyAttributes(ncin, attributes, skip=None):
    """
    Arguments
    ---------
    ncin              : Instance of an object with a createAttribute method
                        (i.e. NcDataset, NcGroup, NcVariable)
    attributes        : Dictionary
                        key   : string
                        value : string/any numeric type
    skip (optional)   : string or list/tuple of strings
                        Name(s) of attribute(s) to skip

    Return
    ------
    None

    Purpose
    -------
    Copy the given attributes to ncin
    """
    for k, v in attributes.items():
        if k not in _tupelize(skip):
            if k == "missing_value":
                v = ncin.dtype.type(v)
            ncin.createAttribute(k, v)
            
           
def copyVariable(ncin, var, data=True):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createCopy method
                      (i.e. NcDataset, NcGroup, NcVariable)
    var             : Instance of NcVariable
    data (optional) : boolean

    Return
    ------
    NcVariable
    
    Purpose
    -------
    Copy the given variables to ncin. Copy the data if data=True
    """
    invardef = var.definition
    invar = ncin.createVariable(
        invardef.pop("name"),invardef.pop("dtype"),**invardef
    )
    invar.copyAttributes(var.attributes)
    if data:
        invar[:] = var[:] 
    return invar

def copyVariables(ncin, variables, skip=None, data=True):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createCopy method
                      (i.e. NcDataset, NcGroup, NcVariable)
    variables       : Dictionary
                      key   : variables name (string)
                      value : instance of NcVariable
    skip (optional) : string or list/tuple of strings
                      Name(s) of variable(s) to skip
    data (optional) : boolean

    Return
    ------
    NcVariable
    
    Purpose
    -------
    Copy the given variables to ncin. Copy the data if data=True
    """
    for v in variables.values():
        if v.name not in _tupelize(skip):
            ncin.copyVariable(v, data)

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

def getDates(ncin, timesteps=None, timevar="time", units=None, calendar=None):
    """
    Arguments
    ---------
    ncin                 : Instance of an object holding variables (NcDataset/NcGroup)
    timesteps (optional) : list/tuple/nd.array of Numerical values.
                           The time_steps to return dates for. If not given the content of the
                           entire time variable will be returned.
    units (optional)     : string
                           time units following the CF Conventions. Needs to be given, if not
                           available as an attribute of the time variable.
    calendar (optional)  : string
                           calendar name following the CF conventions. Needs to be given, if not
                           available as an attribute of the time variable.

    Return
    ------
    List of datetime objects

    Purpose
    -------
    Return datetime objects associated to the time variable of ncin
    """
    var = ncin.variables[timevar]        
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
    """
    Arguments
    ---------
    ncin  : Instance of an object with a _FillValue attribute
            (i.e. NcVariable)
    
    Return
    ------
    Numeric

    Purpose
    -------
    Return the value of the attribute _FillValue
    """
    try:
        return ncin.getncattr("_FillValue")
    except AttributeError:
        return None
        
def setAttribute(ncin, name, value):
    """
    Arguments
    ---------
    ncin  : Instance of an object with a setncatts method
            (i.e. NcDataset/NcGroup/NcVariable)
    name  : string
    value : string or any numeric type
    
    Return
    ------
    None

    Purpose
    -------
    Set/Write the attribute given as name, value
    """
    ncin.setncattr(name, value)

def setAttributes(ncin, attdict):
    """
    Arguments
    ---------
    ncin    : Instance of an object with a setncatts method
              (i.e. NcDataset/NcGroup/NcVariable)
    attdict : dictionary
              key: attribute name (string)
              value: attribute value (string or any numeric type)

    Return
    ------
    None

    Purpose
    -------
    Set/Write the attributes given in attdict
    """
    ncin.setncatts(attdict)
 
def filterVariables(ncin, dims=None, ndim=None):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a variables attribute
                      (i.e. NcDataset/NcGroup)
    dims (optional) : tuple/list of dimension strings
    ndim (optional) : int number of dimensions

    Return
    ------
    OrderedDict:
        key   : variable name (string)
        value : NcVariable instance
    
    Purpose
    -------
    Return all Variables that are based on the dimension(s) given in dims
    and/or have ndims dimensions.
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
    """
    Arguments
    ---------
    lengths : tuple/list of integers

    Return
    ------
    OrderedDict:
        key   : dimension name (string)
        value : NcDimension instance

    Purpose
    -------
    Return all Dimensions with a length given in the argument lengths.
    """
    out = OrderedDict()
    dims = {len(d) : d  for d in ncin.dimensions.values()}
    for l in lengths:
        try:
            d = dims[l]
            out[d.name] = d
        except KeyError:
            pass
    return out

def getGroups(ncin):
    out = OrderedDict()
    for g in getattr(ncin, "groups").values():
        out[g.name] = NcGroup(ncin, g.name, id=g._grpid)
    return out

def getVariables(ncin):
    out = OrderedDict()
    for v in getattr(ncin, "variables").values():
        out[v.name] = NcVariable(ncin, v.name, v.dtype, v.dimensions, id=v._varid)
    return out

def getAttributes(ncin):
    out = OrderedDict()
    for k in ncin.ncattrs():
        if not k.startswith("_"):
            out[k] = ncin.getncattr(k)
    return out

  
class NcDataset(Dataset):
    def __init__(self, *args, **kwargs):
        super(NcDataset, self).__init__(*args, **kwargs)
        for k, v in zip(self.groups, getGroups(self).values()):
            self.groups[k] = v 
        for k, v in zip(self.variables, getVariables(self).values()):
            self.variables[k] = v 
         
    def createGroup(self, name):
        grp = NcGroup(self, name)
        self.groups[name] = grp
        return grp

    def createVariable(self, *args, **kwargs):
        var = NcVariable(self, *args, **kwargs)
        self.variables[var.name] = var
        return var
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()
            
    copyDataset      = copyDataset
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

class NcGroup(Group):
    def __init__(self, *args, **kwargs):
        super(NcGroup,self).__init__(*args, **kwargs)
        for k, v in zip(self.groups, getGroups(self).values()):
            self.groups[k] = v 
        for k,v in zip(self.variables, getVariables(self).values()):
            self.variables[k] = v 
      
    def createGroup(self, name):
        grp = NcGroup(self, name)
        self.groups[name] = grp
        return grp

    def createVariable(self, *args, **kwargs):
        var = NcVariable(self, *args, **kwargs)
        self.variables[var.name] = var
        return var
    
    copyAttributes   = copyAttributes
    copyGroup        = copyGroup
    copyGroups       = copyGroups
    createAttribute  = setAttribute
    createAttributes = setAttributes
    filterVariables  = filterVariables
    filterDimensions = filterDimensions
    attributes       = property(fget=getAttributes)
         
class NcVariable(Variable):
    def __init__(self,*args,**kwargs):
        super(NcVariable,self).__init__(*args,**kwargs)

    copyAttributes   = copyAttributes
    createAttribute  = setAttribute
    createAttributes = setAttributes
    attributes       = property(fget=getAttributes)
    definition       = property(fget=getVariableDefinition)
    fill_value       = property(fget=getFillValue)

# Just to be consistent...
NcDimension = Dimension

if __name__ == "__main__":
    pass
