#!/usr/bin/env python
"""
    UFZ Computational Hydrosystems Python Utilities
    Module offers miscellaneous functions in different categories

    Get help on each function by typing
    >>> help()
    help> ufz.function
    Or
    >>> import ufz
    >>> help(ufz.function)

    Provided functions (alphabetic)
    ------------------
    around         round to the passed power of ten.
    calcvpd        calculates vapour pressure deficit
    cellarea       calc areas of grid cells in m^2
    date2dec	   converts arrays with calendar date to decimal date
    dec2date	   converts arrays with decimal date to calendar date
    fread          reads in float array from ascii file
    gap_filling    gap fills eddy flux data (CO2,LE,H-fluxes) 
    lif            count number of lines in file
    outlier        Rossner''s extreme standardized deviate outlier test
    pack           similar to Fortran pack function with mask
    position       position arrays of subplots to be used with add_axes
    readnetcdf     reads variables or information from netcdf file
    semivariogram  calculates semivariogram from spatial data
    sread          reads in string array from ascii file
    tsym           Raw unicodes for common sybols
    unpack         similar to Fortran unpack function with mask

    
    Provided functions per category
        Array manipulation
        Ascii files
        Data processing
        Date & Time
        Grids
        Miscellaneous
        Plotting
        Special files
    -------------------------------
    Array manipulation
    ------------------
    pack           similar to Fortran pack function with mask
    unpack         similar to Fortran unpack function with mask

    Ascii files
    -----------
    fread          reads in float array from ascii file
    lif            count number of lines in file
    sread          reads in string array from ascii file

    Data processing
    ---------------
    calcvpd        calculates vapour pressure deficit
    gap_filling    gap fills flux data (CO2,LE,H-fluxes) 
    outlier        Rossner''s extreme standardized deviate outlier test
    semivariogram  calculates semivariogram from spatial data

    Date & Time
    -----------
    date2dec	   converts arrays with calendar date to decimal date
    dec2date	   converts arrays with decimal date to calendar date

    Grids
    -----
    cellarea       calc areas of grid cells in m^2

    Miscellaneous
    -------------
    around         round to the passed power of ten.

    Plotting
    --------
    position       position arrays of subplots to be used with add_axes
    tsym           Raw unicodes for common sybols

    Special files
    -------------
    readnetcdf     reads variables or information from netcdf file
        

    History
    -------
    Written,  MC, Jul 2009
    Modified, MC, Jul 2009 - lif, fread, sread
                           - readnetcdf, cellarea
                           - pack, unpack
              MC, Aug 2009 - position
              MG, Jul 2010 - outlier
	      AP, Jan 2011 - date2dec, dec2date
	      AP, Feb 2011 - semivariogram
              TR, May 2011 - gap_filling
              TR, May 2011 - calcvpd
              MC, Jun 2011 - /usr/bin/python to /usr/bin/env python
                           - tsym, around
"""
# Routines provided
from around      import *
from calcvpd     import *
from cellarea    import *
from date2dec    import *
from dec2date    import *
from fread       import *
from gap_filling import *
from lif         import *
from outlier     import *
from pack        import *
from position    import *
from readnetcdf  import *
from sread       import *
from tsym        import *
from unpack      import *

# Information
version = '1.3'
date = '24.03.2011'

# Main
if __name__ == '__main__':
    print 'UFZ Computational Hydrosystems Python Utilities.'
    print "Version %s from %s." % (version,date)
    print ('\nUFZ routines are free software and come with '
           'ABSOLUTELY NO WARRANTY.')
    print 'You are welcome to redistribute it.'
    print ('\nCopyright (C) 2009-2011, Computational Hydrosystems, '
          'Helmholtz Centre for Environmental Research - UFZ, Permoserstr. 15, '
          '04318 Leipzig, Germany.')
    print 'All rights reserved.'
