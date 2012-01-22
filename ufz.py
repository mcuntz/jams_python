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

    Provided functions (alphabetic w/o obsolete)
    ------------------
    around         Round to the passed power of ten.
    autostring     Format number (array) with given decimal precision.
    cellarea       Calc areas of grid cells in m^2
    date2dec	   Converts arrays with calendar date to decimal date
    dec2date	   Converts arrays with decimal date to calendar date
    esat           Calculates the saturation vapour pressure of water/ice.
    fread          Reads in float array from ascii file
    gap_filling    Gapfills eddy flux data (CO2, LE, H) 
    lif            Count number of lines in file
    mad            Median absolute deviation test
    outlier        Rossner''s extreme standardized deviate outlier test
    pack           Similar to Fortran pack function with mask
    position       Position arrays of subplots to be used with add_axes
    readnetcdf     Reads variables or information from netcdf file
    semivariogram  Calculates semivariogram from spatial data
    sread          Reads in string array from ascii file
    tsym           Raw unicodes for common symbols
    unpack         Similar to Fortran unpack function with mask

    
    Provided functions per category
    -------------------------------
        Array manipulation
        Ascii files
        Data processing
        Date & Time
        Grids
        Miscellaneous
        Plotting
        Special files
        Obsolete
    -------------------------------
    Array manipulation
    ------------------
    pack           Similar to Fortran pack function with mask
    unpack         Similar to Fortran unpack function with mask

    Ascii files
    -----------
    fread          Reads in float array from ascii file
    lif            Count number of lines in file
    sread          Reads in string array from ascii file

    Data processing
    ---------------
    gap_filling    Gapfills flux data (CO2, LE, H)
    mad            Median absolute deviation test
    outlier        Rossner''s extreme standardized deviate outlier test
    semivariogram  Calculates semivariogram from spatial data

    Date & Time
    -----------
    date2dec	   Converts arrays with calendar date to decimal date
    dec2date	   Converts arrays with decimal date to calendar date

    Grids
    -----
    cellarea       Calc areas of grid cells in m^2

    Miscellaneous
    -------------
    around         Round to the passed power of ten.
    autostring     Format number (array) with given decimal precision.
    esat           Calculates the saturation vapour pressure of water/ice.

    Plotting
    --------
    position       Position arrays of subplots to be used with add_axes
    tsym           Raw unicodes for common symbols

    Special files
    -------------
    readnetcdf     Reads variables or information from netcdf file

    Obsolete
    -------------
    calcvpd        Calculates vapour pressure deficit
        

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
              MC, Nov 2011 - mad
	      MC, Nov 2011 - try netcdf and stats routines
	      MC, Nov 2011 - autostring
	      MC, Jan 2012 - esat, make calcvpd obsolete
"""
# Routines provided

from around      import *
from autostring  import *
from calcvpd     import *
from cellarea    import *
from esat        import *
try:
    from readnetcdf  import *
    from date2dec    import *
    from dec2date    import *
    from gap_filling import *
except ImportError:
    print "No netcdf support in UFZ library. Disabled functions: readnetcdf, date2dec, dec2date, gap_filling."
from fread       import *
from lif         import *
try:
    from outlier     import *
except ImportError:
    print "No extra statistics in scipy, i.e. in UFZ library. Disabled functions: outlier."
from mad         import *
from pack        import *
from position    import *
from sread       import *
from tsym        import *
from unpack      import *

# Information
version = '1.4'
date = '11.11.2011'

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
