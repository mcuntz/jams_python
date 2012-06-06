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
    abc2plot           Write a, b, c, ... on plots
    around             Round to the passed power of ten.
    astr               Wrapper for autostring.
    autostring         Format number (array) with given decimal precision.
    cellarea           Calc areas of grid cells in m^2
    closest            Get the array index of the element that is closest to a given number.
    const              Provides physical, mathematical, computational, and isotope constants.
    cuntz_gleixner     Cuntz-Gleixner model of 13C discrimination
    dewpoint           Calculates the dew point from ambient humidity
    date2dec	       Converts arrays with calendar date to decimal date
    dec2date	       Converts arrays with decimal date to calendar date
    div                Wrapper for division.
    division           Divide two arrays, return "otherwise" if division by 0.
    esat               Calculates the saturation vapour pressure of water/ice.
    fread              Reads in float array from ascii file
    gapfill            Gapfill Eddy flux data
    gap_filling        Gapfills eddy flux data (CO2, LE, H) 
    heaviside          Heaviside (or unit step) operator
    int2roman          Integer to roman numeral conversion
    lhs                Latin Hypercube Sampling of any distribution without correlations
    lif                Count number of lines in file
    mad                Median absolute deviation test
    nee2gpp            Photosynthesis and ecosystem respiration NEE Eddy flux data
    outlier            Rossner''s extreme standardized deviate outlier test
    pack               Similar to Fortran pack function with mask
    pi                 Parameter importance index PI or alternatively B index calculation
    position           Position arrays of subplots to be used with add_axes
    readnetcdf         Reads variables or information from netcdf file
    roman2int          Roman numeral to integer conversion
    saltelli           Parameter sampling for Sobol indices calculation
    semivariogram      Calculates semivariogram from spatial data
    sobol_index        Calculates the first-order and total variance-based sensitivity indices
    sread              Reads in string array from ascii file
    tcherkez           Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.
    tsym               Raw unicodes for common symbols
    unpack             Similar to Fortran unpack function with mask
    writenetcdf        writing variables and information to a NetCDF file
    yrange             Calculates plot range from input array
    zacharias          Soil water content with van Genuchten and Zacharias et al. (2007)
    zacharias_check    Checks validity of parameter set for Zacharias et al. (2007)

    
    Provided functions per category
    -------------------------------
        Array manipulation
        Ascii files
        Data processing
        Date & Time
        Grids
        Isotopes
        Math
        Meteorology
        Miscellaneous
        Plotting
        Special files
        Obsolete
    -------------------------------

    Array manipulation
    ------------------
    closest        Get the array index of the element that is closest to a given number.
    pack           Similar to Fortran pack function with mask
    unpack         Similar to Fortran unpack function with mask


    Ascii files
    -----------
    fread          Reads in float array from ascii file
    lif            Count number of lines in file
    sread          Reads in string array from ascii file


    Data processing
    ---------------
    gapfill            Gapfill Eddy flux data
    gap_filling        Gapfills flux data (CO2, LE, H)
    mad                Median absolute deviation test
    nee2gpp            Photosynthesis and ecosystem respiration NEE Eddy flux data
    outlier            Rossner''s extreme standardized deviate outlier test
    semivariogram      Calculates semivariogram from spatial data


    Date & Time
    -----------
    date2dec	       Converts arrays with calendar date to decimal date
    dec2date	       Converts arrays with decimal date to calendar date


    Grids
    -----
    cellarea           Calc areas of grid cells in m^2


    Isotopes
    --------
    cuntz_gleixner     Cuntz-Gleixner model of 13C discrimination
    tcherkez           Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.


    Math
    ----
    around             Round to the passed power of ten.
    div                Wrapper for division.
    division           Divide two arrays, return "otherwise" if division by 0.
    heaviside          Heaviside (or unit step) operator
    lhs                Latin Hypercube Sampling of any distribution without correlations
    pi                 Parameter importance index PI or alternatively B index calculation
    saltelli           Parameter sampling for Sobol indices calculation
    sobol_index        Calculates the first-order and total variance-based sensitivity indices


    Meteorology
    -----------
    dewpoint           Calculates the dew point from ambient humidity
    esat               Calculates the saturation vapour pressure of water/ice.


    Miscellaneous
    -------------
    astr               Wrapper for autostring
    autostring         Format number (array) with given decimal precision.
    const              Provides physical, mathematical, computational, and isotope constants.
    int2roman          Integer to roman numeral conversion
    roman2int          Roman numeral to integer conversion
    zacharias          Soil water content with van Genuchten and Zacharias et al. (2007)
    zacharias_check    Checks validity of parameter set for Zacharias et al. (2007)


    Plotting
    --------
    abc2plot           Write a, b, c, ... on plots
    position           Position arrays of subplots to be used with add_axes
    tsym               Raw unicodes for common symbols
    yrange             Calculates plot range from input array


    Special files
    -------------
    readnetcdf         Reads variables or information from netcdf file
    writenetcdf        writing variables and information to a NetCDF file

    Obsolete
    --------
    calcvpd            Calculates vapour pressure deficit
        

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
	      MC, Jan 2012 - esat, closest, dewpoint, division, heaviside, tcherkez, yrange, const
                           - make calcvpd obsolete
                           - cuntz_gleixner
              MC, Mar 2012 - gapfill, nee2gpp
              MC, May 2012 - astr, div, sobol_index, pi, roman, zacharias, saltelli
              MZ, Jun 2012 - writenetcdf
"""
# Routines provided
from abc2plot        import *
from around          import *
from autostring      import *
from calcvpd         import *
from cellarea        import *
from closest         import *
import const
from cuntz_gleixner  import *
try:
    from date2dec    import *
    from dec2date    import *
except ImportError:
    print "No netcdf support in UFZ library. Disabled functions: date2dec, dec2date, gap_filling, readnetcdf."
from dewpoint        import *
from division        import *
from esat            import *
from fread           import *
from gapfill         import *
from nee2gpp         import *
try:
    from gap_filling import *
except ImportError:
    pass
from heaviside       import *
from lhs             import *
from lif             import *
from mad             import *
try:
    from outlier     import *
except ImportError:
    print "No extra statistics in scipy, i.e. in UFZ library. Disabled functions: outlier."
from pack            import *
from pi              import *
from position        import *
try:
    from readnetcdf  import *
except ImportError:
    print "No netcdf support in UFZ library. Disabled function: writenetcdf."    
from roman           import int2roman, roman2int
from saltelli        import *
from semivariogram   import *
from sobol_index     import *
from sread           import *
from tcherkez        import *
from tsym            import *
from unpack          import *
try:
    from writenetcdf import *
except ImportError:
    pass
from yrange          import *
from zacharias       import *

# Information
version = '1.6'
date = '25.05.2012'

# Main
if __name__ == '__main__':
    print '\nUFZ Computational Hydrosystems Python Utilities.'
    print "Version %s from %s." % (version,date)
    print ('\nUFZ routines are free software and come with '
           'ABSOLUTELY NO WARRANTY.')
    print 'You are welcome to redistribute it.'
    print ('\nCopyright (C) 2009-2012, Computational Hydrosystems, '
          'Helmholtz Centre for Environmental Research - UFZ, Permoserstr. 15, '
          '04318 Leipzig, Germany.')
    print 'All rights reserved.'
    print '\nIn case of questions or comments contact matthias.cuntz(at)ufz.de\n'
