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
    abc2plot               Write a, b, c, ... on plots.
    area_poly              Area of a polygon
    around                 Round to the passed power of ten.
    astr                   Wrapper for autostring.
    autostring             Format number (array) with given decimal precision.
    cellarea               Calc areas of grid cells in m^2.
    closest                Get the array index of the element that is closest to a given number.
    colors                 Wrapper for colour.
    colours                Define UFZ colours.
    const                  Provides physical, mathematical, computational, and isotope constants.
    convex_hull            Calculate subset of points that make a convex hull around a set of 2D points. 
    cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
    define_brewer          Defines and registers Brewer colormap.
    dewpoint               Calculates the dew point from ambient humidity.
    date2dec               Converts arrays with calendar date to decimal date.
    dec2date               Converts arrays with decimal date to calendar date.
    div                    Wrapper for division.
    division               Divide two arrays, return 'otherwise' if division by 0.
    esat                   Calculates the saturation vapour pressure of water/ice.
    fill_nonfinite         Fill missing values by linear interpolation.
    find_in_path           Look for file in system path.
    fread                  Reads in float array from ascii file.
    functions              Common functions that are used in curve_fit or fmin parameter estimations.
    gapfill                Gapfill Eddy flux data.
    get_brewer             Defines and returns Brewer colormap.
    heaviside              Heaviside (or unit step) operator.
    in_poly                Determines whether a 2D point falls in a polygon.
    int2roman              Integer to roman numeral conversion.
    interpol               One-dimensional linear interpolation on first dimension.
    jab                    Jackknife-after-Bootstrap error.
    kernel_regression      Multi-dimensional non-parametric regression.
    kernel_regression_h    Optimal bandwidth for kernel regression.
    kriging                Krig a surface from a set of 2D points.
    lhs                    Latin Hypercube Sampling of any distribution without correlations.
    lif                    Count number of lines in file.
    mad                    Median absolute deviation test.
    means                  Calculate daily, monthly, yearly, etc. means of data depending on date stamp.
    nee2gpp                Photosynthesis and ecosystem respiration from NEE Eddy flux data.
    nee2gpp_global         nee2gpp using one fit for whole time period
    nee2gpp_lasslop        nee2gpp using the daytime method of Lasslop et al. (2010)
    nee2gpp_reichstein     nee2gpp using several fits as in Reichstein et al. (2005)
    outlier                Rossner''s extreme standardized deviate outlier test.
    pack                   Similar to Fortran pack function with mask.
    pi                     Parameter importance index PI or alternatively B index calculation.
    plot_brewer            Plots available Brewer color maps in pdf file.
    position               Position arrays of subplots to be used with add_axes.
    print_brewer           Prints available Brewer colormap names.
    readhdf                Reads variables or information from hdf4 and hdf5 files.
    readhdf4               Reads variables or information from hdf4 files.
    readhdf5               Reads variables or information from hdf5 file.
    readnetcdf             Reads variables or information from netcdf file.
    rgb                    Interpolate between colours; make continuous colour maps.
    roman2int              Roman numeral to integer conversion.
    savitzky_golay         Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter.
    savitzky_golay2d       Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter.
    saltelli               Parameter sampling for Sobol indices calculation.
    semivariogram          Calculates semivariogram from spatial data.
    sg                     Wrapper savitzky_golay.
    sg2d                   Wrapper savitzky_golay2d.
    sobol_index            Calculates the first-order and total variance-based sensitivity indices.
    sread                  Reads in string array from ascii file.
    srrasa                 Generates stratified random 2D points within a given rectangular area.
    srrasa_trans           Generates stratified random 2D transects within a given rectangular area.
    tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.
    tsym                   Raw unicodes for common symbols.
    unpack                 Similar to Fortran unpack function with mask.
    volume_poly            Volume of function above a polygon
    writenetcdf            Write netCDF4 file.
    xkcd                   Make plot look handdrawn.
    yrange                 Calculates plot range from input array.
    zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
    zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).

    
    Provided functions per category
    -------------------------------
        Array manipulation
        Ascii files
        Data processing
        Date & Time
        Grids / Polygons
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
    closest                Get the array index of the element that is closest to a given number.
    pack                   Similar to Fortran pack function with mask.
    unpack                 Similar to Fortran unpack function with mask.


    Ascii files
    -----------
    fread                  Reads in float array from ascii file.
    lif                    Count number of lines in file.
    sread                  Reads in string array from ascii file.


    Data processing
    ---------------
    convex_hull            Calculate subset of points that make a convex hull around a set of 2D points.
    fill_nonfinite         Fill missing values by linear interpolation.
    gapfill                Gapfill Eddy flux data.
    interpol               One-dimensional linear interpolation on first dimension.
    kriging                Krig a surface from a set of 2D points.
    kernel_regression      Multi-dimensional non-parametric regression.
    kernel_regression_h    Optimal bandwidth for kernel regression.
    mad                    Median absolute deviation test.
    means                  Calculate daily, monthly, yearly, etc. means of data depending on date stamp.
    nee2gpp                Photosynthesis and ecosystem respiration from NEE Eddy flux data.
    nee2gpp_global         nee2gpp using one fit for whole time period
    nee2gpp_lasslop        nee2gpp using the daytime method of Lasslop et al. (2010)
    nee2gpp_reichstein     nee2gpp using several fits as in Reichstein et al. (2005)
    outlier                Rossner''s extreme standardized deviate outlier test.
    savitzky_golay         Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter.
    savitzky_golay2d       Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter.
    sg                     Wrapper savitzky_golay.
    sg2d                   Wrapper savitzky_golay2d.
    semivariogram          Calculates semivariogram from spatial data.
    srrasa                 Generates stratified random 2D points within a given rectangular area.
    srrasa_trans           Generates stratified random 2D transects within a given rectangular area.


    Date & Time
    -----------
    date2dec               Converts arrays with calendar date to decimal date.
    dec2date               Converts arrays with decimal date to calendar date.


    Grids / Polygons
    ----------------
    area_poly              Area of a polygon
    cellarea               Calc areas of grid cells in m^2.
    in_poly                Determines whether a 2D point falls in a polygon.
    volume_poly            Volume of function above a polygon


    Isotopes
    --------
    cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
    tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.


    Math
    ----
    around                 Round to the passed power of ten.
    div                    Wrapper for division.
    division               Divide two arrays, return 'otherwise' if division by 0.
    functions              Common functions that are used in curve_fit or fmin parameter estimations.
    heaviside              Heaviside (or unit step) operator.
    jab                    Jackknife-after-Bootstrap error.
    lhs                    Latin Hypercube Sampling of any distribution without correlations.
    pi                     Parameter importance index PI or alternatively B index calculation.
    saltelli               Parameter sampling for Sobol indices calculation.
    sobol_index            Calculates the first-order and total variance-based sensitivity indices.


    Meteorology
    -----------
    dewpoint               Calculates the dew point from ambient humidity.
    esat                   Calculates the saturation vapour pressure of water/ice.


    Miscellaneous
    -------------
    astr                   Wrapper for autostring.
    autostring             Format number (array) with given decimal precision.
    const                  Provides physical, mathematical, computational, and isotope constants.
    find_in_path           Look for file in system path.
    int2roman              Integer to roman numeral conversion.
    roman2int              Roman numeral to integer conversion.
    zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
    zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).


    Plotting
    --------
    abc2plot               Write a, b, c, ... on plots.
    colors                 Wrapper for colour.
    colours                Define UFZ colours.
    define_brewer          Defines and registers Brewer colormap.
    get_brewer             Defines and returns Brewer colormap.
    plot_brewer            Plots available Brewer color maps in pdf file.
    position               Position arrays of subplots to be used with add_axes.
    print_brewer           Prints available Brewer colormap names.
    rgb                    Interpolate between colours; make continuous colour maps.
    tsym                   Raw unicodes for common symbols.
    xkcd                   Make plot look handdrawn.
    yrange                 Calculates plot range from input array.


    Special files
    -------------
    readhdf                Reads variables or information from hdf4 and hdf5 files.
    readhdf4               Reads variables or information from hdf4 files.
    readhdf5               Reads variables or information from hdf5 file.
    readnetcdf             Reads variables or information from netcdf file.
    writenetcdf            Write netCDF4 file.


    License
    -------
    This file is part of the UFZ Python library.

    The UFZ Python library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The UFZ Python library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2009-2013 Matthias Cuntz, Arndt Piayda, Matthias Zink, Tino Rau, Maren Goehler,
                        Stephan Thober, Juliane Mai


    History
    -------
    Written,  MC, Jul 2009
    Modified, MC, Jul 2009 - lif, fread, sread, readnetcdf, cellarea, pack, unpack
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
              MC, Jan 2012 - esat, closest, dewpoint, division, heaviside, tcherkez, yrange, const, cuntz_gleixner
                           - calcvpd obsolete
              MC, Mar 2012 - gapfill, nee2gpp
              MC, May 2012 - astr, div, sobol_index, pi, roman, zacharias, saltelli
              MZ, Jun 2012 - writenetcdf
              MC, Jun 2012 - roman -> romanliterals, interpol
              MZ, Jun 2012 - readhdf5
              MC, Jun 2012 - readhdf4, readhdf
              MC, Sep 2012 - brewer
              MC, Oct 2012 - savitzky_golay
              MC, Nov 2012 - added netcdftime but no import so that available w/o netcdf
              AP, Nov 2012 - convex_hull, in_poly, kriging, semivariogram update, srrasa, srrasa_trans
              MC, Nov 2012 - nee2gpp, nee2gpp_global, nee2gpp_lasslop, nee2gpp_reichstein
              MC, Dec 2012 - functions
                           - gap_filling obsolete
              MC, Feb 2013 - area_poly
              MC & JM, Feb 2013 - volume_poly
              MC, Feb 2013 - ported to Python 3
              MC, Mar 2013 - find_in_path, xkcd
              MC, Apr 2013 - rgb
              MC, Jun 2013 - colours
              MC, Jul 2013 - fill_nonfinite, means
"""
from __future__ import print_function

# Routines provided
from abc2plot          import *
from area_poly         import area_poly
from around            import *
from autostring        import *
from brewer            import define_brewer, get_brewer, plot_brewer, print_brewer
try:
    from calcvpd       import *
except ImportError:
    pass
from cellarea          import *
from closest           import *
from colours           import colours, colors
import const
from convex_hull       import convex_hull 
from cuntz_gleixner    import *
from date2dec          import *
from dec2date          import *
from dewpoint          import *
from division          import *
from esat              import *
from fill_nonfinite    import fill_nonfinite
from find_in_path      import *
from fread             import *
import functions
from gapfill           import *
try:
    from gap_filling   import *
except ImportError:
    pass
from heaviside         import *
from in_poly           import in_poly 
from interpol          import *
from kernel_regression import kernel_regression, kernel_regression_h
from kriging           import kriging 
from lhs               import *
from lif               import *
from jab               import *
from mad               import *
from means             import *
from nee2gpp           import nee2gpp, nee2gpp_global, nee2gpp_lasslop, nee2gpp_reichstein
try:
    from outlier       import *
except ImportError:
    print("No extra statistics in scipy, i.e. in UFZ library. Disabled functions: outlier.")
from pack              import *
from pi                import *
from position          import *
try:
    from readhdf       import *
except ImportError:
    print("No hdf4 and/or hdf5 support in UFZ library. Disabled functions: readhdf.")
try:
    from readhdf4      import *
except ImportError:
    print("No hdf4 support in UFZ library. Disabled functions: readhdf4.")
try:
    from readhdf5      import *
except ImportError:
    print("No hdf5 support in UFZ library. Disabled functions: readhdf5.")
try:
    from readnetcdf    import *
except ImportError:
    print("No netcdf support in UFZ library. Disabled functions: readnetcdf, writenetcdf.")
from rgb               import rgb_blend, rgb_range, rgb_gradient
from romanliterals     import int2roman, roman2int
from saltelli          import *
from savitzky_golay    import *
from semivariogram     import semivariogram
from sobol_index       import *
from sread             import *
from srrasa            import srrasa, srrasa_trans 
from tcherkez          import *
from tsym              import *
from unpack            import *
from volume_poly       import volume_poly
try:
    from writenetcdf   import *
except ImportError:
    pass
from xkcd              import xkcd
from yrange            import *
from zacharias         import *

# Information
version = '2.0'
date    = '25.02.2013'

# Main
if __name__ == '__main__':
    print('\nUFZ Computational Hydrosystems Python Library.')
    print("Version %s from %s." % (version,date))
    print('\nThis is the README file. See als the license file gpl.txt.')
    import io
    f = io.open('README','r')
    for line in f:
        print(line,end='')
    f.close()
