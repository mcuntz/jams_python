#!/usr/bin/env python
"""
    UFZ Computational Hydrosystems Python Utilities

    Package offers miscellaneous functions and sub-modules in different categories.

    Get help on each function by typing
    >>> import ufz
    >>> help(ufz.function)

    Provided functions and modules (alphabetic w/o obsolete functions)
    ------------------------------------------------------------------
    abc2plot               Write a, b, c, ... on plots.
    alpha_equ_h2o          Equilibrium fractionation between liquid water and vapour.
    alpha_kin_h2o          Kinetic fractionation of molecular diffusion of water vapour.
    area_poly              Area of a polygon.
    argsort                Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables.
    around                 Round to the passed power of ten.
    astr                   Wrapper for autostring.
    autostring             Format number (array) with given decimal precision.
    cellarea               Calc areas of grid cells in m^2.
    clockplot              The clockplot of mHM.
    closest                Index in array which entry is closest to a given number.
    colors                 Wrapper for colour.
    colours                Define UFZ colours.
    const                  Provides physical, mathematical, computational, and isotope constants.
    convex_hull            Calculate subset of points that make a convex hull around a set of 2D points.
    correlate              Computes the cross-correlation function of two series x and y.
    cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
    directory_from_gui     Open directory selection dialog, returns selected directory
    dewpoint               Calculates the dew point from ambient humidity.
    date2dec               Converts arrays with calendar date to decimal date.
    dec2date               Converts arrays with decimal date to calendar date.
    div                    Wrapper for division.
    division               Divide two arrays, return 'otherwise' if division by 0.
    dumpnetcdf             Convenience function for writenetcdf
    eddybox                Module containing Eddy Covaraince utilities, see eddysuite.py for details
    eddysuite              Example file for processing Eddy data with eddybox and EddySoft
    elementary_effects     Morris measures mu, stddev and mu*
    ellipse_area           Area of ellipse (or circle)
    encrypt                Module to encrypt and decrypt text using a key system as well as a cipher.
    errormeasures          Definition of different error measures.
    esat                   Calculates the saturation vapour pressure of water/ice.
    file                   Module with file list function.
    file_from_gui          Open file selection dialog for one single file, returns selected files
    files_from_gui         Open file selection dialog, returns selected files
    fill_nonfinite         Fill missing values by linear interpolation.
    find_in_path           Look for file in system path.
    fread                  Reads in float array from ascii file.
    fsread                 Simultaneous read of float and string array from ascii file.
    ftp                    Module with functions for interacting with an open FTP connection.
    functions              Module with common functions that are used in curve_fit or fmin parameter estimations.
    gap2lai                Calculation of leaf area index from gap probability observations.
    get_angle              Returns the angle in radiant from each point in xy1 to each point in xy2.
    get_brewer             Registers and returns Brewer colormap.
    get_nearest            Returns a value z for each point in xy near to the xyz field.
    grid_mid2edge          Longitude and latitude grid edges from grid midpoints.
    hdfread                Wrapper for readhdf.
    hdf4read               Wrapper for readhdf4.
    hdf5read               Wrapper for readhdf5.
    head                   Return list with first n lines of file.
    heaviside              Heaviside (or unit step) operator.
    homo_sampling          Generation of homogeneous, randomly distributed points in a given rectangular area.
    in_poly                Determines whether a 2D point falls in a polygon.
    inpoly                 Wrapper for in_poly.
    int2roman              Integer to roman numeral conversion.
    interpol               One-dimensional linear interpolation on first dimension.
    jab                    Jackknife-after-Bootstrap error.
    kernel_regression      Multi-dimensional non-parametric regression.
    kernel_regression_h    Optimal bandwidth for kernel regression.
    kriging                Krig a surface from a set of 2D points.
    lagcorr                Calculate time lag of maximum or minimum correlation of two arrays.
    leafmodel              Model to compute photosynthesis and stomatal conductance of canopies.
    leafprojection         Calculation of leaf projection from leaf angle observations.
    lhs                    Latin Hypercube Sampling of any distribution without correlations.
    lif                    Count number of lines in file.
    line_dev_mask          Maskes elements of an array deviating from a line fit.
    logtools               Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle.
    mad                    Median absolute deviation test.
    means                  Calculate daily, monthly, yearly, etc. means of data depending on date stamp.
    morris_sampling        Sampling of optimised trajectories for Morris measures / elementary effects
    ncread                 Wrapper for readnetcdf.
    netcdfread             Wrapper for readnetcdf.
    outlier                Rossner''s extreme standardized deviate outlier test.
    pack                   Similar to Fortran pack function with mask.
    pca                    Principal component analysis (PCA) upon the first dimension of an 2D-array.
    pi                     Parameter importance index PI or alternatively B index calculation.
    plot_brewer            Plots available Brewer color maps in pdf file.
    position               Position arrays of subplots to be used with add_axes.
    print_brewer           Prints available Brewer colormap names.
    readhdf                Reads variables or information from hdf4 and hdf5 files.
    readhdf4               Reads variables or information from hdf4 files.
    readhdf5               Reads variables or information from hdf5 file.
    readnc                 Wrapper for readnetcdf.
    readnetcdf             Reads variables or information from netcdf file.
    register_brewer        Registers and registers Brewer colormap.
    rgb                    Interpolate between colours; make continuous colour maps.
    roman2int              Roman numeral to integer conversion.
    rossner                Wrapper for outlier.
    t2sap                  Conversion of temperature difference to sap flux density.
    savitzky_golay         Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter.
    savitzky_golay2d       Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter.
    saltelli               Parameter sampling for Sobol indices calculation.
    sce                    Shuffle-Complex-Evolution algorithm for function min(max)imisation
    semivariogram          Calculates semivariogram from spatial data.
    sendmail               Send an e-mail.
    sg                     Wrapper savitzky_golay.
    sg2d                   Wrapper savitzky_golay2d.
    sigma_filter           Mask values deviating more than z standard deviations from a given function.
    signature2plot         Write a copyright notice on a plot.
    tail                   Return list with last n lines of file.
    maskgroup              Masks elements in a 1d array gathered in small groups.
    savez                  Save several numpy arrays into a single file in uncompressed ``.npz`` format.
    savez_compressed       Save several arrays into a single file in compressed ``.npz`` format.
    sobol_index            Calculates the first-order and total variance-based sensitivity indices.
    sread                  Reads in string array from ascii file.
    srrasa                 Generates stratified random 2D points within a given rectangular area.
    srrasa_trans           Generates stratified random 2D transects within a given rectangular area.
    tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.
    tee                    Prints arguments on screen and in file.
    timestepcheck          Fills missing time steps in ascii data files
    tsym                   Raw unicodes for common symbols.
    unpack                 Similar to Fortran unpack function with mask.
    volume_poly            Volume of function above a polygon
    writenetcdf            Write netCDF4 file.
    xkcd                   Make plot look handdrawn.
    yrange                 Calculates plot range from input array.
    zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
    zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).


    Provided functions and modules per category
    -------------------------------------------
        Array manipulation
        Ascii files
        Data processing
        Date & Time
        Grids / Polygons
        Isotopes
        Math
        Meteorology
        Miscellaneous
        Models
        Plotting
        Special files
    -------------------------------------------

    Array manipulation
    ------------------
    argsort                Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables.
    closest                Index in array which entry is closest to a given number.
    pack                   Similar to Fortran pack function with mask.
    maskgroup              Masks elements in a 1d array gathered in small groups.
    unpack                 Similar to Fortran unpack function with mask.


    Ascii files
    -----------
    fread                  Reads in float array from ascii file.
    fsread                 Simultaneous read of float and string array from ascii file.
    head                   Return list with first n lines of file.
    lif                    Count number of lines in file.
    sread                  Reads in string array from ascii file.
    tail                   Return list with last n lines of file.


    Data processing
    ---------------
    convex_hull            Calculate subset of points that make a convex hull around a set of 2D points.
    eddybox                Module containing Eddy Covaraince utilities, see eddybox folder for details
    eddysuite              Example file for processing Eddy data with eddybox and EddySoft
    fill_nonfinite         Fill missing values by linear interpolation.
    gap2lai                Calculation of leaf projection and leaf area index from gap probability observations.
    interpol               One-dimensional linear interpolation on first dimension.
    kriging                Krig a surface from a set of 2D points.
    kernel_regression      Multi-dimensional non-parametric regression.
    kernel_regression_h    Optimal bandwidth for kernel regression.
    leafprojection         Calculation of leaf projection from leaf angle observations.
    line_dev_mask          Maskes elements of an array deviating from a line fit.
    logtools               Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle.
    mad                    Median absolute deviation test.
    means                  Calculate daily, monthly, yearly, etc. means of data depending on date stamp.
    outlier                Rossner''s extreme standardized deviate outlier test.
    pca                    Principal component analysis (PCA) upon the first dimension of an 2D-array.
    rossner                Wrapper for outlier.
    t2sap                  Conversion of temperature difference to sap flux density.
    savitzky_golay         Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter.
    savitzky_golay2d       Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter.
    semivariogram          Calculates semivariogram from spatial data.
    sg                     Wrapper savitzky_golay.
    sg2d                   Wrapper savitzky_golay2d.
    sigma_filter           Mask values deviating more than z standard deviations from a given function.
    srrasa                 Generates stratified random 2D points within a given rectangular area.
    srrasa_trans           Generates stratified random 2D transects within a given rectangular area.
    timestepcheck          Fills missing time steps in ascii data files


    Date & Time
    -----------
    date2dec               Converts arrays with calendar date to decimal date.
    dec2date               Converts arrays with decimal date to calendar date.


    Grids / Polygons
    ----------------
    area_poly              Area of a polygon
    cellarea               Calc areas of grid cells in m^2.
    grid_mid2edge          Longitude and latitude grid edges from grid midpoints.
    homo_sampling          Generation of homogeneous, randomly distributed points in a given rectangular area.
    in_poly                Determines whether a 2D point falls in a polygon.
    inpoly                 Wrapper for in_poly.
    volume_poly            Volume of function above a polygon
    get_angle              Returns the angle in radiant from each point in xy1 to each point in xy2.
    get_nearest            Returns a value z for each point in xy near to the xyz field.


    Isotopes
    --------
    alpha_equ_h2o          Equilibrium fractionation between liquid water and vapour
    alpha_kin_h2o          Kinetic fractionation of molecular diffusion of water vapour
    cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
    tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.


    Math
    ----
    around                 Round to the passed power of ten.
    correlate              Computes the cross-correlation function of two series x and y.
    div                    Wrapper for division.
    division               Divide two arrays, return 'otherwise' if division by 0.
    elementary_effects     Morris measures mu, stddev and mu*
    ellipse_area           Area of ellipse (or circle)
    errormeasures          Definition of different error measures.
    functions              Module with common functions that are used in curve_fit or fmin parameter estimations.
    heaviside              Heaviside (or unit step) operator.
    jab                    Jackknife-after-Bootstrap error.
    lagcorr                Calculate time lag of maximum or minimum correlation of two arrays.
    lhs                    Latin Hypercube Sampling of any distribution without correlations.
    morris_sampling        Sampling of optimised trajectories for Morris measures / elementary effects
    pi                     Parameter importance index PI or alternatively B index calculation.
    saltelli               Parameter sampling for Sobol indices calculation.
    sce                    Shuffle-Complex-Evolution algorithm for function min(max)imisation
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
    directory_from_gui     Open directory selection dialog, returns selected directory
    encrypt                Module to encrypt and decrypt text using a key system as well as a cipher.
    file                   Module with file list function.
    file_from_gui          Open file selection dialog for one single file, returns selected files
    files_from_gui         Open file selection dialog, returns selected files
    find_in_path           Look for file in system path.
    ftp                    Module with functions for interacting with an open FTP connection.
    int2roman              Integer to roman numeral conversion.
    roman2int              Roman numeral to integer conversion.
    sendmail               Send an e-mail.
    tee                    Prints arguments on screen and in file.
    zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
    zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).


    Models
    ------
    leafmodel              Model to compute photosynthesis and stomatal conductance of canopies


    Plotting
    --------
    abc2plot               Write a, b, c, ... on plots.
    clockplot              The clockplot of mHM.
    colors                 Wrapper for colour.
    colours                Define UFZ colours.
    get_brewer             Registers and returns Brewer colormap.
    plot_brewer            Plots available Brewer color maps in pdf file.
    position               Position arrays of subplots to be used with add_axes.
    print_brewer           Prints available Brewer colormap names.
    register_brewer        Registers and registers Brewer colormap.
    rgb                    Interpolate between colours; make continuous colour maps.
    signature2plot         Write a copyright notice on a plot.
    tsym                   Raw unicodes for common symbols.
    xkcd                   Make plot look handdrawn.
    yrange                 Calculates plot range from input array.


    Special files
    -------------
    dumpnetcdf             Convenience function for writenetcdf
    hdfread                Wrapper for readhdf.
    hdf4read               Wrapper for readhdf4.
    hdf5read               Wrapper for readhdf5.
    ncread                 Wrapper for readnetcdf.
    netcdfread             Wrapper for readnetcdf.
    readhdf                Reads variables or information from hdf4 and hdf5 files.
    readhdf4               Reads variables or information from hdf4 files.
    readhdf5               Reads variables or information from hdf5 file.
    readnc                 Wrapper for readnetcdf.
    readnetcdf             Reads variables or information from netcdf file.
    savez                  Save several numpy arrays into a single file in uncompressed ``.npz`` format.
    savez_compressed       Save several arrays into a single file in compressed ``.npz`` format.
    writenetcdf            Write netCDF4 file.


    License
    -------

    This file is part of the UFZ Python package.

    Not all files in the package are free software. The license is given in the 'License' section
    of the docstring of each routine.

    There are 3 possibilities:
    1. The routine is not yet released under the GNU Lesser General Public License.
       This is marked by a text such as
            This file is part of the UFZ Python package.

            It is NOT released under the GNU Lesser General Public License, yet.

            If you use this routine, please contact Matthias Cuntz.

            Copyright 2012-2013 Matthias Cuntz
       If you want to use this routine for publication or similar, please contact the author for possible co-authorship.

    2. The routine is already released under the GNU Lesser General Public License
       but if you use the routine in a publication or similar, you have to cite the respective publication, e.g.
            If you use this routine in your work, you should cite the following reference
            Goehler M, J Mai, and M Cuntz (2013)
                Use of eigendecomposition in a parameter sensitivity analysis of the Community Land Model,
                J Geophys Res 188, 904-921, doi:10.1002/jgrg.20072

    3. The routine is released under the GNU Lesser General Public License. The following applies:
       The UFZ Python package is free software: you can redistribute it and/or modify
       it under the terms of the GNU Lesser General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       The UFZ Python package is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
       GNU Lesser General Public License for more details.

       You should have received a copy of the GNU Lesser General Public License
       along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
       If not, see <http://www.gnu.org/licenses/>.

    Copyright 2009-2014 Matthias Cuntz, Arndt Piayda, Matthias Zink, Tino Rau, Maren Goehler,
                        Stephan Thober, Juliane Mai, Andreas Wiedemann


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
              MC, Nov 2012 - nee2gpp, nee2gpp_falge, nee2gpp_lasslop, nee2gpp_reichstein
              MC, Dec 2012 - functions
                           - gap_filling obsolete
              MC, Feb 2013 - area_poly
              MC & JM, Feb 2013 - volume_poly
              MC, Feb 2013 - ported to Python 3
              MC, Mar 2013 - find_in_path, xkcd
              MC, Apr 2013 - rgb
              MC, Jun 2013 - colours
              MC, Jul 2013 - fill_nonfinite, means
              MC, Oct 2013 - morris, sce, inpoly, rossner, netcdfread, ncread, readnc, hdfread, hdf4read, hdf5read
              AP, Feb 2014 - maskgroup
                           - line_dev_mask
              MC, Feb 2014 - removed all import *
                           - sigma_filter
              AP, Mar 2014 - lagcorr
              MC, Apr 2014 - correlate
              MC, May 2014 - signature2plot
                           - adapted new CHS license scheme
              AP, May 2014 - get_nearest
              AP, Jun 2014 - get_angle
              AW, Jun 2014 - t2sap
              AP, Jul 2014 - errormeasures, homo_sampling, sltclean, meteo4slt, eddycorr, eddyspec
              AP, Aug 2014 - planarfit, timestepcheck, fluxplot, itc, spikeflag, ustarflag
                           - fluxflag, fluxfill
              AP, Sep 2014 - energyclosure, fluxpart
              ST, Sep 2014 - dumpnetcdf
              MC, Sep 2014 - alpha_equ_h2o, alpha_kin_h2o
              AP, Sep 2014 - leafmodel
                           - profile2storage
                           - eddybox -> moved eddycorr, eddyspec, energyclosure,
                                              fluxfill, fluxflag, fluxpart,
                                              fluxplot, gapfill, itc, meteo4slt,
                                              nee2gpp, nee2gpp_falge,
                                              nee2gpp_lasslop, nee2gpp_reichstein,
                                              planarfit, profile2storage,
                                              sltclean, spikeflag, ustarflag into
                                              eddybox module
              AP, Sep 2014 - eddysuite
              MC, Oct 2014 - ufz module -> ufz package
                           - clockplot, ellipse_area, savez, savez_compressed, grid_mid2edge, tee
              MC, Nov 2014 - pca, head
              AP, Nov 2014 - gap2lai, leafprojection
              MC, Dec 2014 - directory_from_gui, file_from_gui, files_from_gui
                           - logtools, 
                           - sendmail, argsort, tail
                           - ftp
                           - file
                           - encrypt
              MC, Feb 2015 - fsread
"""
from __future__ import print_function

# sub-packages without dependencies to all ufz
from ufz.const     import const
from ufz.encrypt   import encrypt
from ufz.functions import functions

# Routines
from .abc2plot          import abc2plot
from .alpha_equ_h2o     import alpha_equ_h2o
from .alpha_kin_h2o     import alpha_kin_h2o
from .area_poly         import area_poly
from .argsort           import argsort
from .around            import around
from .autostring        import autostring, astr
from .brewer            import register_brewer, get_brewer, plot_brewer, print_brewer
try:
    from .calcvpd       import calcvpd
except ImportError:
    pass # obsolete
from .cellarea          import cellarea
from .clockplot         import clockplot
from .closest           import closest
from .colours           import colours, colors
from .convex_hull       import convex_hull
from .correlate         import correlate
from .cuntz_gleixner    import cuntz_gleixner
from .date2dec          import date2dec
from .dec2date          import dec2date
from .dewpoint          import dewpoint
from .division          import division, div
from .ellipse_area      import ellipse_area
from .errormeasures     import bias, mae, mse, rmse, nse, pear2
from .esat              import esat
from .fgui              import directory_from_gui, file_from_gui, files_from_gui
from .fill_nonfinite    import fill_nonfinite
from .find_in_path      import find_in_path
from .fread             import fread
from .fsread            import fsread
try:
    from .gap_filling   import gap_filling
except ImportError:
    pass # obsolete
from .gap2lai           import gap2lai, leafprojection 
from .get_angle         import get_angle
from .get_nearest       import get_nearest
from .grid_mid2edge     import grid_mid2edge
from .head              import head
from .heaviside         import heaviside
from .homo_sampling     import homo_sampling
from .in_poly           import in_poly, inpoly
from .interpol          import interpol
from .kernel_regression import kernel_regression, kernel_regression_h
from .kriging           import kriging
from .lagcorr           import lagcorr
from .lhs               import lhs
from .lif               import lif
from .line_dev_mask     import line_dev_mask
from .jab               import jab
from .mad               import mad
from .maskgroup         import maskgroup
from .means             import means
from .morris            import morris_sampling, elementary_effects
from .npyio             import savez, savez_compressed
try:
    from .outlier       import outlier, rossner
except:
    print("No extra statistics in scipy, i.e. in UFZ package. Disabled functions: outlier, rossner.")
from .pack              import pack
from .pca               import pca, check_pca
from .pi                import pi
from .position          import position
from .readhdf           import readhdf, hdfread
from .readhdf4          import readhdf4, hdf4read
from .readhdf5          import readhdf5, hdf5read
from .readnetcdf        import readnetcdf, netcdfread, ncread, readnc
from .rgb               import rgb_blend, rgb_range, rgb_gradient
from .romanliterals     import int2roman, roman2int
from .sap_app           import t2sap
from .saltelli          import saltelli
from .savitzky_golay    import savitzky_golay, sg, savitzky_golay2d, sg2d
from .sce               import sce
from .semivariogram     import semivariogram
from .sendmail          import sendmail
from .sigma_filter      import sigma_filter
from .signature2plot    import signature2plot
from .sobol_index       import sobol_index
from .sread             import sread
from .srrasa            import srrasa, srrasa_trans
from .tail              import tail
from .tcherkez          import tcherkez
from .tee               import tee
from .timestepcheck     import timestepcheck
from .tsym              import tsym
from .unpack            import unpack
from .volume_poly       import volume_poly
from .writenetcdf       import writenetcdf, dumpnetcdf
from .xkcd              import xkcd
from .yrange            import yrange
from .zacharias         import zacharias, zacharias_check

# sub-packages with dependencies to all ufz have to be loaded separately as in scipy
from . import eddybox
from . import file
from . import ftp
from . import leafmodel
from . import logtools

# Information
__author__   = "Matthias Cuntz"
__version__  = '3.1.0'
__revision__ = "Revision: 1923"
__date__     = 'Date: 05.12.2014'

# Main
if __name__ == '__main__':
    print('\nUFZ Computational Hydrosystems Python Package.')
    print("Version {:s} from {:s}.".format(__version__,__date__))
    print('\nThis is the README file. See als the license file LICENSE.\n\n')
    f = open('README','r')
    for line in f: print(line,end='')
    f.close()
