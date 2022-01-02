#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
JAMS Python Utilities

Package offers miscellaneous functions and sub-modules in different categories.

Get help on each function by typing
>>> import jams
>>> help(jams.function)

Provided functions and modules (alphabetic w/o obsolete functions)
------------------------------------------------------------------
apply_undef            Use a function on masked arguments.
area_poly              Area of a polygon.
around                 Round to the passed power of ten.
astr                   Wrapper for autostring.
autostring             Format number (array) with given decimal precision.
baseflow               Calculate baseflow from discharge timeseries
cellarea               Calc areas of grid cells in m^2.
climate_index_knoben   Determines continuous climate indexes based on Knoben et al. (2018).
clockplot              The clockplot of mHM.
convex_hull            Calculate subset of points that make a convex hull around a set of 2D points.
correlate              Computes the cross-correlation function of two series x and y.
cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
dag                    Generation and plotting of (connected) directed acyclic graphs with one source node.
dielectric_water       Dielectric constant of liquid water.
directories_from_gui   Open directory selection dialogs, returns consecutiveley selected directories.
directory_from_gui     Open directory selection dialog, returns selected directory.
delta_isogsm2          Calculate delta values from downloaded IsoGSM2 data.
dewpoint               Calculates the dew point from ambient humidity.
date2dec               Converts arrays with calendar date to decimal date.
dec2date               Converts arrays with decimal date to calendar date.
dfgui                  A minimalistic GUI for analyzing Pandas DataFrames based on wxPython.
distributions          Module for pdfs of additional distributions.
dumpnetcdf             Convenience function for writenetcdf
eddybox                Module containing Eddy Covaraince utilities, see eddysuite.py for details
eddysuite              Example file for processing Eddy data with eddybox and EddySoft
ellipse_area           Area of ellipse (or circle)
encrypt                Module to encrypt and decrypt text using a key system as well as a cipher.
errormeasures          Definition of different error measures.
esat                   Calculates the saturation vapour pressure of water/ice.
fftngo                 Fast fourier transformation for dummies (like me)
Field                  Generates random hydraulic conductivity fields.
files                  Module with file list function.
file_from_gui          Open file selection dialog for one single file, returns selected files
files_from_gui         Open file selection dialog, returns selected files
fill_nonfinite         Fill missing values by interpolation.
find_in_path           Look for file in system path.
Filtered_Incompr_Field Generates random filtered velocity fields.
ftp                    Module with functions for interacting with an open FTP connection.
fwrite                 Writes an array to ascii file
gap2lai                Calculation of leaf area index from gap probability observations.
geoarray               Pythonic gdal wrapper
get_angle              Returns the angle in radiant from each point in xy1 to each point in xy2.
get_era_interim        Download ERA-Interim data suitable to produce MuSICA input data.
get_era5               Download ERA5 data suitable to produce MuSICA input data.
get_isogsm2            Get IsoGSM2 output.
get_nearest            Returns a value z for each point in xy near to the xyz field.
grid_mid2edge          Longitude and latitude grid edges from grid midpoints.
hdfread                Wrapper for readhdf.
hdf4read               Wrapper for readhdf4.
hdf5read               Wrapper for readhdf5.
head                   Return list with first n lines of file.
heaviside              Heaviside (or unit step) operator.
homo_sampling          Generation of homogeneous, randomly distributed points in a given rectangular area.
Incompr_Field          Generates random velocity fields.
in_poly                Determines whether a 2D point falls in a polygon.
inpoly                 Wrapper for in_poly.
interpol               One-dimensional linear interpolation on first dimension.
intersection           Intersection of two curves from x,y coordinates.
jab                    Jackknife-after-Bootstrap error.
jConfigParser          Extended Python ConfigParser.
kernel_regression      Multi-dimensional non-parametric regression.
kernel_regression_h    Optimal bandwidth for kernel regression.
kriging                Krig a surface from a set of 2D points.
lagcorr                Calculate time lag of maximum or minimum correlation of two arrays.
lat_fmt                Set lat label string (called by Basemap.drawparallels) if LaTeX package clash.
leafmodel              Model to compute photosynthesis and stomatal conductance of canopies.
leafprojection         Calculation of leaf projection from leaf angle observations.
level1                 Module with functions dealing with CHS level1 data files, data and flags.
lhs                    Latin Hypercube Sampling of any distribution without correlations.
lif                    Count number of lines in file.
line_dev_mask          Maskes elements of an array deviating from a line fit.
logtools               Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle.
lon_fmt                Set lon label string (called by Basemap.drawmeridians) if LaTeX package clash.
lowess                 Locally linear regression in n dimensions.
mad                    Median absolute deviation test.
mat2nc                 Converts Matlab file *.mat into NetCDF *.nc.
means                  Calculate daily, monthly, yearly, etc. means of data depending on date stamp.
nc2nc                  Copy netcdf file deleting, renaming, replacing variables and attribues.
ncread                 Wrapper for readnetcdf.
netcdfread             Wrapper for readnetcdf.
netcdf4                Convenience layer around netCDF4
outlier                Rossner''s extreme standardized deviate outlier test.
pack                   Similar to Fortran pack function with mask.
pareto_metrics         Performance metrics to compare Pareto fronts.
pca                    Principal component analysis (PCA) upon the first dimension of an 2D-array.
pet_oudin              Daily potential evapotranspiration following the Oudin formula.
pi                     Parameter importance index PI or alternatively B index calculation.
plot                   Module with code snippets for plotting.
pritay                 Daily reference evapotranspiration after Priestley & Taylor.
pso                    Particle swarm optimization
qa                     Module of quality (error) measures.
readhdf                Reads variables or information from hdf4 and hdf5 files.
readhdf4               Reads variables or information from hdf4 files.
readhdf5               Reads variables or information from hdf5 file.
readnc                 Wrapper for readnetcdf.
readnetcdf             Reads variables or information from netcdf file.
river_network          a class for creating a river network from a DEM including flow direction, flow accumulation and channel order
rolling                Reshape an array in a "rolling window" style.
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
tail                   Return list with last n lines of file.
maskgroup              Masks elements in a 1d array gathered in small groups.
samevalue              Checks if abs. differences of array values within a certain window are smaller than threshold.
savez                  Save several numpy arrays into a single file in uncompressed ``.npz`` format.
savez_compressed       Save several arrays into a single file in compressed ``.npz`` format.
smax                   Calculating smooth maximum of two numbers
smin                   Calculating smooth minimum of two numbers
sobol_index            Calculates the first-order and total variance-based sensitivity indices.
srrasa                 Generates stratified random 2D points within a given rectangular area.
srrasa_trans           Generates stratified random 2D transects within a given rectangular area.
tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.
timestepcheck          Fills missing time steps in ascii data files
tsym                   Raw unicodes for common symbols.
unpack                 Similar to Fortran unpack function with mask.
volume_poly            Volume of function above a polygon
writenetcdf            Write netCDF4 file.
xkcd                   Make plot look handdrawn.
xlsread                Wrapper for xread.
xlsxread               Wrapper for xread.
xread                  Simultaneous read of float and string array from Excel file.
yrange                 Calculates plot range from input array.
zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).


Deprecated
----------
abc2plot               Write a, b, c, ... on plots.
alpha_equ_h2o          Equilibrium fractionation between liquid water and vapour.
alpha_kin_h2o          Kinetic fractionation of molecular diffusion of water vapour.
argmax                 Wrapper for numpy.argmax, numpy.ma.argmax, and using max for Python iterables.
argmin                 Wrapper for numpy.argmin, numpy.ma.argmin, and using min for Python iterables.
argsort                Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables.
ascii2ascii            Convert date notations between to ascii date format DD.MM.YYYY hh:mm:ss.
ascii2en               Convert date notations to English date format YYYY-MM-DD hh:mm:ss.
ascii2fr               Convert date notations to French date format DD/MM/YYYY hh:mm:ss.
ascii2us               Convert date notations to American date format MM/DD/YYYY hh:mm:ss.
closest                Index in array which entry is closest to a given number.
color                  Module with color functions for plotting.
const                  Provides physical, mathematical, computational, and isotope constants.
div                    Wrapper for division.
division               Divide two arrays, return 'otherwise' if division by 0.
elementary_effects     Morris measures mu, stddev and mu*
en2ascii               Convert date notations from English YYYY-MM-DD to ascii date format DD.MM.YYYY hh:mm:ss.
fr2ascii               Convert date notations from French DD/MM/YYYT to ascii date format DD.MM.YYYY hh:mm:ss.
fread                  Reads in float array from ascii file.
fsread                 Simultaneous read of float and string array from ascii file.
functions              Module with common functions that are used in curve_fit or fmin parameter estimations.
get_brewer             Registers and returns Brewer colormap.
int2roman              Integer to roman numeral conversion.
mcPlot                 Matthias Cuntz' standard plotting class.
morris_sampling        Sampling of optimised trajectories for Morris measures / elementary effects
plot_brewer            Plots available Brewer color maps in pdf file.
position               Position arrays of subplots to be used with add_axes.
print_brewer           Prints available Brewer colormap names.
register_brewer        Registers and registers Brewer colormap.
roman2int              Roman numeral to integer conversion.
screening              Samples trajectories, runs model and returns measures of Morris Elemenary Effects
signature2plot         Write a copyright notice on a plot.
sread                  Reads in string array from ascii file.
str2tex                Convert strings to LaTeX strings in math environement used by matplotlib's usetex
tee                    Prints arguments on screen and in file.
us2ascii               Convert date notations from American MM/DD/YYYY to ascii format DD.MM.YYYY hh:mm:ss.


Provided functions and modules per category
-------------------------------------------
    Array manipulation
    Ascii files
    Data processing
    Date & Time
    Grids / Polygons
    Hydrology
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
pack                   Similar to Fortran pack function with mask.
samevalue              Checks if abs. differences of array values within a certain window are smaller than threshold.
maskgroup              Masks elements in a 1d array gathered in small groups.
rolling                Reshape an array in a "rolling window" style.
smax                   Calculating smooth maximum of two numbers
smin                   Calculating smooth minimum of two numbers
unpack                 Similar to Fortran unpack function with mask.


Ascii files
-----------
fwrite                 Writes an array to ascii file
head                   Return list with first n lines of file.
lif                    Count number of lines in file.
tail                   Return list with last n lines of file.


Data processing
---------------
convex_hull            Calculate subset of points that make a convex hull around a set of 2D points.
eddybox                Module containing Eddy Covaraince utilities, see eddybox folder for details
eddysuite              Example file for processing Eddy data with eddybox and EddySoft
fill_nonfinite         Fill missing values by interpolation.
gap2lai                Calculation of leaf projection and leaf area index from gap probability observations.
interpol               One-dimensional linear interpolation on first dimension.
kriging                Krig a surface from a set of 2D points.
kernel_regression      Multi-dimensional non-parametric regression.
kernel_regression_h    Optimal bandwidth for kernel regression.
leafprojection         Calculation of leaf projection from leaf angle observations.
level1                 Module with functions dealing with CHS level1 data files, data and flags.
line_dev_mask          Mask elements of an array deviating from a line fit.
logtools               Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle.
lowess                 Locally linear regression in n dimensions.
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


Hydrology
---------
baseflow               Calculate baseflow from discharge timeseries
river_network          a class for creating a river network from a DEM including flow direction, flow accumulation and channel order


Isotopes
--------
cuntz_gleixner         Cuntz-Gleixner model of 13C discrimination.
delta_isogsm2          Calculate delta values from downloaded IsoGSM2 data.
get_isogsm2            Get IsoGSM2 output.
tcherkez               Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.


Math
----
around                 Round to the passed power of ten.
correlate              Computes the cross-correlation function of two series x and y.
dag                    Generation and plotting of (connected) directed acyclic graphs with one source node.
distributions          Module for pdfs of additional distributions.
ellipse_area           Area of ellipse (or circle)
errormeasures          Definition of different error measures.
fftngo                 Fast fourier transformation for dummies (like me)
heaviside              Heaviside (or unit step) operator.
intersection           Intersection of two curves from x,y coordinates.
jab                    Jackknife-after-Bootstrap error.
lagcorr                Calculate time lag of maximum or minimum correlation of two arrays.
lhs                    Latin Hypercube Sampling of any distribution without correlations.
pareto_metrics         Performance metrics to compare Pareto fronts.
pi                     Parameter importance index PI or alternatively B index calculation.
pso                    Particle swarm optimization
qa                     Module of quality assessment (error) measures.
saltelli               Parameter sampling for Sobol indices calculation.
sce                    Shuffle-Complex-Evolution algorithm for function min(max)imisation
sobol                  Generates Sobol sequences
sobol_index            Calculates the first-order and total variance-based sensitivity indices.


Meteorology
-----------
climate_index_knoben   Determines continuous climate indexes based on Knoben et al. (2018).
dewpoint               Calculates the dew point from ambient humidity.
dielectric_water       Dielectric constant of liquid water.
esat                   Calculates the saturation vapour pressure of water/ice.
get_era_interim        Download ERA-Interim data suitable to produce MuSICA input data.
get_era5               Download ERA5 data suitable to produce MuSICA input data.
pet_oudin              Daily potential evapotranspiration following the Oudin formula.
pritay                 Daily reference evapotranspiration after Priestley & Taylor


Miscellaneous
-------------
apply_undef            Use a function on masked arguments.
astr                   Wrapper for autostring.
autostring             Format number (array) with given decimal precision.
directories_from_gui   Open directory selection dialogs, returns consecutiveley selected directories
directory_from_gui     Open directory selection dialog, returns selected directory
encrypt                Module to encrypt and decrypt text using a key system as well as a cipher.
files                  Module with file list function.
file_from_gui          Open file selection dialog for one single file, returns selected files
files_from_gui         Open file selection dialog, returns selected files
find_in_path           Look for file in system path.
ftp                    Module with functions for interacting with an open FTP connection.
sendmail               Send an e-mail.
zacharias              Soil water content with van Genuchten and Zacharias et al. (2007).
zacharias_check        Checks validity of parameter set for Zacharias et al. (2007).


Models
------
Field                  Generates random hydraulic conductivity fields.
Filtered_Incompr_Field Generates random filtered velocity fields.
Incompr_Field          Generates random velocity fields.
leafmodel              Model to compute photosynthesis and stomatal conductance of canopies


Plotting
--------
clockplot              The clockplot of mHM.
dfgui                  A minimalistic GUI for analyzing Pandas DataFrames based on wxPython.
lat_fmt                Set lat label string (called by Basemap.drawparallels) if LaTeX package clash.
lon_fmt                Set lon label string (called by Basemap.drawmeridians) if LaTeX package clash.
plot                   Module with code snippets for plotting.
tsym                   Raw unicodes for common symbols.
xkcd                   Make plot look handdrawn.
yrange                 Calculates plot range from input array.


Special files
-------------
dumpnetcdf             Convenience function for writenetcdf
geoarray               Pythonic gdal wrapper
hdfread                Wrapper for readhdf.
hdf4read               Wrapper for readhdf4.
hdf5read               Wrapper for readhdf5.
jConfigParser          Extended Python ConfigParser.
mat2nc                 Converts Matlab file *.mat into NetCDF *.nc.
nc2nc                  Copy netcdf file deleting, renaming, replacing variables and attribues.
ncread                 Wrapper for readnetcdf.
netcdfread             Wrapper for readnetcdf.
netcdf4                Convenience layer around netCDF4
readhdf                Reads variables or information from hdf4 and hdf5 files.
readhdf4               Reads variables or information from hdf4 files.
readhdf5               Reads variables or information from hdf5 file.
readnc                 Wrapper for readnetcdf.
readnetcdf             Reads variables or information from netcdf file.
savez                  Save several numpy arrays into a single file in uncompressed ``.npz`` format.
savez_compressed       Save several arrays into a single file in compressed ``.npz`` format.
writenetcdf            Write netCDF4 file.
xlsread                Wrapper for xread.
xlsxread               Wrapper for xread.
xread                  Simultaneous read of float and string array from Excel file.


License
-------
This file is part of the JAMS Python package, distributed under the MIT
License. The JAMS Python package originates from the former UFZ Python library,
Department of Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany.

Copyright (c) 2009-2021 Matthias Cuntz, Juliane Mai, Stephan Thober, Arndt
Piayda

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History
-------
Written,  Matthias Cuntz, Jul 2009
Modified, Matthias Cuntz, Jul 2009
              - lif, fread, sread, readnetcdf, cellarea, pack, unpack
          Matthias Cuntz, Aug 2009 - position
          Maren Goehler, Jul 2010  - outlier
          Arndt Piayda, Jan 2011   - date2dec, dec2date
          Arndt Piayda, Feb 2011   - semivariogram
          Tino Rau, May 2011       - gap_filling
          Tino Rau, May 2011       - calcvpd
          Matthias Cuntz, Jun 2011
              - /usr/bin/python to /usr/bin/env python
              - tsym, around
          Matthias Cuntz, Nov 2011 - mad
          Matthias Cuntz, Nov 2011 - try netcdf and stats routines
          Matthias Cuntz, Nov 2011 - autostring
          Matthias Cuntz, Jan 2012
              - esat, closest, dewpoint, division, heaviside, tcherkez, yrange,
              - const, cuntz_gleixner
              - calcvpd obsolete
          Matthias Cuntz, Mar 2012 - gapfill, nee2gpp
          Matthias Cuntz, May 2012
              - astr, div, sobol_index, pi, roman, zacharias, saltelli
          Matthias Zink, Jun 2012  - writenetcdf
          Matthias Cuntz, Jun 2012 - roman -> romanliterals, interpol
          Matthias Zink, Jun 2012  - readhdf5
          Matthias Cuntz, Jun 2012 - readhdf4, readhdf
          Matthias Cuntz, Sep 2012 - brewer
          Matthias Cuntz, Oct 2012 - savitzky_golay
          Matthias Cuntz, Nov 2012
              - added netcdftime but no import so that available w/o netcdf
          Arndt Piayda, Nov 2012
              - convex_hull, in_poly, kriging, semivariogram update
              - srrasa, srrasa_trans
          Matthias Cuntz, Nov 2012
              - nee2gpp, nee2gpp_falge, nee2gpp_lasslop, nee2gpp_reichstein
          Matthias Cuntz, Dec 2012
              - functions
              - gap_filling obsolete
          Matthias Cuntz, Feb 2013 - area_poly
          Matthias Cuntz & Juliane Mai, Feb 2013 - volume_poly
          Matthias Cuntz, Feb 2013 - ported to Python 3
          Matthias Cuntz, Mar 2013 - find_in_path, xkcd
          Matthias Cuntz, Apr 2013 - rgb
          Matthias Cuntz, Jun 2013 - colours
          Matthias Cuntz, Jul 2013 - fill_nonfinite, means
          Matthias Cuntz, Oct 2013
              - morris, sce, inpoly, rossner, netcdfread, ncread, readnc
              - hdfread, hdf4read, hdf5read
          Arndt Piayda, Feb 2014
              - maskgroup
              - line_dev_mask
          Matthias Cuntz, Feb 2014
              - removed all import *
              - sigma_filter
          Arndt Piayda, Mar 2014   - lagcorr
          Matthias Cuntz, Apr 2014 - correlate
          Matthias Cuntz, May 2014
              - signature2plot
              - adapted new CHS license scheme
          Arndt Piayda, May 2014   - get_nearest
          Arndt Piayda, Jun 2014   - get_angle
          Andreas Wiedemann, Jun 2014 - t2sap
          Arndt Piayda, Jul 2014
              - errormeasures, homo_sampling, sltclean, meteo4slt
              - eddycorr, eddyspec
          Arndt Piayda, Aug 2014
              - planarfit, timestepcheck, fluxplot, itc, spikeflag, ustarflag
              - fluxflag, fluxfill
          Arndt Piayda, Sep 2014   - energyclosure, fluxpart
          Stephan Thober, Sep 2014 - dumpnetcdf
          Matthias Cuntz, Sep 2014 - alpha_equ_h2o, alpha_kin_h2o
          Arndt Piayda, Sep 2014
              - leafmodel
              - profile2storage
              - eddybox -> moved eddycorr, eddyspec, energyclosure,
                                 fluxfill, fluxflag, fluxpart,
                                 fluxplot, gapfill, itc, meteo4slt,
                                 nee2gpp, nee2gpp_falge,
                                 nee2gpp_lasslop, nee2gpp_reichstein,
                                 planarfit, profile2storage,
                                 sltclean, spikeflag, ustarflag
                           into eddybox module
          Arndt Piayda, Sep 2014   - eddysuite
          Matthias Cuntz, Oct 2014
              - ufz module -> ufz package
              - clockplot, ellipse_area, savez, savez_compressed
              - grid_mid2edge, tee
          Matthias Cuntz, Nov 2014 - pca, head
          Arndt Piayda, Nov 2014   - gap2lai, leafprojection
          Matthias Cuntz, Dec 2014
              - directory_from_gui, file_from_gui, files_from_gui
              - logtools,
              - sendmail, argsort, tail
              - ftp
              - file
              - encrypt
          Matthias Cuntz, Feb 2015
              - fsread
              - ascii2ascii, ascii2eng, eng2ascii
          Matthias Cuntz, Mar 2015
              - module level1 with get_flag, set_flag, read_data, write_data
              - rename file to files
              - dielectric_water
              - color
              - redone all __init__.py
          David Schaefer, Sep 2015 - hollickLyneFilter
          Arndt Piayda, Sep 2015   - confidence intervals to errormeasures
          Andreas Wiedemann, Sep 2015 - samevalue
          Matthias Cuntz, Oct 2015 - str2tex, lat_fmt, lon_fmt
          Matthias Cuntz, Oct 2015 - directories_from_gui
          Stephan Thober, Nov 2015 - kge
          Stephan Thober, Dec 2015 - river_network
          Stephan Thober, Feb 2016
              - function for writing 2d arrays to ascii file
          Juliane Mai, Feb 2016    - pareto_metrics
          Stephan Thober, Mar 2016 - smax, smin
          David Schaefer, Mar 2016 - netcdf4
          Matthias Cuntz, May 2016 - qa
          Matthias Cuntz, May 2016 - distributions
          Matthias Cuntz, Oct 2016 - rm colours and rgb from main directory
          Arndt Piayda, Oct 2016   - fftngo
          Juliane Mai, Oct 2016    - mat2nc
          Juliane Mai, Oct 2016    - dag
          Arndt Piayda, Oct 2016   - pritay
          David Schaefer, Oct 2016 - added geoarray
          Matthias Cuntz, Nov 2016 - ported to Python 3
          Matthias Cuntz, Nov 2016 - pso
          Arndt Piayda, Dec 2016   - rolling
          Stephan Thober, Aug 2017 - added fwrite
          Matthias Cuntz, Nov 2017 - xread
          Juliane Mai, Dec 2017    - pawn_index
          Matthias Cuntz, Dec 2017 - screening
          Matthias Cuntz, Jan 2018 - lowess
          Matthias Cuntz, Jan 2018 - apply_undef
          Matthias Cuntz, Mar 2018
              - ascii2en, en2ascii, ascii2fr, fr2ascii, ascii2us, us2ascii
          Matthias Cuntz, Jul 2018 - plot
          Matthias Cuntz, Nov 2018 - intersection, jConfigParser
          Matthias Cuntz, Jan 2019
              - dfgui, delta_isogsm2, get_era5, get_era_interim, get_isogsm2
          Matthias Cuntz, Feb 2019 - xlsread, xlsxread
          Matthias Cuntz, Apr 2019 - nc2nc
          Matthias Cuntz, Jul 2019 - argmax, argmin
          Juliane Mai, Feb 2020    - pet_oudin
          Juliane Mai, Feb 2020    - climate_index_knoben
          Matthias Cuntz, Dec 2020 - mcPlot
          Matthias Cuntz, Oct 2021 - started deprecation

"""

# sub-packages without dependencies to rest of jams
from . import const
from . import encrypt
from . import functions
from . import plot
from . import qa

# Routines
from .abc2plot             import abc2plot
from .alpha_equ_h2o        import alpha_equ_h2o
from .alpha_kin_h2o        import alpha_kin_h2o
from .apply_undef          import apply_undef
from .area_poly            import area_poly
from .argsort              import argmax, argmin, argsort
from .around               import around
from .ascii2ascii          import ascii2ascii, ascii2en, ascii2fr, ascii2us, ascii2eng, en2ascii, fr2ascii, us2ascii, eng2ascii
from .autostring           import autostring, astr
from .baseflow             import hollickLyneFilter
from .brewer               import register_brewer, get_brewer, plot_brewer, print_brewer
try:
    from .calcvpd          import calcvpd
except ImportError:
    pass # obsolete
from .cellarea             import cellarea
from .climate_index_knoben import climate_index_knoben
from .clockplot            import clockplot
from .closest              import closest
from .convex_hull          import convex_hull
from .correlate            import correlate
from .cuntz_gleixner       import cuntz_gleixner
try:
    from .dag              import create_network, source_nodes, sink_nodes, plot_network
except ImportError:
    pass # networkx not installed
from .date2dec             import date2dec
from .dec2date             import dec2date
from .delta_isogsm2        import delta_isogsm2
from .dewpoint             import dewpoint
#                          import dfgui
from .dielectric_water     import dielectric_water
from .division             import division, div
from .ellipse_area         import ellipse_area
from .errormeasures        import bias, mae, mse, rmse, nse, kge, pear2
from .esat                 import esat
from .fftngo               import fftngo
from .fgui                 import directories_from_gui, directory_from_gui, file_from_gui, files_from_gui
# from .field_gen          import Field, Incompr_Field, Filtered_Incompr_Field
from .fill_nonfinite       import fill_nonfinite
from .find_in_path         import find_in_path
from .fread                import fread
from .fsread               import fsread
from .fwrite               import fwrite
try:
    from .gap_filling      import gap_filling
except                     ImportError:
    pass # obsolete
from .gap2lai              import gap2lai, leafprojection
try:
    from .geoarray         import geoarray
except                     ImportError:
    pass
from .get_angle            import get_angle
try:
    from .get_era_interim  import get_era_interim
except                     ImportError:
    pass
from .get_era5             import get_era5
try:
    from .get_isogsm2      import get_isogsm2
except                     ImportError:
    pass
from .get_nearest          import get_nearest
from .grid_mid2edge        import grid_mid2edge
from .head                 import head
from .heaviside            import heaviside
from .homo_sampling        import homo_sampling
from .in_poly              import in_poly, inpoly
from .interpol             import interpol
from .intersection         import intersection
from .jab                  import jab
from .jconfigparser        import jConfigParser
from .kernel_regression    import kernel_regression, kernel_regression_h
from .kriging              import kriging
from .lagcorr              import lagcorr
from .latlon_fmt           import lat_fmt, lon_fmt
from .lhs                  import lhs
from .lif                  import lif
from .line_dev_mask        import line_dev_mask
from .lowess               import lowess
from .mad                  import mad
from .maskgroup            import maskgroup
from .mat2nc               import mat2nc
from .mcplot               import mcPlot
from .means                import means
from .morris               import morris_sampling, elementary_effects
from .nc2nc                import nc2nc
try:
    from .npyio            import savez, savez_compressed
except ImportError:
    pass # old numpy version
from .netcdf4              import netcdf4
try:
    from .outlier          import outlier, rossner
except:
    pass # No extra statistics in scipy and hence in JAMS. Disabled functions: outlier, rossner.
from .pack                 import pack
from .pareto_metrics       import sn, cz, hi, ef, aed, is_dominated, point_to_front
try:
    from .pawn_index       import pawn_index
except:
    pass # No statsmodels installed.
from .pca                  import pca, check_pca
from .pet_oudin            import pet_oudin
from .pi                   import pi
from .position             import position
from .pritay               import pritay
from .pso                  import pso
try:
    from .readhdf          import readhdf,  hdfread
except ImportError:
    pass # not installed
try:
    from .readhdf4         import readhdf4, hdf4read
except ImportError:
    pass # not installed
from .readhdf5             import readhdf5, hdf5read
from .readnetcdf           import readnetcdf, netcdfread, ncread, readnc
from .river_network        import river_network, upscale_fdir
from .rolling              import rolling
from .romanliterals        import int2roman, roman2int
from .saltelli             import saltelli
from .samevalue            import samevalue
from .sap_app              import t2sap
from .savitzky_golay       import savitzky_golay, sg, savitzky_golay2d, sg2d
from .sce                  import sce
from .screening            import screening
from .semivariogram        import semivariogram
from .sendmail             import sendmail
from .sigma_filter         import sigma_filter
from .signature2plot       import signature2plot
from .smooth_minmax        import smin, smax
from .sobol_index          import sobol_index
from .sread                import sread
from .srrasa               import srrasa, srrasa_trans
from .str2tex              import str2tex
from .tail                 import tail
from .tcherkez             import tcherkez
from .tee                  import tee
from .timestepcheck        import timestepcheck
from .tsym                 import tsym
from .unpack               import unpack
from .volume_poly          import volume_poly
from .writenetcdf          import writenetcdf, dumpnetcdf
from .xkcd                 import xkcd
try:
    from .xread            import xread, xlsread, xlsxread
except ImportError:
    pass # not installed
from .yrange               import yrange
from .zacharias            import zacharias, zacharias_check

# sub-packages with dependencies to rest jams have to be loaded separately as in scipy
from . import color
# ToDo: from here on redo __init__.py
from . import distributions
from . import eddybox
from . import files
from . import ftp
from . import leafmodel
from . import level1
from . import logtools


# Information
__author__   = "Matthias Cuntz"
__version__  = '4.5.0'
__date__     = 'Date: 13.12.2020'

# Main
if __name__ == '__main__':
    print('\nJAMS Python Package.')
    print("Version {:s} from {:s}.".format(__version__,__date__))
    print('\nThis is the README file. See als the license file LICENSE.\n\n')
    f = open('README','r')
    for line in f: print(line,end='')
    f.close()
