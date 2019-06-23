
# This is the Python package of JAMS.

JAMS is a general Python package offering miscellaneous functions in
different categories, such as reading different file formats, julian
date routines, or meteorological functions.

It has several subpackages offering constants, working with Eddy
covariance data and software such as EddySoft, offering special
functions, or objective functions be used with scipy.optimize.fmin or
scipy.optimize.curvefit, and much more.

Created June 2009 by Matthias Cuntz  
while at the Department Computational Hydrosystems, Helmholtz Centre
for Environmental Research - UFZ, Permoserstr. 15, 04318 Leipzig, Germany

It is distributed under the MIT License (see LICENSE file and below).

Copyright (c) 2012-2019 Matthias Cuntz, Juliane Mai, Stephan Thober, Arndt Piayda

Contact Matthias Cuntz - mc (at) macu (dot) de


---------------------------------------------------------------

### Installation

The library is maintained with a git repository at:

    https://github.com/mcuntz/jams_python/

To use it, checkout the git repository

    git clone https://github.com/mcuntz/jams_python.git

and either add it to your Python path, for example in bash:

    export PYTHONPATH=/path/to/the/jams/package

or install it with setup.py after changing into the downloaded directory:

    python setup.py install

or install it with pip:

    pip install ./

One can also install it directly with pip from the git repository:

    pip install git+https://github.com/mcuntz/jams_python.git

Append --user on pip commands if you have no root access.


---------------------------------------------------------------

### Documentation

The documentation of the package is in the docstring of __init__.py so
that one can get help on the Python prompt:  
\>\>\> import jams  
\>\>\> help(jams)

The individual functions also provide their help as doctrings.  
Getting help, for example, on fread for reading numbers from an ascii file:  
\>\>\> import jams  
\>\>\> help(jams.fread)

or  
\>\>\> from jams import fread  
\>\>\> help(fread)

One can produce html versions of the documentation in a directory
called html by calling the script:

    ./bin/makehtml jams
    ./bin/makehtml jams/eddybox

This script works only for the one directory and not on subdirectories, i.e. subpackages.  
A more complete documentation is therefore given by opening a documentation server with

    pydoc -p 1024

and then call the url

    http://localhost:1024/jams.html

in a web browser.


---------------------------------------------------------------

### Dependencies

The package is compatible with Python 2 (> 2.6) and 3 (> 3.2).  
Note that all packages but one used by the package are already Python 3 ready (except pyhdf).  
Reading HDF4 files is thus disabled in python 3.

The package uses different third-party packages if they are
installed. Otherwise the functions are disabled. For example, reading
of special files needs netCDF4 for netcdf and pyhdf/h5py for
HDF4/HDF5, resp. If the latter are not installed, the functions
readhdf, readhdf4, readhdf5 are disabled.

Essential third-party packages, which are given as dependencies in setup.py, are numpy, scipy, netcdf4, and matplotlib.

The full list of all third-party packages used is:  
bottleneck, bs4, cartopy, cdsapi, ecmwfapi, h5py, matplotlib, mpi4py, netCDF4, networkx, numpy, pandas, pygraphviz, pyhdf, schwimmbad, scipy, seaborn, tkinter, wx, xlrd.


---------------------------------------------------------------

###  Content

#### Provided functions and modules (alphabetic w/o obsolete functions)
| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | abc2plot | Write a, b, c, ... on plots. |
| | alpha\_equ\_h2o | Equilibrium fractionation between liquid water and vapour. |
| | alpha\_kin\_h2o | Kinetic fractionation of molecular diffusion of water vapour. |
| | apply\_undef | Use a function on masked arguments. |
| | area\_poly | Area of a polygon. |
| | argsort | Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables. |
| | around | Round to the passed power of ten. |
| | ascii2ascii | Convert date notations between to ascii date format DD.MM.YYYY hh:mm:ss. |
| | ascii2en | Convert date notations to English date format YYYY-MM-DD hh:mm:ss. |
| | ascii2fr | Convert date notations to French date format DD/MM/YYYY hh:mm:ss. |
| | ascii2us | Convert date notations to American date format MM/DD/YYYY hh:mm:ss. |
| | astr | Wrapper for autostring. |
| | autostring | Format number (array) with given decimal precision. |
| | baseflow | Calculate baseflow from discharge timeseries |
| | cellarea | Calc areas of grid cells in m^2. |
| | clockplot | The clockplot of mHM. |
| | closest | Index in array which entry is closest to a given number. |
| | color | Module with color functions for plotting. |
| | const | Provides physical, mathematical, computational, and isotope constants. |
| | convex\_hull | Calculate subset of points that make a convex hull around a set of 2D points. |
| | correlate | Computes the cross-correlation function of two series x and y. |
| | cuntz\_gleixner | Cuntz-Gleixner model of 13C discrimination. |
| | dag | Generation and plotting of (connected) directed acyclic graphs with one source node. |
| | dielectric\_water | Dielectric constant of liquid water. |
| | directories\_from\_gui | Open directory selection dialogs, returns consecutiveley selected directories. |
| | directory\_from\_gui | Open directory selection dialog, returns selected directory. |
| | delta\_isogsm2 | Calculate delta values from downloaded IsoGSM2 data. |
| | dewpoint | Calculates the dew point from ambient humidity. |
| | date2dec | Converts arrays with calendar date to decimal date. |
| | dec2date | Converts arrays with decimal date to calendar date. |
| | dfgui | A minimalistic GUI for analyzing Pandas DataFrames based on wxPython. |
| | distributions | Module for pdfs of additional distributions. |
| | div | Wrapper for division. |
| | division | Divide two arrays, return 'otherwise' if division by 0. |
| | dumpnetcdf | Convenience function for writenetcdf |
| | eddybox | Module containing Eddy Covaraince utilities, see eddysuite.py for details |
| | eddysuite | Example file for processing Eddy data with eddybox and EddySoft |
| | elementary\_effects | Morris measures mu, stddev and mu* |
| | ellipse\_area | Area of ellipse (or circle) |
| | encrypt | Module to encrypt and decrypt text using a key system as well as a cipher. |
| | en2ascii | Convert date notations from English YYYY-MM-DD to ascii date format DD.MM.YYYY hh:mm:ss. |
| | errormeasures | Definition of different error measures. |
| | esat | Calculates the saturation vapour pressure of water/ice. |
| | fftngo | Fast fourier transformation for dummies (like me) |
| | Field | Generates random hydraulic conductivity fields. |
| | files | Module with file list function. |
| | file\_from\_gui | Open file selection dialog for one single file, returns selected files |
| | files\_from\_gui | Open file selection dialog, returns selected files |
| | fill\_nonfinite | Fill missing values by interpolation. |
| | find\_in\_path | Look for file in system path. |
| | Filtered\_Incompr\_Field | Generates random filtered velocity fields. |
| | fr2ascii | Convert date notations from French DD/MM/YYYT to ascii date format DD.MM.YYYY hh:mm:ss. |
| | fread | Reads in float array from ascii file. |
| | fsread | Simultaneous read of float and string array from ascii file. |
| | ftp | Module with functions for interacting with an open FTP connection. |
| | functions | Module with common functions that are used in curve\_fit or fmin parameter estimations. |
| | fwrite | Writes an array to ascii file |
| | gap2lai | Calculation of leaf area index from gap probability observations. |
| | geoarray | Pythonic gdal wrapper |
| | get\_angle | Returns the angle in radiant from each point in xy1 to each point in xy2. |
| | get\_brewer | Registers and returns Brewer colormap. |
| | get\_era\_interim | Download ERA-Interim data suitable to produce MuSICA input data. |
| | get\_era5 | Download ERA5 data suitable to produce MuSICA input data. |
| | get\_isogsm2 | Get IsoGSM2 output. |
| | get\_nearest | Returns a value z for each point in xy near to the xyz field. |
| | grid\_mid2edge | Longitude and latitude grid edges from grid midpoints. |
| | hdfread | Wrapper for readhdf. |
| | hdf4read | Wrapper for readhdf4. |
| | hdf5read | Wrapper for readhdf5. |
| | head | Return list with first n lines of file. |
| | heaviside | Heaviside (or unit step) operator. |
| | homo\_sampling | Generation of homogeneous, randomly distributed points in a given rectangular area. |
| | Incompr\_Field | Generates random velocity fields. |
| | in\_poly | Determines whether a 2D point falls in a polygon. |
| | inpoly | Wrapper for in\_poly. |
| | int2roman | Integer to roman numeral conversion. |
| | interpol | One-dimensional linear interpolation on first dimension. |
| | intersection | Intersection of two curves from x,y coordinates. |
| | jab | Jackknife-after-Bootstrap error. |
| | jConfigParser | Extended Python ConfigParser. |
| | kernel\_regression | Multi-dimensional non-parametric regression. |
| | kernel\_regression\_h | Optimal bandwidth for kernel regression. |
| | kriging | Krig a surface from a set of 2D points. |
| | lagcorr | Calculate time lag of maximum or minimum correlation of two arrays. |
| | lat\_fmt | Set lat label string (called by Basemap.drawparallels) if LaTeX package clash. |
| | leafmodel | Model to compute photosynthesis and stomatal conductance of canopies. |
| | leafprojection | Calculation of leaf projection from leaf angle observations. |
| | level1 | Module with functions dealing with CHS level1 data files, data and flags. |
| | lhs | Latin Hypercube Sampling of any distribution without correlations. |
| | lif | Count number of lines in file. |
| | line\_dev\_mask | Maskes elements of an array deviating from a line fit. |
| | logtools | Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle. |
| | lon\_fmt | Set lon label string (called by Basemap.drawmeridians) if LaTeX package clash. |
| | lowess | Locally linear regression in n dimensions. |
| | mad | Median absolute deviation test. |
| | mat2nc | Converts Matlab file *.mat into NetCDF *.nc. |
| | means | Calculate daily, monthly, yearly, etc. means of data depending on date stamp. |
| | morris\_sampling | Sampling of optimised trajectories for Morris measures / elementary effects |
| | nc2nc | Copy netcdf file deleting, renaming, replacing variables and attribues. |
| | ncread | Wrapper for readnetcdf. |
| | netcdfread | Wrapper for readnetcdf. |
| | netcdf4 | Convenience layer around netCDF4 |
| | outlier | Rossner''s extreme standardized deviate outlier test. |
| | pack | Similar to Fortran pack function with mask. |
| | pareto\_metrics | Performance metrics to compare Pareto fronts. |
| | pca | Principal component analysis (PCA) upon the first dimension of an 2D-array. |
| | pi | Parameter importance index PI or alternatively B index calculation. |
| | plot | Module with code snippets for plotting. |
| | plot\_brewer | Plots available Brewer color maps in pdf file. |
| | position | Position arrays of subplots to be used with add\_axes. |
| | print\_brewer | Prints available Brewer colormap names. |
| | pritay | Daily reference evapotranspiration after Priestley & Taylor. |
| | pso | Particle swarm optimization |
| | qa | Module of quality (error) measures. |
| | readhdf | Reads variables or information from hdf4 and hdf5 files. |
| | readhdf4 | Reads variables or information from hdf4 files. |
| | readhdf5 | Reads variables or information from hdf5 file. |
| | readnc | Wrapper for readnetcdf. |
| | readnetcdf | Reads variables or information from netcdf file. |
| | register\_brewer | Registers and registers Brewer colormap. |
| | river\_network | a class for creating a river network from a DEM including flow direction, flow accumulation and channel order |
| | rolling | Reshape an array in a "rolling window" style. |
| | roman2int | Roman numeral to integer conversion. |
| | rossner | Wrapper for outlier. |
| | t2sap | Conversion of temperature difference to sap flux density. |
| | savitzky\_golay | Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter. |
| | savitzky\_golay2d | Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter. |
| | saltelli | Parameter sampling for Sobol indices calculation. |
| | sce | Shuffle-Complex-Evolution algorithm for function min(max)imisation |
| | screening | Samples trajectories, runs model and returns measures of Morris Elemenary Effects |
| | semivariogram | Calculates semivariogram from spatial data. |
| | sendmail | Send an e-mail. |
| | sg | Wrapper savitzky\_golay. |
| | sg2d | Wrapper savitzky\_golay2d. |
| | sigma\_filter | Mask values deviating more than z standard deviations from a given function. |
| | signature2plot | Write a copyright notice on a plot. |
| | str2tex | Convert strings to LaTeX strings in math environement used by matplotlib's usetex |
| | us2ascii | Convert date notations from American MM/DD/YYYY to ascii format DD.MM.YYYY hh:mm:ss. |
| | tail | Return list with last n lines of file. |
| | maskgroup | Masks elements in a 1d array gathered in small groups. |
| | samevalue | Checks if abs. differences of array values within a certain window are smaller than threshold. |
| | savez | Save several numpy arrays into a single file in uncompressed ``.npz`` format. |
| | savez\_compressed | Save several arrays into a single file in compressed ``.npz`` format. |
| | smax | Calculating smooth maximum of two numbers |
| | smin | Calculating smooth minimum of two numbers |
| | sobol\_index | Calculates the first-order and total variance-based sensitivity indices. |
| | sread | Reads in string array from ascii file. |
| | srrasa | Generates stratified random 2D points within a given rectangular area. |
| | srrasa\_trans | Generates stratified random 2D transects within a given rectangular area. |
| | tcherkez | Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle. |
| | tee | Prints arguments on screen and in file. |
| | timestepcheck | Fills missing time steps in ascii data files |
| | tsym | Raw unicodes for common symbols. |
| | unpack | Similar to Fortran unpack function with mask. |
| | volume\_poly | Volume of function above a polygon |
| | writenetcdf | Write netCDF4 file. |
| | xkcd | Make plot look handdrawn. |
| | xlsread | Wrapper for xread. |
| | xlsxread | Wrapper for xread. |
| | xread | Simultaneous read of float and string array from Excel file. |
| | yrange | Calculates plot range from input array. |
| | zacharias | Soil water content with van Genuchten and Zacharias et al. (2007). |
| | zacharias\_check | Checks validity of parameter set for Zacharias et al. (2007). |

#### Provided functions and modules per category
 * Array manipulation
 * Ascii files
 * Data processing
 * Date & Time
 * Grids / Polygons
 * Hydrology
 * Isotopes
 * Math
 * Meteorology
 * Miscellaneous
 * Models
 * Plotting
 * Special files

**Array manipulation**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | argsort | Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables. |
| | closest | Index in array which entry is closest to a given number. |
| | pack | Similar to Fortran pack function with mask. |
| | samevalue | Checks if abs. differences of array values within a certain window are smaller than threshold. |
| | maskgroup | Masks elements in a 1d array gathered in small groups. |
| | rolling | Reshape an array in a "rolling window" style. |
| | smax | Calculating smooth maximum of two numbers |
| | smin | Calculating smooth minimum of two numbers |
| | unpack | Similar to Fortran unpack function with mask. |

**Ascii files**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | fread | Reads in float array from ascii file. |
| | fsread | Simultaneous read of float and string array from ascii file. |
| | fwrite | Writes an array to ascii file |
| | head | Return list with first n lines of file. |
| | lif | Count number of lines in file. |
| | sread | Reads in string array from ascii file. |
| | tail | Return list with last n lines of file. |

**Data processing**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | convex\_hull | Calculate subset of points that make a convex hull around a set of 2D points. |
| | eddybox | Module containing Eddy Covaraince utilities, see eddybox folder for details |
| | eddysuite | Example file for processing Eddy data with eddybox and EddySoft |
| | fill\_nonfinite | Fill missing values by interpolation. |
| | gap2lai | Calculation of leaf projection and leaf area index from gap probability observations. |
| | interpol | One-dimensional linear interpolation on first dimension. |
| | kriging | Krig a surface from a set of 2D points. |
| | kernel\_regression | Multi-dimensional non-parametric regression. |
| | kernel\_regression\_h | Optimal bandwidth for kernel regression. |
| | leafprojection | Calculation of leaf projection from leaf angle observations. |
| | level1 | Module with functions dealing with CHS level1 data files, data and flags. |
| | line\_dev\_mask | Maskes elements of an array deviating from a line fit. |
| | logtools | Module with control file functions of Logtools, the Logger Tools Software of Olaf Kolle. |
| | lowess | Locally linear regression in n dimensions. |
| | mad | Median absolute deviation test. |
| | means | Calculate daily, monthly, yearly, etc. means of data depending on date stamp. |
| | outlier | Rossner''s extreme standardized deviate outlier test. |
| | pca | Principal component analysis (PCA) upon the first dimension of an 2D-array. |
| | rossner | Wrapper for outlier. |
| | t2sap | Conversion of temperature difference to sap flux density. |
| | savitzky\_golay | Smooth (and optionally differentiate) 1D data with a Savitzky-Golay filter. |
| | savitzky\_golay2d | Smooth (and optionally differentiate) 2D data with a Savitzky-Golay filter. |
| | semivariogram | Calculates semivariogram from spatial data. |
| | sg | Wrapper savitzky\_golay. |
| | sg2d | Wrapper savitzky\_golay2d. |
| | sigma\_filter | Mask values deviating more than z standard deviations from a given function. |
| | srrasa | Generates stratified random 2D points within a given rectangular area. |
| | srrasa\_trans | Generates stratified random 2D transects within a given rectangular area. |
| | timestepcheck | Fills missing time steps in ascii data files |

**Date & Time**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | ascii2ascii | Convert date notations to ascii date format DD.MM.YYYY hh:mm:ss |
| | ascii2en | Convert date notations to English date format YYYY-MM-DD hh:mm:ss. |
| | ascii2fr | Convert date notations to French date format DD/MM/YYYY hh:mm:ss. |
| | ascii2us | Convert date notations to American date format MM/DD/YYYY hh:mm:ss. |
| | date2dec | Converts arrays with calendar date to decimal date. |
| | dec2date | Converts arrays with decimal date to calendar date. |
| | en2ascii | Convert date notations from English YYYY-MM-DD to ascii date format DD.MM.YYYY hh:mm:ss. |
| | fr2ascii | Convert date notations from French DD/MM/YYYY to ascii date format DD.MM.YYYY hh:mm:ss. |
| | us2ascii | Convert date notations from American MM/DD/YYYY to ascii format DD.MM.YYYY hh:mm:ss. |

**Grids / Polygons**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | area\_poly | Area of a polygon |
| | cellarea | Calc areas of grid cells in m^2. |
| | grid\_mid2edge | Longitude and latitude grid edges from grid midpoints. |
| | homo\_sampling | Generation of homogeneous, randomly distributed points in a given rectangular area. |
| | in\_poly | Determines whether a 2D point falls in a polygon. |
| | inpoly | Wrapper for in\_poly. |
| | volume\_poly | Volume of function above a polygon |
| | get\_angle | Returns the angle in radiant from each point in xy1 to each point in xy2. |
| | get\_nearest | Returns a value z for each point in xy near to the xyz field. |

**Hydrology**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | baseflow | Calculate baseflow from discharge timeseries |
| | river\_network | a class for creating a river network from a DEM including flow direction, flow accumulation and channel order |
    
**Isotopes**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | alpha\_equ\_h2o | Equilibrium fractionation between liquid water and vapour |
| | alpha\_kin\_h2o | Kinetic fractionation of molecular diffusion of water vapour |
| | cuntz\_gleixner | Cuntz-Gleixner model of 13C discrimination. |
| | delta\_isogsm2 | Calculate delta values from downloaded IsoGSM2 data. |
| | get\_isogsm2 | Get IsoGSM2 output. |
| | tcherkez | Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle. |

**Math**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | around | Round to the passed power of ten. |
| | correlate | Computes the cross-correlation function of two series x and y. |
| | dag | Generation and plotting of (connected) directed acyclic graphs with one source node. |
| | distributions | Module for pdfs of additional distributions. |
| | div | Wrapper for division. |
| | division | Divide two arrays, return 'otherwise' if division by 0. |
| | elementary\_effects | Morris measures mu, stddev and mu* |
| | ellipse\_area | Area of ellipse (or circle) |
| | errormeasures | Definition of different error measures. |
| | fftngo | Fast fourier transformation for dummies (like me)     |
| | functions | Module with common functions that are used in curve\_fit or fmin parameter estimations. |
| | heaviside | Heaviside (or unit step) operator. |
| | intersection | Intersection of two curves from x,y coordinates. |
| | jab | Jackknife-after-Bootstrap error. |
| | lagcorr | Calculate time lag of maximum or minimum correlation of two arrays. |
| | lhs | Latin Hypercube Sampling of any distribution without correlations. |
| | morris\_sampling | Sampling of optimised trajectories for Morris measures / elementary effects |
| | pareto\_metrics | Performance metrics to compare Pareto fronts. |
| | pi | Parameter importance index PI or alternatively B index calculation. |
| | pso | Particle swarm optimization |
| | qa | Module of quality assessment (error) measures. |
| | saltelli | Parameter sampling for Sobol indices calculation. |
| | sce | Shuffle-Complex-Evolution algorithm for function min(max)imisation |
| | screening | Samples trajectories, runs model and returns measures of Morris Elemenary Effects |
| | sobol | Generates Sobol sequences |
| | sobol\_index | Calculates the first-order and total variance-based sensitivity indices. |

**Meteorology**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | dewpoint | Calculates the dew point from ambient humidity. |
| | dielectric\_water | Dielectric constant of liquid water. |
| | esat | Calculates the saturation vapour pressure of water/ice. |
| | get\_era\_interim | Download ERA-Interim data suitable to produce MuSICA input data. |
| | get\_era5 | Download ERA5 data suitable to produce MuSICA input data. |
| | pritay | Daily reference evapotranspiration after Priestley & Taylor |    

**Miscellaneous**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | apply\_undef | Use a function on masked arguments. |
| | astr | Wrapper for autostring. |
| | autostring | Format number (array) with given decimal precision. |
| | const | Provides physical, mathematical, computational, and isotope constants. |
| | directories\_from\_gui | Open directory selection dialogs, returns consecutiveley selected directories |
| | directory\_from\_gui | Open directory selection dialog, returns selected directory |
| | encrypt | Module to encrypt and decrypt text using a key system as well as a cipher. |
| | files | Module with file list function. |
| | file\_from\_gui | Open file selection dialog for one single file, returns selected files |
| | files\_from\_gui | Open file selection dialog, returns selected files |
| | find\_in\_path | Look for file in system path. |
| | ftp | Module with functions for interacting with an open FTP connection. |
| | int2roman | Integer to roman numeral conversion. |
| | roman2int | Roman numeral to integer conversion. |
| | sendmail | Send an e-mail. |
| | tee | Prints arguments on screen and in file. |
| | zacharias | Soil water content with van Genuchten and Zacharias et al. (2007). |
| | zacharias\_check | Checks validity of parameter set for Zacharias et al. (2007). |

**Models**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | Field | Generates random hydraulic conductivity fields. |
| | Filtered\_Incompr\_Field | Generates random filtered velocity fields. |
| | Incompr\_Field | Generates random velocity fields. |
| | leafmodel | Model to compute photosynthesis and stomatal conductance of canopies |

**Plotting**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | abc2plot | Write a, b, c, ... on plots. |
| | clockplot | The clockplot of mHM. |
| | color | Module with color functions for plotting. |
| | dfgui | A minimalistic GUI for analyzing Pandas DataFrames based on wxPython. |
| | get\_brewer | Registers and returns Brewer colormap. |
| | lat\_fmt | Set lat label string (called by Basemap.drawparallels) if LaTeX package clash. |
| | lon\_fmt | Set lon label string (called by Basemap.drawmeridians) if LaTeX package clash. |
| | plot | Module with code snippets for plotting. |
| | plot\_brewer | Plots available Brewer color maps in pdf file. |
| | position | Position arrays of subplots to be used with add\_axes. |
| | print\_brewer | Prints available Brewer colormap names. |
| | register\_brewer | Registers and registers Brewer colormap. |
| | signature2plot | Write a copyright notice on a plot. |
| | str2tex | Convert strings to LaTeX strings in math environement used by matplotlib's usetex |
| | tsym | Raw unicodes for common symbols. |
| | xkcd | Make plot look handdrawn. |
| | yrange | Calculates plot range from input array. |

**Special files**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | dumpnetcdf | Convenience function for writenetcdf |
| | geoarray | Pythonic gdal wrapper |
| | hdfread | Wrapper for readhdf. |
| | hdf4read | Wrapper for readhdf4. |
| | hdf5read | Wrapper for readhdf5. |
| | jConfigParser | Extended Python ConfigParser. |
| | mat2nc | Converts Matlab file *.mat into NetCDF *.nc. |
| | nc2nc | Copy netcdf file deleting, renaming, replacing variables and attribues. |
| | ncread | Wrapper for readnetcdf. |
| | netcdfread | Wrapper for readnetcdf. |
| | netcdf4 | Convenience layer around netCDF4 |
| | readhdf | Reads variables or information from hdf4 and hdf5 files. |
| | readhdf4 | Reads variables or information from hdf4 files. |
| | readhdf5 | Reads variables or information from hdf5 file. |
| | readnc | Wrapper for readnetcdf. |
| | readnetcdf | Reads variables or information from netcdf file. |
| | savez | Save several numpy arrays into a single file in uncompressed ``.npz`` format. |
| | savez\_compressed | Save several arrays into a single file in compressed ``.npz`` format. |
| | writenetcdf | Write netCDF4 file. |
| | xlsread | Wrapper for xread. |
| | xlsxread | Wrapper for xread. |
| | xread | Simultaneous read of float and string array from Excel file. |


---------------------------------------------------------------

###  License

This file is part of the JAMS Python package, distributed under the MIT License.

Copyright (c) 2012-2019 Matthias Cuntz, Juliane Mai, Stephan Thober, Arndt Piayda

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
