
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

To use it, either checkout the git repository

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
>>> import jams  
>>> help(jams)

The individual functions also provide their help as doctrings.  
Getting help, for example, on fread for reading numbers from an ascii file:  
>>> import jams  
>>> help(jams.fread)

or  
>>> from jams import fread  
>>> help(fread)

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

##  License

This file is part of the JAMS Python package.

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
