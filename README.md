
# This is the Python package of JAMS.

Created June 2009 by Matthias Cuntz at the  
Department Computational Hydrosystems  
Helmholtz Centre for Environmental Research - UFZ  
Permoserstr. 15, 04318 Leipzig, Germany

Copyright 2009-2019 JAMS  
Contact Matthias Cuntz - mc (at) macu (dot) de

---------------------------------------------------------------

The library is maintained with a git repository at:

    https://bitbucket.org/mcuntz/jams_python

The package has to be in your Python path. For example in bash:

    export PYTHONPATH=/path/to/the/jams/package

It can also be installed with the usual setup.py commands using distutils:

    python setup.py install

If one wants to use the development capabilities of setuptools, you can use something like

    python -c "import setuptools; execfile('setup.py')" develop

This basically creates an .egg-link file and updates an easy-install.pth file so that the project
is on sys.path by default.

Distutils also allows to make Windows installers with
    python setup.py bdist_wininst


The documentation of the package is in the docstring of __init__.py so that one can get help on the
Python prompt:  
>>> import jams  
>>> help(jams)

The individual functions also provide their help as doctrings.  
Getting, for example, help on fread for reading floats from ascii files:  
>>> import jams  
>>> help(jams.fread)

One can produce html versions of the documentation in a directory html by calling the script:

    ./jams/bin/makehtml

This script works only for the basic package and not the subpackages.  
A more complete documentation is therefore given with by opening a documentation server with

    pydoc -p 1024

and then call the url

    http://localhost:1024/jams.html

in a web browser.


The package is compatible with Python 2 (> 2.6) and 3 (> 3.2).  
Note that all packages but one used by the package are already Python 3 ready (except pyhdf).  
Reading HDF4 files is thus disabled in python 3.

Essential third-party packages are numpy and scipy.  
Reading of special files needs pyhdf, h5py and netCDF4 for HDF4, HDF5 and NetCDF files, resp.  
Some functions provide visual checks using matplotlib for plotting.  
mad.py uses bottleneck, if available.  
dag.py uses graphviz, if available.  

---------------------------------------------------------------

##  License

This file is part of the JAMS Python package.

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

Copyright 2009-2019 Matthias Cuntz
