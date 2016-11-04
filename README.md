
# This is the Python package of JAMS.

Created June 2009 by Matthias Cuntz at the  
Department Computational Hydrosystems  
Helmholtz Centre for Environmental Research - UFZ  
Permoserstr. 15, 04318 Leipzig, Germany

Copyright 2009-2016 JAMS  
Contact Matthias Cuntz - mc (at) macu.de

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

Not all files in the package are free software. The license is given in the 'License' section
of the docstring of each routine.

There are 3 possibilities:

1. The routine is not yet released under the GNU Lesser General Public License. This is marked by a
   text such as  
        This file is part of the JAMS Python package.  
        It is NOT released under the GNU Lesser General Public License, yet.  
        If you use this routine, please contact Matthias Cuntz.  
        Copyright 2012-2013 Matthias Cuntz  
    If you want to use this routine for publication or similar, please contact the author for possible co-authorship.

2. The routine is already released under the GNU Lesser General Public License
   but if you use the routine in a publication or similar, you have to cite the respective
   publication, e.g.  
        If you use this routine in your work, you should cite the following reference  
        Goehler M, J Mai, and M Cuntz (2013)  
            Use of eigendecomposition in a parameter sensitivity analysis of the Community Land
            Model,  
            J Geophys Res 188, 904-921, doi:10.1002/jgrg.20072

3. The routine is released under the GNU Lesser General Public License. The following applies:  
   The JAMS Python package is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The JAMS Python package is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
   If not, see <http://www.gnu.org/licenses/>.

   Copyright 2009-2016 Matthias Cuntz
