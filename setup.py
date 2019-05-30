#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
jams: JAMS Python Utilities

JAMS is a general package offering miscellaneous functions in different categories,
such as reading different file formats, julian date routines, or meteorological functions.

It has several subpackages offering constants, working with Eddy covariance data and EddySoft,
offering special functions, or objective functions be used with scipy.optimize.fmin or scipy.optimize.curvefit,
and much more.

Copyright (c) 2012-2019 Matthias Cuntz, Juliane Mai, Stephan Thober, Arndt Piayda
"""
DOCLINES = __doc__.split("\n")

readme = open('README').read()

import sys
if sys.version_info[:2] < (2, 6) or (3, 0) <= sys.version_info[0:2] < (3, 2):
    raise RuntimeError("Python version 2.6, 2.7 or >= 3.2 required.")

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: End Users/Desktop
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Natural Language :: English
Operating System :: MacOS
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 3
Topic :: Scientific/Engineering
Topic :: Software Development
Topic :: Utilities
"""

MAJOR               = 5
MINOR               = 0
MICRO               = 0
ISRELEASED          = True
VERSION             = '{:d}.{:d}.{:d}'.format(MAJOR, MINOR, MICRO)

from setuptools import setup, find_packages

metadata = dict(
    name = 'jams',
    version=VERSION,
    maintainer = "Matthias Cuntz",
    maintainer_email = "mc (at) macu (dot) de",
    description = DOCLINES[0],
    # long_description = "\n".join(DOCLINES[2:]),
    long_description = readme,
    url = "https://github.com/mcuntz/jams_python.git",
    author = "JAMS = Matthias Cuntz, Juliane Mai, Stephan Thober, Arndt Piayda",
    author_email = "mc (at) macu (dot) de",
    license = 'MIT -  see LICENSE',
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
    include_package_data=True,
    scripts=['jams/bin/makehtml'],
    install_requires=['numpy, scipy'],
    packages=find_packages(exclude=['templates', 'tests*']),
    )

setup(**metadata)
