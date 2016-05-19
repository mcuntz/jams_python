#!/usr/bin/env python
"""
    Provide pdfs of additional distributions.


    Provided distributions
    ----------------------
    ep          Exponential Power
    exponential Exponential
    gauss       Gauss (Normal)
    laplace     Laplace
    normal      Normal (Gauss)
    sep         Skew Exponential Power
    sstudentt   Skewed Student t
    studentt    Student t


    Example
    -------
    >>> import numpy as np
    >>> print(str(laplace(0.)))
    0.5
    >>> print(str(laplace(0., 2., 2.) - 0.25/np.e))
    0.0


    License
    -------
    This file is part of the UFZ Python package.

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

    Copyright 2016 Matthias Cuntz


    History
    -------
    Written,  MC, May 2016
"""

# UFZ colours
from .distributions import ep, exponential, gauss, laplace, normal, sep, sstudentt, studentt

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 2681"
__date__     = 'Date: 12.05.2016'
