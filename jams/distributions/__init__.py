#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    distributions.py        -> pdfs of additional distributions.
    sample_distributions.py -> sampling from additional distributions.


    Provided distributions
    ----------------------
    ep          Exponential Power of Box and Tiao (1992)
    exponential Exponential
    gauss       Gauss (Normal)
    laplace     Laplace
    multigauss  Multivariate Gauss (Normal) probability density function (pdf).
    multinorm   Multivariate normal (Gauss) probability density function (pdf).
    multinormal Multivariate normal (Gauss) probability density function (pdf).
    norm        Normal (Gauss)
    normal      Normal (Gauss)
    sep         Skew Exponential Power
    sep_fs      Skew Exponential Power after Fernandez & Steel (1998)
    st          Skew Student t
    st_fs       Skew Student t after Fernandez & Steel (1998)
    t           Student t


    Provided sampling from distributions
    ----------------------
    sample_ep      sample from exponential Power of Box and Tiao (1992)
    sample_sep     sample from skew exponential power distribution
    sample_sep_fs  sample from skew exponential power distribution after Fernandez and Steel (1998)
    sample_st      Sample from skew Student t
    sample_st_fs   Sample from skew Student t after Fernandez & Steel (1998)
    sample_t       Sample from Student t


    Examples
    --------
    >>> import numpy as np
    >>> print(str(laplace(0.)))
    0.5
    >>> print(str(laplace(0., 2., 2.) - 0.25/np.e))
    0.0


    License
    -------
    This file is part of the JAMS Python package.

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

    Copyright 2016 Matthias Cuntz, Juliane Mai, Dmitri Kavetski


    History
    -------
    Written,  MC,       May 2016
    Modified, JM+DK+MC, May 2016  - sampling from distributions
    Modified, MC,       Dec 2017  - multinormal
"""

from .distributions        import exponential, laplace
from .distributions        import gauss, normal, norm
from .distributions        import multigauss, multinormal, multinorm
from .distributions        import ep, sep, sep_fs, sep_fs_mean, sep_fs_std
from .distributions        import st, st_fs, st_fs_mean, st_fs_std, t

from .sample_distributions import sample_ep, sample_sep, sample_sep_fs
from .sample_distributions import sample_st, sample_st_fs, sample_t

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.2'
__revision__ = "Revision: 20171222"
__date__     = 'Date: 22.12.2017'
