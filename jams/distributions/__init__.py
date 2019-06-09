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
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2016-2017 Matthias Cuntz, Juliane Mai, Dmitri Kavetski - mc (at) macu (dot) de

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
