#!/usr/bin/env python
"""
    distributions.py        --> Provide pdfs of additional distributions.
    sample_distributions.py --> Provide sampling from additional distributions.

    Provided distributions
    ----------------------
    ep          Exponential Power
    exponential Exponential
    gauss       Gauss (Normal)
    laplace     Laplace
    normal      Normal (Gauss)
    sep         Skew Exponential Power
    ssep        Standardised skew Exponential Power
    sstudentt   Skewed Student t
    ssstudentt  Standardised skewed Student t
    studentt    Student t


    Provided sampling from distributions
    ----------------------
    sample_ep01    sample from exponential power function EP(0,1,beta)
    sample_ep      sample from (general) exponential power function EP(loc,sca,beta)
    sample_sep_fs  sample from skew exponential power distribution obtained by using the approach of Fernandez and Steel
    sample_ssep    sample from standardized skew exponential power distribution (mean=zero, std. dev.=1)

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

    Copyright 2016 Matthias Cuntz, Juliane Mai, Dmitri Kavetski


    History
    -------
    Written,  MC,    May 2016
    Modified, JM+DK, May 2016  - adding sampling from distributions
"""

from .distributions        import gauss, normal
from .distributions        import exponential, laplace
from .distributions        import ep, sep, sep_mean, sep_std, ssep
from .distributions        import ssstudentt, sstudentt, sstudentt_mean, sstudentt_std, studentt
from .sample_distributions import sample_ep01, sample_ep, sample_sep_fs, sample_ssep



# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 2681"
__date__     = 'Date: 12.05.2016'
