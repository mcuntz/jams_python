#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Provides physical, mathematical, computational, and isotope constants.


    Definition
    ----------
    Pi = 3.141592653589793238462643383279502884197
    ...
    Define the following constants:
        Mathematical
            Pi, Pi2, Pi3, TwoPi, Sqrt2
        Physical
            Gravity, T0, P0, T25, sigma, R, Na, REarth
        Isotope
            R13VPDB, R18VSMOW, R2VSMOW
        Computational
            tiny, huge, eps
        Material
            mmol_co2, mmol_h2o, mmol_air
            density_quartz, cheat_quartz, cheat_water, cheat_air, latentheat_vaporization


    Examples
    --------
    area = jams.const.Pi * radius**2


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014-2019 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Oct 2014
"""
from .const import Pi, Pi2, Pi3, TwoPi, Sqrt2
from .const import Gravity, T0, P0, T25, sigma, R, Na, REarth
from .const import mmol_co2, mmol_h2o, mmol_air
from .const import density_quartz, cheat_quartz, cheat_water, cheat_air, latentheat_vaporization
from .const import R13VPDB, R18VSMOW, R2VSMOW
from .const import tiny, huge, eps

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.1'
__revision__ = "Revision: 2071"
__date__     = 'Date: 24.03.2015'
