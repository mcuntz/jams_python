#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
'''
    leafmodel: contains routines to model photosynthesis and stomatal conductance
    of canopies or leaves.
    
    
    major functions contained:
    --------------------------
    opt_leafmodel:  optimization against observed transpiration and
                    gross primary production
    cost_leafmodel: cost function for optimization
    leafmodel:      model for photosynthesis and stomatal conductance
    twoleafmodel:   uses leafmodel in a sunlit and shaded scheme
    farquhar:       photosysnthesis model
    leuning:        stomatal conductance modeling
    ball_berry:     stomatal conductance modeling


    Examples
    --------
    ci, gcmol_c, gcmol_h, gsmol_c, gsmol_h, no_conv, Jc, Je, Rd, GPPmod, Emod = jams.leafmodel.leafmodel(
        ci_ini, Tl, PAR, ea, ca, ga, gb, P, par1, par2, gs_method=gs_method, temp_method=temp_method, n=n, eta=eta)


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014 Arndt Piayda, Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  AP & MC, Jul 2014
'''
from .leafmodel import *

# Information
__author__   = "Arndt Piayda, Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1796"
__date__     = 'Date: 05.10.2014'
