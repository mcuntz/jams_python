#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
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
    This file is part of the JAMS Python package.

    It is NOT released under the GNU Lesser General Public License, yet.

    If you use this routine, please contact Arndt Piayda.

    Copyright 2014 Arndt Piayda, Matthias Cuntz


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
