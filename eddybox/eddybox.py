#!/usr/bin/env python

'''
    UFZ Eddy Covariance utilities 
    
    Get help on each function by typing
    >>> help()
    help> eddybox.function
    Or
    >>> from ufz import eddybox
    >>> help(eddybox.function)
    
    
    Provided functions (alphabetic w/o obsolete)
    ------------------
    eddycorr               Calculate time lags between wind and concentrations for EddyFlux.
    eddyspec               Performs spectrum analysis with EddySpec and SpecMean and determines inductances.
    energyclosure          Computes energy closure and correction for Eddy covaraince data
    fluxfill               Wrapper function for gapfill with file management and plotting.
    fluxflag               Quality flag calculation for Eddy Covariance data
    fluxpart               Wrapper function for nee2gpp including file management and plotting
    fluxplot               Plotting routine for Eddy Covariance or other ascii data file
    gapfill                Gapfill Eddy flux data.
    itc                    Calculation of integral turbulence characteristics after Thomas & Foken (2002)
    meteo4slt              EddyFlux supply with meteorological data.
    nee2gpp                Photosynthesis and ecosystem respiration from NEE Eddy flux data.
    nee2gpp_falge          nee2gpp using one fit for whole time period
    nee2gpp_lasslop        nee2gpp using the daytime method of Lasslop et al. (2010)
    nee2gpp_reichstein     nee2gpp using several fits as in Reichstein et al. (2005)
    planarfit              Planar fit of Eddy Covariance wind components
    profile2storage        Calculate storage fluxes from profile data to correct eddy data
    sltclean               Moves *.slt files in a deleted folder to exclude from processing (EddySoft files).
    spikeflag              Spike detection for Eddy Covariance data (and basically all other data)
    ustarflag              Friction velocity flagging for Eddy Covariance data
    
    
    Example
    -------
    see eddysuite.py
    
    
    License
    -------
    This file is part of the UFZ Python library.
    
    It is NOT released under the GNU Lesser General Public License, yet.
    
    If you use this routine, please contact Arndt Piayda.
    
    Copyright 2014 Arndt Piayda, Matthias Cuntz
    
    
    History
    -------
    Written  AP, Sep 2014
'''

from eddycorr          import eddycorr
from eddyspec          import eddyspec
from energyclosure     import energyclosure
from fluxfill          import fluxfill
from fluxflag          import fluxflag
from fluxpart          import fluxpart
from fluxplot          import fluxplot
from gapfill           import gapfill
from itc               import itc
from meteo4slt         import meteo4slt
from nee2gpp           import nee2gpp, nee2gpp_falge, nee2gpp_lasslop, nee2gpp_reichstein
from planarfit         import planarfit
from profile2storage   import profile2storage 
from sltclean          import sltclean
from spikeflag         import spikeflag 
from ustarflag         import ustarflag 

if __name__ == '__main__':
    import doctest
    doctest.testmod()