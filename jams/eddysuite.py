#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
from jams import eddybox as eb
from jams.timestepcheck import timestepcheck

'''
    Example file for processing Eddy Covariance data with eddybox and EddySoft.


    License
    -------
    This file is part of the JAMS Python package.

    Copyright (c) 2014 Arndt Piayda

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
    Written  AP, Sep 2014
'''

'''
0. Folder structure for eddy data processing:
e.g. |-2013                             :containing everything
       |-config                 :containg log and config files
       |-flux                   :for final fluxes
         |-itc                  :for itc plots
         |-spike                        :for spike plots
         |-ustar                        :for ustar plots and thresholds
       |-lags                   :for lags from eddycorr and plots
       |-meteo                  :containing meteo files
       |-online                 :containing online calculated files (unused)
       |-profile         :containing profile data
       |-pfit                   :for planar fit plots and files
       |-raw                            :for uncorrected raw fluxes
       |-slt                            :containing slt files
         |-deleted              :for slts to small
       |-spec                   :for inductance files and plots
       eddysuite.py             :this file
'''

'''
1. clean slt files
'''
#eb.sltclean('slt')

'''
2. calculate raw fluxes
Use EddyFlux to calculate raw, uncorrected fluxes
- no meteo files, no time lag, no filter, no inductances, no coord rotation
- no extra output, no spike detection, p=1000, set sonic direction correctly
- save results in raw/raw.csv
- save settings in raw/raw.cal
'''

'''
3. sync meteorological data with available *.slt data, requires the meteo file to be timestep checked
Use only full year meteo data with no missing time stamps
'''
#eb.meteo4slt('slt', 'meteo/meteo.csv', [1,2,6], 'meteo/meteo4slt.csv')
#eb.fluxplot('meteo/meteo.csv', 'meteo/meteo.pdf', units=False, plot=True)

'''
4. calculate time lags
Use EddyCorr to calculate time lags between w+c and w+h
- transform to meteorological data, every 1. value, scan rate 20 Hz
- lag: e.g. 300, the bigger, the longer it takes
- substract mean, extreme values
'''

'''
5. modell missing lags and attach lags to meteo file
You will be asked to give maximum, minimum and median lags. For moisture depdendend
modelling of the H2O lag, user is asked to give limits to used moisture values and
lags for fitting as well as initial guess for the model
'''
#eb.eddycorr('lags', 'slt', '35_corr.csv', '36_corr.csv', 'meteo/meteo4slt.csv', 'lags.csv', plot=True)

'''
6. Calculate inductances for carbon and water fluxes
Calculate inductances for multiple times within the year, e.g. ones per month
Use EddySpec in between:
- Co- and Quadrature spectra between w+t, w+c and w+h with the files shown
- transform to meteorological data, normalize to 1, NO inductances
- scan rate 20, dynamic smooth, multiply by freq., log.equid. freq., Welch
- for w+t lag=0, for w+c lag=c lag determined before (pos.), for w+h lag=h lag determined before (pos.)
Use SpecMean to average each spectrum for w+t, w+c and w+h
'''
#eb.eddyspec('spec', 'lags/35_corr.csv', 'lags/36_corr.csv', 'raw/raw.csv', 'slt', plot=True)

'''
7. Do the planar fit rotation for entire year or vegetations seasons
Use EddyPFit in between with:
- file format csv, minimum 100/sector, start at 0, width 30 deg, upper limit 1 m/s
- lower limit .1 m/s
'''
#eb.planarfit('pfit', 'raw/raw.csv', '2013.csv', plot=True)

'''
8. Use EddyFlux to calculate fluxes with all corrections
- Enable spike detection in EddyConf
- use inductances from spec folder, use PF, sect. with pfitmatrix file from pfit folder
- extra output, save the settings in flux
Do the timestepcheck and correct time stamps not in regular interval. Plot afterwards.
'''
#timestepcheck('flux', 'flux.csv', 'flux/flux_checked.csv', '01.01.2013 00:30', '01.01.2014 00:00')
#eb.fluxplot('flux/flux_checked.csv', 'flux/flux_checked.pdf', units=True, plot=True)

'''
9. Calculation of quality flags
- ustar flagging only possible with one full year data set
'''
swdr = 21       # column number of short wave downward radiation in meteo.csv
T    = 1        # column number of air temperature in meteo.csv
lat  = 50.450127 # latitude of tower position in decimal degrees
#eb.fluxflag('flux/flux_checked.csv', 'meteo/meteo.csv', 'flux', swdr, T, lat, spike_f=True, itc_f=True, ustar_f=True, plot=True)

'''
10. Calculation of storage fluxes
'''
#eb.profile2storage('flux/fluxflags.csv', 'flux/flux_checked.csv', 'profile/profile.csv', 'flux', heights=[0.1,0.3,1.0,2.0,5.0,9.0,15.0,23.8,30.0],
#                   CO2=[0,1,2,3,4,5,6,7,8], H2O=[9,10,11,12,13,14,15,16,17], delimiter=[',', ',', ';'], plot=True)
#eb.profile2storage('flux/flux+stor.csv', 'flux/flux_checked.csv', 'meteo/meteo.csv', 'flux', heights=[2.0,30.0], T=[2,1], plot=True)

'''
11. Gap filling
'''
rg = 19
tair = 1
rh = 5
#eb.fluxfill('flux/flux+stor.csv', 'meteo/meteo.csv', 'flux', rg, tair, rh, undef=-9999, plot=True)

'''
12. Partitioning
'''
swdr = 21
tair = 1
rh = 5
#eb.fluxpart('flux/fluxfilled.csv', 'meteo/meteo.csv', 'flux', swdr, tair, rh, method='reichstein', nogppnight=True, undef=-9999, plot=True)

'''
13. Energy balance closure
'''
swdr = 21
Rn = 32
G = 33
#eb.energyclosure('flux/fluxpart.csv', 'meteo/meteo.csv', 'flux', Rn, G, swdr, method='year', undef=-9999, plot=True)

'''
14. You're done! :-)
'''



