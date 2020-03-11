#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def climate_index_knoben(time, precip, tave, pet, snow=None, color=True, indicators=False):
    """
    Calculates climate index (color and/or indicators) based on Knoben et al. (2018).

    Addition: The fraction of precipitation as snow is based on snow time series given rather 
    than using Eq. 3 of Knoben et al. (2018).

    The daily precipitation, avergage daily temperature, and potential evapotranspiration 
    is converted into climate indicators aridity I_m, Seasonality I_m,r, and Precipitation 
    as snow f_S. These three are used to derive a color which is used to identify the 
    continuous climatic zones.

    Knoben, W. J. M., Woods, R. A., & Freer, J. E. (2018). 
    A Quantitative Hydrological Climate Classification Evaluated With Independent Streamflow Data. 
    Water Resources Research, 54(7), 5088-5109. 
    http://doi.org/10.1029/2018WR022913


    Definition
    ----------
    climate_index_knoben(time, precip, tave, pet, snow=None, color=True, indicators=False):


    Input
    -----
    time         array of datetime objects of time stamps
    precip       array of daily precipitation (snow+rain)
    tave         array of daily average temperature
    pet          array of daily potential evapotranspiration


    Options
    -------
    snow         array of daily average snow (uses this estimate rather 
                 than Eq. 3 of Knoben et al. (2018).
    color        triplet of RGB values for 
                 - aridity (red), 
                 - seasonality (green), and 
                 - precipitation as snow (blue) (range: [0,1]) 
    indictators  triplet of indicator values for 
                 - aridity (arid to wet) [-1,1], 
                 - seasonality (constant to seasonal) [0,2], and 
                 - precipitation as snow (no snow to all snow) [0,1] 

    Output
    ------
    Dictionary of outputs selected:
    color = True:
       {color: {red: <red>, green: <green>, blue: <blue>}}
    indicators = True:
       {indicators: {aridity: <aridity>, seasonality: <seasonality>, precip_as_snow: <precip_as_snow>}}


    Restrictions
    ------------
    None.


    Examples
    --------
    # Climate index (color and/or indicators) based on Knoben et al. (2018).
    >>> import numpy as np
    >>> import datetime as datetime
    >>>
    >>> np.random.seed(1)
    >>> time   = np.array([ datetime.datetime(1950,1,1,0,0)+datetime.timedelta(itime) for itime in range(1000) ])
    >>> precip = np.random.random(1000)*20   
    >>> tave   = np.random.random(1000)*20-10.   
    >>> pet    = np.random.random(1000)*100   
    >>>
    >>> climate_index = climate_index_knoben(time, precip, tave, pet, color=True, indicators=True)
    >>>
    >>> from autostring import astr
    >>> print(astr([climate_index['color']['red'],climate_index['color']['green'],climate_index['color']['blue']],3,pp=True))
    ['0.896' '0.039' '0.165']
    >>> print(astr([climate_index['indicators']['aridity'],climate_index['indicators']['seasonality'],climate_index['indicators']['precip_as_snow']],3,pp=True))
    ['-0.792' ' 0.077' ' 0.165']


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2020 Juliane Mai - juliane.mai@uwaterloo.ca 

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
    Written,  JM, Feb 2020
    """

    import numpy  as np
    import pandas as pd

    # put all in dataframe to make averaging easier
    data = np.transpose( np.array([ precip, tave, pet]) )
    df = pd.DataFrame(data, columns = ['precip', 'tave', 'pet'], index=pd.DatetimeIndex(time), dtype=float)

    # monthly sum of precip
    precip_total_monthly     = df.precip.resample("M").agg(['sum'])                                  # 12 x 61 values
    precip_total_monthly_ave = precip_total_monthly.groupby(precip_total_monthly.index.month).mean() # 12 values

    # average monthly temperature
    tave_ave_monthly     = df.tave.resample("M").agg(['mean'])                                       # 12 x 61 values
    tave_ave_monthly_ave = tave_ave_monthly.groupby(tave_ave_monthly.index.month).mean()             # 12 values

    # average monthly pet
    pet_total_monthly     = df.pet.resample("M").agg(['sum'])                                        # 12 x 61 values
    pet_total_monthly_ave = pet_total_monthly.groupby(pet_total_monthly.index.month).mean()          # 12 values

    # >>>>>>>>>>>>>> TODO: write these three to file

    # Equation 1 in Knoben et al. (2018): moisture index MI
    nmonths = 12
    MI = np.zeros(nmonths)
    for itime in range(nmonths):
        if precip_total_monthly_ave.values[itime] > pet_total_monthly_ave.values[itime]:
            MI[itime] = 1.0 - pet_total_monthly_ave.values[itime] / precip_total_monthly_ave.values[itime]
        elif precip_total_monthly_ave.values[itime] < pet_total_monthly_ave.values[itime]:
            MI[itime] = precip_total_monthly_ave.values[itime] / pet_total_monthly_ave.values[itime] - 1.0
        else:
            MI[itime] = 0.0

    # Equation 2 in Knoben et al. (2018)
    #  ---> Red   = Aridity I_m
    #       -1 = arid
    #        1 = wet
    I_m = np.mean(MI)

    # Equation 3 in Knoben et al. (2018)
    #  ---> Green = Seasonality I_m,r
    #       0 = constant
    #       2 = seasonal
    I_mr = np.max(MI) - np.min(MI)

    if (snow is None):
        # Equation 4 in Knoben et al. (2018)
        #  ---> Blue  = Fraction of precipitation as snow f_S
        #       0 = no snow
        #       1 = all snow
        f_S = np.sum( precip_total_monthly_ave.values[tave_ave_monthly_ave.values <= 0.0] ) / np.sum( precip_total_monthly_ave.values )
    else:
        # if snow farction is given, ratio of total snow amount and total precipitation amount is used to derive snow fraction
        f_S = np.sum( snow ) / np.sum( precip_total_monthly.values )
        # print('snow amount: ',np.sum( snow ))
        # print('prec amount: ',np.sum( precip_total_monthly.values ))

    # print("Aridity I_m                           = ", I_m)
    # print("Seasonality I_m,r                     = ", I_mr)
    # print("Fraction of precipitation as snow f_S = ", f_S)


    # color mapping
    # R = 1 - (aridity_Im + 1)./2       --> transforms [-1,1] to [1,0] --> [most arid, least arid]
    # G = (aridity_seasonality_Imr)./2  --> transforms [ 0,2] to [0,1] --> [least seasonal, most seasonal]
    # B = precipitation_as_snow_fs      --> transforms [ 0,1] to [0,1] --> [least snow, most snow]
    red   = 1.0 - (I_m + 1.0)/ 2.0
    green = I_mr/ 2.0
    blue  = f_S

    # print("Red   = ", red)
    # print("Green = ", green)
    # print("Blue  = ", blue)

    output = {}
    if color:
        output['color'] = {'red': red, 'green': green, 'blue': blue}
    if indicators:
        output['indicators'] = { 'aridity': I_m, 'seasonality': I_mr, 'precip_as_snow': f_S}

    return output

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
