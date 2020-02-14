#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def pet_oudin(temperature, lat, doy):
    """
        Calculates the daily potential evapotranspiration in [mm/d] following the Oudin formula.

        For the mathematical details of the PE formulation, see:
        Oudin, L., Hervieu, F., Michel, C., Perrin, C., Andreassian, V., Anctil, F. and Loumagne, C., 2005. 
        Which potential evapotranspiration input for a rainfall-runoff model? 
        Part 2 - Towards a simple and efficient PE model for rainfall-runoff modelling. 
        Journal of Hydrology 303(1-4), 290-306.

        For the calculation of extra-atmospheric global radiation, see
        Appendix C of the article by Morton, F.I., 1983. 
        Operational estimates of areal evapotranspiration and their significance to the science and 
        practice of hydrology. 
        Journal of Hydrology 66 (1/4), 1-76.


        Definition
        ----------
        def pet_oudin(temperature, lat, doy):


        Input
        -----
        temperature  (array of) temperature in degrees Celsius
        lat          (array of) latitudes in degrees N
        doy          (array of) day of the year (1 ... 366)


        Options
        -------
        none

        Output
        ------
        array with potential evapotranspiration (PET) based on Oudin formula in [mm/d]


        Restrictions
        ------------
        Shape of input arrays must be same


        Examples
        --------
        # PET based on Oudin
        >>> import numpy as np
        >>> temp = np.array([ 32.0,  32.0,  32.0,  32.0,       32.0,       32.0])
        >>> lat  = np.array([ 48.73, 48.73, 48.73, 9.27785325, 5.56671363, 5.56671363])
        >>> doy  = np.array([ 1,     2,     3,     364,        365,        366])
        >>> from autostring import astr
        >>> print(astr(pet_oudin(temp, lat, doy),3,pp=True))
        ['1.301' '1.309' '1.317' '4.844' '5.123' '5.127']


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 20020 Juliane Mai - juliane.mai@uwaterloo.ca 

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

    # check that inputs have all same shape
    ishape = np.shape(temperature)
    if not(np.shape(lat) == ishape):
        raise ValueError("pet_oudin: Array of latitudes must have same shape as temperature array")
    if not(np.shape(doy) == ishape):
        raise ValueError("pet_oudin: Array of day-of-year must have same shape as temperature array")

    # initialize
    pet_oudin = np.ones(ishape) * -9999.0

    teta  = 0.4093 * np.sin(doy/58.1 - 1.405)
    cosGz = np.maximum(0.001, np.cos(lat/57.3-teta))
    Gz    = np.arccos(cosGz)
    cosOM = np.maximum(-1.0, np.minimum(1.0, 1.0-cosGz/np.cos(lat/57.3)/np.cos(teta)))
    OM    = np.arccos(cosOM)
    Eta   = 1.0 + np.cos(doy/58.1)/30.0
    cosPz = cosGz + np.cos(lat/57.3)*np.cos(teta)*(np.sin(OM)/OM-1.0)

    # Global Radiation [W/m2]
    globrad = 446.0 * OM * Eta * cosPz

    # Potential Evapotranspiration [mm/d]
    pet_oudin = np.maximum(0.0, globrad/28.5 * (temperature+5.0)/100.)
    
    return pet_oudin

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
