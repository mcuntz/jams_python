#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Almost all control file functions of Logtools, the Logger Tools Software of
    Olaf Kolle, MPI-BGC Jena, (c) 2012.

    The following is a copy of section 4 of the logtools manual.

    Some functions are renamed compared to the original logger tools
        chs -> varchs
        add -> varadd
        sub -> varsub
        mul -> varmul
        div -> vardiv
        sqr -> varsqr/varsqrt
        exp -> varexp
        log -> varlog
        pot -> varpot

    Some functions are slightly enhanced, which is reflected in the documentation below.

    Not all functions are implemented (yet). Missing functions are:
        varset
        met_torad
        met_psy_rh
        met_dpt_rh
        write

    All functions have an additional keyword undef, which defaults to -9999.:
        elements are excluded from the calculations if any of the inputs equals undef.
    Only bit_test and the if-statements ifeq, ifne, ifle, ifge, iflt, igt do not have the undef keyword.

    The Looger Tools control functions are:
        1. Assignment # not implemented


        2. Change sign
            x = varchs(a) means x = -a, where a is a variable or a number.

            Definition
            ----------
            def varchs(var1, undef=-9999.):


        3. Addition
            x = varadd(a, b) means x = a + b, where a and b are variables or numbers.

            Definition
            ----------
            def varadd(var1, var2, undef=-9999.):


        4. Subtraction
            x = varsub(a, b) means x = a - b, where a and b are variables or numbers.

            Definition
            ----------
            def varsub(var1, var2, undef=-9999.):


        5. Multiplication
            x = varmul(a, b) means x = a * b, where a and b are variables or numbers.

            Definition
            ----------
            def varmul(var1, var2, undef=-9999.):


        6. Division
            x = vardiv(a, b) means x = a/b, where a and b are variables or numbers.

            Definition
            ----------
            def vardiv(var1, var2, undef=-9999.):


        7. Square root
            x = varsqr(a) means x = sqrt(a), where a is a variable or a number.
            x = varsqrt(a) means x = sqrt(a), where a is a variable or a number.

            Definition
            ----------
            def varsqr(var1, undef=-9999.):
            def varsqrt(var1, undef=-9999.):


        8. Exponentiation of e
            x = varexp(a) means x = exp(a), where a is a variable or a number.

            Definition
            ----------
            def varexp(var1, undef=-9999.):


        9. Natural logarithm
            x = varlog(a) means x = ln(a), where a is a variable or a number.

            Definition
            ----------
            def varlog(var1, undef=-9999.):


        10. Exponentiation
            x = varpot(a, b) means x = a**b, where a and b are variables or numbers.

            Definition
            ----------
            def varpot(var1, var2, undef=-9999.):


        11. Apply linear function
            x = lin(y, a0, a1) means x = a0 + a1 * y,
            where a0 and a1 are variables or numbers.

            Definition
            ----------
            def lin(var1, a, b, undef=-9999.):


        12. Apply 2nd order function
            x=quad(y,a0,a1,a2) means x = a0 +a1*y + a2*y**2,
            where a0, a1 and a2 are variables or numbers.

            Definition
            ----------
            def quad(var1, a, b, c, undef=-9999.):


        13. Apply 3rd order function
            x=cubic(y,a0,a1,a2,a3) means x = a0 +a1*y+a2*y**2+a3*y**3,
            where a0, a1, a2 and a3 are variables or numbers.

            Definition
            ----------
            def cubic(var1, a, b, c, d, undef=-9999.):


        14. Calculate fraction of day from hours, minutes and seconds
            x = hms(h, m, s) means x = (h + m/60 + s/3600)/24,
            where h, m and s (hours, minutes and seconds) are variables or numbers.

            Definition
            ----------
            def hms(h, m, s, undef=-9999.):


        15. Bitwise test
            x = bit_test(y, b, start=0) means x = 1 if bit b ist set in y otherwise x = 0.
            Returns a list of b is an array.
            Counting of b starts at start.
            For the behaviour of the original logger tools, set start=1.
            Negative b's is not implemented.

            Definition
            ----------
            def bit_test(var1, var2, start=0):


        16. Replacement of underflows by new value
            x = setlow(y,lo,ln=None) means IF (y > lo) THEN x = ln ELSE x = y,
            where lo and ln are variables or numbers.
            ln is optional. If not given lo will be used.
            This function may be used to adjust small negative values of short wave radiation
            during nighttime to zero values.

            Definition
            ----------
            def setlow(dat, low, islow=None, undef=-9999.):


        17. Replacement of overflows by new value
            x = sethigh(y,lo,ln=None) means IF (y < lo) THEN x = ln ELSE x = y,
            where lo and ln are variables or numbers.
            ln is optional. If not given lo will be used.
            This function may be used to adjust relative humidity values of a little bit more than 100 % to 100 %.

            Definition
            ----------
            def sethigh(dat, high, ishigh=None, undef=-9999.):


        18. Replacement of underflows or overflows by the missing value
            x = limits(y, ll, lh) means
            IF (y > ll) OR (y < lh) THEN x = [missing value] ELSE x = y,
            where ll and lh are variables or numbers.
            This function may be used to check values lying in between certain limits.
            If one of the limits is exceeded the value is set to the missing value defined in the control file.

            Definition
            ----------
            def limits(dat, mini, maxi, undef=-9999.):


        19. Calculation of mean value
            x = mean(y1, y2, ..., yn) means x = (y1 + y2 + ... + yn)/n,
            where y1, y2, ..., yn are variables or numbers.

            Definition
            ----------
            def mean(var1, axis=None, undef=-9999.):


        20. Calculation of minimum value
            x = mini(y1,y2,...,yn) means x = min(y1,y2,...,yn),
            where y1, y2, ..., yn are variables or numbers.

            Definition
            ----------
            def mini(var1, axis=None, undef=-9999.):


        21. Calculation of maximum value
            x = maxi(y1,y2,...,yn) means x = max(y1,y2,...,yn),
            where y1, y2, ..., yn are variables or numbers.

            Definition
            ----------
            def maxi(var1, axis=None, undef=-9999.):


        22. Calculation of total radiation from net radiometer # no implemented


        23. Calculation of long wave radiation from net radiometer
            x = met_lwrad(y, Tp) where
            y is the output voltage of the net radiometer in mV,
            Tp is the temperature of the net radiometer body in oC.
            The total radiation in W m-2 is calculated according to the following formula:
            x=y*fl +sigma*(Tp +273.16)**4
            where sigma = 5.67051 * 10**8 W m-2 K-4 is the Stephan-Boltzmann-Constant and
            fl is the factor for long wave radiation (reciprocal value of sensitivity) in W m-2 per mV.
            The function assumes that fl was already applied before.
            All parameters may be variables or numbers.

            Definition
            ----------
            def met_lwrad(dat, tpyr, undef=-9999.): # assumes that dat was already multiplied with calibration factor


        24. Calculation of radiation temperature from long wave radiation
            x = met_trad(Rl, epsilon) where
            Rl is the long wave radiation in W m-2,
            epsilon is the long wave emissivity of the surface (between 0 and 1).
            The radiation temperature in oC is calculated according to the following formula:
            x= sqrt4(Rl/(sigma*epsilon)) - 273.16 
            where sigma = 5.67051 * 108 W m-2 K-4 is the Stephan-Boltzmann-Constant.
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_trad(dat, eps, undef=-9999.):


        25. Calculation of albedo from short wave downward and upward radiation
            x = met_alb(Rsd, Rsu) where
            Rsd is the short wave downward radiation in Wm-2, Rsu is the short wave upward radiation in Wm-2,
            The albedo in % is calculated according to the following formula:
            x = 100 * ( Rsu / Rsd )
            If Rsd > 50 W m-2 or Rsu > 10 W m-2 the result is [missing value].
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_alb(swd, swu, swdmin=50., swumin=10., undef=-9999.):


        26. Calculation of albedo from short wave downward and upward radiation with limits
            x=met_albl(Rsd,Rsu,Rsd_limit,Rsu_limit)where
            Rsd is the short wave downward radiation in Wm-2,
            Rsu is the short wave upward radiation in Wm-2,
            Rsd_limit is the short wave downward radiation limit in Wm-2,
            Rsu_limit is the short wave upward radiation limit in Wm-2,
            The albedo in % is calculated according to the following formula:
            x = 100 * ( Rsu / Rsd )
            If Rsd > Rsd_limit or Rsu > Rsu_limit the result is [missing value].
            All four parameters may be variables or numbers.

            Definition
            ----------
            def met_albl(swd, swu, swdmin, swumin, undef=-9999.):


        27. Calculation of saturation water vapour pressure
            x = met_vpmax(T) where
            T is the air temperature in oC.
            The saturation water vapour pressure in mbar (hPa) is calculated according to the following formula:
            x = 6.1078 * exp(17.08085 * T / (234.175 + T))
            The parameter may be a variable or a number.

            Definition
            ----------
            def met_vpmax(temp, undef=-9999.):


        28. Calculation of actual water vapour pressure
            x = met_vpact(T,rh) where T is the air temperature in oC, rh is the relative humidity in %.
            The actual water vapour pressure in mbar (hPa) is calculated according to the following for- mulas:
            Es = 6.1078*exp(17.08085*T/ (234.175 + T))
            x = Es * rh/100
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_vpact(temp, rh, undef=-9999.):


        29. Calculation of water vapour pressure deficit
            x = met_vpdef(T, rh) where T is the air temperature in oC, rh is the relative humidity in %.
            The water vapour pressure deficit in mbar (hPa) is calculated according to the following for- mulas:
            Es = 6.1078*exp(17.08085*T/ (234.175 + T))
            E = Es * rh/100
            x = Es - E
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_vpdef(temp, rh, undef=-9999.):


        30. Calculation of specific humidity
            x = met_sh(T, rh, p) where
            T is the air temperature in oC,
            rh is the relative humidity in %,
            p is the air pressure in mbar (hPa).
            The specific humidity in g kg-1 is calculated according to the following formulas:
            Es = 6.1078*exp(17.08085*T/ (234.175 + T))
            E = Es * rh/100
            x = 622 * E/(p-0.378*E)
            All parameters may be variables or numbers.

            Definition
            ----------
            def met_sh(temp, rh, p, undef=-9999.):


        31. Calculation of potential temperature
            x = met_tpot(T, p) where
            T is the air temperature in oC,
            p is the air pressure in mbar (hPa).
            The potential temperature in K is calculated according to the following formula:
            x = (T + 273.16) * (1000/p)**0.286
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_tpot(temp, p, undef=-9999.):


        32. Calculation of air density
            x = met_rho(T, rh, p) where
            T is the air temperature in oC,
            rh is the relative humidity in %,
            p is the air pressure in mbar (hPa).
            The air density in kg m-3 is calculated according to the following formulas:
            Es = 6.1078*exp(17.08085*T/ (234.175 + T))
            E = Es * rh/100
            sh = 622 * E/(p-0.378*E)
            Tv = ((T + 273.16) * (1 + 0.000608 * sh)) - 273.16
            x = p * 100 / (287.05 * (Tv + 273.16))
            All parameters may be variables or numbers.

            Definition
            ----------
            def met_rho(temp, rh, p, undef=-9999.):


        33. Calculation of dew point temperature
            x = met_dpt(T, rh) where
            T is the air temperature in oC, rh is the relative humidity in %.
            The dew point temperature in oC is calculated according to the following formulas:
            Es = 6.1078*exp(17.08085*T/(234.175 + T))
            E = Es * rh/100
            x = 234.175 * ln(E/6.1078)/(17.08085 - ln(E/6.1078))
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_dpt(temp, rh, undef=-9999.):


        34. Calculation of water vapour concentration
            x = met_h2oc(T, rh, p) where T is the air temperature in oC,
            rh is the relative humidity in %,
            p is the air pressure in mbar (hPa).
            The water vapour concentration in mmol mol-1 is calculated according to the following formu- las:
            Es = 6.1078*exp(17.08085*T/ (234.175 + T))
            E = Es * rh/100
            x = 0.1 * E /(0.001*p*100*0.001)
            All parameters may be variables or numbers.

            Definition
            ----------
            def met_h2oc(temp, rh, p, undef=-9999.):


        35. Calculation of relative humidity from dry and wet bulb temperature # not implemented


        36. Calculation of relative humidity from dew point temperature # not implemented


        37. Calculation of relative humidity from water vapour concentration
            x = met_h2oc_rh(T, [H2O], p) where
            T is the air temperature in oC,
            [H2O] is the water vapour concentration in mmolmol-1, p is the air pressure in mbar (hPa).
            The relative humidity in % is calculated according to the following formulas:
            Es = 6.1078*exp(17.08085*T/(234.175 + T))
            E = 10 * [H2O] * 0.001 * p * 100 * 0.001
            x = 100 * E / Es
            All parameters may be variables or numbers.

            Definition
            ----------
            def met_h2oc_rh(temp, h, p, undef=-9999.):


        38. Rotation of wind direction
            x = met_wdrot(wd, a) where
            wd is the wind direction in degree,
            a is the rotation angle in degree (positive is clockwise).
            The rotated wind direction is calculated according to the following formulas:
            x = wd + a
            IF x > 0 THEN x = x + 360
            IF x >= 360 THEN x = x - 360
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_wdrot(wd, a, undef=-9999.):


        39. Rotation of u-component of wind vector
            x = met_urot(u, v, a) where
            u is the u-component of the wind vector,
            v is the v-component of the wind vector,
            a is the rotation angle in degree (positive is clockwise).
            The rotated u-component is calculated according to the following formula:
            x = u * cos (a) + v * sin (a)
            All three parameters may be variables or numbers.

            Definition
            ----------
            def met_urot(u, v, a, undef=-9999.):


        40. Rotation of v-component of wind vector
            x = met_vrot(u, v, a) where
            u is the u-component of the wind vector,
            v is the v-component of the wind vector,
            a is the rotation angle in degree (positive is clockwise).
            The rotated v-component is calculated according to the following formula:
            x = -u * sin (a) + v * cos (a)
            All three parameters may be variables or numbers.

            Definition
            ----------
            def met_vrot(u, v, a, undef=-9999.):


        41. Calculation of wind velocity from u- and v-component of wind vector
            x = met_uv_wv(u, v) where
            u is the u-component of the wind vector, v is the v-component of the wind vector.
            The horizontal wind velocity is calculated according to the following formula:
            x = sqrt(u**2 + v**2)
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_uv_wv(u, v, undef=-9999.):


        42. Calculation of wind direction from u- and v-component of wind vector
            x = met_uv_wd(u, v) where
            u is the u-component of the wind vector, v is the v-component of the wind vector.
            The horizontal wind velocity is calculated according to the following formulas:
            IF u = 0 AND v = 0 THEN x = 0
            IF u = 0 AND v > 0 THEN x = 360
            IF u = 0 AND v < 0 THEN x = 180
            IF u < 0 THEN x = 270 - arctan(v/u)
            IF u > 0 THEN x = 90 - arctan(v/u)
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_uv_wd(u, v, undef=-9999.):


        43. Calculation of u-component of wind vector from wind velocity and wind direction
            x = met_wvwd_u(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
            The u-component of the wind vector is calculated according to the following formula:
            x = -wv * sin (wd)
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_wvwd_u(wv, wd, undef=-9999.):


        44. Calculation of v-component of wind vector from wind velocity and wind direction
            x = met_wvwd_v(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
            The v-component of the wind vector is calculated according to the following formula:
            x = -wv * cos (wd)
            Both parameters may be variables or numbers.

            Definition
            ----------
            def met_wvwd_v(wv, wd, undef=-9999.):


        45. If-statements
            x = ifeq(y,a0,a1,a2) means IF y == a0 THEN x = a1 ELSE x = a2
            x = ifne(y,a0,a1,a2) means IF y != a0 THEN x = a1 ELSE x = a2
            x = ifle(y,a0,a1,a2) means IF y <= a0 THEN x = a1 ELSE x = a2
            x = ifge(y,a0,a1,a2) means IF y >= a0 THEN x = a1 ELSE x = a2
            x = iflt(y,a0,a1,a2) means IF y > a0 THEN x = a1 ELSE x = a2
            x = ifgt(y,a0,a1,a2) means IF y < a0 THEN x = a1 ELSE x = a2
            All parameters may be variables or numbers.

            Definition
            ----------
            def ifeq(var1, iif, ithen, ielse):
            def ifne(var1, iif, ithen, ielse):
            def ifle(var1, iif, ithen, ielse):
            def ifge(var1, iif, ithen, ielse):
            def iflt(var1, iif, ithen, ielse):
            def ifgt(var1, iif, ithen, ielse):


        46. Write variables to a file # not implemented


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Jun-Dec 2014
    Modified, CM, Jun 2014 - corrected type in met_tpot
"""
from .logtools import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1911"
__date__     = 'Date: 02.12.2014'
