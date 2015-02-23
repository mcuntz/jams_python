#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.division import division
from ufz.esat     import esat
from ufz.const    import sigma

__all__ = ['varchs', 'varadd', 'varsub', 'varmul', 'vardiv', 'varsqr',
           'varsqrt', 'varexp', 'varlog', 'varpot', 'lin', 'quad',
           'cubic', 'hms', 'bit_test', 'setlow', 'sethigh', 'limits',
           'mean', 'mini', 'maxi', 'met_lwrad', 'met_trad', 'met_alb',
           'met_albl', 'met_vpmax', 'met_vpact', 'met_vpdef', 'met_sh',
           'met_tpot', 'met_rho', 'met_dpt', 'met_h2oc', 'met_h2oc_rh',
           'met_wdrot', 'met_urot', 'met_vrot', 'met_uv_wv', 'met_uv_wd',
           'met_wvwd_u', 'met_wvwd_v', 'ifeq', 'ifne', 'ifle', 'ifge',
           'iflt', 'ifgt']


# Not implemented: varset


def varchs(var1, undef=-9999.):
    """
        Change sign
        x = varchs(a) means x = -a, where a is a variable or a number.


        Definition
        ----------
        def varchs(var1, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Changed sign

        History
        -------
        Written,  MC, Dec 2014
    """
    return np.where(var1==undef, undef, -var1)


def varadd(var1, var2, undef=-9999.):
    """
        Addition
        x = varadd(a, b) means x = a + b, where a and b are variables or numbers.


        Definition
        ----------
        def varadd(var1, var2, undef=-9999.):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Addition


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef) | (var2==undef), undef, var1 + var2)


def varsub(var1, var2, undef=-9999.):
    """
        Subtraction
        x = varsub(a, b) means x = a - b, where a and b are variables or numbers.


        Definition
        ----------
        def varsub(var1, var2, undef=-9999.):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Subtraction


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef) | (var2==undef), undef, var1 - var2)


def varmul(var1, var2, undef=-9999.):
    """
        Multiplication
        x = varmul(a, b) means x = a * b, where a and b are variables or numbers.


        Definition
        ----------
        def varmul(var1, var2, undef=-9999.):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Multiplication


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef) | (var2==undef), undef, var1 * var2)


def vardiv(var1, var2, undef=-9999.):
    """
        Division
        x = vardiv(a, b) means x = a/b, where a and b are variables or numbers.


        Definition
        ----------
        def vardiv(var1, var2, undef=-9999.):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Division


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef) | (var2==undef), undef, division(var1, var2, undef))


def varsqr(var1, undef=-9999.):
    """
        Square root
        x = varsqr(a)  means x = sqrt(a), where a is a variable or a number.


        Definition
        ----------
        def varsqr(var1, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Square root


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, np.sqrt(var1))


def varsqrt(var1, undef=-9999.):
    """
        Square root
        x = varsqrt(a) means x = sqrt(a), where a is a variable or a number.


        Definition
        ----------
        def varsqrt(var1, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Square root


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, np.sqrt(var1))


def varexp(var1, undef=-9999.):
    """
        Exponentiation of e
        x = varexp(a) means x = exp(a), where a is a variable or a number.


        Definition
        ----------
        def varexp(var1, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Exponentiation


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, np.exp(var1))


def varlog(var1, undef=-9999.):
    """
        Natural logarithm
        x = varlog(a) means x = ln(a), where a is a variable or a number.


        Definition
        ----------
        def varlog(var1, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Natural logarithm


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, np.log(var1))


def varpot(var1, var2, undef=-9999.):
    """
        Exponentiation
        x = varpot(a, b) means x = a**b, where a and b are variables or numbers.


        Definition
        ----------
        def varpot(var1, var2, undef=-9999.):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Exponentiation


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef) | (var2==undef), undef, var1**var2)


def lin(var1, a, b, undef=-9999.):
    """
        Apply linear function
        x = lin(y, a0, a1) means x = a0 + a1 * y,
        where a0 and a1 are variables or numbers.


        Definition
        ----------
        def lin(var1, a, b, undef=-9999.):


        Input
        -----
        var1    ND-array
        a0      scalar
        a1      scalar


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        linear function


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, a + b*var1)


def quad(var1, a, b, c, undef=-9999.):
    """
        Apply 2nd order function
        x=quad(y,a0,a1,a2) means x = a0 +a1*y + a2*y**2,
        where a0, a1 and a2 are variables or numbers.


        Definition
        ----------
        def quad(var1, a, b, c, undef=-9999.):


        Input
        -----
        var1    ND-array
        a0      scalar
        a1      scalar
        a2      scalar


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        2nd order function


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, a + b*var1 + c*var1*var1)


def cubic(var1, a, b, c, d, undef=-9999.):
    """
        Apply 3rd order function
        x=cubic(y,a0,a1,a2,a3) means x = a0 +a1*y+a2*y**2+a3*y**3,
        where a0, a1, a2 and a3 are variables or numbers.


        Definition
        ----------
        def cubic(var1, a, b, c, d, undef=-9999.):


        Input
        -----
        var1    ND-array
        a0      scalar
        a1      scalar
        a2      scalar
        a3      scalar


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        3rd order function


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((var1==undef), undef, a + b*var1 + c*var1*var1 + d*var1*var1*var1)


def hms(h, m, s, undef=-9999.):
    """
        Calculate fraction of day from hours, minutes and seconds
        x = hms(h, m, s) means x = (h + m/60 + s/3600)/24,
        where h, m and s (hours, minutes and seconds) are variables or numbers.


        Definition
        ----------
        def hms(h, m, s, undef=-9999.):


        Input
        -----
        h    ND-array
        m    ND-array
        s    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        fraction of day


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((h==undef) | (m==undef) | (s==undef), undef, (h+m/60.+s/3600.)/24.)


def bit_test(var1, var2, start=0):
    """
        Bitwise test
        x = bit_test(y, b, start=0) means x = 1 if bit b ist set in y otherwise x = 0.
        Returns a list of b is an array.
        Counting of b starts at start.
        For the behaviour of the original logger tools, set start=1.
        Negative b's is not implemented.


        Definition
        ----------
        def bit_test(var1, var2, start=0):


        Input
        -----
        var1    ND-array
        var2    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        Bitwise test


        History
        -------
        Written,  MC, Jun 2014
    """
    if np.size(var2) > 1:
        return [ (var1 >> i+start)%2 for i in var2 ]
    else:
        return (var1 >> var2+start)%2


def setlow(dat, low, islow=None, undef=-9999.):
    """
        Replacement of underflows by new value
        x = setlow(y,lo,ln=None) means IF (y < lo) THEN x = ln ELSE x = y,
        where lo and ln are variables or numbers.
        ln is optional. If not given lo will be used.
        This function may be used to adjust small negative values of short wave radiation
        during nighttime to zero values.


        Definition
        ----------
        def setlow(dat, low, islow=None, undef=-9999.):


        Input
        -----
        dat    ND-array
        low    ND-array


        Optional Input
        --------------
        islow   if not None, use islow in case of y < low
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        underflows replaced by new value


        History
        -------
        Written,  MC, Jun 2014
    """
    if islow is None:
        out = np.maximum(dat, low)
    else:
        out = np.where(dat < low, islow, dat)
    return np.where(dat == undef, undef, out)


def sethigh(dat, high, ishigh=None, undef=-9999.):
    """
        Replacement of overflows by new value
        x = sethigh(y,lo,ln=None) means IF (y > lo) THEN x = ln ELSE x = y,
        where lo and ln are variables or numbers.
        ln is optional. If not given lo will be used.
        This function may be used to adjust relative humidity values of a little bit more than 100 % to 100 %.


        Definition
        ----------
        def sethigh(dat, high, ishigh=None, undef=-9999.):


        Input
        -----
        dat     ND-array
        high    ND-array


        Optional Input
        --------------
        ishigh  if not None, use ishigh in case of y > high
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        overflows replaced by new value


        History
        -------
        Written,  MC, Jun 2014
    """
    if ishigh is None:
        out = np.minimum(dat, high)
    else:
        out = np.where(dat > high, ishigh, dat)
    return np.where(dat == undef, undef, out)


def limits(dat, mini, maxi, undef=-9999.):
    """
        Replacement of underflows or overflows by the missing value
        x = limits(y, ll, lh) means
        IF (y < ll) OR (y > lh) THEN x = [missing value] ELSE x = y,
        where ll and lh are variables or numbers.
        This function may be used to check values lying in between certain limits.
        If one of the limits is exceeded the value is set to the missing value defined in the control file.


        Definition
        ----------
        def limits(dat, mini, maxi, undef=-9999.):


        Input
        -----
        dat     ND-array
        mini    ND-array
        maxi    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        underflows or overflows replaced by the missing value


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((dat >= mini) & (dat <= maxi), dat, undef)


def mean(var1, axis=None, undef=-9999.):
    """
        Calculation of mean value
        x = mean(y1, y2, ..., yn) means x = (y1 + y2 + ... + yn)/n,
        where y1, y2, ..., yn are variables or numbers.


        Definition
        ----------
        def mean(var1, axis=None, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        axis    if not None, average over given axis
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        mean value


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.ma.mean(np.ma.array(var1, mask=(var1==undef)), axis=axis).filled(undef)


def mini(var1, axis=None, undef=-9999.):
    """
        Calculation of minimum value
        x = mini(y1,y2,...,yn) means x = min(y1,y2,...,yn),
        where y1, y2, ..., yn are variables or numbers.


        Definition
        ----------
        def mini(var1, axis=None, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        axis    if not None, average over given axis
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        minimum value


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.ma.amin(np.ma.array(var1, mask=(var1==undef)), axis=axis).filled(undef)


def maxi(var1, axis=None, undef=-9999.):
    """
        Calculation of maximum value
        x = maxi(y1,y2,...,yn) means x = max(y1,y2,...,yn),
        where y1, y2, ..., yn are variables or numbers.


        Definition
        ----------
        def maxi(var1, axis=None, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        axis    if not None, average over given axis
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        maximum value


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.ma.amax(np.ma.array(var1, mask=(var1==undef)), axis=axis).filled(undef)


# Not implemented: met_torad


def met_lwrad(dat, tpyr, undef=-9999.): # assumes that dat was already multiplied with calibration factor
    """
        Calculation of long wave radiation from net radiometer
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


        Input
        -----
        dat     ND-array
        tpyr    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        total radiation in W m-2


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((dat==undef) | (tpyr==undef), undef,
                    dat + sigma * (tpyr+273.15)**4)


def met_trad(dat, eps, undef=-9999.):
    """
        Calculation of radiation temperature from long wave radiation
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


        Input
        -----
        dat     ND-array
        eps     ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        radiation temperature in oC


        History
        -------
        Written,  MC, Jun 2014
    """
    const = 1./(eps*sigma)
    trad  = np.ma.sqrt(np.ma.sqrt(const*np.ma.array(dat, mask=(dat==undef)))) - 273.15
    return trad.filled(undef)


def met_alb(swd, swu, swdmin=50., swumin=10., undef=-9999.):
    """
        Calculation of albedo from short wave downward and upward radiation
        x = met_alb(Rsd, Rsu) where
        Rsd is the short wave downward radiation in Wm-2, Rsu is the short wave upward radiation in Wm-2,
        The albedo in % is calculated according to the following formula:
        x = 100 * ( Rsu / Rsd )
        If Rsd < 50 W m-2 or Rsu < 10 W m-2 the result is [missing value].
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_alb(swd, swu, swdmin=50., swumin=10., undef=-9999.):


        Input
        -----
        swd    ND-array
        swu    ND-array


        Optional Input
        --------------
        swdmin  If Rsd < 50 W m-2 the result is undef.
        swumin  If Rsu < 10 W m-2 the result is undef.
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        albedo in %


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((swd==undef) | (swu==undef) | (swd<swdmin) | (swu<swumin),
                    undef, division(swu*100., swd, undef))


def met_albl(swd, swu, swdmin, swumin, undef=-9999.):
    """
        Calculation of albedo from short wave downward and upward radiation with limits
        x=met_albl(Rsd,Rsu,Rsd_limit,Rsu_limit)where
        Rsd is the short wave downward radiation in Wm-2,
        Rsu is the short wave upward radiation in Wm-2,
        Rsd_limit is the short wave downward radiation limit in Wm-2,
        Rsu_limit is the short wave upward radiation limit in Wm-2,
        The albedo in % is calculated according to the following formula:
        x = 100 * ( Rsu / Rsd )
        If Rsd < Rsd_limit or Rsu < Rsu_limit the result is [missing value].
        All four parameters may be variables or numbers.


        Definition
        ----------
        def met_albl(swd, swu, swdmin, swumin, undef=-9999.):


        Input
        -----
        swd      ND-array
        swu      ND-array
        swdmin   short wave downward radiation limit in Wm-2
        swumin   short wave upward radiation limit in Wm-2


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        albedo in %


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((swd==undef) | (swu==undef) | (swd<swdmin) | (swu<swumin),
                    undef, division(swu*100., swd, undef))


def met_vpmax(temp, undef=-9999.):
    """
        Calculation of saturation water vapour pressure
        x = met_vpmax(T) where
        T is the air temperature in oC.
        The saturation water vapour pressure in mbar (hPa) is calculated according to the following formula:
        x = 6.1078 * exp(17.08085 * T / (234.175 + T))
        The parameter may be a variable or a number.


        Definition
        ----------
        def met_vpmax(temp, undef=-9999.):


        Input
        -----
        var1    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        saturation water vapour pressure in mbar (hPa)


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=(temp==undef)))*0.01
    return es.filled(undef)


def met_vpact(temp, rh, undef=-9999.):
    """
        Calculation of actual water vapour pressure
        x = met_vpact(T,rh) where T is the air temperature in oC, rh is the relative humidity in %.
        The actual water vapour pressure in mbar (hPa) is calculated according to the following for- mulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        x = Es * rh/100
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_vpact(temp, rh, undef=-9999.):


        Input
        -----
        temp    ND-array
        rh      ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        actual water vapour pressure in mbar (hPa)


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef))))*0.01
    ea = es*rh*0.01
    return ea.filled(undef)


def met_vpdef(temp, rh, undef=-9999.):
    """
        Calculation of water vapour pressure deficit
        x = met_vpdef(T, rh) where T is the air temperature in oC, rh is the relative humidity in %.
        The water vapour pressure deficit in mbar (hPa) is calculated according to the following for- mulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        E = Es * rh/100
        x = Es - E
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_vpdef(temp, rh, undef=-9999.):


        Input
        -----
        temp    ND-array
        rh      ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        water vapour pressure deficit in mbar (hPa)


        History
        -------
        Written,  MC, Jun 2014
    """
    es  = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef))))*0.01
    ea  = es*rh*0.01
    vpd = es - ea
    return vpd.filled(undef)


def met_sh(temp, rh, p, undef=-9999.):
    """
        Calculation of specific humidity
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


        Input
        -----
        temp    ND-array
        rh      ND-array
        p       ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        specific humidity in g kg-1


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea = es*rh*0.01
    sh =  division(622.*ea, (p-0.378*ea), undef)
    return sh.filled(undef)


def met_tpot(temp, p, undef=-9999.):
    """
        Calculation of potential temperature
        x = met_tpot(T, p) where
        T is the air temperature in oC,
        p is the air pressure in mbar (hPa).
        The potential temperature in K is calculated according to the following formula:
        x = (T + 273.16) * (1000/p)**0.286
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_tpot(temp, p, undef=-9999.):


        Input
        -----
        temp    ND-array
        p       ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        potential temperature in K


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((temp==undef) | (p==undef), undef, (temp+273.15)*division(1000.,p)**0.286)


def met_rho(temp, rh, p, undef=-9999.):
    """
        Calculation of air density
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


        Input
        -----
        temp    ND-array
        rh      ND-array
        p       ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        air density in kg m-3


        History
        -------
        Written,  MC, Jun 2014
    """
    es  = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea  = es*rh*0.01
    sh  = division(622.*ea, (p-0.378*ea), undef)
    Tv  = ((temp+273.15)*(1+0.000608*sh)) - 273.15
    rho = division(p*100., (287.05*(Tv+273.15)), undef)
    return rho.filled(undef)


def met_dpt(temp, rh, undef=-9999.):
    """
        Calculation of dew point temperature
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


        Input
        -----
        temp    ND-array
        rh      ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        dew point temperature in oC


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef))))*0.01
    ea = es*rh*0.01
    dpt = 234.175 * np.ma.log(ea/6.1078) / (17.08085 - np.ma.log(ea/6.1078))
    return dpt.filled(undef)


def met_h2oc(temp, rh, p, undef=-9999.):
    """
        Calculation of water vapour concentration
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


        Input
        -----
        temp    ND-array
        rh      ND-array
        p       ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        water vapour concentration in mmol mol-1


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea = es*rh*0.01
    c  = division(1000.*ea, p, undef)
    return c.filled(undef)


# Not implemented: met_psy_rh


# Not implemented: met_dpt_rh


def met_h2oc_rh(temp, h, p, undef=-9999.):
    """
        Calculation of relative humidity from water vapour concentration
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


        Input
        -----
        temp    ND-array
        rh      ND-array
        p       ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        relative humidity in %


        History
        -------
        Written,  MC, Jun 2014
    """
    es = esat(np.ma.array(temp+273.15, mask=((temp==undef)|(h==undef)|(p==undef))))*0.01
    ea = 0.001 * h * p
    c  = 100.*ea/es
    return c.filled(undef)


def met_wdrot(wd, a, undef=-9999.):
    """
        Rotation of wind direction
        x = met_wdrot(wd, a) where
        wd is the wind direction in degree,
        a is the rotation angle in degree (positive is clockwise).
        The rotated wind direction is calculated according to the following formulas:
        x = wd + a
        IF x < 0 THEN x = x + 360
        IF x <= 360 THEN x = x - 360
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_wdrot(wd, a, undef=-9999.):


        Input
        -----
        wd    ND-array
        a    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        rotated wind direction


        History
        -------
        Written,  MC, Jun 2014
    """
    rot = np.ma.array(wd+a, mask=((wd==undef)|(a==undef)))
    rot = np.ma.where(rot < 0., rot+360., rot)
    rot = np.ma.where(rot >= 360., rot-360., rot)
    return rot.filled(undef)


def met_urot(u, v, a, undef=-9999.):
    """
        Rotation of u-component of wind vector
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


        Input
        -----
        u    ND-array
        v    ND-array
        a    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        rotated u-component


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((u==undef) | (v==undef) | (a==undef), undef, u*np.cos(np.deg2rad(a)) + v*np.sin(np.deg2rad(a)))


def met_vrot(u, v, a, undef=-9999.):
    """
        Rotation of v-component of wind vector
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


        Input
        -----
        u    ND-array
        v    ND-array
        a    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        rotated v-component


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((u==undef) | (v==undef) | (a==undef), undef, -u*np.sin(np.deg2rad(a)) + v*np.cos(np.deg2rad(a)))


def met_uv_wv(u, v, undef=-9999.):
    """
        Calculation of wind velocity from u- and v-component of wind vector
        x = met_uv_wv(u, v) where
        u is the u-component of the wind vector, v is the v-component of the wind vector.
        The horizontal wind velocity is calculated according to the following formula:
        x = sqrt(u**2 + v**2)
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_uv_wv(u, v, undef=-9999.):


        Input
        -----
        u    ND-array
        v    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        horizontal wind velocity


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((u==undef) | (v==undef), undef, np.sqrt(u*u + v*v))


def met_uv_wd(u, v, undef=-9999.):
    """
        Calculation of wind direction from u- and v-component of wind vector
        x = met_uv_wd(u, v) where
        u is the u-component of the wind vector, v is the v-component of the wind vector.
        The horizontal wind velocity is calculated according to the following formulas:
        IF u = 0 AND v = 0 THEN x = 0
        IF u = 0 AND v < 0 THEN x = 360
        IF u = 0 AND v > 0 THEN x = 180
        IF u > 0 THEN x = 270 - arctan(v/u)
        IF u < 0 THEN x = 90 - arctan(v/u)
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_uv_wd(u, v, undef=-9999.):


        Input
        -----
        u    ND-array
        v    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        horizontal wind velocity


        History
        -------
        Written,  MC, Jun 2014
    """
    wd = np.ma.zeros(u.shape)
    wd.mask = (u==undef) | (v==undef)
    wd = np.ma.where((u==0.) & (v==0.), 0., wd)
    wd = np.ma.where((u==0.) & (v<0.), 360., wd)
    wd = np.ma.where((u==0.) & (v>0.), 180., wd)
    wd = np.ma.where((u>0.), 270.-np.rad2deg(np.ma.arctan(v/u)), wd)
    wd = np.ma.where((u<0.), 90.-np.rad2deg(np.ma.arctan(v/u)), wd)
    return wd.filled(undef)


def met_wvwd_u(wv, wd, undef=-9999.):
    """
        Calculation of u-component of wind vector from wind velocity and wind direction
        x = met_wvwd_u(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
        The u-component of the wind vector is calculated according to the following formula:
        x = -wv * sin (wd)
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_wvwd_u(wv, wd, undef=-9999.):


        Input
        -----
        wv    ND-array
        wd    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        u-component of the wind vector


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((wv==undef) | (wd==undef), undef, -wv*np.sin(np.deg2rad(wd)))


def met_wvwd_v(wv, wd, undef=-9999.):
    """
        Calculation of v-component of wind vector from wind velocity and wind direction
        x = met_wvwd_v(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
        The v-component of the wind vector is calculated according to the following formula:
        x = -wv * cos (wd)
        Both parameters may be variables or numbers.


        Definition
        ----------
        def met_wvwd_v(wv, wd, undef=-9999.):


        Input
        -----
        wv    ND-array
        wd    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        v-component of the wind vector


        History
        -------
        Written,  MC, Jun 2014
    """
    return np.where((wv==undef) | (wd==undef), undef, -wv*np.cos(np.deg2rad(wd)))


def ifeq(var1, iif, ithen, ielse):
    """
        If-statements
        x = ifeq(y,a0,a1,a2) means IF y == a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def ifeq(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y == a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 == iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


def ifne(var1, iif, ithen, ielse):
    """
        If-statements
        x = ifne(y,a0,a1,a2) means IF y != a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def ifne(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y != a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 != iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


def ifle(var1, iif, ithen, ielse):
    """
        If-statements
        x = ifle(y,a0,a1,a2) means IF y >= a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def ifle(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y >= a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 <= iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


def ifge(var1, iif, ithen, ielse):
    """
        If-statements
        x = ifge(y,a0,a1,a2) means IF y <= a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def ifge(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y <= a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 <= iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


def iflt(var1, iif, ithen, ielse):
    """
        If-statements
        x = iflt(y,a0,a1,a2) means IF y < a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def iflt(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y < a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 < iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


def ifgt(var1, iif, ithen, ielse):
    """
        If-statements
        x = ifgt(y,a0,a1,a2) means IF y > a0 THEN x = a1 ELSE x = a2
        All parameters may be variables or numbers.


        Definition
        ----------
        def ifgt(var1, iif, ithen, ielse):


        Input
        -----
        var1     ND-array
        iif      ND-array
        ithen    ND-array
        ielse    ND-array


        Optional Input
        --------------
        undef   elements are excluded from the calculations if any of the inputs equals undef (default: -9999.)


        Output
        ------
        IF y > a0 THEN x = a1 ELSE x = a2


        History
        -------
        Written,  MC, Jun 2014
    """
    out = np.where(var1 > iif, ithen, ielse)
    return np.where((var1==undef) | (iif==undef) | (ithen==undef) | (ielse==undef), undef, out)


# Not implemented: write


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
