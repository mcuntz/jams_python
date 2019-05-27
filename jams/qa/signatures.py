#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
"""
    Defines signatures of discharge time series.

    They are:
        autocorrelation
        flow duration curves
        rising and declining limb densities
        maximum monthly flow
        moments
        peak distribution

    Translation of the module mo_mrm_signatures.f90 by mainly Juliane Mai, Jun 2015
    RunoffRatio and ZeroFlowRatio were not ported.


    License
    -------
    This file is part of the JAMS Python package.

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

    Copyright 2016 Matthias Cuntz


    History
    -------
    Written, MC, May 2016
"""

all = ['autocorrelation', 'flowdurationcurve', 'limbdensities', 'maximummonthlyflow', 'moments', 'peakdistribution']

# --------------------------------------------------------------------

def autocorrelation(dat, lags):
    """
        NAME
        autocorrelation

        PURPOSE
        Autocorrelation of a given data series.

        Calculates the autocorrelation of a data series at given lags.
        The function is basically a wrapper of the function jams.correlate.

        DEFINITION
        def autocorrelation(dat, lags):

        INTENT(IN)
        data   1D array_like of data
        lags   scalar or 1D array_like of lags where autocorrelation is requested
               positive and negative lags are possible

        RETURN
        autocorrelation of data at given lags

        RESTRICTIONS
        Works only with 1D arrays

        EXAMPLE
            None

        LITERATURE
            Used as hydrologic signature with lag 1 in
               Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
               A framework to assess the realism of model structures using hydrological signatures.
               Hydrology and Earth System Sciences, 17(5), 1893-1912. doi:10.5194/hess-17-1893-2013

        HISTORY
        Written, Juliane Mai, Jun 2015 in Fortran
        Modified, Matthias Cuntz, May 2016 - Python port
    """
    from jams.correlate import correlate
    autocorrelation = correlate(dat, dat)
    return autocorrelation[lags]

# --------------------------------------------------------------------

def flowdurationcurve(dat, quantiles=None, concavity_index=False,
                      mid_segment_slope=False, mhigh_segment_volume=False,
                      high_segment_volume=False, low_segment_volume=False):
    """
        NAME
            flowdurationcurve

        PURPOSE
        Flow duration curves.

        Calculates the flow duration curves for a given data vector. The Flow duration curve at a
        certain quantile x is the data point p where x% of the data points are above the value p.
        Hence the function percentile of the module mo_percentile is used. But percentile is
        determining the point p where x% of the data points are below that value. Therfore, the
        given quantiles are transformed by (1.0-quantile) to get the percentiles of exceedance probabilities.

        Optionally, the concavity index CI can be calculated [Zhang2014]. CI is defined by
                \f[ CI = \frac{q_{10\%}-q_{99\%}}{q_{1\%}-q_{99\%}} \f]
        where \f$ q_{x} \f$ is the data point where x% of the data points are above that value.
        Hence, exceedance probabilities are used.\n

        Optionally, the FDC mid-segment slope \f$FDC_{MSS}\f$ as used by Shafii et. al (2014) can be returned.
        The \f$FDC_{MSS}\f$ is defined as
                \f[ FDC_{MSS} = \log(q_{m_1})-\log(q_{m_2}) \f]
        where \f$ m_1 \f$ and \f$ m_2 \f$ are the lowest and highest flow exceedance probabilities within the
        midsegment of FDC. The settings \f$m_1=0.2\f$ and \f$0.7\f$ are used by Shafii et. al (2014) and are
        implemented like that.\n

        Optionally, the FDC medium high-segment volume \f$FDC_{MHSV}\f$ as used by Shafii et. al (2014) can be
        returned. The \f$FDC_{MHSV}\f$ is defined as
                \f[ FDC_{MHSV} = \sum_{h=1}^{H} q_h \f]
        where \f$ h=1,2,...,H \f$ are flow indeces located within the high-flow segment (exceedance probabilities
        lower than \f$m_1\f$). \f$H\f$ is the index of the maximum flow. The settings \f$m_1=0.2\f$ is used here
        to be consistent with the definitions of the low-segment (0.7-1.0) and the mid-segment (0.2-0.7).\n

        Optionally, the FDC high-segment volume \f$FDC_{HSV}\f$ as used by Shafii et. al (2014) can be returned.
        The \f$FDC_{HSV}\f$ is defined as
                \f[ FDC_{HSV} = \sum_{h=1}^{H} q_h \f]
        where \f$ h=1,2,...,H \f$ are flow indeces located within the high-flow segment (exceedance probabilities
        lower than \f$m_1\f$). \f$H\f$ is the index of the maximum flow. The settings \f$m_1=0.02\f$ is used by
        Shafii et. al (2014) and is implemented like that.\n

        Optionally, the FDC low-segment volume \f$FDC_{LSV}\f$ as used by Shafii et. al (2014) can be returned.
        The \f$FDC_{LSV}\f$ is defined as
                \f[ FDC_{LSV} = -\sum_{l=1}^{L} (\log(q_l) - \log(q_L)) \f]
        where \f$ l=1,2,...,L \f$ are flow indeces located within the low-flow segment (exceedance probabilities
        larger than \f$m_1\f$). \f$L\f$ is the index of the minimum flow. The settings \f$m_1=0.7\f$ is used by
        Shafii et. al (2014) and is implemented like that.\n

        CALLING SEQUENCE
        def flowdurationcurve(dat, quantiles=None, concavity_index=False,
                              mid_segment_slope=False, mhigh_segment_volume=False,
                              high_segment_volume=False, low_segment_volume=False):

        INTENT(IN)
        data              1D array_like data series

        INTENT(IN), OPTIONAL
        quantiles              Scalar or 1D array_like percentages of exceedance
        concavity_index        True: concavity index as defined by Sauquet et al. (2011)
        mid_segment_slope      True: mid-segment slope as defined by Shafii et al. (2014)
        mhigh_segment_volume   True: medium high-segment volume
        high_segment_volume    True: high-segment volume as defined by Shafii et al. (2014)
        low_segment_volume     True: low-segment volume as defined by Shafii et al. (2014)

        RETURN
        Flow duration curve at resp. quantiles or optional signatures of the flow duration curve

        RESTRICTIONS
        Thresholds in mid_segment_slope, mhigh_segment_volume, high_segment_volume, low_segment_volume are hard coded.

        EXAMPLE
            None

        LITERATURE
            FDC is used as hydrologic signature (quantiles not specified) in
               Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
               A framework to assess the realism of model structures using hydrological signatures.
               Hydrology and Earth System Sciences, 17(5), 1893-1912. doi:10.5194/hess-17-1893-2013
            Concavity Index used as hydrologic signature in
               Zhang, Y., Vaze, J., Chiew, F. H. S., Teng, J., & Li, M. (2014).
               Predicting hydrological signatures in ungauged catchments using spatial interpolation, index model, and
               rainfall-runoff modelling.
               Journal of Hydrology, 517(C), 936-948. doi:10.1016/j.jhydrol.2014.06.032
            Concavity index is defined using exceedance probabilities by
               Sauquet, E., & Catalogne, C. (2011).
               Comparison of catchment grouping methods for flow duration curve estimation at ungauged sites in France.
               Hydrology and Earth System Sciences, 15(8), 2421-2435. doi:10.5194/hess-15-2421-2011
            mid_segment_slope, high_segment_volume, low_segment_volume used as hydrologic signature in
               Shafii, M., & Tolson, B. A. (2015).
               Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
               Water Resources Research, 51(5), 3796-3814. doi:10.1002/2014WR016520

        HISTORY
        Written Remko Nijzink and Juliane Mai, March 2014 - June 2015 in Fortran
        Modified, Matthias Cuntz, May 2016 - Python port
    """
    if isinstance(dat, np.ma.MaskedArray):
        mdat = dat.data[~dat.mask]
    else:
        mdat = dat

    allind = concavity_index+mid_segment_slope+mhigh_segment_volume+high_segment_volume+low_segment_volume
    if quantiles is not None:
        if allind > 0:
            raise IOError('Either quantiles or indexes can be returned.')
    else:
        if allind > 1:
            raise IOError('Only one index can be returned.')

    if quantiles is not None:
        fdc = np.percentile(mdat, (1.-np.array(quantiles))*100.)
        return fdc
    elif concavity_index:
        p10 = np.percentile(mdat, (1.-0.10)*100.)
        p99 = np.percentile(mdat, (1.-0.99)*100.)
        p01 = np.percentile(mdat, (1.-0.01)*100.)
        concavity_index = (p10-p99) / (p01-p99)
        return concavity_index
    elif mid_segment_slope:
        p20 = np.percentile(mdat, (1.-0.20)*100.)
        p70 = np.percentile(mdat, (1.-0.70)*100.)
        mid_segment_slope = np.log(p20) - np.log(p70)
        return mid_segment_slope
    elif mhigh_segment_volume:
        # medium high-flows are defined to be between 0.0 and 0.2 as to be constistent
        # with the mid-segment (0.2-0.7) and low-segment (0.7-1.0) definitions
        p20 = np.percentile(mdat, (1.-0.20)*100.)
        ii = np.where(mdat>=p20)[0]
        mhigh_segment_volume = np.sum(mdat[ii])
        return mhigh_segment_volume
    elif high_segment_volume:
        # high-flows are defined to be between 0.0 and 0.02 by Shafii et. al (2014)
        p02 = np.percentile(mdat, (1.-0.02)*100.)
        ii = np.where(mdat>=p02)[0]
        high_segment_volume = np.sum(mdat[ii])
        return high_segment_volume
    elif low_segment_volume:
        # low-flows are defined to be between 0.7 and 1.0 by Shafii et. al (2014)
        min_flow_value = np.amin(mdat)
        p70 = np.percentile(mdat, (1.-0.7)*100.)
        ii = np.where(mdat<=p70)[0]
        low_segment_volume = -1.0 * np.sum(np.log(mdat[ii]) - np.log(min_flow_value))
        return low_segment_volume
    else:
        raise ValueError('Should not be here.')

# --------------------------------------------------------------------

def limbdensities(dat):
    """
        NAME
        limbdensities

        PURPOSE
        Calculates limb densities

        Calculates rising and declinging limb densities. The peaks of the given series are
        first determined by looking for points where preceding and subsequent datapoint are lower.
        Second, the number of datapoints with rising values (nrise) and declining values (ndecline)
        are counted basically by comparing neighbors.

        The duration the data increase (nrise) divided by the number of peaks (npeaks)
        gives the rising limb density RLD
            \f[ RLD=t_{rise}/n_{peak} \f]
        whereas the duration the data decrease (ndecline) divided by the number of peaks (npeaks)
        gives the declining limb density DLD
            \f[ DLD=t_{fall}/n_{peak}. \f]

        CALLING SEQUENCE
        def limbdensities(dat, rld=False, dld=False):

        INTENT(IN)
        dat    1D array_like data series

        RETURN
        rising limb density, declining limb density

        RESTRICTIONS
            None

        EXAMPLE
            None

        LITERATURE
            Rising limb density used as hydrologic signature in
               Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
               A framework to assess the realism of model structures using hydrological signatures.
               Hydrology and Earth System Sciences, 17(5), 1893-1912. doi:10.5194/hess-17-1893-2013

        HISTORY
            Written Remko Nijzink and Juliane Mai, March 2014 - Jun 2015 in Fortran
            Modified Matthias Cuntz, May 2016 - Python port
    """
    if isinstance(dat, np.ma.MaskedArray):
        ndat = dat.count()
    else:
        ndat = dat.size

    thres_rise = 1.0
    dd = np.ma.diff(dat)
    nrise = np.sum((dd-thres_rise) > 0.)
    ndecline = ndat - nrise
    # dd[0:-1]*dd[1:]<0 gives min and max: multiply with dd[0:-1] -> <0 are maxima
    dd1 = dd[0:-1]*dd[1:]
    dextreme = np.where(dd1 < 0., dd1, 0.)
    npeak = np.sum(dextreme*dd[0:-1] < 0.)

    if npeak>0:
        return float(nrise)/float(npeak), float(ndecline)/float(npeak)
    else:
        return 0., 0.

# --------------------------------------------------------------------

def maximummonthlyflow(date, dat):
    """
    NAME
    maximummonthlyflow

    PURPOSE
    Maximum of average flows per months.

    Maximum of average flow per month is defined as
       \f[ max_{monthly flow} = Max( F(i), i=1,..12 ) \f]
    where \$f F(i) $\f is the average flow of month i.

    CALLING SEQUENCE
    def maximummonthlyflow(date, dat):

    INTENT(IN)
    date          1d array_like julian date
    data          1d array_like of data

    RETURN
    Maximum of average flow per month

    RESTRICTIONS
    None

    EXAMPLE
        None

    LITERATURE
        used as hydrologic signature in
           Shafii, M., & Tolson, B. A. (2015).
           Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
           Water Resources Research, 51(5), 3796-3814. doi:10.1002/2014WR016520

    HISTORY
    Written Juliane Mai, Jun 2015 in Fortran
    Modified Matthias Cuntz, May 2016 - Python port
    """
    from jams.means import means
    mdate, mdat = means(date, dat, month=True)
    return mdat.max()

# --------------------------------------------------------------------

def moments(dat, mean_data=False, stddev_data=False, median_data=False,
            max_data=False, mean_log=False, stddev_log=False, median_log=False, max_log=False):
    """
    NAME
    moments

    PURPOSE
    Moments of data and log-transformed data, e.g. mean and standard deviation.

    Returns several moments of data series given, i.e.
        * mean               of data
        * standard deviation of data
        * median             of data
        * maximum/ peak      of data
        * mean               of log-transformed data
        * standard deviation of log-transformed data
        * median             of log-transformed data
        * maximum/ peak      of log-transformed data

    CALLING SEQUENCE
    def moments(dat, mean_data=False, stddev_data=False, median_data=False,
                max_data=False, mean_log=False, stddev_log=False, median_log=False, max_log=False):

    INTENT(IN)
    dat        1D array_like of data

    INTENT(OUT), OPTIONAL
    mean_data     mean               of data
    stddev_data   standard deviation of data
    median_data   median             of data
    max_data      maximum/ peak      of data
    mean_log      mean               of log-transformed data
    stddev_log    standard deviation of log-transformed data
    median_log    median             of log-transformed data
    max_log       maximum/ peak      of log-transformed data

    RETURN
    moment(s) of data depending on optional arguments
    order is given by the above order of the parameters

    RESTRICTIONS
    None

    EXAMPLE
    None

    LITERATURE
        mean_log and stddev_log used as hydrologic signature in
           Zhang, Y., Vaze, J., Chiew, F. H. S., Teng, J., & Li, M. (2014).
           Predicting hydrological signatures in ungauged catchments using spatial interpolation, index model, and
           rainfall-runoff modelling.
           Journal of Hydrology, 517(C), 936-948. doi:10.1016/j.jhydrol.2014.06.032
        mean_data, stddev_data, median_data, max_data, mean_log, and stddev_log used as hydrologic signature in
           Shafii, M., & Tolson, B. A. (2015).
           Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
           Water Resources Research, 51(5), 3796-3814. doi:10.1002/2014WR016520

    HISTORY
    Written Juliane Mai, Jun 2015 in Fortran
    Modified Matthias Cuntz, May 2016 - Python port
    """
    allind = mean_data+stddev_data+median_data+max_data+mean_log+stddev_log+median_log+max_log
    if allind == 0:
        raise IOError('At least one moment must be returned.')

    if isinstance(dat, np.ma.MaskedArray):
        mdat = dat.data[~dat.mask]
    else:
        mdat = dat

    out = []
    if mean_data:   out.append(np.mean(mdat))
    if stddev_data: out.append(np.std(mdat, ddof=1))
    if median_data: out.append(np.median(mdat))
    if max_data:    out.append(np.amax(mdat))

    allind = mean_log+stddev_log+median_log+max_log
    if allind > 0:
        undef = -9999
        ii = np.where(mdat > 0.)[0]
        if ii.size > 0:
            logdat = np.log(mdat[ii])
            if mean_log:   out.append(np.mean(logdat))
            if stddev_log: out.append(np.std(logdat, ddof=1))
            if median_log: out.append(np.median(logdat))
            if max_log:    out.append(np.amax(logdat))
        else:
            if mean_log:   out.append(undef)
            if stddev_log: out.append(undef)
            if median_log: out.append(undef)
            if max_log:    out.append(undef)

    return out

# --------------------------------------------------------------------

def peakdistribution(dat, quantiles=None, slope=False):
    """
    NAME
    peakdistribution

    PURPOSE
    Calculates the peak distribution.

    First, the peaks of the time series given are identified. For the peak distribution
    only this subset of data points are considered. Second, the peak distribution at the
    quantiles given is calculated. Calculates the peak distribution at the quantiles given
    using mo_percentile. Since the exceedance probabilities are usually used in
    hydrology the function percentile is used with (1.0-quantiles). \n

    Optionally, the slope of the peak distribution between 10th and 50th percentile, i.e.
       \f[ slope = \frac{\mathrm{peak\_{data}}_{0.1}-\mathrm{peak\_{data}}_{0.5}}{0.9-0.5} \f]
    can be returned.\n

    CALLING SEQUENCE
    def peakdistribution(dat, quantiles=None, slope=False):

    INTENT(IN)
    dat                       1D array_like data array

    INTENT(IN), OPTIONAL
    quantiles                 requested quantiles for distribution
    slope                     True: slope of the peak distribution between 10th and 50th percentile

    RETURN
    Distribution of peak values at resp. quantiles   or   slope of the peak distribution

    RESTRICTIONS
    None

    EXAMPLE
    None

    LITERATURE
        slope used as hydrologic signature in
           Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
           A framework to assess the realism of model structures using hydrological signatures.
           Hydrology and Earth System Sciences, 17(5), 1893-1912. doi:10.5194/hess-17-1893-2013

    HISTORY
    Written Remko Nijzink and Juliane Mai, March 2014 - Jun 2015 in Fortran
    Modified Matthias Cuntz, May 2016 - Python port
    """
    if quantiles is not None:
        if slope:
            raise IOError('Either quantiles or slope of peak distribution can be returned.')
    else:
        if not slope:
            raise IOError('Either quantiles or slope of peak distribution must be returned.')

    # find peaks
    data_peak = []
    for i in range(1,len(dat)-1):
        if ((dat[i-1] <= dat[i]) and (dat[i+1] <= dat[i])): data_peak.append(dat[i])
    data_peak = np.array(data_peak)

    if quantiles is not None:
        return np.percentile(data_peak, (1.-np.array(quantiles))*100.)
    else:
        # calculate slope between 10% and 50% quantiles, per definition
        quantiles = [0.1, 0.5]
        p10, p50 = np.percentile(data_peak, (1.-np.array(quantiles))*100.)
        return (p10-p50)/(0.9-0.5)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
