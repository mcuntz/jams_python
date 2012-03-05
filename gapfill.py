#!/usr/bin/env python
import numpy as np

def gapfill(datein, datain, rgin, tairin, vpdin,
            data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
            rg_dev=50., tair_dev=2.5, vpd_dev=5,
            longestmarginalgap=60, undef=-9999., ddof=1,
            err=False, shape=False):
    """
        Fills gaps of flux data from Eddy covariance measurements according to
        Reichstein et al. (2005).
        If there is a gap in the data, look for similar meteorological conditions
        (defined as maximum possible deviations) in a certain time window and fill
        with the average of these 'similar' values.

        The routine can also do the same search for similar meteorological conditions
        for every data point and calculate its standard deviation as a measure of uncertainty.


        Definition
        ----------
        def gapfill(datein, datain, rgin, tairin, vpdin,
                    data_flag=None, rg_flag=None, tair_flag=None, vpd_flag=None,
                    rg_dev=50., tair_dev=2.5, vpd_dev=5,
                    longestmarginalgap=60, undef=-9999.,
                    err=False):


        Input
        -----
        datein        julian days
        datain        fluxes to fill
        rgin          global radiation [W m-2]
        tairin        air Temperature [deg C]
        vpdin         vapour pressure deficit [hPa]


        Optional Input
        -------------
        data_flag     flags of fluxes: 0=good quality; >0=bad data (default: 0)
        rg_flag       flags of global radiation: 0=good quality; >0=bad data (default: 0)
        tair_flag     flags of air temperature: 0=good quality; >0=bad data (default: 0)
        vpd_flag      flags of vapour pressure deficit: 0=good quality; >0=bad data (default: 0)


        Parameters
        ----------
        rg_dev               threshold for maximum deviation of global radiation (default: 50)
        tair_dev             threshold for maximum deviation of air Temperature (default: 2.5)
        vpd_dev              threshold for maximum deviation of vpd (default: 5)
        longestmarginalgap   avoid extraploation into a gap longer than longestmarginalgap days (default: 60)
        undef                undefined values in data  (default: -9999.)
        ddof                 Degrees of freedom tu use in calculation of standard deviation
                             for error estimate (default: 1)
        err                  if True, fill every data point, i.e. used for error generation (default: False)
        shape                if False then outputs are 1D arrays;
                             if True, output have the same shape as datain
                             if a shape tuple is given, then this tuple is used to reshape


        Ouput
        -----
        if err=False:
            filled_data, quality_class
        else:
            err_data


        Restrictions
        ------------
        If err=True, there is no error estimate if there are no meteorological
          conditions in the vicinity of the data point.


        Literature
        ----------
        Reichstein et al. (2005) On the separation of net ecosystem exchange into
          assimilation and ecosystem respiration: review and improved algorithm.
          Global Change Biology,11,9, p. 1424-1439.


        Examples
        --------
        >>> import numpy as np
        >>> dates = np.array([ 2454258.76042,  2454258.78125,  2454258.80208,  2454258.82292, \
        2454258.84375, \
        2454258.86458,  2454258.88542,  2454258.90625,  2454258.92708,  2454258.94792, \
        2454258.96875,  2454258.98958,  2454259.01042,  2454259.03125,  2454259.05208, \
        2454259.07292,  2454259.09375,  2454259.11458,  2454259.13542,  2454259.15625, \
        2454259.17708,  2454259.19792,  2454259.21875,  2454259.23958,  2454259.26042, \
        2454259.28125,  2454259.30208,  2454259.32292,  2454259.34375,  2454259.36458, \
        2454259.38542,  2454259.40625,  2454259.42708,  2454259.44792,  2454259.46875, \
        2454259.48958,  2454259.51042,  2454259.53125,  2454259.55208,  2454259.57292, \
        2454259.59375,  2454259.61458,  2454259.63542,  2454259.65625,  2454259.67708, \
        2454259.69792,  2454259.71875,  2454259.73958,  2454259.76042,  2454259.78125, \
        2454259.80208,  2454259.82292,  2454259.84375,  2454259.86458,  2454259.88542, \
        2454259.90625,  2454259.92708,  2454259.94792,  2454259.96875,  2454259.98958, \
        2454260.01042,  2454260.03125,  2454260.05208,  2454260.07292,  2454260.09375, \
        2454260.11458,  2454260.13542,  2454260.15625,  2454260.17708,  2454260.19792, \
        2454260.21875,  2454260.23958,  2454260.26042,  2454260.28125,  2454260.30208, \
        2454260.32292,  2454260.34375,  2454260.36458,  2454260.38542,  2454260.40625, \
        2454260.42708,  2454260.44792,  2454260.46875,  2454260.48958,  2454260.51042, \
        2454260.53125,  2454260.55208,  2454260.57292,  2454260.59375,  2454260.61458, \
        2454260.63542,  2454260.65625,  2454260.67708,  2454260.69792,  2454260.71875, \
        2454260.73958,  2454260.76042,  2454260.78125,  2454260.80208,  2454260.82292, \
        2454260.84375,  2454260.86458,  2454260.88542,  2454260.90625,  2454260.92708, \
        2454260.94792,  2454260.96875,  2454260.98958,  2454261.01042,  2454261.03125, \
        2454261.05208,  2454261.07292,  2454261.09375,  2454261.11458,  2454261.13542, \
        2454261.15625,  2454261.17708,  2454261.19792,  2454261.21875,  2454261.23958, \
        2454261.26042,  2454261.28125,  2454261.30208,  2454261.32292,  2454261.34375, \
        2454261.36458,  2454261.38542,  2454261.40625,  2454261.42708,  2454261.44792, \
        2454261.46875,  2454261.48958,  2454261.51042,  2454261.53125,  2454261.55208, \
        2454261.57292,  2454261.59375,  2454261.61458,  2454261.63542,  2454261.65625, \
        2454261.67708,  2454261.69792,  2454261.71875,  2454261.73958,  2454261.76042, \
        2454261.78125,  2454261.80208,  2454261.82292,  2454261.84375,  2454261.86458])
        >>> data = np.array([  7.64160000e-01,   2.49872000e+00,   3.13540000e-01,   1.57725000e+00, \
         1.09717000e+00, \
         1.48275000e+00,   2.96011000e+02,  -8.91389000e+00,  -1.75670000e+01,  -1.90672000e+01, \
        -1.61016000e+01,  -1.50144000e+01,   1.19916000e+01,  -3.49755000e+01,  -4.86368000e+01, \
        -2.03383000e+01,  -1.36902000e+01,  -1.59123000e+01,  -9.46660000e+00,  -1.70566000e+01, \
        -1.99620000e+00,  -4.01141000e+00,  -6.49928000e+00,  -3.67885000e+00,  -1.27924000e+00, \
        -6.94000000e-03,   4.73680000e-01,   3.43840000e-01,   2.22053000e+00,   2.53193000e+00, \
         2.96486000e+00,   3.55745000e+00,   4.49445000e+00,   8.03341000e+00,   8.64464000e+00, \
         5.99874000e+00,   3.36288000e+00,   4.29297000e+00,   5.97186000e+00,   4.81254000e+00, \
         3.44986000e+00,   2.54239000e+00,   4.15090000e+00,   3.80287000e+00,   8.28172000e+00, \
         5.34056000e+00,   3.35551000e+00,   1.71703000e+00,  -2.26169000e+00,  -1.99390000e-01, \
        -5.21457000e+00,  -5.31815000e+00,  -8.80319000e+00,  -1.93282000e+01,  -1.60800000e+01, \
        -1.18457000e+01,  -1.33145000e+01,  -1.54858000e+01,  -9.25691000e+00,  -1.36931000e+01, \
        -7.62091000e+00,  -1.32067000e+01,  -1.10190000e+01,  -8.03985000e+00,  -9.39667000e+00, \
        -6.81524000e+00,  -9.77187000e+00,  -8.31190000e+00,  -7.15069000e+00,  -6.68666000e+00, \
        -5.38734000e+00,  -3.87582000e+00,  -1.18202000e+00,   5.16840000e-01,   1.48799000e+00, \
         1.95226000e+00,   3.42510000e+00,   3.73472000e+00,   3.78379000e+00,   3.02995000e+00, \
         4.86265000e+00,   4.28708000e+00,   6.66179000e+00,   4.81036000e+00,   4.09262000e+00, \
         4.56051000e+00,   4.68438000e+00,   5.23003000e+00,   8.04232000e+00,   4.66788000e+00, \
         5.43978000e+00,   9.47543000e+00,   9.70409000e+00,   5.95810000e+00,   2.90853000e+00, \
         2.57788000e+00,   4.55150000e-01,  -1.25265000e+00,  -8.71061000e+00,  -6.58599000e+00, \
        -5.95496000e+00,  -1.77837000e+01,  -1.21494000e+01,  -1.97989000e+01,  -1.66717000e+01, \
        -7.98233000e+00,  -3.30503000e+01,  -2.20931000e+01,   5.61759000e+00,  -2.60516000e+01, \
         3.29603000e+01,  -8.65377000e+00,  -9.92088000e+00,  -9.48236000e+00,  -6.41284000e+00, \
        -8.03108000e+00,  -7.10870000e+00,  -2.72479000e+00,  -5.17970000e-01,   1.86730000e-01, \
         1.71264000e+00,   9.51260000e-01,   2.44303000e+00,   1.22984000e+01,  -9.88986000e+00, \
        -5.05748000e+01,  -1.20842000e+01,   2.05030000e+02,   2.15621000e+02,  -2.05614000e+00, \
        -3.50786000e+00,  -2.90432000e+00,  -2.36699000e+00,  -2.77420000e-01,  -2.80950000e-01, \
        -1.44877000e+01,   1.09837000e+01,   9.89526000e+00,  -1.59730000e-01,  -2.57200000e-02, \
         2.16320000e-01,  -1.19154000e+00,   7.58938000e+00,   4.13524000e+00,  -1.01411000e+01, \
         1.60020000e+01,  -1.31191000e+01,  -7.97806000e+00,  -3.40220000e+00,  -8.60500000e-02])
        >>> rg = np.array([  73.92,   95.72,  111.88,  144.01,  187.81,  203.01,  359.17,  387.19, \
        403.73,  448.68, \
        547.46,  829.08,  766.92,  870.2 ,  779.35,  681.81,  541.53,  702.85,  618.89,  595.75, \
        325.41,  187.52,  310.31,  203.59,  117.52,   68.28,   25.29,    0.  ,    0.  ,    0.  , \
          0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  , \
          0.  ,    0.  ,    0.  ,    0.  ,    0.  ,   22.97,   68.95,  127.95,  201.1 ,  290.61, \
        385.85,  483.1 ,  579.68,  667.09,  740.14,  806.13,  861.59,  908.45,  931.4 ,  953.39, \
        950.52,  933.31,  909.4 ,  864.46,  815.69,  739.19,  614.3 ,  571.94,  510.64,  412.53, \
        331.34,  235.34,  136.17,   60.47,   28.69,    0.  ,    0.  ,    0.  ,    0.  ,    0.  , \
          0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  , \
          0.  ,    0.  ,    4.93,   24.38,   35.16,   78.13,  198.9 ,  296.82,  351.62,  484.06, \
        570.7 ,  609.04,  733.45,  773.61,  759.27,  890.28,  896.97,  917.05,  930.44,  918.01, \
        885.5 ,  833.86,  781.26,  715.28,  639.35,  568.21,  497.25,  246.43,  149.65,  104.04, \
         54.49,   22.57,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  , \
          0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  , \
          0.  ,   12.34,   31.65,   84.72,  178.15,  286.59,  273.2 ,  277.32,  136.65,  117.52])
        >>> tair = np.array([ 15.71,  15.86,  15.84,  15.83,  16.22,  16.56,  16.78,  17.35,  18.8 , \
        19.27,  19.53, \
        20.37,  21.73,  22.87,  23.53,  23.94,  24.45,  24.81,  25.33,  25.66,  25.17,  24.39, \
        24.47,  24.17,  23.71,  23.29,  22.93,  22.67,  22.27,  22.02,  21.58,  21.22,  20.54, \
        19.88,  19.35,  18.84,  18.46,  17.95,  17.54,  17.15,  16.79,  16.56,  16.27,  16.03, \
        15.83,  15.69,  15.82,  16.15,  16.84,  18.21,  19.41,  20.31,  21.19,  22.03,  22.54, \
        23.01,  23.95,  24.38,  24.83,  25.14,  25.79,  26.18,  26.56,  26.76,  27.04,  27.34, \
        27.56,  27.73,  27.97,  27.69,  27.39,  27.07,  26.42,  25.32,  24.29,  23.38,  22.7 , \
        22.14,  21.66,  21.13,  20.58,  20.26,  19.81,  19.25,  18.6 ,  18.  ,  17.64,  17.31, \
        16.99,  16.71,  16.48,  16.25,  16.16,  16.25,  16.45,  16.44,  17.07,  18.71,  19.5 , \
        20.83,  22.12,  23.56,  24.  ,  24.96,  25.5 ,  26.17,  26.68,  27.38,  27.26,  27.72, \
        28.22,  28.47,  28.48,  28.35,  28.48,  28.46,  28.64,  27.99,  27.4 ,  26.78,  26.33, \
        25.9 ,  25.68,  25.4 ,  23.44,  20.66,  19.49,  19.08,  18.82,  18.44,  17.78,  17.34, \
        17.36,  17.53,  16.96,  17.  ,  16.85,  17.06,  17.03,  17.41,  17.84,  17.72,  16.79, \
        16.71,  17.63,  18.47,  19.48,  19.94,  20.09,  20.05])
        >>> vpd = np.array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, \
         0.00000000e+00, \
         0.00000000e+00,   0.00000000e+00,   1.98020999e-02,   1.01950016e+00,   2.59108448e+00, \
         2.79224405e+00,   3.61097769e+00,   5.84943300e+00,   7.99721893e+00,   9.51159857e+00, \
         1.03734359e+01,   1.13086021e+01,   1.19930192e+01,   1.28226314e+01,   1.34387706e+01, \
         1.14214143e+01,   9.31369610e+00,   8.80611609e+00,   8.28757853e+00,   7.00627076e+00, \
         5.54488177e+00,   4.72630946e+00,   4.51475873e+00,   4.46007333e+00,   5.39821005e+00, \
         5.12626351e+00,   4.76245510e+00,   3.81822052e+00,   3.03916126e+00,   2.51422732e+00, \
         2.00061176e+00,   1.65631656e+00,   1.31620085e+00,   1.10226324e+00,   9.38552406e-01, \
         7.26256912e-01,   6.40382295e-01,   4.80746063e-01,   3.09557999e-01,   2.69670728e-01, \
         1.95994866e-01,   2.15598651e-01,   3.85335122e-01,   9.01120083e-01,   1.83959118e+00, \
         2.70389046e+00,   3.71674716e+00,   5.43279632e+00,   7.09609170e+00,   8.43959755e+00, \
         1.00043485e+01,   1.31456014e+01,   1.57169883e+01,   1.63964561e+01,   1.82352581e+01, \
         1.90860008e+01,   2.06865461e+01,   2.16767294e+01,   2.26013587e+01,   2.31193279e+01, \
         2.39297487e+01,   2.45714542e+01,   2.48168287e+01,   2.49404901e+01,   2.44248239e+01, \
         2.43646736e+01,   2.40191893e+01,   2.33245459e+01,   1.97874046e+01,   1.64823803e+01, \
         1.34491912e+01,   1.16935674e+01,   1.03958655e+01,   8.98260645e+00,   7.76840720e+00, \
         6.61359691e+00,   6.05671080e+00,   5.26662608e+00,   4.10487999e+00,   2.87052587e+00, \
         1.91862877e+00,   1.33108665e+00,   8.88844215e-01,   5.80680504e-01,   3.99318894e-01, \
         2.62348372e-01,   2.03133448e-01,   1.83610005e-01,   2.58533479e-01,   5.23695983e-01, \
         6.35511933e-01,   1.20617136e+00,   3.27860665e+00,   4.16923392e+00,   5.83056602e+00, \
         7.58772705e+00,   1.00807880e+01,   1.22007393e+01,   1.60818029e+01,   1.79773729e+01, \
         2.01990530e+01,   2.17255275e+01,   2.32932803e+01,   2.36730370e+01,   2.55088493e+01, \
         2.72206272e+01,   2.79299836e+01,   2.65858273e+01,   2.48812733e+01,   2.49144961e+01, \
         2.60890796e+01,   2.65593102e+01,   2.27786127e+01,   2.06929075e+01,   1.91088445e+01, \
         1.80261873e+01,   1.64711586e+01,   1.61588594e+01,   1.55680917e+01,   1.08733644e+01, \
         4.79603812e+00,   2.37770382e+00,   1.81006901e+00,   1.36827315e+00,   1.46336785e+00, \
         1.62773793e+00,   1.46442951e+00,   1.68424260e+00,   2.06293605e+00,   9.46642981e-01, \
         9.29677819e-01,   5.75548103e-01,   5.83261735e-01,   4.07507984e-01,   4.77055621e-01, \
         1.24584895e+00,   1.58105608e+00,   5.54248696e-01,   4.37349264e-01,   8.46521127e-01, \
         1.55111406e+00,   2.17254953e+00,   3.07376831e+00,   3.17295993e+00,   3.04788792e+00])
        >>> data_flag = np.array([128, 128, 128, 128, 128, 128, 137, 129, 128,   1,   1,   1,   1, \
         1,   1,  32,   1,   1, \
         1,   0,   0,   1,   1,  32,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0, \
         0,   0,  32, 128, 128, 128, 128, 129, 128, 128, 128, 128, 128,  32,   0,   0,   0,   0, \
         0,   0,   0,  32,   1,   1,   1,   1,   1,  32,   1,   0,   0,   0,   0,   0,   0,   0, \
         1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  32, 128, 128, 128, \
         128, 128, 128, 128, 128, 129,   1,   1,  32,   0,   1,   0,   0,   0,  32,   1,   1,   1, \
         1,   1,   1,  32,   1,   0,   0,   0,   0,   0,   1,   0,   0,  32,   1,   1,   1, 137, \
         1,   9,   9, 129, 128, 128, 128, 128, 128, 129, 128, 128, 128, 128, 128, 128, 129, 128, \
         129,   1,   1,   1,  32,   1])
        >>> rg_flag   = np.zeros(150)
        >>> tair_flag = np.zeros(150)
        >>> vpd_flag  = np.zeros(150)
        >>> rg_dev   = 50.
        >>> tair_dev = 2.5
        >>> vpd_dev  = 5.
        >>> out, qout = gapfill(dates, data, rg, tair, vpd, \
                                data_flag, rg_flag, tair_flag, vpd_flag, \
                                rg_dev, tair_dev, vpd_dev)
        >>> print out
        [  7.28490000e-01   6.02317500e-01   8.98950000e-02  -9.66602750e+00
          -1.10983017e+01  -3.30030500e+00  -4.19937000e+00  -5.95061500e+00
          -5.95061500e+00  -5.95207000e+00  -9.87069200e+00  -1.25801000e+01
          -1.39628500e+01  -9.99900000e+03  -1.33583667e+01  -1.17411333e+01
          -7.88111200e+00  -1.41147000e+01  -1.74201500e+01  -1.70566000e+01
          -1.99620000e+00  -7.52963500e+00  -3.69177000e+00  -3.30030500e+00
           8.98950000e-02  -6.94000000e-03   4.73680000e-01   3.52585444e+00
           2.22053000e+00   2.53193000e+00   2.96486000e+00   3.55745000e+00
           4.49445000e+00   8.03341000e+00   8.64464000e+00   5.99874000e+00
           3.36288000e+00   4.29297000e+00   5.60643556e+00   5.10896000e+00
           4.51968000e+00   4.46154400e+00   4.07724500e+00   4.07212000e+00
           4.42674000e+00   9.80835000e-01   7.28490000e-01  -9.99900000e+03
          -3.30030500e+00  -3.36944333e+00  -5.21457000e+00  -5.31815000e+00
          -8.80319000e+00  -1.93282000e+01  -1.60800000e+01  -1.18457000e+01
          -1.33145000e+01  -1.46377000e+01  -9.99900000e+03  -9.99900000e+03
          -9.99900000e+03  -9.99900000e+03  -8.14880000e+00  -8.14880000e+00
          -1.49863667e+01  -6.81524000e+00  -9.77187000e+00  -8.31190000e+00
          -7.15069000e+00  -6.68666000e+00  -5.38734000e+00  -3.87582000e+00
          -3.26662222e-01   5.16840000e-01   1.48799000e+00   1.95226000e+00
           3.42510000e+00   3.73472000e+00   3.78379000e+00   3.02995000e+00
           4.86265000e+00   4.28708000e+00   6.66179000e+00   4.81036000e+00
           4.09262000e+00   4.56051000e+00   5.60643556e+00   5.10896000e+00
           5.10896000e+00   4.46154400e+00   4.46154400e+00   4.07724500e+00
           4.07724500e+00   9.80835000e-01   8.36842000e-01   7.79452000e-01
          -3.30030500e+00  -3.69177000e+00  -4.19937000e+00  -6.58599000e+00
          -1.32934450e+01  -1.77837000e+01  -1.21494000e+01  -1.97989000e+01
          -1.26952667e+01  -1.46377000e+01  -9.99900000e+03  -9.99900000e+03
          -9.99900000e+03  -9.99900000e+03  -8.14880000e+00  -1.25801000e+01
          -1.33378480e+01  -9.48236000e+00  -6.41284000e+00  -8.03108000e+00
          -7.10870000e+00  -2.72479000e+00  -3.39086556e+00   1.86730000e-01
           1.71264000e+00   9.80835000e-01   4.34774714e+00   4.34774714e+00
           3.18516400e+00   4.58357063e+00   5.04460071e+00   5.20458077e+00
           5.23307500e+00   5.38540455e+00   5.47450000e+00   5.30306375e+00
           5.30306375e+00   5.60643556e+00   5.10896000e+00   5.10896000e+00
           4.51968000e+00   5.10896000e+00   5.10896000e+00   5.60643556e+00
           5.47450000e+00   5.60643556e+00   1.04778750e+00   6.02317500e-01
          -5.70623667e+00  -3.36944333e+00  -3.30030500e+00  -2.86560333e+00
          -1.40794975e+01   8.98950000e-02]
        >>> print np.array(qout,dtype=np.int)
        [1 1 1 2 2 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 0 2 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0
         0 1 1 1 1 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0 2 0 0 0 0 2 2 1 0 0 0 0 0 0 0 2 0
         0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 1 2 0 0 0 0 2
         1 1 0 0 0 0 0 2 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1
         2 1]
        >>> sout = gapfill(dates, data, rg, tair, vpd, data_flag, err=True)
        >>> print np.where(sout != -9999., np.array(np.abs(sout/out)*100.,dtype=np.int), -1)
        [  -1   -1   -1   -1   -1   -1   -1   -1   -1   15   -1   -1   21   -1   17
           -1   -1   19    2    3   -1   -1   -1   -1   -1 4896   -1   25   79   67
           68   55   43   22   19   27   48   41   32   34   19   21   12   15    4
           -1   -1   -1   -1   -1   -1   16   72   -1   17   -1   -1   -1   -1   -1
           -1   -1   -1   -1   -1   27   14   11    0   -1   -1   20   -1  142   10
           44   21   25   33   52   37   40   26   34   40   39   32   34   34   21
           21   12   12   -1   -1   -1   -1   -1   -1   13   47   28   31   27   -1
           -1   -1   -1   -1   -1   -1   -1   -1   19   37   11    0   29   -1  429
           43   -1   -1   -1   23   41   33   31   32   31   32   32   32   32   34
           34   19   34   34   32   32   32   -1   -1   -1   -1   -1   -1   -1   -1]


        History
        -------
        Written  MC, Mar 2012 - modified gap_filling.py
    """

    # -------------------------------------------------------------
    # Checks

    # remember shape is any
    inshape = np.shape(datain)
    date = np.squeeze(datein)
    data = np.squeeze(datain)
    rg   = np.squeeze(rgin)
    tair = np.squeeze(tairin)
    vpd  = np.squeeze(vpdin)
    if np.size(np.shape(date)) != 1:
        raise ValueError('Error gapfill: squeezed dates must be 1D array.')
    if np.size(np.shape(data)) != 1:
        raise ValueError('Error gapfill: squeezed data must be 1D array.')
    if np.size(np.shape(rg)) != 1:
        raise ValueError('Error gapfill: squeezed rg must be 1D array.')
    if np.size(np.shape(tair)) != 1:
        raise ValueError('Error gapfill: squeezed tair must be 1D array.')
    if np.size(np.shape(vpd)) != 1:
        raise ValueError('Error gapfill: squeezed vpd must be 1D array.')

    # check flags
    ndata = np.size(data)
    if (data_flag != None):
        data_flg = np.squeeze(data_flag)
    else:
        data_flg = np.where(np.squeeze(data) == undef, 1, 0)
    if (rg_flag != None):
        rg_flg = np.squeeze(rg_flag)
    else:
        rg_flg = np.zeros(ndata)
    if (tair_flag != None):
        tair_flg = np.squeeze(tair_flag)
    else:
        tair_flg = np.zeros(ndata)
    if (vpd_flag != None):
        vpd_flg = np.squeeze(vpd_flag)
    else:
        vpd_flg = np.zeros(ndata)

    if ((np.size(date) != ndata) | (np.size(rg) != ndata) | (np.size(tair) != ndata) |
        (np.size(vpd) != ndata) | (np.size(data_flg) != ndata) | (np.size(rg_flg) != ndata) |
        (np.size(tair_flg) != ndata) | (np.size(vpd_flg) != ndata)):
        raise ValueError('Error gapfill: inputs must have the same size.')
        
    # -------------------------------------------------------------
    # Parameters

    # number of data points per week; basic factor of the time
    # window
    ddate    = date-np.roll(date,1)
    ddate[0] = ddate[1]
    if np.any((ddate-ddate[0]) > 1e-4): # ca 2.secs
        raise ValueError('Error gapfill: dates must be equally spaced.')
    week    = np.int(np.around(7./ddate[0]))
    nperday = week / 7
    hour    = (np.array(np.floor((date-np.trunc(date))*24.), dtype=np.int) + 12) % 24

    ndata = np.size(data)
    if err:
        # error estimate
        data_std  = np.ones(ndata)*undef
    else:
        # array for filled values
        data_fill = np.where(data_flg == 0, data, undef)
        # gap quality classes
        quality   = np.zeros(ndata)

    #------------------------------------------------------------------
    # Gap filling

    # flag for all meteorological conditions
    meteo_flg  = (tair_flg == 0) & (vpd_flg == 0) & (rg_flg == 0)
    # flag for all meteorological conditions and data
    total_flag = meteo_flg & (data_flg == 0)

    # Check for large margins at beginning
    largemargin = np.zeros(ndata)
    firstvalid  = np.amin(np.squeeze(np.where(data_flg==0)))
    lastvalid   = np.amax(np.squeeze(np.where(data_flg==0)))
    nn          = nperday*longestmarginalgap
    if firstvalid > nn:
        largemargin[0:(firstvalid-nn)] = 1
    if lastvalid < (ndata-nn):
        largemargin[(lastvalid+nn):]   = 1

    # Fill loop over all data points
    for j in xrange(ndata):
        if not err:
            # no reason to go further, no gap -> continue
            if data_flg[j] == 0 or largemargin[j] == 1: continue
        # 3 Methods
        #   1. tair, vpd and global radiation;
        #   2. just global radiation;
        #   3. none meteorolgical conditions: take the mean of +- hour

        # for better overview: dynamic calculation of radiation threshold
        # minimum 20; maximum 50 [Wm-2] according to private correspondence
        # with Markus Reichstein
        rg_devmax = np.maximum(20,np.minimum(rg[j],rg_dev))

        # Method 1: all met conditions
        if meteo_flg[j]:
            # search for values around the met-conditions in a window of time
            # (one week in the first iteration and odd weeks in the next)
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,ndata-1))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                           (np.abs(tair[win]-tair[j]) < tair_dev) &
                           (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    # assign also quality category of gap filling
                    quality[j]   = 1
                continue
            else: # --> extend time window to two weeks
                j1  = j - np.arange(1,2*week+1) + 1
                j2  = j + np.arange(1,2*week)
                jj  = np.append(j1,j2)
                win = np.sort(np.clip(jj,0,ndata-1))
                conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                               (np.abs(tair[win]-tair[j]) < tair_dev) &
                               (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        quality[j]   = 1
                    continue

        # if you come here, no error estimate
        if err: continue

        # If nothing is found under similar meteo within two weeks,
        # look for global radiation within one week ->

        # Method 2: just global radiation available
        if rg_flg[j] == 0:
            j1  = j - np.arange(1,week+1) + 1
            j2  = j + np.arange(1,week)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,ndata-1))
            # get boolean array where meteo-conditions are in a given width
            conditions = ( (np.abs(rg[win]  -rg[j]) < rg_devmax) &
                           total_flag[win] )
            num4avg = np.sum(conditions)
            # we need at least two samples with similar conditions
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j]   = 1
                continue

        # If still nothing is found under similar rg within one week,
        # take the same hour within 1-7 days

        # Method 3: same hour
        for i in xrange(2):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,ndata-1))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (data_flg[win] == 0)
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    if i == 0:
                        quality[j] = 1
                    else:
                        quality[j] = 2
                continue

        # sanity check
        if err:
            if data_std[j]  != undef: continue
        else:
            if data_fill[j] != undef: continue

        # If still nothing is found, start a new cycle with increased window size
        # Method 4: same as 1 but for 3-12 weeks
        if meteo_flg[j]:
            for multi in xrange(3,12):
                j1  = j - np.arange(1,multi*week+1) + 1
                j2  = j + np.arange(1,multi*week)
                jj  = np.append(j1,j2)
                win = np.sort(np.clip(jj,0,ndata-1))
                conditions = ( (np.abs(rg[win]  -rg[j])   < rg_devmax) &
                               (np.abs(tair[win]-tair[j]) < tair_dev) &
                               (np.abs(vpd[win] -vpd[j])  < vpd_dev) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        # assign also quality category of gap filling
                        if multi <= 2:
                            quality[j] = 1
                        elif multi > 4:
                            quality[j] = 3
                        else:
                            quality[j] = 2
                        break
            # Check because continue does not support to jump out of two loops
            if err:
                if data_std[j]  != undef: continue
            else:
                if data_fill[j] != undef: continue

        # Method 5: same as 2 but for 2-12 weeks
        if rg_flg[j] == 0:
            for multi in xrange(2,12):
                j1  = j - np.arange(1,multi*week+1) + 1
                j2  = j + np.arange(1,multi*week)
                jj  = np.append(j1,j2)
                win = np.sort(np.clip(jj,0,ndata-1))
                # get boolean array where meteo-conditions are in a given width
                conditions = ( (np.abs(rg[win]  -rg[j]) < rg_devmax) &
                               total_flag[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=(~conditions))
                    if err:
                        data_std[j] = np.ma.std(dat, ddof=ddof)
                    else:
                        data_fill[j] = np.ma.mean(dat)
                        if multi ==0:
                            quality[j] = 1
                        elif multi <= 2:
                            quality[j] = 2
                        else:
                            quality[j] = 3
                        break
            if err:
                if data_std[j]  != undef: continue
            else:
                if data_fill[j] != undef: continue

        # Method 6: same as 3 but for 3-120 days
        for i in xrange(3,120):
            t_win = nperday * (2*i+1)/2
            j1  = j - np.arange(1,t_win+1) + 1
            j2  = j + np.arange(1,t_win)
            jj  = np.append(j1,j2)
            win = np.sort(np.clip(jj,0,ndata-1))
            conditions = (np.abs(hour[win]-hour[j]) < 1.1) & (data_flg[win] == 0)
            num4avg = np.sum(conditions)
            if num4avg >= 2:
                dat = np.ma.array(data[win], mask=(~conditions))
                if err:
                    data_std[j] = np.ma.std(dat, ddof=ddof)
                else:
                    data_fill[j] = np.ma.mean(dat)
                    quality[j] = 3
                break

    if shape != False:
        if shape != True:
            if err:
                return np.reshape(data_std,shape)
            else:
                return np.reshape(data_fill,shape), np.reshape(quality,shape)
        else
            if err:
                return np.reshape(data_std,inshape)
            else:
                return np.reshape(data_fill,inshape), np.reshape(quality,inshape)
    else:
        if err:
            return data_std
        else:
            return data_fill, quality


if __name__ == '__main__':
    import doctest
    doctest.testmod()
