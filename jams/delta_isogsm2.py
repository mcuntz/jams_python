#!/usr/bin/env python
"""
Delta values from downloaded IsoGSM2 data.

Calculates delta values either from the original yearly IsoGSM2 files or
from the (multi-year) files produced get_isogsm2.


--------------------------------------------------------
usage: delta_isogsm2.py [-h] [-o output_file] [-t timestep]
                        [isogsm_output_file]

Calculate delta values from IsoGSM2 output.

positional arguments:
  isogsm_output_file    Output file of IsoGSM2.

optional arguments:
  -h, --help            show this help message and exit
  -o output_file, --output output_file
                        Output file name (default: ifile-delta.csv).
  -t timestep, --timestep timestep
                        Output timestep in seconds (default: IsoGSM2 timestep
                        of 6 hours = 21600).


Example
-------
# Hesse
python get_isogsm2.py -l 48.6742166667,7.06461666667
python delta_isogsm2.py -o delta_hesse_isogsm2_1979-2017.dat \
    lat48.571lon7.500_isogsm2_6hrly_1979-2017.dat

for i in $(ls lat*_isogsm2_*.dat) ; do python delta_isogsm2.py ${i} ; done


History
-------
Written,  Matthias Cuntz, Oct 2018
Modified, Matthias Cuntz, Nov 2018
              - make it callable function as well as script.
          Matthias Cuntz, Nov 2018
              - flak8 compatible
              - missing value of IsoGSM2 is 9.999e20
                -> precip, shm > 1e+15 = -1000
         Matthias Cuntz, Sep 2021
              - no error but only message when called without IsoGSM2 output
         Matthias Cuntz, Sep 2022
              - add timestep

"""
from datetime import datetime
import numpy as np
import netCDF4 as nc


# -------------------------------------------------------------------------
# Command line arguments
#

def delta_isogsm2(ifile, output=None, timestep=21600):
    '''
    Delta values from downloaded IsoGSM2 data

    Calculate delta values either from the original yearly IsoGSM2 files
    or from the (multi-year) files produced by get_isogsm2.
    Values can get interpolated to a desired output timestep.

    The script expects that the year(s) is(are) given at the end of the
    filename before the suffix, separated by _ from the rest of the
    filename, because there is no date info in the IsoGSM2 files.
    A year range is given by -, e.g.
        x001y192_isogsm2_6hrly_1999.dat
    or
        lat48.571lon7.500_isogsm2_6hrly_1979-2017.dat

    Parameters
    ----------
    ifile : string
        filename of IsoGSM2 output
    output : string, optional
        output filename. Default: ifile-delta.csv
    timestep : int, optional
        output time step in seconds. Default: 21600 (= 6 hours)
        IsoGSM2 output is 6-hourly. Interpolation will be linear for
        everything but precipitation. The latter will be distributed
        equally within the 6 hours.

    Returns
    -------
    filename of output file

    Notes
    -----
    The output file has the following columns:
       Date - timestep in the middle of the time period
       precip (mm) - cumulative precipitation during time step
       d18oprecip (permil) - delta-18O of precipitation
       d2hprecip (permil) - delta-2H of precipitation
       sh2m (kg/kg) - specific humidity at 2 m
       d18ovap (permil) - delta-18O of atmospheric vapour
       d2hvap (permil) - delta-2H of atmospheric vapour
       tair (K) - air temperature
       rh (percent) - relative humidity
       p (Pa) - surface pressure

    All delta values are relative to V-SMOW.

    Examples
    --------
    # All available IsoGSM2 output for FR-Hes
    latlon = '48.6742166667,7.06461666667'
    fisogsm2 = get_isogsm2(latlon)
    # delta values at half-hourly resolution
    disogsm2 = delta_isogsm2(fisogsm2, timestep=30*60)

    '''
    if timestep > 21600:
        raise ValueError('Timestep must be less or equal than 6 hours,'
                         ' i.e. 21600 seconds')
    if (86400 % timestep) > 0:
        raise ValueError(f'Seconds per day (86400) must be a multiple'
                         f' of timestep: {timestep}.')
    if (21600 % timestep) > 0:
        raise ValueError(f'IsoGSM2 timestep, 6 hours = 21600 seconds, must be'
                         f' a multiple of timestep: {timestep}.')

    if (output is None):
        output = ifile[0:ifile.rfind('.')] + '-delta.csv'

    # -------------------------------------------------------------------------
    # Data
    #

    # Lon= 7.5 Lat= 48.571
    # PWAT PWAT18O PWATHDO PRATE PRATE18O PRATEHDO CPRAT CPRAT18O CPRATHDO
    # LHTFL LHTFL18O LHTFLHDO SPFH2m SPFH2m18O SPFH2mHDO TMP2m RH2m PRESsfc
    # HGT500
    # 0.6678350E+01 0.6431830E+01 0.4832670E+01 0.3923000E-04 0.3868800E-04
    # 0.3539200E-04 0.0000000E+00 0.0000000E+00 0.0000000E+00 -0.1431808E+02
    # -0.1416832E+02 -0.1319024E+02 0.2095527E-02 0.2035712E-02 0.1645568E-02
    # 0.2643900E+03 0.1000000E+03 0.9299900E+05 0.5445400E+04

    # print('Read ', ifile)
    # Determine years from filename
    ff1 = ifile[0:ifile.rfind('.')]
    yrs = ff1[ff1.rfind('_') + 1:]
    if '-' in yrs:
        try:
            yr1, yr2 = [ int(yy) for yy in yrs.split('-') ]
        except:
            raise ValueError('Could not determine years from file name. (1)')
    else:
        try:
            yr1 = yr2 = int(yrs)
        except:
            raise ValueError('Could not determine years from file name. (2)')
    # nyr = yr2 - yr1 + 1

    # Read file
    ff = open(ifile, 'r')
    lon, lat = [ float(ll) for ll in ff.readline().strip().split()[1::2] ]
    head = ff.readline().split()
    ff.close()
    dat = np.loadtxt(ifile, skiprows=2)
    ndat = dat.shape[0]

    # Make dates
    units = 'days since ' + str(yr1) + '-01-01 00:00:00'
    # date0 is central in time step
    gtimestep = 21600
    second = 0
    minute = 0
    hour = 3
    date0 = nc.date2num(datetime(yr1, 1, 1, hour, minute, second),
                        units)
    # fractional days since 01.01.Year1
    fdate = gtimestep / 86400.
    jdate = date0 + np.arange(ndat) * fdate
    # datetime objects
    ddate = nc.num2date(jdate, units)
    # ascii date (eng)
    adate = [ d.isoformat(sep=' ') for d in ddate ]

    # Get deltas
    vprecip = 'PRATE'
    vvap    = 'SPFH2m'
    vtair   = 'TMP2m'
    vrh     = 'RH2m'
    vp      = 'PRESsfc'

    precip     = dat[:, head.index(vprecip)]
    # permil
    d18oprecip = (dat[:, head.index(vprecip + '18O')] / precip - 1.) * 1000.
    d2hprecip  = (dat[:, head.index(vprecip + 'HDO')] / precip - 1.) * 1000.
    # current undef of IsoGSM is 9.999e20
    ii = np.where((precip < 1.e-15) | (precip > 1.e+15))[0]
    if ii.size > 0:
        d18oprecip[ii] = -1000.
        d2hprecip[ii]  = -1000.
    sh2m    = dat[:, head.index(vvap)]  # specific humidity (kg/kg)
    d18ovap = (dat[:, head.index(vvap + '18O')] / sh2m - 1.) * 1000.
    d2hvap  = (dat[:, head.index(vvap + 'HDO')] / sh2m - 1.) * 1000.
    ii = np.where((sh2m < 1.e-15) | (sh2m > 1.e+15))[0]
    if ii.size > 0:
        d18ovap[ii] = -1000.
        d2hvap[ii]  = -1000.
    tair = dat[:, head.index(vtair)]  # K
    rh   = dat[:, head.index(vrh)]    # %
    p    = dat[:, head.index(vp)]     # Pa

    # interpolate to timestep
    if timestep == 21600:
        odate = adate
    else:
        # date0 is central in time step
        ihoy = int(timestep // 2)
        second = ihoy % 60
        ihoy   = ihoy // 60
        minute = ihoy % 60
        ihoy   = ihoy // 60
        hour = ihoy % 24
        ihoy = ihoy // 24
        # hour0 = np.rint(fdate * 24.)
        odate0 = nc.date2num(datetime(yr1, 1, 1, hour, minute, second),
                             units)
        fodate = timestep / 86400.
        nodat = (ndat // 4) * (86400 // timestep)
        jodate = odate0 + np.arange(nodat) * fodate
        dodate = nc.num2date(jodate, units)
        odate = [ d.isoformat(sep=' ') for d in dodate ]
        # continues variables
        sh2m = np.interp(jodate, jdate, sh2m)
        d18ovap = np.interp(jodate, jdate, d18ovap)
        d2hvap = np.interp(jodate, jdate, d2hvap)
        tair = np.interp(jodate, jdate, tair)
        rh = np.interp(jodate, jdate, rh)
        p = np.interp(jodate, jdate, p)
        # rain isotopes
        # interpolate over gaps and set to undef later
        ii = np.where(d18oprecip != -1000.)[0]
        if ii.size > 0.:
            d18oprecip = np.interp(jodate, jdate[ii], d18oprecip[ii])
        else:
            d18oprecip = np.full(nodat, -1000.)
        ii = np.where(d2hprecip != -1000.)[0]
        if ii.size > 0.:
            d2hprecip = np.interp(jodate, jdate[ii], d2hprecip[ii])
        else:
            d2hprecip = np.full(nodat, -1000.)
        # rain
        # distribute equally over 6 hours
        precip = precip * (timestep / gtimestep)
        oprecip = np.full(nodat, 0.)
        n = gtimestep // timestep
        for i in range(ndat):
            oprecip[i * n:(i + 1) * n] = precip[i]
        precip = oprecip
        ii = np.where(precip == 0.)[0]
        if ii.size > 0:
            d18oprecip[ii] = -1000.
            d2hprecip[ii] = -1000.

    # Write file
    # print('Write ', output)
    with open(output, 'w') as ff:
        # print('Lon={:.2f},Lat={:.3f}'.format(lon,lat), file=ff)
        astr = ('Date'
                ',precip (mm),d18oprecip (permil),d2hprecip (permil)'
                ',sh2m (kg/kg),d18ovap (permil),d2hvap (permil)'
                ',tair (K),rh (percent),p (Pa)')
        print(astr, file=ff)
        for i in range(tair.size):
            astr = odate[i]
            astr = astr + ',{:.6e},{:.4f},{:.4f}'.format(
                precip[i], d18oprecip[i], d2hprecip[i])
            astr = astr + ',{:.6e},{:.4f},{:.4f}'.format(
                sh2m[i], d18ovap[i], d2hvap[i])
            astr = astr + ',{:.2f},{:.4f},{:.0f}'.format(
                tair[i], rh[i], p[i])
            print(astr, file=ff)

    return output


if __name__ == '__main__':

    import argparse

    timestep = 21600  # 6 * 3600
    ofile  = None
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''Calculate delta values from IsoGSM2 output.''')
    parser.add_argument('-o', '--output', action='store',
                        default=ofile, dest='ofile', metavar='output_file',
                        help='Output file name (default: ifile-delta.csv).')
    parser.add_argument('-t', '--timestep', action='store',
                        default=timestep, dest='timestep',
                        metavar='timestep',
                        help=('Output timestep in seconds (default:'
                              ' IsoGSM2 timestep of 6 hours = 21600).'))
    parser.add_argument('ifile', nargs='?', default=None,
                        metavar='isogsm_output_file',
                        help='Output file of IsoGSM2.')

    args = parser.parse_args()
    ofile    = args.ofile
    timestep = int(args.timestep)
    ifile    = args.ifile

    del parser, args

    # Check input
    if (ifile is None):
        print('No IsoGSM2 output file given.')
        import sys
        sys.exit()

    import time as ptime
    t1 = ptime.time()

    dfile = delta_isogsm2(ifile, output=ofile, timestep=timestep)
    print('File: ', dfile)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.) if (t2 - t1) > 60. else
             '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed', strin)
