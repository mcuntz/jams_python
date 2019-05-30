#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Delta values from downloaded IsoGSM2 data.

    Calculates delta values either from the original yearly IsoGSM2 files or
    from the (multi-year) files produced get_isogsm2.


    --------------------------------------------------------
    usage: delta_isogsm2.py [-h] [-o output_file] [isogsm_output_file]

    Calculate delta values from IsoGSM2 output.

    positional arguments:
      isogsm_output_file    Output file of IsoGSM2.

    optional arguments:
      -h, --help            show this help message and exit
      -o output_file, --output output_file
                            Output file name (default: ifile-delta.csv).


    Example
    -------
    # Hesse
    python get_isogsm2.py -l 48.6742166667,7.06461666667
    python delta_isogsm2.py -o delta_hesse_isogsm2_1979-2017.dat lat48.571lon7.500_isogsm2_6hrly_1979-2017.dat

    for i in $(\ls lat*_isogsm2_*.dat) ; do python delta_isogsm2.py ${i} ; done


    History
    -------
    Written,  MC, Oct 2018
    Modified, MC, Nov 2018 - make it callable function as well as script.
"""

# -------------------------------------------------------------------------
# Command line arguments
#

def delta_isogsm2(ifile, output=None):
    '''
        Delta values from downloaded IsoGSM2 data.
    
        Calculates delta values either from the original yearly IsoGSM2 files or
        from the (multi-year) files produced get_isogsm2.

        The script expects that the year(s) are given at the end of the file name before the suffix,
        separated by _ from the rest of the filename. A year range is given by -, e.g.
            x001y192_isogsm2_6hrly_1999.dat
            lat48.571lon7.500_isogsm2_6hrly_1979-2017.dat


        Definition
        ----------
        def delta_isogsm2(ifile, output=None)


        Input
        -----
        ifile         string
                      filename of IsoGSM2 output


        Parameters
        ----------
        output        string (default: ifile-delta.csv')
                      output filename


        Ouput
        -----
        Returns filename of the output file


        Examples
        --------
        >>> latlon = '48.6742166667,7.06461666667'
        >>> ff = get_isogsm2(latlon)
        >>> df = [ delta_isogsm2(ifile) for ifile in ff ]
        >>> import os
        >>> file1 = 'lat48.571lon7.500_isogsm2_6hrly_1979-2017-delta.csv'
        >>> if not os.path.exists(file1): print('No file: ', file1)


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2018 Matthias Cuntz - mc (at) macu (dot) de

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
        Written  MC, Oct 2018
        Modified MC, Nov 2018 - make it callable function as well as script.
    '''
    if (output is None):
        output = ifile[0:ifile.rfind('.')]+'-delta.csv'

    from datetime import datetime
    import numpy as np
    import netCDF4 as nc


    # -------------------------------------------------------------------------
    # Data
    #

    # Lon= 7.5 Lat= 48.571
    # PWAT PWAT18O PWATHDO PRATE PRATE18O PRATEHDO CPRAT CPRAT18O CPRATHDO LHTFL LHTFL18O LHTFLHDO SPFH2m SPFH2m18O SPFH2mHDO TMP2m RH2m PRESsfc HGT500
    # 0.6678350E+01 0.6431830E+01 0.4832670E+01 0.3923000E-04 0.3868800E-04 0.3539200E-04 0.0000000E+00 0.0000000E+00 0.0000000E+00 -0.1431808E+02 -0.1416832E+02 -0.1319024E+02 0.2095527E-02 0.2035712E-02 0.1645568E-02 0.2643900E+03 0.1000000E+03 0.9299900E+05 0.5445400E+04

    # print('Read ', ifile)
    # Determine years from filename
    ff1 = ifile[0:ifile.rfind('.')]
    yrs = ff1[ff1.rfind('_')+1:]
    if '-' in yrs:
        try:
            yr1, yr2 = [ int(yy) for yy in yrs.split('-') ]
        except:
            raise ValueError('\nCould not determine years from file name.')
    else:
        try:
            yr1 = y2 = int(yrs)
        except:
            raise ValueError('\nCould not determine years from file name.')
    nyr = yr2 - yr1 + 1

    # Read file
    ff = open(ifile, 'r')
    lon, lat = [ float(ll) for ll in ff.readline().strip().split()[1::2] ]
    head = ff.readline().split()
    ff.close()
    dat = np.loadtxt(ifile, skiprows=2)
    ndat = dat.shape[0]

    # Make dates
    units = 'days since '+str(yr1)+'-01-01 00:00:00'
    date0 = nc.date2num(datetime(yr1,1,1,3,0,0), units) # 0.125 # central in time step
    jdate = date0 + np.arange(ndat)*0.25                # fractional days since 01.01.Year1
    ddate = nc.num2date(jdate, units)                   # datetime objects
    adate = [ d.isoformat(sep=' ') for d in ddate ]     # ascii date (eng)

    # Get deltas
    vprecip = 'PRATE'
    vvap    = 'SPFH2m'
    vtair   = 'TMP2m'
    vrh     = 'RH2m'
    vp      = 'PRESsfc'

    precip     = dat[:,head.index(vprecip)] # mm
    d18oprecip = (dat[:,head.index(vprecip+'18O')] / precip - 1.) * 1000. # permil
    d2hprecip  = (dat[:,head.index(vprecip+'HDO')] / precip - 1.) * 1000.
    ii = np.where(precip < 1.e-15)[0]
    if ii.size > 0:
        d18oprecip[ii] = -1000.
        d2hprecip[ii]  = -1000.
    sh2m    = dat[:,head.index(vvap)] # specific humidity (kg/kg)
    d18ovap = (dat[:,head.index(vvap+'18O')] / sh2m - 1.) * 1000.
    d2hvap  = (dat[:,head.index(vvap+'HDO')] / sh2m - 1.) * 1000.
    ii = np.where(sh2m < 1.e-15)[0]
    if ii.size > 0:
        d18ovap[ii] = -1000.
        d2hvap[ii]  = -1000.
    tair = dat[:,head.index(vtair)] # K
    rh   = dat[:,head.index(vrh)]   # %
    p    = dat[:,head.index(vp)]    # Pa

    # Write file
    # print('Write ', output)
    ff = open(output, 'w')
    # print('Lon={:.2f},Lat={:.3f}'.format(lon,lat), file=ff)
    astr = 'Date'
    astr = astr + ',precip (mm),d18oprecip (permil),d2hprecip (permil)'
    astr = astr + ',sh2m (kg/kg),d18ovap (permil),d2hvap (permil)'
    astr = astr + ',tair (K),rh (percent),p (Pa)'
    print(astr, file=ff)
    for i in range(tair.size):
        astr = adate[i]
        astr = astr + ',{:.6e},{:.4f},{:.4f}'.format(precip[i], d18oprecip[i], d2hprecip[i])
        astr = astr + ',{:.6e},{:.4f},{:.4f}'.format(sh2m[i], d18ovap[i], d2hvap[i])
        astr = astr + ',{:.2f},{:.4f},{:.0f}'.format(tair[i], rh[i], p[i])
        print(astr, file=ff)
    ff.close()

    return output


if __name__ == '__main__':

    import argparse

    ofile  = None
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''Calculate delta values from IsoGSM2 output.''')
    parser.add_argument('-o', '--output', action='store',
                        default=ofile, dest='ofile', metavar='output_file',
                        help='Output file name (default: ifile-delta.csv).')
    parser.add_argument('ifile', nargs='?', default=None, metavar='isogsm_output_file',
                        help='Output file of IsoGSM2.')

    args  = parser.parse_args()
    ofile = args.ofile
    ifile = args.ifile

    del parser, args

    # Check input
    if (ifile is None):
        raise IOError('\nIsoGSM2 output file has to be given.\n')

    import time as ptime
    t1 = ptime.time()

    dfile = delta_isogsm2(ifile, output=ofile)
    print('File: ', dfile)

    t2    = ptime.time()
    strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    print('Time elapsed', strin)
