#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
'''
Get IsoGSM2 output for one or several lat,lon.

--------------------------------------------------------
usage: get_isogsm2.py [-h] [-l lat,lon] [-o] [-p path] [-u url] [latlon_file]

Get IsoGSM2 output for one or several lat,lon.

positional arguments:
  latlon_file           File with lat lon info (only these two header columns
                        must begin with lat and lon (case-insensitive). File
                        only used if no -l lat,lon given.

optional arguments:
  -h, --help            show this help message and exit
  -l lat,lon, --latlon lat,lon
                        latitude,longitude to extract from IsoGSM2 (default:
                        read from input file). If latitude or longitude is
                        negative, use the notation -l=lat,lon instead of -l
                        lat,lon.
  -o, --override        Do not check that output file already exists. Override
                        existing output file (default: False).
  -p path, --path path  Output directory (default: current directory ".").
  -u url, --url url     URL of IsoGSM2 base directory (default:
                        http://isotope.iis.u-tokyo.ac.jp/~kei/tmp/isogsm2).


Example
-------
# Hesse
python get_isogsm2.py -l 48.6742166667,7.06461666667

# Lucas site collection
python get_isogsm2.py Leaf_water_data_combined_v4.csv


History
-------
Written,  Matthias Cuntz, Oct 2018
Modified, Matthias Cuntz, Nov 2018
              - make it callable function as well as script.
          Matthias Cuntz, Jan 2021
              - (almost) flake8 compatible
'''
import os
import glob
import fileinput
import numpy as np
# Reading from HTTP
try:
    isrequests = True
    import requests_not_used
except ImportError:  # no requests module
    isrequests = False
    try:                 # Python 3
        import urllib.request as liburl
    except ImportError:  # Python 2
        import urllib as liburl
# Parsing HTML pages
try:
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup
try:
    bparser = 'lxml'
    import lxml
except:
    bparser = None
# Reading user input
try:               # Python 2
    input = raw_input
except NameError:  # Python 3
    pass


__all__ = ['get_isogsm2']


isogsm2base = 'http://isotope.iis.u-tokyo.ac.jp/~kei/tmp/isogsm2'
# IsoGSM2 T62 grid
isogsm2lons = np.arange(192) * 1.875
isogsm2lats = np.array(
    [-88.542, -86.653, -84.753, -82.851, -80.947, -79.043, -77.139, -75.235,
     -73.331, -71.426, -69.522, -67.617, -65.713, -63.808, -61.903, -59.999,
     -58.094, -56.189, -54.285, -52.380, -50.475, -48.571, -46.666, -44.761,
     -42.856, -40.952, -39.047, -37.142, -35.238, -33.333, -31.428, -29.523,
     -27.619, -25.714, -23.809, -21.904, -20.000, -18.095, -16.190, -14.286,
     -12.381, -10.476, -8.571, -6.667, -4.762, -2.857, -0.952,
     0.952, 2.857, 4.762, 6.667, 8.571, 10.476, 12.381, 14.286, 16.190, 18.095,
     20.000, 21.904, 23.809, 25.714, 27.619, 29.523, 31.428, 33.333, 35.238,
     37.142, 39.047, 40.952, 42.856, 44.761, 46.666, 48.571, 50.475, 52.380,
     54.285, 56.189, 58.094, 59.999, 61.903, 63.808, 65.713, 67.617, 69.522,
     71.426, 73.331, 75.235, 77.139, 79.043, 80.947, 82.851, 84.753, 86.653,
     88.542])


# -------------------------------------------------------------------------
# Functions
#

# Read URL
def get_url(url, **kwargs):
    ''' Returns the content of website or data file as string '''
    if isrequests:
        rr = requests.get(url, **kwargs)
        return rr.text
    else:
        if 'proxies' in kwargs:
            if kwargs['proxies'] is not None:
                proxy  = liburl.ProxyHandler(proxies)
                auth   = liburl.HTTPBasicAuthHandler()
                opener = liburl.build_opener(proxy, auth, liburl.HTTPHandler)
                liburl.install_opener(opener)
        rr = liburl.urlopen(url)
        return rr.read().decode('utf-8')


# Read URL, write to file
def url2filename(url, ifile, **kwargs):
    ''' Returns file with the local name "ifile" from an
        url such as http://website/file.dat '''
    rr = get_url(url, **kwargs)

    ff = open(ifile, 'w')
    print(rr, file=ff)
    ff.close()

    return


# Read URL, write to file with filename from URL
def url2file(url, **kwargs):
    ''' Returns file from an url such as http://website/file.dat
        having the same local name as the remote name '''
    if url[-1] == '/':
        raise ValueError('url does not end with valid filename.')

    ifile = url[url.rfind('/') + 1:]
    url2filename(url, ifile, **kwargs)

    return


# Get IsoGSM2 available years
def get_years(url, **kwargs):
    ''' Get all the years from IsoGSM2 base website '''
    page = get_url(url, **kwargs)

    bpage = BeautifulSoup(page, features=bparser)

    # '/~kei/tmp/', '1979/', ...
    href = [ i.get('href') for i in bpage.findAll('a') ]
    # '/~kei/tmp', '1979', ...
    hdir = [ i[:-1] for i in href if i.endswith('/') ]
    # 1979, 1980, ...
    yrs = []
    for hh in hdir:
        try:
            yrs.append(int(hh))
        except:
            pass
    return yrs


# Get IsoGSM2 files in year
def get_filenames(url, **kwargs):
    ''' Get all files in yearly directory of IsoGSM2 '''
    page = get_url(url, **kwargs)

    bpage = BeautifulSoup(page, features=bparser)

    # '/~kei/tmp/isogsm2/', 'x001y001_isogsm2_6hrly_2013.dat', ...
    href = [ i.get('href') for i in bpage.findAll('a') ]
    # 'x001y001_isogsm2_6hrly_2013.dat', ...
    # files = [ i for i in href if i.startswith('x') ]
    # 'x001y001_isogsm2_6hrly_2013.dat', ...
    files = [ i for i in href
              if (not i.startswith('?')) and (not i.endswith('/')) ]
    return files


# Get lat lon in IsoGSM2 grid
def closest(vec, num, value=False):
    ''' Index in array at which the entry is closest to a given number. '''
    out = np.argmin(np.abs(np.array(vec) - num))
    if value:
        return vec.flat[out]
    else:
        return out


def lonlat2xy(ilon, ilat):
    ''' Get x and y coordinates in IsoGSM2 grid of given lat and lon. '''
    if ilon < 0:
        ilon += 360.  # IsoGSM2 grid is 0-360.
    x = closest(isogsm2lons, ilon) + 1
    y = closest(isogsm2lats, ilat) + 1
    return (x, y)


# -------------------------------------------------------------------------
# Main routine
#

def get_isogsm2(latlon, baseurl=isogsm2base, path='.', override=False):
    '''Get IsoGSM2 output.

    Downloads all available years of IsoGSM2 output for a specific lat,lon,
    a list of lats and lons, or lat/lon read from columns of an input file.

    If override=False, checks if the data is already available in the
    download directory.
    It expects files with the same naming convention than its own, i.e.
        x, y = lonlat2xy(lon, lat)
        ofile = 'lat{:.3f}lon{:.3f}_isogsm2_6hrly_{:s}-{:s}.dat'.format(
                 isogsm2lats[y-1], isogsm2lons[x-1],
                 str(yrs[0]), str(yrs[-1]))
    It checks if the same lat/lon is present in the download directory by
    checking ONLY filenames.


    Definition
    ----------
    def get_isogsm2(latlon, baseurl=isogsm2base, path='.', override=False):


    Input
    -----
    latlon    string or list/tuple/ND-array
              latlon can be
                  string of one 'lat,lon', e.g. 48.5,7,
                  list, tuple or nd-array with iterables of lats and lons, e.g.
                      [[48.5,49],[7,15]], ([48.5,49],[7,15]),
                      (3,2)-array: np.array([[48.5,49,50],[7,15,-10]])
                  filename of file wherein two columns in the first header line
                  start with strings 'lat' and 'lon' (case-insensitive).
                  Only these two header columns must begin with lat and lon.


    Parameters
    ----------
    baseurl     string
                URL of IsoGSM2 base directory
                default: http://isotope.iis.u-tokyo.ac.jp/~kei/tmp/isogsm2
    path        string (default: '.')
                Output path
    override    boolean (default: False)
                If True, check if data is already available in path.
                Expect files with the naming convention
                x, y = lonlat2xy(lon, lat)
                ofile = 'lat{:.3f}lon{:.3f}_isogsm2_6hrly_{:s}-{:s}.dat'.format(
                    isogsm2lats[y-1], isogsm2lons[x-1],
                    str(yrs[0]), str(yrs[-1]))
                Check if lat/lon was already downloaded in path by checking
                ONLY filenames.


    Ouput
    -----
    Returns filenames of the output files that were newly written or already
    present in path:
       'lat{:.3f}lon{:.3f}_isogsm2_6hrly_{:s}-{:s}.dat'.format(
           isogsm2lats[y-1], isogsm2lons[x-1], str(yrs[0]), str(yrs[-1]))


    Examples
    --------
    >>> latlon = '48.6742166667,7.06461666667'
    >>> ff = get_isogsm2(latlon)
    >>> import os
    >>> file1 = 'lat48.571lon7.500_isogsm2_6hrly_1979-2017.dat'
    >>> if not os.path.exists(file1): print('No file: ', file1)


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2018-2021 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Matthias Cuntz, Oct 2018
    Modified, Matthias Cuntz, Nov 2018
                  - make it callable function as well as script.
              Matthias Cuntz, Jan 2021
                  - (almost) flake8 compatible
    '''
    # get lats and lons
    if isinstance(latlon, str):
        if ',' in latlon:
            ilat, ilon = [ float(i) for i in latlon.split(',') ]
            lats = [ilat]
            lons = [ilon]
        else:
            if not os.path.exists(latlon):
                raise IOError(f'Did not find latlon: {latlon}. latlon must'
                              f' be either a file name or a string "lat,lon".')
            lats = []
            lons = []
            ff = open(latlon, 'r')
            lhead = ff.readline().rstrip()
            for delim in [',', ';', None]:
                head  = lhead.split(delim)
                nhead = len(head)
                if nhead > 1:
                    break
            if nhead == 1:
                ff.close()
                raise IOError(f'Cannot detect delimter in latlon file:'
                              f' {latlon}')
            ii = [ i for i in range(nhead)
                   if head[i].lower().startswith('lat') ]
            if len(ii) != 1:
                ff.close()
                raise IOError(f'Found no or more than one column starting'
                              f' with "lat" in latlon file: {latlon}')
            iilat = ii[0]
            ii = [ i for i in range(nhead)
                   if head[i].lower().startswith('lon') ]
            if len(ii) != 1:
                ff.close()
                raise IOError(f'Found no or more than one column starting'
                              f' with "lon" in latlon file: {latlon}')
            iilon = ii[0]
            for line in ff:
                s  = line.rstrip()
                ss = s.split(delim)
                lats.append(ss[iilat])
                lons.append(ss[iilon])
            ff.close()
            if (len(lats) == 0) or (len(lons) == 0):
                raise IOError(
                    'Could not get lat,lon list from latlon file: ' +
                    latlon)
    elif isinstance(latlon, (tuple, list, np.ndarray)):
        if isinstance(latlon, (tuple, list)):
            lats = latlon[0]
            lons = latlon[1]
        else:
            lats = latlon[0, :]
            lons = latlon[1, :]
    else:
        raise IOError('latlon must be either a string "lat,lon", a filename,'
                      ' or an iterable with lats and lons.')
    nlat = len(lats)

    # make output directory
    if not os.path.exists(path):
        os.makedirs(path)

    # Check existing files
    hasfile  = [False] * nlat
    isofiles = [''] * nlat
    files = glob.glob(path + '/lat*.dat')
    if (len(files) > 0) and (not override):
        for ll in range(nlat):
            x, y = lonlat2xy(float(lons[ll]), float(lats[ll]))
            ilat = '{:.3f}'.format(isogsm2lats[y - 1])
            ilon = '{:.3f}'.format(isogsm2lons[x - 1])
            for ff in files:
                fs = os.path.basename(ff).split('_')[0][3:]
                lat, lon = fs.split('lon')
                if ((lat == ilat) and (lon == ilon)):
                    hasfile[ll]  = True
                    isofiles[ll] = ff
    # Return if all exist
    if np.all(hasfile):
        return isofiles

    # Detect IsoGSM2 years
    try:
        yrs = get_years(baseurl)
    except Exception as e:
        if str(e) == 'HTTP Error 407: Proxy Authentication Required':
            print('Enter your http proxy details:'
                  ' http://user:password@proxy_server:port')
            print('For example: http://mcuntz:password@revelec.inra.fr:3128')
            user = input('User (mcuntz):')
            if user == '':
                user = 'mcuntz'
            print('Password: ', end='', flush=True)
            import getpass
            password = getpass.getpass()
            proxy_server = input('Proxy server (revelec.inra.fr):')
            if proxy_server == '':
                proxy_server = 'revelec.inra.fr'
            port = input('Port (3128):')
            if port == '':
                port = '3128'
            proxy_string  = 'http://' + str(user) + ':' + str(password)
            proxy_string += '@' + str(proxy_server) + ':' + str(int(port))
            proxies = {"http": proxy_string}
            # needs proxies only once during session
            yrs = get_years(baseurl, proxies=proxies)
        else:
            raise IOError(str(e))

    for ll in range(nlat):
        if not hasfile[ll]:
            print('  Get lat lon: ', lats[ll], lons[ll])
            x, y = lonlat2xy(float(lons[ll]), float(lats[ll]))
            llfiles = []
            for yy in yrs:
                iifile = 'x{:03d}y{:03d}_isogsm2_6hrly_{:s}.dat'.format(
                    x, y, str(yy))
                url = baseurl + '/' + str(yy) + '/' + iifile
                print('    ', url)
                try:
                    url2file(url)
                    llfiles.append(path + '/' + iifile)
                except:
                    print('Could not find URL: ', url)

            # concatenate files
            isofiles[ll] = (
                path + '/' +
                'lat{:.3f}lon{:.3f}_isogsm2_6hrly_{:s}-{:s}.dat'.format(
                    isogsm2lats[y - 1], isogsm2lons[x - 1], str(yrs[0]),
                    str(yrs[-1])))
            print('  Write ', isofiles[ll])
            with open(isofiles[ll], 'w') as fout:
                fin = fileinput.input(llfiles)
                for line in fin:
                    # write header only of first file and
                    # strip emtpy lines at end of files
                    if ( ((fileinput.filelineno() > 2) or
                          (fileinput.lineno() < 3)) and
                         line.strip() != '' ):
                        fout.write(line)
                fin.close()

            # clean downloaded files
            for ll in llfiles:
                os.remove(ll)

    return isofiles


# -------------------------------------------------------------------------
# Command line arguments
#

if __name__ == '__main__':

    import argparse

    latlon   = None
    override = False
    path     = '.'
    baseurl  = isogsm2base
    parser  = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''Get IsoGSM2 output for one or several lat,lon.''')
    parser.add_argument(
        '-l', '--latlon', action='store',
        default=latlon, dest='latlon', metavar='lat,lon',
        help=('latitude,longitude to extract from IsoGSM2 (default: read'
              ' from input file). If latitude or longitude is negative,'
              ' use the notation -l=lat,lon instead of -l lat,lon.'))
    parser.add_argument(
        '-o', '--override', action='store_true', default=override,
        dest='override',
        help=('Do not check that output file already exists. '
              'Override existing output file (default: False).'))
    parser.add_argument(
        '-p', '--path', action='store', default=path, dest='path',
        metavar='path',
        help='Output directory (default: current directory ".").')
    parser.add_argument(
        '-u', '--url', action='store',
        default=baseurl, dest='baseurl', metavar='url',
        help='URL of IsoGSM2 base directory (default: ' + baseurl + ').')
    parser.add_argument(
        'ifile', nargs='?', default=None, metavar='latlon_file',
        help=('File with lat lon info (only these two header columns must'
              ' begin with lat and lon (case-insensitive). File only used'
              ' if no -l lat,lon given.'))

    args     = parser.parse_args()
    latlon   = args.latlon
    override = args.override
    path     = args.path
    baseurl  = args.baseurl
    ifile    = args.ifile

    del parser, args

    # Check input
    if (latlon is None) and (ifile is None):
        raise IOError('-l lat,lon or input file must be given.')

    if (latlon is not None):
        ilatlon = latlon
    else:
        ilatlon = ifile

    import time as ptime
    t1 = ptime.time()

    ofiles = get_isogsm2(ilatlon, baseurl=baseurl, path=path,
                         override=override)
    print('Files: ', ofiles)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.) if (t2 - t1) > 60. else
             '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed', strin)
