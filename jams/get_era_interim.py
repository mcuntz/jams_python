#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
'''
    Download ERA-Interim data suitable to produce MuSICA input data.
    
    If override=False (default), the script checks if the data is already available
    in the local download directory (path).
    It expects files with the same naming convention than its own, i.e.
         path + '/' + 'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
    and
         path + '/' + 'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
    Filenames can also end on .nc?.
    It checks if the area and years are included in a file in the download directory by checking ONLY filenames.

    Written  Matthias Cuntz, Feb 2018
    Modified Jerome Ogee,    May 2018 - figured out precip download, with the help of Sebastian Lafont
             Matthias Cuntz, Nov 2018 - retrieve all years at once because of long wait times per request
             Matthias Cuntz, Nov 2018 - make it callable function as well as script.


    --------------------------------------------------------
    usage: get_era_interim.py [-h] [-a area] [-y years] [-p path] [-o]

    Download ERA-Interim data suitable to produce MuSICA input data.

    optional arguments:
      -h, --help            show this help message and exit
      -a area, --area area  area format as either lat,lon or
                            NorthLat/WestLon/SouthLat/EastLon (default: global
                            90/-180/-90/180.
      -y years, --years years
                            years format is startyear,endyear (default:
                            1979,current-1).
      -p path, --path path  Output directory (default: current directory ".").
      -o, --override        Do not check that output file already exists that
                            includes request. Override existing output file
                            (default: False).


    Example
        # Hesse
        python get_era_interim.py -a 48.6742166667,7.06461666667 -y 1995,2017 -p era


    --------------------------------------------------------
    Script was originally adapted from extract-ERA-Interim-locations.py from
        https://confluence.ecmwf.int/display/CKB/How+to+extract+ERA-Interim+data+for+specific+locations
    It ignores efficiency suggestions from
        https://confluence.ecmwf.int/display/WEBAPI/ERA-Interim+daily+retrieval+efficiency
        "Iterate efficiently over several years and months for ERA-Interim request.
         Data is organised as follows:
         type of data (analysis, forecast, 4D analysis increments)
             year
                 month
                     type of level (model level, pressure level, surface)
                         dates, times, steps, levels, parameters (same tape file)"
    because the web-api server is very busy so it tries to get as much data as possible
    once the queued request gets treated. The web-api server also seems to take care of the efficiency
    by splitting the request by months.


    --------------------------------------------------------
    Original disclaimer of extract-ERA-Interim-locations.py:

    (C) Copyright 1996-2012 ECMWF.

    This software is licensed under the terms of the Apache Licence Version 2.0
    which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
    In applying this licence, ECMWF does not waive the privileges and immunities
    granted to it by virtue of its status as an intergovernmental organisation nor
    does it submit to any jurisdiction.

    Author: Karl Hennermann, ECMWF, copernicus-support@ecmwf.int
    Purpose: retrieve meteorological data from ECMWF's ERA-Interim dataset for given locations and time period
    Date/Time: 2017-02-26-1300

    Requires:
       ECMWF Web API key, https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
       ECcodes, https://software.ecmwf.int/wiki/display/ecc
       pandas, http://pandas.pydata.org/

    Credits:
       Christopher Blunck, http://pydoc.net/Python/weather/0.9.1/weather.units.temp/
       Brian McNoldy, http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
'''
import os
import glob

__all__ = ['get_era_interim']


# --------------------------------------------------------------------
# Retrieval function
#

def get_sfc_data(typ, date, time, step, param, area, target):

    from ecmwfapi import ECMWFDataServer
    server = ECMWFDataServer()
    
    rdict = {
        "class":   "ei",
        "expver":  "1",
        "dataset": "interim",
        "stream":  "oper",
        "levtype": "sfc",       # "sfc", "ml", "pl"
        "type":    typ,         # "an", "fc"
        "date":    date,        # "2010-01-01", "2010-01-01/2010-01-15/2010-01-30", "2010-01-01/to/2015-12-31"
        "time":    time,        # "00:00:00/06:00:00/12:00:00/18:00:00"
        "step":    step,        # 0/3/6/9/12
        "param":   param,       # 129/134/151/...
        "area":    area,        # within 90/-180/-90/180
        "format": "netcdf",     # file format
        "target":  target,      # filename
        }
    if area == "90/-180/-90/180": # full resolution for whole grid, i.e. "resol": "av"
        rdict["resol"] = "av"
    else:                         # otherwise "grid": "0.75/0.75"
        rdict["grid"] = "0.75/0.75"

    server.retrieve(rdict)

    return target


# --------------------------------------------------------------------
# Main routine
#

def get_era_interim(area=None, years=None, path='.', override=False):
    '''
        Download ERA-Interim data suitable to produce MuSICA input data.
    
        If override=False, checks if the data is already available in the download directory.
        It expects files with the same naming convention than its own, i.e.
             path+'/'+'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
             path+'/'+'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
        Filenames can also end on .nc?.
        It checks if the area and years are included in a file by checking ONLY filenames.


        Definition
        ----------
        def get_era_interim(area=None, years=None, path='.', override=False):


        Parameters
        ----------
        area          string (default: '90/-180/-90/180')
                      Area. Can be one point given as lat,lon
                      or box as NorthLat/WestLon/SouthLat/EastLon. Default is the globe.
        years         Iterable of len(2) [default: (1979,current_year-1)]
                      Year range.
        path          string (default: '.')
                      Output path
        override      boolean (default: False)
                      If True, check if data is already available in path.
                      Expect files with the naming convention
                          path+'/'+'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
                          path+'/'+'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
                      Check if area and years are included in a file in path by checking ONLY filenames.


        Ouput
        -----
        Returns file names of the two output files in path, either the newly written files
        or the files that contain the output requested.
        *_an_* for analysis variables and *_fc_* for forecast variables:
            path+'/'+'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
            path+'/'+'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)


        Restrictions
        ------------
        Existing files will only be checked by filename not by content.


        Examples
        --------
        >>> area  = '48/7/47/8'
        >>> years = (1995,2017)
        >>> get_era_interim(area=area, years=years, path='.')
        >>> import os
        >>> file1 = 'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
        >>> if not os.path.exists(file1): print('No file: ', file1)

        >>> file2 = 'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(*years)
        >>> if not os.path.exists(file2): print('No file: ', file2)


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

        Copyright 2018 Matthias Cuntz


        History
        -------
        Written  Matthias Cuntz, Feb 2018
        Modified Jerome Ogee,    May 2018 - figured out precip download, with the help of Sebastian Lafont
                 Matthias Cuntz, Nov 2018 - retrieve all years at once because of long wait times per request
                 Matthias Cuntz, Nov 2018 - make it callable function as well as script
    '''
    # check parameters
    if area is None:
        area == '90/-180/-90/180'
    else:
        if '/' in area:
            sarea = area.split('/')
            assert len(sarea) == 4, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            assert float(sarea[0]) <= 90.,   'area must be in 90/-180/-90/180'
            assert float(sarea[1]) >= -180., 'area must be in 90/-180/-90/180'
            assert float(sarea[2]) >= -90.,  'area must be in 90/-180/-90/180'
            assert float(sarea[3]) <= 180.,  'area must be in 90/-180/-90/180'
            assert float(sarea[0]) >= float(sarea[2]), 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            assert float(sarea[1]) <= float(sarea[3]), 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        else:
            assert ',' in area, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            sarea = area.split(',')
            assert len(sarea) == 2, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            lat, lon = sarea
            area = lat+'/'+lon+'/'+lat+'/'+lon

    if years is None:
        import time as ptime
        yearstart = 1979
        curryear  = int(ptime.asctime().split()[-1])
        yearend   = curryear - 1
    else:
        assert len(years) == 2, 'years must be iterable with two elements.'
        yearstart, yearend = [ int(i) for i in years ]
    # The time period to analyse. Valid formats:
    #    A single date as "2010-01-01"
    #    Multiple dates as "2010-01-01/2010-01-15/2010-01-30"
    #    A time period as "2010-01-01/to/2015-12-31"
    startdate = '{:04d}-{:02d}-{:02d}'.format(yearstart, 1, 1)
    lastdate  = '{:04d}-{:02d}-{:02d}'.format(yearend, 12, 31)
    dates = (startdate+"/to/"+lastdate)
    
    # make output directory
    if not os.path.exists(path): os.makedirs(path)

    # output files
    targetan = path + '/' + 'era-interim_an_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
    targetfc = path + '/' + 'era-interim_fc_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)

    # Check existing files
    files = glob.glob(path+'/*.nc*')
    files = [ ff for ff in files if ('_an_' in ff) or ('_fc_' in ff) ]
    hasan = False
    hasfc = False
    if (len(files) > 0) and (not override):
        ilat1, ilon1, ilat2, ilon2 = [ float(i) for i in area.split('/') ]
        for ff in files:
            fs = ff.split('_')
            typ  = fs[-6]
            lat1 = float(fs[-5])
            lon1 = float(fs[-4])
            lat2 = float(fs[-3])
            lon2 = float(fs[-2])
            fs = fs[-1]
            yrs = fs[:fs.rfind('.')]
            yr1, yr2 = [ int(i) for i in yrs.split('-') ]
            if ((ilat1 <= lat1) and (ilat2 >= lat2) and (ilon1 >= lon1) and (ilon2 <= lon2) and
                (yearstart >= yr1 ) and (yearend <= yr2)):
                if typ == 'an':
                    hasan = True
                    targetan = ff
                elif typ == 'fc':
                    hasfc = True
                    targetfc = ff

    # List of meteorological parameter(s) to analyse.
    #    Examples:
    #        param = "10u/10v/2t/tp"
    #        param = "165/166/167/228"
    #    The first example uses the parameters' short names, the second uses IDs. Both are valid and yield the same result.
    #    The parameters in the example are:
    #        10 metre U wind component, 10u, ID 165
    #        10 metre V wind component, 10v, ID 166
    #        2 metre temperature, 2t, ID 167
    #        Total precipitation, tp, ID 228
    # Look up parameter names, short names and IDs:
    #     http://apps.ecmwf.int/codes/grib/param-db
    # See which parameters are available in ERA-Interim:
    #     http://www.ecmwf.int/en/elibrary/8174-era-interim-archive-version-20
    #
    # MuSICA needs the following forcing data:
    #     Name      | Description                | unit
    #     windspeed | Wind speed                 | m/s
    #        precip | Precipitation              | mm/dt
    #          tair | Air temperature            | C
    #         rhair | Air relative humidity      | %
    #         psurf | Atmospheric pressure       | hPa
    #        swdown | Incoming solar radiation   | W/m2
    #        lwdown | Incoming thermal radiation | W/m2
    #        co2air | CO2 mixing ratio           | ppmv
    # Optional
    #       winddir | Wind direction             | deg
    #
    # ERA variables to calculate MuSICA variable
    #     ----------------------------------------------------------------------------------------------
    #     ERA short name | ERA description                                        | ERA unit | GRIB code
    #     ----------------------------------------------------------------------------------------------
    #                par | Surface photosynthetically active radiation            | W/m2/s   | 58
    #                  z | Surface geopotential                                   | m2/s2    | 129
    #                 sp | Surface pressure                                       | Pa       | 134
    #                lsp | Large-scale precipiation                               | m        | 142
    #                 cp | Convective precipiation                                | m        | 143
    #                 sf | Snowfall                                               | m        | 144
    #                msl | Mean sea level pressure                                | Pa       | 151
    #                tcc | Total cloud cover                                      | (0-1)    | 164
    #                10u | 10 metre U wind component                              | m/s      | 165
    #                10v | 10 metre V wind component                              | m/s      | 166
    #                 2t | 2 metre temperature                                    | K        | 167
    #                 2d | 2 metre dewpoint temperature                           | K        | 168
    #               ssrd | Surface solar radiation downwards                      | W/m2 s   | 169
    #                lsm | Land/sea mask                                          | (0-1)    | 172
    #               strd | Surface thermal radiation downwards                    | W/m2 s   | 175
    #               mx2t | Maximum 2m temperature since last post-processing step | K        | 201
    #               mn2t | Minimum 2m temperature since last post-processing step | K        | 202
    #                 tp | Total precipitation                                    | m        | 228
    #     ----------------------------------------------------------------------------------------------
    # # Jerome
    # paramsan = "134/164/165/166/167/168"
    # paramsfc = "228"
    # # all
    # paramsan = "129/134/151/164/165/166/167/168/172"
    # paramsfc = "58/142/143/144/169/175/201/202/228"
    # minimal
    paramsan = "129/134/164/165/166/167/168"
    paramsfc = "169/175/228"

    # get analysis
    if not hasan:
        data_file_an = get_sfc_data("an", dates, "00:00:00/06:00:00/12:00:00/18:00:00", "0",
                                    paramsan, area, targetan)
    # get forecast
    if not hasfc:
        data_file_fc = get_sfc_data("fc", dates, "00:00:00/12:00:00", "3/6/9/12",
                                    paramsfc, area, targetfc)

    return targetan, targetfc

            
# --------------------------------------------------------------------
# Script
#

if __name__ == "__main__":

    import argparse

    area     = None
    years    = None
    path     = '.'
    override = False
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                       description='''Download ERA-Interim data suitable to produce MuSICA input data.''')
    parser.add_argument('-a', '--area', action='store', default=area, dest='area', metavar='area',
                        help='area format as either lat,lon or NorthLat/WestLon/SouthLat/EastLon '
                        '(default: global 90/-180/-90/180.')
    parser.add_argument('-y', '--years', action='store', default=years, dest='years', metavar='years',
                        help='years format is startyear,endyear (default: 1979,current-1).')
    parser.add_argument('-p', '--path', action='store', default=path, dest='path', metavar='path',
                        help='Output directory (default: current directory ".").')
    parser.add_argument('-o', '--override', action='store_true', default=override, dest='override',
                        help='Do not check that output file already exists that includes request. '
                        'Override existing output file (default: False).')

    args     = parser.parse_args()
    area     = args.area
    years    = args.years
    path     = args.path
    override = args.override

    del parser, args

    if years is not None:
        if ',' in years:
            syears = years.split(',')
            assert len(syears) == 2, 'years format is startyear,endyear.'
            lyears = [ int(i) for i in syears ]
        else:
            year = int(years)
            lyears = (year, year)
    else:
        lyears = years

    import time as ptime
    t1 = ptime.time()

    anfile, fcfile = get_era_interim(area=area, years=lyears, path=path, override=override)
    print('Files: ', anfile, fcfile)

    t2    = ptime.time()
    strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    print('Time elapsed', strin)
