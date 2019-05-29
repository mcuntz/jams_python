#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
'''
    Download ERA5 data suitable to produce MuSICA input data.

    If override=False (default), the script checks if the data is already available
    in the local download directory (path).
    It expects files with the same naming convention than its own, i.e.
         path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(year)
    or
         path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
    Filenames can also end on .nc?.
    It checks if the area and years are included in a file in the download directory by checking ONLY filenames.

    Be aware that the request is processed once it is queued even if you abort this script.
    Queued request can be deleted after login at:
       https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
    You can (re-)download the data from this site later as well.

    Written  Matthias Cuntz, Jan 2019 - from get_era_interim.py


    --------------------------------------------------------
    usage: get_era5.py [-h] [-a area] [-y years] [-p path] [-o]

    Download ERA5 data suitable to produce MuSICA input data.

    optional arguments:
      -h, --help            show this help message and exit
      -a area, --area area  area format as either lat,lon or
                            NorthLat/WestLon/SouthLat/EastLon (default: global 90/-180/-90/180).
                            Min and Max of Lat and Lon must be different by at least 0.25 degree.
      -y years, --years years
                            years format is either single year or startyear,endyear
                            (default: 1979,current-1).
      -p path, --path path  Output directory (default: current directory ".").
      -o, --override        Do not check that output file already exists that
                            includes request. Override existing output file
                            (default: False).


    Example
        # Hesse
        python get_era5.py -a 48.6742166667,7.06461666667 -y 1995,2017 -p era


    --------------------------------------------------------

    Script was originally adapted from CDS Web API
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
    Further help comes from
        https://confluence.ecmwf.int/display/CKB/C3S+ERA5%3A+Web+API+to+CDS+API
'''
import os
import glob

__all__ = ['get_era5']


# --------------------------------------------------------------------
# Retrieval function
#

def get_era5_single_level5(variables, date, time, area, target, grid=None):

    import cdsapi
    # Disable: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate
    #     verification is strongly advised.
    try:
        import urllib3
        urllib3.disable_warnings()
    except:
        import requests
        from requests.packages.urllib3.exceptions import InsecureRequestWarning
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

    server = cdsapi.Client()

    request = {
        'product_type': 'reanalysis',
        'format':       'netcdf',     # grib or netcdf
        'variable':     variables,    # ['10m_v_component_of_wind', '2m_dewpoint_temperature', ...]
        'date':         date,         # "2010-01-01", "2010-01-01/2015-12-31"
        'time':         time,         # "[00:00, 01:00, 02:00, ..., 23:00]
        'area':         area,         # # North, West, South, East within "90/-180/-90/180" or [90, -180, -90, 180]
        }

    if grid is not None:
        # 'grid': [1.0,1.0]
        request['grid'] = grid        # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude). Default: 0.25 x 0.25

    # print('Request: ', request)
    server.retrieve('reanalysis-era5-single-levels', request, target)

    return target

def get_era5_single_level5_datestarget(variables, time, area, datestarget, grid=None):
    date, target = datestarget
    return get_era5_single_level5(variables, date, time, area, target, grid)

# --------------------------------------------------------------------
# Main routine
#

def get_era5(area=None, years=None, path='.', override=False):
    '''
        Download ERA5 data suitable to produce MuSICA input data.
    
        If override=False, checks if the data is already available in the download directory.
        It expects files with the same naming convention than its own, i.e.
            path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(year)
        or
            path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
        Filenames can also end on .nc?.
        It checks if the area and years are included in a file by checking ONLY filenames.

        Be aware that the request is processed once it is queued even if you abort this function.
        Queued request can be deleted after login at:
           https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
        You can (re-)download the data from this site later as well.


        Definition
        ----------
        def get_era5(area=None, years=None, path='.', override=False):


        Parameters
        ----------
        area          string (default: '90/-180/-90/180')
                      Area. Can be one point given as lat,lon
                      or box as NorthLat/WestLon/SouthLat/EastLon. Default is the globe.
                      Min and Max of Lat and Lon of box must be different by at least 0.25 degree.
        years         Scalar or iterable of len(2) [default: (1979,current_year-1)]
                      Year or year range.
        path          string (default: '.')
                      Output path
        override      boolean (default: False)
                      If True, check if data is already available in path.
                      Expect file with the naming convention
                          path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(year)
                      or
                          path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
                      Check if area and years are included in a file in path by checking ONLY filenames.


        Ouput
        -----
        Returns file names of the output files in path, either the newly written files
        or the files that contains the output requested:
            for yy in years: path+'/'+'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(yy)


        Restrictions
        ------------
        Existing files will only be checked by filename not by content.


        Examples
        --------
        >>> area  = '48/7/47/8'
        >>> years = (1995,2017)
        >>> get_era5(area=area, years=years, path='.')
        >>> import os
        >>> file1 = 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(years[0])
        >>> if not os.path.exists(file1): print('No file: ', file1)


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de

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

        Copyright 2019 Matthias Cuntz


        History
        -------
        Written  Matthias Cuntz, Jan 2019 - from get_era_interim.py
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
            assert float(sarea[0]) >= float(sarea[2])+0.25, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            assert float(sarea[1])+0.25 <= float(sarea[3]), 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        else:
            assert ',' in area, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            sarea = area.split(',')
            assert len(sarea) == 2, 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
            lat, lon = sarea
            # single point does not work. Needs to encompass an actual grid point.
            area = str(float(lat)-0.125)+'/'+str(float(lon)-0.125)+'/'+str(float(lat)+0.125)+'/'+str(float(lon)+0.125)

    if years is None:
        import time as ptime
        yearstart = 1979
        curryear  = int(ptime.asctime().split()[-1])
        yearend   = curryear - 1
    else:
        assert len(years) <= 2, 'years must be scaler or iterable with two elements.'
        if len(years) == 1:
            yearstart, yearend = years, years
        else: 
            yearstart, yearend = [ int(i) for i in years ]
    
    # make output directory
    if not os.path.exists(path): os.makedirs(path)

    # output files, single year files
    targetera5 = [ path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(yy) for yy in range(yearstart,yearend+1) ]
    targetyrs = targetera5[:]

    # Check existing files
    files = glob.glob(path+'/*.nc*')
    files = [ ff for ff in files if os.path.basename(ff).startswith('era5_') ]
    hasyrs = list(range(yearstart,yearend+1))
    if (len(files) > 0) and (not override):
        ilat1, ilon1, ilat2, ilon2 = [ float(i) for i in area.split('/') ]
        for ff in files:
            fs = ff.split('_')
            lat1 = float(fs[-5])
            lon1 = float(fs[-4])
            lat2 = float(fs[-3])
            lon2 = float(fs[-2])
            if ((ilat1 <= lat1) and (ilat2 >= lat2) and (ilon1 >= lon1) and (ilon2 <= lon2)):
                fs = fs[-1]
                yrs = fs[:fs.rfind('.')]
                if '-' in yrs:
                    # check merged file
                    yr1, yr2 = [ int(i) for i in yrs.split('-') ]
                else:
                    # check single year file
                    yr2 = yr1 = int(yrs)
                for yr in range(yr1,yr2+1):
                    if (yr >= yearstart) and (yr <= yearend):
                        if yr in hasyrs:
                            ii = hasyrs.index(yr)
                            hasyrs.remove(yr)
                            _ = targetyrs.pop(ii)

    # List of meteorological parameter(s) to analyse.
    #    Examples:
    #        variables = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', ...]
    #        param = "165/166/167/228"
    #    The first example uses the variable names, the second uses parameter IDs.
    #    Both are valid and yield the same result.
    #    The parameters in the example are:
    #        10 metre U wind component, 10u, ID 165
    #        10 metre V wind component, 10v, ID 166
    #        2 metre temperature, 2t, ID 167
    #        Total precipitation, tp, ID 228
    #
    # Look up variable names by selecting variables and examining the API:
    #     https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
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
    # Standard
    variables  = ['10m_u_component_of_wind', '10m_v_component_of_wind',
                  'total_precipitation', 'snowfall',
                  '2m_temperature',
                  '2m_dewpoint_temperature',
                  'surface_pressure',
                  'surface_solar_radiation_downwards',
                  'surface_thermal_radiation_downwards']
    # # Test alternative variables
    # variables  = ['mean_total_precipitation_rate', 'mean_snowfall_rate',
    #               'mean_surface_downward_short_wave_radiation_flux',
    #               'mean_surface_downward_long_wave_radiation_flux',
    #               'mean_top_downward_short_wave_radiation_flux', 'total_cloud_cover']
    # # Extra information
    # variables += ['land_sea_mask', 'orography', 'slope_of_sub_gridscale_orography',
    #               'soil_type',
    #               'type_of_high_vegetation', 'type_of_low_vegetation'
    #               'high_vegetation_cover', 'low_vegetation_cover', 
    #               'leaf_area_index_high_vegetation', 'leaf_area_index_low_vegetation']

    times = ['00:00', '01:00', '02:00',
             '03:00', '04:00', '05:00',
             '06:00', '07:00', '08:00',
             '09:00', '10:00', '11:00',
             '12:00', '13:00', '14:00',
             '15:00', '16:00', '17:00',
             '18:00', '19:00', '20:00',
             '21:00', '22:00', '23:00']
        
    # get reanalysis
    if len(hasyrs) > 0:
        if len(hasyrs) == 1:
            yy = 0
            target = targetyrs[yy]
            # The time period to analyse. Valid formats:
            #    A single date as "2010-01-01"
            #    A time period as "2010-01-01/2015-12-31"
            startdate = '{:04d}-{:02d}-{:02d}'.format(hasyrs[yy], 1, 1)
            lastdate  = '{:04d}-{:02d}-{:02d}'.format(hasyrs[yy], 12, 31)
            dates = (startdate+"/"+lastdate)
            print('Retrieve ', dates, target)
            data_file_era5 = get_era5_single_level5(variables, dates, times, area, target)
        else:
            # Do one by year because CDS allows only download of 100000 items at a time.
            # Items are time steps and variables but a map is one item.
            # So 10 variables, 24 hours a day per year are: 10*24*365 = 87600.
            datestarget = []
            for yy in range(len(hasyrs)):
                target = targetyrs[yy]
                # The time period to analyse. Valid formats:
                #    A single date as "2010-01-01"
                #    A time period as "2010-01-01/2015-12-31"
                startdate = '{:04d}-{:02d}-{:02d}'.format(hasyrs[yy], 1, 1)
                lastdate  = '{:04d}-{:02d}-{:02d}'.format(hasyrs[yy], 12, 31)
                dates = (startdate+"/"+lastdate)
                datestarget += [(dates, target)]
                print('Pool ', dates, target)
            from functools import partial
            getit = partial(get_era5_single_level5_datestarget, variables, times, area)
            from multiprocessing import Pool
            pool = Pool(processes=len(hasyrs))
            pool.map(getit, datestarget)

    return targetera5

            
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
                                       description='''Download ERA5 data suitable to produce MuSICA input data.''')
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

    era5file = get_era5(area=area, years=lyears, path=path, override=override)
    print('File: ', era5file)

    t2    = ptime.time()
    strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    print('Time elapsed', strin)

'''
Examples from different sources:


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[
            '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','high_vegetation_cover','land_sea_mask',
            'leaf_area_index_high_vegetation','leaf_area_index_low_vegetation','low_vegetation_cover',
            'mean_snowfall_rate','mean_surface_downward_long_wave_radiation_flux',
            'mean_surface_downward_short_wave_radiation_flux',
            'mean_top_downward_short_wave_radiation_flux','mean_total_precipitation_rate','orography',
            'slope_of_sub_gridscale_orography','snowfall','soil_type',
            'surface_pressure','surface_solar_radiation_downwards','surface_thermal_radiation_downwards',
            'total_cloud_cover','total_precipitation','type_of_high_vegetation',
            'type_of_low_vegetation'
        ],
        'year':[
            '1979','1980','1981',
            '1982','1983','1984',
            '1985','1986','1987',
            '1988','1989','1990',
            '1991','1992','1993',
            '1994','1995','1996',
            '1997','1998','1999',
            '2000','2001','2002',
            '2003','2004','2005',
            '2006','2007','2008',
            '2009','2010','2011',
            '2012','2013','2014',
            '2015','2016','2017',
            '2018'
        ],
        'month':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ]
    },
    'download.nc')


    c.retrieve('reanalysis-era5-pressure-levels', {
           "variable": "temperature",
           "pressure_level": "1000",
           "product_type": "reanalysis",
           "date": "2017-12-01/2017-12-31",
           "time": "12:00",
           "area": "75/-15/28/45",
           "format": "grib",
       }, 'download.grib')


    c.retrieve('reanalysis-era5-pressure-levels', {
        'variable'      : 'temperature',
        'pressure_level': '1000',
        'product_type'  : 'reanalysis',
        'date'          : '2008-01-01',
        'area'          : [60, -10, 50, 2], # North, West, South, East. Default: global
        'grid'          : [1.0, 1.0], # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude). Default: 0.25 x 0.25
        'time'          : '12:00',
        'format'        : 'netcdf' # Supported format: grib and netcdf. Default: grib
        }, 'test.nc')

'''
