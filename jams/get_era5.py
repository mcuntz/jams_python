#!/usr/bin/env python
'''
Download ERA5 or ERA5-Land data from Copernicus Climate Data Store.

If override=False (default), the script checks if the data is already available
in the local download directory (path).
It expects files with the same naming convention than its own, i.e.
     path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(year)
or
     path + '/' + 'era5_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
Filenames can also end on .nc?.
It checks if the area and years are included in a file in the download
directory by checking ONLY filenames.

Be aware that the request is processed once it is queued even if you abort
this script. Queued request can be deleted after login at:
   https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
You can (re-)download the data from this site later as well.

Written  Matthias Cuntz, Jan 2019 - from get_era_interim.py
Modified Stephan Thober, Mar 2020 - added era5-land capability
         Matthias Cuntz, Jun 2020 - return correct file list and not the list
                                    of projected filenames for chosen area
                                  - input variable list with MuSICA variables
                                    as default
                                  - finalised era5-land capability
         Matthias Cuntz, Feb 2021 - bug in single point -> lat1 > lat2
         Matthias Cuntz, Feb 2021 - added grib format
                                  - default download era5-land in grib format
         Matthias Cuntz, Sep 2021 - no default area: -a must be given.


--------------------------------------------------------
usage: get_era5.py [-h] [-a area] [-f format] [-o] [-p path]
                   [-r reanalysis_model] [-v variables] [-y years]

        Download ERA5 or ERA5-Land data from Copernicus Climate Data Store
        https://climate.copernicus.eu/climate-data-store.


optional arguments:
  -h, --help            show this help message and exit
  -a area, --area area  area format as either lat,lon or
                        NorthLat/WestLon/SouthLat/EastLon
                        (default: global 90/-180/-90/180).
  -f format, --format format
                        Output format netcdf or grib.
                        (default: netcdf if era5, grib if era5-land).
  -o, --override        Do not check that output file already exists that
                        includes request. Override existing output file
                        (default: False).
  -p path, --path path  Output directory (default: current directory ".").
  -r reanalysis_model, --reanalyis-model reanalysis_model
                        Reanalyis model to download, either: era5, era5-land or
                        era5land (default: era5)
  -v variables, --variables variables
                        Comma-separated variable list var1,var2,...
                        (default: forcing variables of ecosystem model MuSICA:
                        10m_u_component_of_wind,10m_v_component_of_wind,
                        2m_temperature,2m_dewpoint_temperature,total_precipitation,
                        snowfall,surface_pressure,
                        surface_solar_radiation_downwards,
                        surface_thermal_radiation_downwards).
  -y years, --years years
                        years format is startyear,endyear
                        (default: 1979,current-1 for ERA5 and
                        1981,current-1 for ERA5-Land.)


Examples
--------
    # Hesse
    python get_era5.py -r ERA5-Land -v 2m_temperature,2m_dewpoint_temperature \
                       -a 48.6742166667,7.06461666667 -y 1995,2017 -p era
    # Tumbarumba
    python get_era5.py -r ERA5-Land -v 10m_u_component_of_wind,10m_v_component_of_wind,2m_temperature,2m_dewpoint_temperature,total_precipitation,snowfall,surface_pressure,surface_solar_radiation_downwards,surface_thermal_radiation_downwards --area=-35.7833302,148.0166666 -y 2016,2018 -p era


--------------------------------------------------------

Script was originally adapted from CDS Web API
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
Further help comes from
    https://confluence.ecmwf.int/display/CKB/C3S+ERA5%3A+Web+API+to+CDS+API
'''
from __future__ import division, absolute_import, print_function
import os
import glob


__all__ = ['get_era5']


# --------------------------------------------------------------------
# Retrieval function
#

def _get_era5_single_level5(variables, date, time, area, target,
                            grid=None, reanalysis_model='era5',
                            output_format='netcdf'):
    """
    ToDo
    """
    import cdsapi

    # check reanalysis mode
    retrieve_name = {'era5':      'reanalysis-era5-single-levels',
                     'era5land':  'reanalysis-era5-land',
                     'era5-land': 'reanalysis-era5-land'}
    if reanalysis_model not in retrieve_name.keys():
        estr = 'Value of variable reanalysis_model'
        estr = estr + ' {:s} is not in retrieve_name.keys()'.format(
            reanalysis_model)
        raise ValueError(estr)

    # Disable: InsecureRequestWarning: Unverified HTTPS request is being made.
    #          Adding certificate verification is strongly advised.
    try:
        import urllib3
        urllib3.disable_warnings()
    except:
        import requests
        from requests.packages.urllib3.exceptions import InsecureRequestWarning
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

    server = cdsapi.Client()

    request = {
        'format':   output_format,  # grib or netcdf
        'variable': variables,      # ['10m_v_component_of_wind', ...]
        'date':     date,           # "2010-01-01", "2010-01-01/2015-12-31"
        'time':     time,           # "[00:00, 01:00, 02:00, ..., 23:00]
        'area':     area,           # North, West, South, East within
                                    # "90/-180/-90/180" or [90, -180, -90, 180]
    }

    if reanalysis_model == 'era5':
        request['product_type'] = 'reanalysis'

    if grid:
        # 'grid': [1.0,1.0]
        # Latitude/longitude grid: east-west (longitude)
        # and north-south resolution (latitude). Default: 0.25 x 0.25
        request['grid'] = grid

    # print('Request: ', request)
    server.retrieve(retrieve_name[reanalysis_model], request, target)

    return target


def _get_era5_single_level5_datestarget(variables, time, area, datestarget,
                                        grid=None, reanalysis_model='era5',
                                        output_format='netcdf'):
    """
    Wrapper function for `_get_era5_single_level5` to parallelise with
    Python's `multiprocessing.Pool`.
    """
    date, target = datestarget
    return _get_era5_single_level5(variables, date, time, area, target,
                                   grid, reanalysis_model=reanalysis_model,
                                   output_format=output_format)


# --------------------------------------------------------------------
# Main routine
#

def get_era5(vars=['10m_u_component_of_wind', '10m_v_component_of_wind',
                   '2m_temperature', '2m_dewpoint_temperature',
                   'total_precipitation', 'snowfall', 'surface_pressure',
                   'surface_solar_radiation_downwards',
                   'surface_thermal_radiation_downwards'],
             area='90/-180/-90/180', years=None, path='.',
             override=False, reanalysis_model='era5',
             output_format=''):
    """
    Download ERA5 or ERA5-Land data from Copernicus Climate Data Store.

    If override=False, checks if the data is already available in
    the download directory.
    It expects files with the same naming convention than its own, i.e.
        path + '/' + 'era5<land>_'+area.replace('/','_')+'_{:04d}.nc'.format(year)
    or
        path + '/' + 'era5<land>_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)
    Filenames can also end on .nc?.
    It checks if the area and years are included in a file by checking
    ONLY filenames.

    Be aware that the request is processed once it is queued even if you
    abort this function.
    Queued request can be deleted after login at:
       https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
    You can (re-)download the data from this site later as well.

    Parameters
    ----------
    vars : list of str, optional
        List of names of variables to retrieve.
        Default are the input variables for the ecosystem model MuSICA:

        `vars=['10m_u_component_of_wind', '10m_v_component_of_wind',
               'total_precipitation',
               'snowfall', '2m_temperature', '2m_dewpoint_temperature',
               'surface_pressure',
               'surface_solar_radiation_downwards',
               'surface_thermal_radiation_downwards']`

        Look up variable names for ERA5 by selecting variables and examining
        the API (Show API request):
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form

        and for ERA5-Land:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
    area : string, optional
        Area. Can be one point given as lat,lon
        or box as NorthLat/WestLon/SouthLat/EastLon.
        Default is the globe: '90/-180/-90/180'.

        Min and Max of Lat and Lon of box must be different by
        at least 0.25 degree for ERA5 and 0.1 degree for ERA5-Land.
    years : scalar or iterable of len(2), optional
        Year or year range (default: [1979,current_year-1] for ERA5
        and [1981,current_year-1] for ERA5-Land)
    path : string, optional
        Output path (default: '.')
    override : bool, optional
        If True, check if data is already available in path (default: False).

        Expects file with the naming convention

            path + '/' + 'era5<land>_'+area.replace('/','_')+'_{:04d}.nc'.format(year)

        or

            path + '/' + 'era5<land>_'+area.replace('/','_')+'_{:04d}-{:04d}.nc'.format(yearstart, yearend)

        Checks if `area` and `years` are included in a file in path by
        checking ONLY filenames.
    reanalysis_model : string, optional
        Reanalyis model to download. Can be era5, era5-land or era5land
        (default: 'era5').
        Output filenames will be adapted accordingly.
    output_format : string, optional
        File format of output file.
        Default is 'netcdf' if `reanalysis_model=='era5'` and
        'grib' if `reanalysis_model=='era5-land'.
        Output filenames will be suffixed by '.nc' or '.grb', respectively.

        Model output at ECMWF is stored in grib format. There are limitations
        on the conversion to netCDF using the current ECMWF infrastructure.
        One gets errors like 'One or more variable sizes violate format
        constraints.':
        https://confluence.ecmwf.int/display/CKB/Common+Error+Messages+for+CDS+Requests
        Download in grib format in this case and use the climate data operators
        'cdo' to convert to netCDF format:

        Check the file content to know latitude, longitude and year:

            `cdo -t ecmwf sinfov gribfile`

        Construct the netcdf name, e.g. era5-land_lat1_lon1_lat2_lon2_year.nc
        Convert to compressed netCDF4 format:

            `cdo -f nc4 -z zip -t ecmwf copy gribfile netcdffile`

        If you want to have lowercase variable names, this is in bash:

            `gribfile=ecmwf_output.grb`
            `# cdo -t ecmwf sinfov ${gribfile}`
            `year=$(echo $(cdo -s showyear ${gribfile}))`
            `netcdffile=era5-land_75_-15_28_45_${year}.nc`
            `vars=$(cdo -s -t ecmwf showname ${gribfile})`
            `copt="" ; for i in ${vars} ; do copt="${copt} -chname,${i},$(echo ${i} | tr A-Z a-z)" ; done`
            `cdo -L -f nc4 -z zip -t ecmwf ${copt} ${gribfile} ${netcdffile}`

    Returns
    -------
    list
        Returns filenames of the output files in path, either the newly written
        files or the files that contain the output requested (netcdf example):

            [ path + '/' + 'era5<land>_' + area.replace('/','_')
            + '_{:04d}.nc'.format(yy) for yy in years ]

    Warnings
    --------
    Existing files will only be checked by filename not by content.

    Examples
    --------
    area  = '48/7/47/8'
    years = (1995,2017)
    ofile = get_era5(area=area, years=years, path='.')
    file1 = 'era5_'+area.replace('/','_')+'_{:04d}.nc'.format(years[0])
    if file1 != ofile[0]: print(
        'Returned filename not recognised: ', ofile, ' Expected: ', file1)
    import os
    if not os.path.exists(ofile): print('No ofile: ', ofile)
    if not os.path.exists(file1): print('No file1: ', file1)

    License
    -------
    This file is part of the JAMS Python package,
    distributed under the MIT License.

    Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de


    History
    -------
    Written  Matthias Cuntz, Jan 2019 - from get_era_interim.py
    Modified Matthias Cuntz, Dec 2019 - default area (global) was not working:
                                        used == instead of =
             Stephan Thober, Mar 2020 - added optional reanalysis_model
                                        argument to download era5land
             Matthias Cuntz, Jun 2020 - allow name era5-land
                                      - return correct file list and not the
                                        list of projected filenames for chosen
                                        area
                                      - input variable list with MuSICA
                                        variables as default
                                      - finalised era5-land capability
                                      - use numpydoc format
             Matthias Cuntz, Feb 2021 - bug in single point -> lat1 > lat2
             Matthias Cuntz, Feb 2021 - added grib format
                                      - download era5-land in grib format
    """
    # Check parameters
    # reanalysis model
    rmodel = reanalysis_model.lower()
    estr = 'Reanalysis model must be era5, era5-land or era5land.'
    assert rmodel in ['era5', 'era5-land', 'era5land'], estr
    if rmodel == 'era5':
        resolution = 0.25
        minyear    = 1979
        if not output_format:
            output_format = 'netcdf'
    else:
        rmodel     = 'era5land'
        resolution = 0.1
        minyear    = 1981
        if not output_format:
            output_format = 'grib'
    output_format = output_format.lower()
    if output_format == 'netcdf':
        suffix = '.nc'
    elif output_format == 'grib':
        suffix = '.grb'
    else:
        estr  = 'Output format must be netcdf or grib. Given: '
        estr += output_format
        raise ValueError(estr)

    # area
    if '/' in area:
        sarea = area.split('/')
        estr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        assert len(sarea) == 4, estr
        estr = 'area must be in 90/-180/-90/180'
        assert float(sarea[0]) <= 90.,   estr
        assert float(sarea[1]) >= -180., estr
        assert float(sarea[2]) >= -90.,  estr
        assert float(sarea[3]) <= 180.,  estr
        estr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        assert float(sarea[0]) >= float(sarea[2])+resolution, estr
        estr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        assert float(sarea[1])+resolution <= float(sarea[3]), estr
    else:
        estr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        assert ',' in area, estr
        sarea = area.split(',')
        estr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
        assert len(sarea) == 2, estr
        lat, lon = sarea
        # single point does not work. Needs to encompass an actual grid point.
        area = (str(float(lat)+resolution/2.) + '/'
                + str(float(lon)-resolution/2.) + '/'
                + str(float(lat)-resolution/2.) + '/'
                + str(float(lon)+resolution/2.))

    # years
    if years:
        estr = 'years must be scaler or iterable with two elements.'
        assert len(years) <= 2, estr
        if len(years) == 1:
            yearstart = yearend = years
        else:
            yearstart, yearend = [ int(i) for i in years ]
        estr = 'Start year must be greater'
        estr = estr + ' {:d} for reanalysis model {:s}.'.format(
            minyear, reanalysis_model)
        assert yearstart >= minyear, estr
    else:
        import time as ptime
        yearstart = minyear
        curryear  = int(ptime.asctime().split()[-1])
        yearend   = curryear - 1

    # Make output directory
    if not os.path.exists(path):
        os.makedirs(path)

    # Single year output files
    hasyrs     = list(range(yearstart, yearend+1))
    targetera5 = [ path + '/' + rmodel + '_' + area.replace('/', '_')
                   + '_{:04d}'.format(yy) + suffix for yy in hasyrs ]
    targetyrs  = targetera5[:]

    # Check existing files
    files = glob.glob(path+'/*.nc*') + glob.glob(path+'/*.gr*b*')
    files = [ ff for ff in files
              if os.path.basename(ff).startswith(rmodel+'_') ]
    if (len(files) > 0) and (not override):
        ilat1, ilon1, ilat2, ilon2 = [ float(i) for i in area.split('/') ]
        for ff in files:
            fs = ff.split('_')
            lat1 = float(fs[-5])
            lon1 = float(fs[-4])
            lat2 = float(fs[-3])
            lon2 = float(fs[-2])
            if ((ilat1 <= lat1) and (ilat2 >= lat2) and
                (ilon1 >= lon1) and (ilon2 <= lon2)):
                fs  = fs[-1]
                yrs = fs[:fs.rfind('.')]
                if '-' in yrs:
                    # check merged file
                    yr1, yr2 = [ int(i) for i in yrs.split('-') ]
                else:
                    # check single year file
                    yr2 = yr1 = int(yrs)
                for yr in range(yr1, yr2+1):
                    if (yr >= yearstart) and (yr <= yearend):
                        if yr in hasyrs:
                            ii = hasyrs.index(yr)
                            hasyrs.remove(yr)
                            fyr = targetyrs.pop(ii)
                            targetera5[targetera5.index(fyr)] = ff

    # Select all times
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
            data_file_era5 = _get_era5_single_level5(
                vars, dates, times, area, target, reanalysis_model=rmodel,
                output_format=output_format)
        else:
            # Do one by year because CDS allows only download of 100000 items
            # at a time.
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
            getit = partial(_get_era5_single_level5_datestarget, vars, times,
                            area, reanalysis_model=rmodel,
                            output_format=output_format)
            from multiprocessing import Pool
            pool = Pool(processes=len(hasyrs))
            pool.map(getit, datestarget, None)

    return targetera5


# --------------------------------------------------------------------
# Script
#

if __name__ == "__main__":

    import argparse

    area             = ''
    oformat          = ''
    override         = False
    path             = '.'
    reanalysis_model = 'era5'
    varis            = '10m_u_component_of_wind'
    varis            = varis + ',' + '10m_v_component_of_wind'
    varis            = varis + ',' + '2m_temperature,2m_dewpoint_temperature'
    varis            = varis + ',' + 'total_precipitation,snowfall'
    varis            = varis + ',' + 'surface_pressure'
    varis            = varis + ',' + 'surface_solar_radiation_downwards'
    varis            = varis + ',' + 'surface_thermal_radiation_downwards'
    years            = None

    parser           = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        Download ERA5 or ERA5-Land data from Copernicus Climate Data Store
        https://climate.copernicus.eu/climate-data-store.
        ''')
    hstr = 'area format as either lat,lon or NorthLat/WestLon/SouthLat/EastLon'
    hstr = hstr + ', e.g. global 90/-180/-90/180, mandatory.'
    parser.add_argument('-a', '--area', action='store', default=area,
                        dest='area', metavar='area', help=hstr)
    hstr  = 'Output format netcdf or grib. (default: netcdf if era5,'
    hstr += ' grib if era5-land).'
    parser.add_argument('-f', '--format', action='store', default=oformat,
                        dest='oformat', metavar='format', help=hstr)
    hstr = 'Do not check that output file already exists that includes'
    hstr = hstr + ' request. Override existing output file (default: False).'
    parser.add_argument('-o', '--override', action='store_true',
                        default=override, dest='override',
                        help=hstr)
    hstr = 'Output directory (default: current directory ".").'
    parser.add_argument('-p', '--path', action='store', default=path,
                        dest='path', metavar='path', help=hstr)
    hstr = 'Reanalyis model to download, either: era5, era5-land or era5land'
    hstr = hstr + ' (default: era5)'
    parser.add_argument('-r', '--reanalyis-model', action='store',
                        default=reanalysis_model, dest='reanalyis_model',
                        metavar='reanalysis_model', help=hstr)
    hstr = 'Comma-separated variable list var1,var2,... (default: forcing'
    hstr = hstr + ' variables of ecosystem model MuSICA: '+varis+').'
    parser.add_argument('-v', '--variables', action='store', default=varis,
                        dest='varis', metavar='variables', help=hstr)
    hstr = 'years format is startyear,endyear (default: 1979,current-1 for'
    hstr = hstr + ' ERA5 and 1981,current-1 for ERA5-Land.)'
    parser.add_argument('-y', '--years', action='store', default=years,
                        dest='years', metavar='years', help=hstr)

    args             = parser.parse_args()
    area             = args.area
    oformat          = args.oformat
    override         = args.override
    path             = args.path
    reanalysis_model = args.reanalyis_model
    varis            = args.varis
    years            = args.years

    del parser, args

    if not area:
        print('area as either lat,lon or NorthLat/WestLon/SouthLat/EastLon')
        print('must be given with -a option, e.g. --area="90/-180/-90/180"')
        print('or a specific location like FR-Hes: -a 48.6742167,7.0646167.')
        import sys
        sys.exit()

    # years
    if years:
        if ',' in years:
            syears = years.split(',')
            assert len(syears) == 2, 'years format is startyear,endyear.'
            lyears = [ int(i) for i in syears ]
        else:
            year   = int(years)
            lyears = (year, year)
    else:
        lyears = None

    # variables
    if ',' in varis:
        svars = varis.split(',')
    else:
        svars = [varis]

    import time as ptime
    t1 = ptime.time()

    era5files = get_era5(vars=svars, area=area, years=lyears, path=path,
                         override=override, reanalysis_model=reanalysis_model,
                         output_format=oformat)
    print('Files: ', era5files)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1) > 60. else
             '[s]: {:d}'.format(int(t2-t1)))
    print('Time elapsed', strin)
