#!/usr/bin/env python

"""
EVENTUALLY MOVE THIS TO GEOCLAW!

GeoClaw util Module  `$CLAW/geoclaw/src/python/geoclaw/tidetools.py`

Module provides provides tide prediction functions.

:Functions:

 - retrieve_constituents - retrieves harmonic constituents from NOAA gauge station
 - retrieve_water_levels - retrieves observed water levels from NOAA's API
 - retrieve_predicted_tide - retrieves predicted tide from NOAA's API
 - datum_value - retrieves datum value for desired datum reference
 - predict_tide - predicts tide with Pytides
 - datetimes - prepares a collection of datetimes from beginning to end dates
 - detide - detides observed water levels
 - surge - predicts surge at NOAA gauge station
"""

from __future__ import absolute_import
from __future__ import print_function
from collections.abc import Iterable
from collections import OrderedDict, namedtuple
from scipy.optimize import leastsq, fsolve
from itertools import takewhile, count
from datetime import datetime, timedelta
from functools import reduce
from six.moves.urllib.parse import urlencode
from six.moves.urllib.request import urlopen

try:
    from itertools import izip, ifilter
except ImportError: #Python3
    izip = zip
    ifilter = filter
    
try:
    import requests
    import json
    import string
    import lxml.html as lh
    import pandas as pd
    import operator as op
    import numpy as np
    import io
    import os
    import os.path
except ImportError as e:
    print(e)

# Here's what we'd like to do:  (?)
#from clawpack.pytides.tide import Tide
#import clawpack.pytides.constituent as cons
#

# For now, hardwire in the path...
import sys, os
CLAW = os.environ['CLAW']  # path to Clawpack files
pathstr = os.path.join(CLAW, 'tidal-examples/pytides')
assert os.path.isdir(pathstr), '*** Need clawpack/tidal-examples/pytides ***'
print('Using Pytides from: %s' % pathstr)
if pathstr not in sys.path:
    sys.path.insert(0,pathstr)

from pytides.tide import Tide



d2r, r2d = np.pi/180.0, 180.0/np.pi

NOAA_API_URL = 'https://tidesandcurrents.noaa.gov/api/datagetter'
NOAA_home = 'https://tidesandcurrents.noaa.gov/harcon.html'

######################## Tide Prediction Functions ########################

def retrieve_constituents(station, time_zone='GMT', units='meters', cache_dir=None,
                         verbose=True):
                         
    """Fetch constituent data for given NOAA tide station.
    By default, retrieved data is cached in the geoclaw scratch directory
    located at:
        $CLAW/geoclaw/scratch
    :Required Arguments:
      - station (string): 7 character station ID
    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages
    :Returns:
      - constituent_dict (dictionary): dictionary of tidal constituents for NOAA gauge station
    """
                        
    def units_num(units):
        if (units == 'meters'):
            return 0
        elif (time_zone == 'feet'):
            return 1
    
    def time_zone_num(time_zone):
        if (time_zone == 'GMT'):
            return 0
        elif (time_zone == 'Local'):
            return 1
        
    def get_noaa_params(station, time_zone, units):
        noaa_params = {
            'unit': units_num(units),
            'timezone': time_zone_num(time_zone),
            'id': station
        }
        return noaa_params
    
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')    #### cache_dir
 
    def get_cache_path(station, time_zone, units):
        filename = '{}_{}_constituents'.format(time_zone, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, 'constituents', station, filename)
    
    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)
            
        component_array = pd.DataFrame(data)
        component_array.to_csv(cache_path, index=False)
            
    def parse(cache_path):
        # read data into structured array, skipping header row if present
        data = pd.read_csv(cache_path)
        component_array = pd.DataFrame(data)
        component_dict = component_array.to_dict(orient='list')
        return component_dict

    #Requests URL
    def fetch_data(station, time_zone, units):
        noaa_params = get_noaa_params(station, time_zone, units)
        cache_path = get_cache_path(station, time_zone, units)
        
        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached constituent data for station {}'.format(station))
            return parse(cache_path)
            

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching constituent data from NOAA for station {}'.format(station))
  
        #Forms URL
        url = '{}?{}'.format(NOAA_home, urlencode(noaa_params))

        page = requests.get(url)
        doc = lh.fromstring(page.content)
        tr_elements = doc.xpath('//tr')
        col = [((t.text_content(),[])) for t in tr_elements[0]]
        for j in range(1, len(tr_elements)):
            T, i = tr_elements[j], 0
            for t in T.iterchildren():
                col[i][1].append(t.text_content())
                i+=1
                
        constituent_dict = {title:column for (title,column) in col}

        # if there were no errors, then cache response
        save_to_cache(cache_path, constituent_dict)
                
        return constituent_dict


    try:
        constituents = fetch_data(station, time_zone, units)
         
    except:
        print('*** Fetching NOAA Constituents failed, returning None')
        constituents = None
    
    return constituents


def fetch_noaa_tide_data(station, begin_date, end_date, datum='MTL', time_zone='GMT', units='metric', cache_dir=None, verbose=True):
    """Fetch water levels and tide predictions at given NOAA tide station.
    The data is returned in 6 minute intervals between the specified begin and
    end dates/times.  A complete specification of the NOAA CO-OPS API for Data
    Retrieval used to fetch the data can be found at:
        https://tidesandcurrents.noaa.gov/api/
    By default, retrieved data is cached in the geoclaw scratch directory
    located at:
        $CLAW/geoclaw/scratch
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of retrieval
      - end_date (datetime): end of date/time range of retrieval
    :Optional Arguments:
      - datum (string): see NOAA API documentation for possible values
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages
    :Returns:
      - date_time (numpy.ndarray): times corresponding to retrieved data
      - water_level (numpy.ndarray): preliminary or verified water levels
      - prediction (numpy.ndarray): tide predictions
    """
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')

    def fetch(product, expected_header, col_idx, col_types):
        noaa_params = get_noaa_params(product)
        cache_path = get_cache_path(product)

        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached {} data for station {}'.format(
                    product, station))
            return parse(cache_path, col_idx, col_types, header=True)

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching {} data from NOAA for station {}'.format(
                product, station))
        full_url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        with urlopen(full_url) as response:
            text = response.read().decode('utf-8')
            with io.StringIO(text) as data:
                # ensure that received header is correct
                header = data.readline().strip()
                if header != expected_header or 'Error' in text:
                    # if not, response contains error message
                    raise ValueError(text)

                # if there were no errors, then cache response
                save_to_cache(cache_path, text)

                return parse(data, col_idx, col_types, header=False)

    def get_noaa_params(product):
        noaa_date_fmt = '%Y%m%d %H:%M'
        noaa_params = {
            'product': product,
            'application': 'Clawpack',
            'format': 'csv',
            'station': station,
            'begin_date': begin_date.strftime(noaa_date_fmt),
            'end_date': end_date.strftime(noaa_date_fmt),
            'time_zone': time_zone,
            'datum': datum,
            'units': units
        }
        return noaa_params

    def get_cache_path(product):
        cache_date_fmt = '%Y%m%d%H%M'
        dates = '{}_{}'.format(begin_date.strftime(cache_date_fmt),
                               end_date.strftime(cache_date_fmt))
        filename = '{}_{}_{}'.format(time_zone, datum, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, product, station, dates, filename)

    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        # write data to cache file
        with open(cache_path, 'w') as cache_file:
            cache_file.write(data)

    def parse(data, col_idx, col_types, header):
        # read data into structured array, skipping header row if present
        a = np.genfromtxt(data, usecols=col_idx, dtype=col_types,
                             skip_header=int(header), delimiter=',',
                             missing_values='')

        # return tuple of columns
        return tuple(a[col] for col in a.dtype.names)

    # only need first two columns of data; first column contains date/time,
    # and second column contains corresponding value
    col_idx = (0, 1)
    col_types = 'datetime64[m], float'

    # fetch water levels and tide predictions
    try:
        date_time, water_level = fetch(
            'water_level', 'Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality',
            col_idx, col_types)
    except:
        print('*** Fetching water_level failed, returning None')
        date_time = None
        water_level = None

    try:
        date_time2, prediction = fetch('predictions', 'Date Time, Prediction',
                                       col_idx, col_types)
        if date_time is None:
            date_time = date_time2

    except:
        print('*** Fetching prediction failed, returning None')
        date_time2 = None
        prediction = None

    # ensure that date/time ranges are the same
    if (date_time is not None) and (date_time2 is not None):
        if not np.array_equal(date_time, date_time2):
            raise ValueError('Received data for different times')

    return date_time, water_level, prediction
    
    
def datum_value(station, datum, time_zone='GMT', units='metric'):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - datum (string): MSL, MTL
    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
    :Returns:
      - datum_value (float): value for requested datum reference
    """
    def get_noaa_params(station, time_zone, units):
        noaa_params = {
            'product': 'datums',
            'units': units,
            'time_zone': time_zone,
            'station': station,
            'application': 'Clawpack',
            'format':'json'
        }
        return noaa_params
    
    #Scrapes MTL/MSL Datum Value
    def get_datum(station, time_zone, units):
        noaa_params = get_noaa_params(station, time_zone, units)
        url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        page_data = requests.get(url)
        data = page_data.json()['datums']
        datum_value = [d['v'] for d in data if d['n'] == datum]
        return datum_value
    
    try:
        datum_value = float(get_datum(station, time_zone, units)[0])
    except:
        print('*** Fetching datum value failed, returning None')
        datum_value = None
        
    return datum_value

     
def predict_tide(station, begin_date, end_date, datum='MTL', time_zone='GMT', units='meters'):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of prediction
      - end_date (datetime): end of date/time range of prediction
    :Optional Arguments:
      - datum (string): MTL for tide prediction
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
    :Returns:
      - heights (float): tide heights
    """
    #These are the NOAA constituents, in the order presented on NOAA's website.
    from pytides.constituent import noaa, _Z0
    
    constituents = [c for c in noaa if c != _Z0]
    noaa_values = retrieve_constituents(station, time_zone, units)
    noaa_amplitudes = [float(amplitude) for amplitude in noaa_values['Amplitude']]
    noaa_phases = [float(phases) for phases in noaa_values['Phase']] 
    
    #We can add a constant offset - set to MTL
    #    MTL = datum_value(args[0], 'MTL')
    desired_datum = datum_value(station, datum)
    MSL = datum_value(station, 'MSL')
    offset = MSL - desired_datum
    constituents.append(_Z0)
    noaa_phases.append(0)
    noaa_amplitudes.append(offset)
       
    #Build the model
    assert(len(constituents) == len(noaa_phases) == len(noaa_amplitudes))
    model = np.zeros(len(constituents), dtype = Tide.dtype)
    model['constituent'] = constituents
    model['amplitude'] = noaa_amplitudes
    model['phase'] = noaa_phases
    tide = Tide(model = model, radians = False)
    
    #Time Calculations
    delta = (end_date-begin_date)/timedelta(hours=1) + .1
    times = Tide._times(begin_date, np.arange(0, delta, .1))
    
    #Height Calculations
    heights_arrays = [tide.at([times[i]]) for i in range(len(times))]
    heights = [val for sublist in heights_arrays for val in sublist]
 
    return heights


def datetimes(begin_date, end_date):
    #Time Calculations
    delta = (end_date-begin_date)/timedelta(hours=1) + .1
    times = Tide._times(begin_date, np.arange(0, delta, .1))
    return times


def detide(NOAA_observed_water_level, predicted_tide):
    # NOAA observed water level - predicted tide 
    return [(NOAA_observed_water_level[i] - predicted_tide[i]) for i in range(len(NOAA_observed_water_level))]

#Surge Implementation
def surge(station, beg_date, end_date, landfall_date):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of prediction
      - end_date (datetime): end of date/time range of prediction
      - landfall_date (datetime): approximate time of landfall for reference
    :Optional Arguments:
      - datum (string): MTL for tide prediction and retrieval
      - time_zone (string): see NOAA API documentation for possible values
    :Returns:
      - times (float): times with landfall event as reference
      - surge (float): surge heights
    """
    predicted_tide = predict_tide(station, beg_date, end_date)
    NOAA_times, NOAA_observed_water_level, NOAA_predicted_tide = fetch_noaa_tide_data(station, beg_date, end_date)

    #detides NOAA observed water levels with predicted tide
    surge = detide(NOAA_observed_water_level, predicted_tide)
    #modifies NOAA times to datetimes
    times = [((pd.to_datetime(time).to_pydatetime())-landfall_date)/timedelta(days=1) for time in NOAA_times]
    
    return times, surge



# ======================
# From fetch_noaa_gauge_info.py


url_base = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations'

url_base_web = 'https://tidesandcurrents.noaa.gov/stationhome.html'

station = '9441102' # Westport

# Datums:

def fetch_datums(station, units='meters', verbose=False):
    
    """
    Fetch datums for CO-OPS station.
    Convert units if required.
    If verbose is True, print out the station info and datums.
    
    Returns:
        datums = list of dictionaries, one for each datum, with keys
                    'name', 'description', 'value'
                 To print the datum names:
                     print([d['name'] for d in datums])
                 values have been converted to specified units.
        datums_info = dictionary with these keys:
                 'accepted', 'superseded', 'epoch', 'units', 'OrthometricDatum',
                 'datums', 'LAT', 'LATdate', 'LATtime', 'HAT', 'HATdate',
                 'HATtime', 'min', 'mindate', 'mintime', 'max', 'maxdate', 
                 'maxtime', 'disclaimers', 'DatumAnalysisPeriod', 'NGSLink', 
                 'ctrlStation', 'self'
        datums_info['datums'] is the list from which
                 datums was created, after possibly changing units.
    """

    assert units in ['meters','feet'], '*** unrecognized units: %s' % units
    
    try:
        station_id = str(int(station))   # make sure it's a string of an int
    except:
        raise ValueError('Station cannot be converted to int')

    url_datums = '%s/%s/datums.json'  % (url_base, station_id)

    with urlopen(url_datums) as response:
        text = response.read().decode('utf-8')
        
    # convert to a dictionary:
    p = json.loads(text)  # dictionary
    datums = p['datums']  # list of dictionaries, with keys 'name', 'value'
    
    noaa_units = p['units']  # probably 'feet'
    
    if noaa_units == 'feet' and units == 'meters':
        for d in datums:
            d['value'] = d['value'] * 0.3048
            
    if noaa_units == 'meters' and units == 'feet':
        for d in datums:
            d['value'] = d['value'] / 0.3048

    if verbose:
        print('Datums info for station %s' % station_id)
        print('Should agree with info at \n    %s?id=%s' \
                % (url_base_web,station_id))
        for k in p.keys(): 
            if k != 'datums':
                print('%s:  %s' % (k.rjust(20), p[k])) 
        print('Datums will be returned in units = %s:' % units)
        for d in datums:
            print('%s: %11.3f %s' % (d['name'].rjust(7),d['value'],units))
            
    datums_info = p   # also return full dictionary
    return datums, datums_info


def fetch_harcon(station, units='meters', verbose=False):

    """
    Fetch harmonic constituents for CO-OPS station.
    Convert units of the amplitudes if required.
    If verbose is True, print out the info and constinuents.
    
    Returns:
        harcon = list of dictionaries, one for each constituent, with keys:
                    'number', 'name', 'description', 'amplitude',
                    'phase_GMT', 'phase_local', 'speed'
                 To print the constituent names:
                     print([h['name'] for h in harcon])
                 amplitudes have been converted to specified units.
        harcon_info = dictionary with keys:
                    'units', 'HarmonicConstiuents', 'self'
            harcon_info['HarmonicConstiuents'] is the list from which
                 harcon was created, after possibly changing units.
    """
    
    assert units in ['meters','feet'], '*** unrecognized units: %s' % units
    
    try:
        station_id = str(int(station))   # make sure it's a string of an int
    except:
        raise ValueError('Station cannot be converted to int')

    url_harcon = '%s/%s/harcon.json'  % (url_base, station_id)

    with urlopen(url_harcon) as response:
        text = response.read().decode('utf-8')
        
    # convert to a dictionary:
    p = json.loads(text)
    hc = p['HarmonicConstituents']  # list of dictionaries

    noaa_units = p['units']  # probably 'feet'
    
    if noaa_units == 'feet' and units == 'meters':
        for h in hc:
            h['amplitude'] = h['amplitude'] * 0.3048
            
    if noaa_units == 'meters' and units == 'feet':
        for h in hc:
            h['amplitude'] = h['amplitude'] / 0.3048

    if verbose:
        print('Harmonic constituent info for station %s' % station_id)
        print('Should agree with info at \n    %s?id=%s' \
                % (url_base_web,station_id))
        for k in p.keys(): 
            if k != 'HarmonicConstituents':
                print('%s:  %s' % (k.rjust(20), p[k])) 
                
        print('Harmonic Constituents will be returned in units = %s:' % units)

        print('  Name    amplitude (%s)     phase (GMT)        speed' % units)
        for h in hc:
            print('%s: %11.3f  %20.3f %18.6f' % (h['name'].rjust(5),h['amplitude'],
                    h['phase_GMT'],h['speed']))
            
    harcon = hc
    harcon_info = p   # also return full dictionary
    
    return harcon, harcon_info
    
def make_pytides_model(station):
    """
    Fetch harmonic constituents for station and return lists of 37
    amplitudes (meters) and phases (GMT) as required for a pytides model.
    """
    
    print('Fetching harmonic constituents...')
    harcon, harcon_info = fetch_harcon(station, units='meters', verbose=False)
    
    numbers = list(range(1,38))  # assume standard 37 constituents
    harcon_numbers = [h['number'] for h in harcon]
    # make sure there are the expected number and in the right order:
    assert harcon_numbers == numbers, '*** unexpected harcon_numbers = %s' \
        % harcon_numbers
        
    amplitudes = [h['amplitude'] for h in harcon]
    phases = [h['phase_GMT'] for h in harcon]
    
    return amplitudes, phases