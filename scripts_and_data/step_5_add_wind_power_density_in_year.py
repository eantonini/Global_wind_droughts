#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script adds the wind power density time series of a given year to the climatological wind power density time series.
'''

import numpy as np
import xarray as xr
import sys
import os

import settings as settings
import utilities_for_data_processing as utilities


# Read the year to be processed
year = int(sys.argv[1])

if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_5_'+str(year)
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_5_'+str(year)

utilities.write_to_log_file(log_file_name, 'Year '+str(year)+'\n\n', new_file=True, write_time=False)
utilities.write_to_log_file(log_file_name, 'Task started\n\n')

# Load the wind power density time series of the year
wind_power_density_time_series = xr.open_dataarray(settings.climate_data_directory+'/{:d}-wind_power_density_time_series.nc'.format(year), chunks=settings.chunks)

hour_of_leap_year = list(np.arange(366*24))
hour_of_non_leap_year = list(np.arange((31+28)*24)) + list(np.arange((31+29)*24,366*24))
if year%4 == 0:
    hour_of_year = hour_of_leap_year
else:
    hour_of_year = hour_of_non_leap_year

# Change the wind power density time series to make it consistent with the climatological wind power density time series
wind_power_density_time_series['time'] = hour_of_year
wind_power_density_time_series = wind_power_density_time_series.rename({'time': 'hour of year'})

# Rename the climatological wind power density time series file to avoid overwriting it and load it
os.rename(settings.result_directory+'/climatological_wind_power_density_time_series.nc', settings.result_directory+'/original_climatological_wind_power_density_time_series.nc')
original_climatological_wind_power_density_time_series = xr.open_dataarray(settings.result_directory+'/original_climatological_wind_power_density_time_series.nc', chunks=settings.chunks_climatological)

# Add the wind power density time series of the year to the climatological wind power density time series
with xr.set_options(arithmetic_join='outer'):
    climatological_wind_power_density_time_series = original_climatological_wind_power_density_time_series + wind_power_density_time_series
climatological_wind_power_density_time_series = climatological_wind_power_density_time_series.fillna(original_climatological_wind_power_density_time_series).rename('Climatological 100m wind power density')
utilities.save(settings.result_directory+'/climatological_wind_power_density_time_series.nc', climatological_wind_power_density_time_series)
os.remove(settings.result_directory+'/original_climatological_wind_power_density_time_series.nc')
utilities.write_to_log_file(log_file_name, 'Climatological wind power calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
