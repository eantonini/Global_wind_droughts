#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script divides each hour of the climatological wind power density time series by the number of years to obtain the final climatological wind power density time series.
'''

import numpy as np
import xarray as xr
import os

import settings as settings
import utilities_for_data_processing as utilities


if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_6'
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_6'

utilities.write_to_log_file(log_file_name, 'Task started\n\n', new_file=True)

# Load the climate data for all years
annual_mean_wind_power_density = xr.open_dataarray(settings.result_directory+'/annual_mean_wind_power_density.nc', chunks=settings.chunks)

# Calculate the start and end year
start_year = int(annual_mean_wind_power_density['year'][0])
end_year = int(annual_mean_wind_power_density['year'][-1])

# Calculate a series of weights for each hour where each weight is the reciprocal of the number of years
hour_of_leap_year = list(np.arange(366*24))
weights_of_non_leap_year = np.array([1]*((31+28)*24) + [0]*24 + [1]*((366-31-29)*24))
weights_of_leap_year = np.array([1]*(366*24))
weights = np.array([0]*(366*24))
for year in range(1979,2022+1):
    if year%4 == 0:
        weights = weights + weights_of_leap_year
    else:
        weights = weights + weights_of_non_leap_year
weights = xr.DataArray(data=1/weights, dims=['hour of year'], coords={'hour of year': hour_of_leap_year}, name='weights')

# Rename the climatological wind power density time series file to avoid overwriting it and load it
os.rename(settings.result_directory+'/climatological_wind_power_density_time_series.nc', settings.result_directory+'/original_climatological_wind_power_density_time_series.nc')
climatological_wind_power_density_time_series = xr.open_dataarray(settings.result_directory+'/original_climatological_wind_power_density_time_series.nc', chunks=settings.chunks)

# Calculate the final climatological wind power density time series by dividing it by the number of years
climatological_wind_power_density_time_series = (climatological_wind_power_density_time_series*weights).rename('Climatological 100m wind power density')
utilities.save(settings.result_directory+'/climatological_wind_power_density_time_series.nc', climatological_wind_power_density_time_series)
os.remove(settings.result_directory+'/original_climatological_wind_power_density_time_series.nc')
utilities.write_to_log_file(log_file_name, 'Climatological wind power calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
