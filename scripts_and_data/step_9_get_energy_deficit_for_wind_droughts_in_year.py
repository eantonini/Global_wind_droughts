#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script calculates the energy deficit for the wind droughts in a given year.
'''

import numpy as np
import xarray as xr
import pandas as pd
import sys

import settings as settings
import utilities_for_data_processing as utilities


# Define the year to be processed
year = int(sys.argv[1])

if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_9_'+str(year)
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_9_'+str(year)

utilities.write_to_log_file(log_file_name, 'Year '+str(year)+'\n\n', new_file=True, write_time=False)
utilities.write_to_log_file(log_file_name, 'Task started\n\n')

# Load the climate data for the year
wind_power_density_time_series = xr.open_dataarray(settings.climate_data_directory+'/{:d}-wind_power_density_time_series.nc'.format(year), chunks=settings.chunks)
min_annual_mean_wind_power_density = xr.open_dataarray(settings.result_directory+'/min_annual_mean_wind_power_density.nc')
climatological_wind_power_density_time_series = xr.open_dataarray(settings.result_directory+'/climatological_wind_power_density_time_series.nc', chunks=settings.chunks_climatological)

# Define the time series of climatological wind power density, accounting for leap and non-leap years
hour_of_non_leap_year = list(np.arange((31+28)*24)) + list(np.arange((31+29)*24,366*24))
if year%4 != 0: # If the year is not a leap year
    climatological_wind_power_density_time_series = climatological_wind_power_density_time_series.loc[climatological_wind_power_density_time_series['hour of year'].isin(hour_of_non_leap_year)]
climatological_wind_power_density_time_series['hour of year'] = pd.date_range(str(year), str(year+1), freq='H')[:-1]
climatological_wind_power_density_time_series = climatological_wind_power_density_time_series.rename({'hour of year': 'time'})

# Calculate the wind droughts
energy_deficit_for_wind_droughts = utilities.calc_energy_deficit(wind_power_density_time_series, min_annual_mean_wind_power_density, 'climatological_target', year, climatological_time_series=climatological_wind_power_density_time_series)
utilities.save(settings.result_directory+'/{:d}-energy_deficit_for_wind_droughts.nc'.format(year), energy_deficit_for_wind_droughts)
utilities.write_to_log_file(log_file_name, 'Wind droughts calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
